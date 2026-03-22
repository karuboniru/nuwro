# Metropolis-Hastings Sampling in NuWro

## Motivation

### Standard method: test run + rejection sampling

NuWro's default event generation is a two-phase procedure. First, `test_events()` runs `number_of_test_events` weighted events through `Chooser` (which selects channels proportional to their integrated cross sections). Each generated weight is fed into `_procesy.add(k, weight, bias)`, which tracks the running maximum and average per channel. After the test run, `_procesy.set_weights_to_avg()` fixes the per-channel max weights to be used as rejection thresholds.

Then `real_events()` generates the actual output. For each channel k it loops, calling `makeevent` and testing `_procesy.accept(k, weight, bias)`. The acceptance condition is a standard von Neumann rejection: `weight / bias / max_weight > Uniform(0,1)`. Accepted events are stored with `e->weight = _procesy.total()` (the total cross section), so every stored event carries the same weight — the output is unweighted.

#### Max-weight exceedance and retroactive discard

`Mxw` (the running max weight per channel) is updated continuously even during `real_events`, not frozen after the test run. If a new event's `weight/bias` exceeds the current `Mxw`, `accept()` enters the overflow branch:

```
prevmax = Mxw
add(i, x, bias)          // Mxw updated to new higher value
Ready *= prevmax / Mxw   // scale down: retroactively "un-accept" some earlier events
Ready++                   // accept the current event unconditionally
return true
```

Scaling `Ready` down by `prevmax/Mxw` (which is < 1) corrects the acceptance count as if all previously accepted events had been drawn under the new, higher max weight. Because `Ready` drops below `desired`, the generation loop continues and produces additional replacement events — all generated under the correct updated `Mxw`.

The earlier events already written to the per-channel `.part` file are not physically deleted mid-run. They are discarded at the final merge step: `real_events` copies only the **last** `desired(k)` entries from each `.part` file into the output tree (for `mixed_order=0`) or draws randomly from the last `desired(k)` entries per channel (for `mixed_order=1`, the default). Events written before the max-weight update are therefore silently dropped from the output.

The limitation is that the test run must be long enough to find a reliable max weight for each channel. An underestimated max weight does not cause permanent bias — the overflow mechanism will self-correct — but it does waste CPU generating events that will ultimately be discarded.

### M-H method: no test run, adaptive

The Metropolis-Hastings mode (`real_events_mh`) eliminates the test run entirely. It learns channel cross sections on-the-fly via an adaptive proposal distribution, and uses M-H accept/reject to produce unweighted events from a Markov chain. The trade-off is that successive output events are correlated (the chain may stay in the same state across multiple output steps), whereas the standard method's output events are independent.

---

## Algorithm

### Core sampler: `Metropolis<T, Args...>` (`src/metropolis.h`)

A generic M-H sampler templated on the state type `T` and extra discriminator arguments `Args...`. For NuWro, `T = event` and `Args = size_t` (channel index).

**State**: the sampler holds a single current state (a heap-allocated `event`) and its weight.

**Proposal and acceptance** (`update_state`):

```
given proposed state x' with weight w(x', i'):
  ratio = w(x', i') / w(x_current, i_current)
  if ratio > 1 or ratio > Uniform(0,1):
      accept: replace current state with x'
  else:
      reject: discard x', keep current state
```

The weight function is set at construction time in `NuWro::NuWro()`:

```
w(e, index) = e->weight / bias / channel_sampling_weight[index]
```

where:
- `e->weight` — raw differential cross-section weight produced by the dynamics code
- `bias = e->in[0].t` — incoming neutrino energy; this corrects for the flux-weighted beam sampling (see below)
- `channel_sampling_weight[index]` — the current proposal probability for channel `index`

The stationary distribution of this M-H chain is proportional to `e->weight / bias / channel_sampling_weight[index]`. When the adaptive channel weights (described below) have converged, `channel_sampling_weight[index]` is proportional to the flux-averaged cross section of that channel, making the ratio `e->weight / bias / channel_sampling_weight[index]` approximately constant across events from the same channel — yielding a nearly uniform (flat-weight) sample.

### Bias and `dismode`

In `real_events_mh`, `dismode` is forced to `true`. This instructs the beam sampler (`_beam->shoot`) to draw neutrino energies from the flux-weighted distribution (the same path used for DIS events in the standard mode). As a consequence, the raw event weight `e->weight` already encodes a factor of the neutrino energy. Dividing by `bias = e->in[0].t` removes this factor and recovers the true cross-section weight at that energy.

### Per-output-event loop: `NuWro::get_event()`

For each event stored to the output file, the inner M-H loop runs `mh_sample_interval` (default 80) proposal steps:

```
for j in 0 .. mh_sample_interval:
    channel_index = sample_channel(frandom())   // draw from channel_sampling_weight
    e_proposed = makeevent(channel_index)        // generate kinematics for that channel
    if e_proposed->weight or bias is NaN:
        retry (don't count this step)
    accumulate channel_weight_sum[channel_index]    += e_proposed->weight / bias
    accumulate channel_weight_sum_fraction[channel_index] += 1 / bias
    sampler.update_state(e_proposed, channel_index) // M-H accept/reject

if use_weighted_channel:
    update_channel_sampling_weight()

return sampler.get_state()   // current accepted event (may be same as previous output)
```

`sample_channel` selects a channel index with probability `channel_sampling_weight[i]` (a normalized discrete distribution). Initially all channels have equal proposal weight `1 / N_channels`.

### Adaptive channel weights (`use_weighted_channel`)

After each inner loop, `update_channel_sampling_weight` recomputes:

```
channel_sampling_weight[i] = channel_weight_sum[i] / channel_weight_sum_fraction[i]
                             (= running estimate of <weight/bias> for channel i,
                                i.e. the flux-averaged cross section sigma_i)
then normalize so they sum to 1.
```

If any estimate is zero, NaN, or Inf, the update is skipped and equal weights are kept.

This is an **adaptive importance sampling** scheme: the proposal distribution is steered toward channels with larger cross sections, reducing variance in the acceptance ratio and improving mixing across channels.

### Output and cross-section reporting

After all events are generated, `real_events_mh` writes:

- A `.xsec` log file alongside the output ROOT file with per-channel average cross section, final event counts, and overall acceptance rate.
- A `TParameter<double>` named `"xsec"` stored in the ROOT `TTree`'s user info, holding the sum of flux-averaged cross sections across all enabled channels.

The overall acceptance rate (fraction of inner-loop proposals that changed the sampler state at least once) is printed and logged. Low acceptance indicates poor mixing — consider reducing `mh_sample_interval` or enabling adaptive channel weights.

---

## Channel index mapping

`initialize_dynamics_list()` builds the list of active channels from `params.txt` flags:

| Channel index (`dyn`) | Physics process | Parameter flag |
|---|---|---|
| 0 | QEL CC | `dyn_qel_cc` |
| 1 | QEL NC | `dyn_qel_nc` |
| 2 | RES CC | `dyn_res_cc` |
| 3 | RES NC | `dyn_res_nc` |
| 4 | DIS CC | `dyn_dis_cc` |
| 5 | DIS NC | `dyn_dis_nc` |
| 6 | COH CC | `dyn_coh_cc` |
| 7 | COH NC | `dyn_coh_nc` |
| 8 | MEC CC | `dyn_mec_cc` |
| 9 | MEC NC | `dyn_mec_nc` |
| 10 | HYP CC | `dyn_hyp_cc` |
| 12 | Leptonic | `dyn_lep` |
| 20 | QEL elastic | `dyn_qel_el` |
| 21 | RES elastic | `dyn_res_el` |

Only channels with their flag set to 1 are included. The M-H proposal and weight normalization operate solely over this active set.

---

## Parameters

All parameters are set in `params.txt` (or equivalent config file).

### `use_mh` (bool, default: `0`)

Set to `1` to activate M-H sampling mode. When enabled:
- `real_events_mh()` is called instead of `real_events()`.
- The standard `Chooser`-based channel selection (`refresh_dyn`) is skipped during `set()`.
- `dismode` is forced to `true` for all channels.
- Output events are unweighted (all carry the same implicit weight equal to the total cross section).

### `mh_sample_interval` (int, default: `80`)

Number of M-H proposal steps taken per output event. A larger value:
- Increases decorrelation between successive output events (better effective sample).
- Increases runtime proportionally.
- Allows the adaptive channel weights more data before the next output event.

A value of 1 means each output event is the immediate result of a single proposal — maximum correlation, minimum cost.

### `use_weighted_channel` (bool, default: `1`)

Controls adaptive channel proposals.

- `1` (enabled): after each inner loop, `channel_sampling_weight` is updated to match the running estimate of the flux-averaged cross section per channel. The M-H chain mixes across channels proportionally to their physical rates. This is the recommended setting.
- `0` (disabled): channel proposals remain uniform (`1/N_channels`) throughout the run. Channels with very small cross sections will be over-proposed, leading to many rejections. Only useful for debugging or when all channels have similar integrated cross sections.

---

## `nuwro2rootracker` integration

When converting M-H output to rootracker format, `nuwro2rootracker` reads the `"xsec"` `TParameter<double>` from the tree's user info. If present (nonzero), it overrides the per-event weight fields:

```
fEvtWght = xsec * 1e38      // total cross section, same for all events
fEvtXSec = xsec * coef      // scaled by flux/target normalization
```

This is correct for unweighted M-H output where every event represents the same physical rate. The code also writes pre-FSI particles (`e->out`) into the StdHep array with status code `2`, in addition to the post-FSI particles (`e->post`).

---

## Practical notes

- **Starting state**: the sampler's `current_weight` initializes to zero, so the first proposal is always accepted regardless of its weight.
- **NaN guard**: events with NaN weight or bias are silently retried (the step counter `j` is decremented), preventing the chain from getting stuck on pathological phase-space points.
- **Autocorrelation and thinning**: `mh_sample_interval` is a thinning interval — the chain advances that many proposal steps between consecutive output events, discarding the intermediate states. This reduces autocorrelation between stored events (successive output events are less likely to be identical). It is not a burn-in: there is no explicit warm-up phase. The chain starts cold (initial equal-probability channel weights, `current_weight = 0`), so the very first output events may reflect transient bias before the adaptive weights have converged. For precision studies, discarding the first O(100) output events is advisable.
- **Reproducibility**: the `Metropolis` sampler seeds its `mt19937` from `std::random_device`, independently of the NuWro `frandom` seed set by `random_seed` in `params.txt`. M-H runs are therefore not reproducible across invocations even with a fixed seed.
