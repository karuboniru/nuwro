#include "Analyser1.h"
#include "Interaction.h"
#include "TFile.h"
#include "TH1.h"
#include "TROOT.h"
#include "TTree.h"
#include "TVectorDfwd.h"
#include "args.h"
#include "beam.h"
#include "beamHist.h"
#include "coh.h"
#include "cohevent2.h"
#include "dis/disevent.h"
#include "dis/resevent2.h"
#include "dis/singlepion.h"
#include "e_el_event.h"
#include "e_spp_event.h"
#include "event1.h"
#include "ff.h"
#include "geomy.h"
#include "hist.h"
#include "hypevent.h"
#include "kaskada7.h"
#include "lepevent.h"
#include "mecevent.h"
#include "nucleusmaker.h"
#include "output.h"
#include "pauli.h"
#include "pdg.h"
#include "qelevent.h"
#include "rew/rewparams.h"
#include "sfevent.h"
#include "stdlib.h"
#include "target_mixer.h"
#include <cmath>
#include <iomanip>
#include <memory>
#include <sstream>
#include <unordered_map>
#include <vector>
// #include <TVectorD.h>
#include <TParameter.h>

extern double SPP[2][2][2][3][40];
// extern double sppweight;
extern "C" {
void shhpythiaitokay_(void);
void youcanspeaknowpythia_(void);
}

params *p1 = NULL;
string data_dir;
#include "nuwro_metropolis.h"

NuWro_metropolis::~NuWro_metropolis() {
  delete _mixer;
  delete _detector;
  delete _beam;
  delete _nucleus;
}

NuWro_metropolis::NuWro_metropolis()
    : sampler([this](event *e) {
        double bias = 1.;
        if (this->dismode) {
          bias = e->in[0].t;
        }
        return e->weight / bias;
      }) {
  _mixer = NULL;
  _detector = NULL;
  _beam = NULL;
  _nucleus = NULL;
  std::random_device rd;
  gen = std::mt19937(rd());
}

NuWro_metropolis::NuWro_metropolis(const char *filename) : NuWro_metropolis() {
  set_dir_by_env();
  init(filename);
}

void NuWro_metropolis ::set(params &par) {
  p = par;

  frandom_init(par.random_seed);

  // dismode = false;

  if (par.target_type == 1)
    _mixer = new target_mixer(par);
  _detector = make_detector(par);

  _beam = create_beam(par, _detector);

  _nucleus = make_nucleus(par);

  ff_configure(par);
  // refresh_dyn(par);
}

void NuWro_metropolis ::refresh_target(params &par) {
  delete _nucleus;
  _nucleus = make_nucleus(par);
}

// void NuWro_metropolis ::refresh_dyn(params &par) { _procesy.reset(par); }

geomy *NuWro_metropolis::make_detector(params &p) {
  if (p.target_type != 2)
    return NULL;

  if (p.geo_file.length()) {
    try {
      if (p.geo_d.norm2() == 0)
        return new geomy(p.geo_file, p.geo_name, p.geom_length_units,
                         p.geom_density_convert);
      else
        return new geomy(p.geo_file, p.geo_name, p.geom_length_units,
                         p.geom_density_convert, p.geo_volume, p.geo_d,
                         p.geo_o);
    } catch (...) {
      cerr << "Failed to make detector." << endl;
      exit(3);
    }
  } else {
    cerr << "Failed to make detector. Parameter geo_file must not be empty if "
            "target_type=2."
         << endl;
    exit(4);
    return NULL;
  }
}

void NuWro_metropolis::init(int argc, char **argv) {
  frame_top("Simulation parameters");

  // dismode=false;
  // dismode = true;
  set_dirs(argv[0]);
  a.read(argc, argv);
  p.read(a.input);
  p.read(a.params, "command line");
  p.list(cout);
  p.list(string(a.output) + ".par");
  p1 = &p;
  rew.init(p);
  _progress.open(a.progress);
  frandom_init(p.random_seed);

  frame_bottom();

  frame_top("Initialize the simulation");

  if (p.beam_test_only == 0 && p.kaskada_redo == 0)
    if (p.dyn_dis_nc or p.dyn_res_nc or p.dyn_dis_cc or p.dyn_res_cc) {
      cout << "     -> Calculating the single-pion functions..." << endl;
      singlepion(p);
    }
  if (p.kaskada_redo == 0) {
    cout << "     -> Building the target nuclei..." << endl;
    _nucleus = make_nucleus(p);
    if (p.target_type == 1)
      _mixer = new target_mixer(p);
    else
      _mixer = NULL;
    cout << "     -> Constructing the detector..." << endl;
    _detector = make_detector(p);

    cout << "     -> Creating the beam..." << endl;
    _beam = create_beam(p, _detector);
    if (_beam == NULL) {
      cerr << "No beam defined." << endl;
      exit(5);
    }
  }

  // load the input data
  cout << "     -> Loading external physical data..." << endl;
  input.initialize(p);
  input.load_data();

  cout << "     -> Configuring form factors..." << endl;
  ff_configure(p);

  cout << "     -> Extablishing the choice of dynamics..." << endl;
  // refresh_dyn(p);

  frame_bottom();
  {
    if (p.dyn_qel_cc)
      enabled_dyns.push_back(0);
    if (p.dyn_qel_nc)
      enabled_dyns.push_back(1);
    if (p.dyn_res_nc)
      enabled_dyns.push_back(3);
    if (p.dyn_res_cc)
      enabled_dyns.push_back(2);
    if (p.dyn_dis_cc)
      enabled_dyns.push_back(4);
    if (p.dyn_dis_nc)
      enabled_dyns.push_back(5);
    if (p.dyn_coh_cc)
      enabled_dyns.push_back(6);
    if (p.dyn_coh_nc)
      enabled_dyns.push_back(7);
    if (p.dyn_mec_cc)
      enabled_dyns.push_back(8);
    if (p.dyn_mec_nc)
      enabled_dyns.push_back(9);
    if (p.dyn_hyp_cc)
      enabled_dyns.push_back(10);
    if (p.dyn_lep)
      enabled_dyns.push_back(12);
    if (p.dyn_qel_el)
      enabled_dyns.push_back(20);
    if (p.dyn_res_el)
      enabled_dyns.push_back(21);
  }
}

void NuWro_metropolis::init(const char *filename) {
  frame_top("Simulation parameters");

  // dismode=false;
  // dismode = true;
  p.read(filename);
  // p.read(a.params, "command line");
  p.list(cout);
  p.list(string(a.output) + ".par");
  p1 = &p;
  rew.init(p);
  _progress.open(a.progress);
  frandom_init(p.random_seed);

  frame_bottom();

  frame_top("Initialize the simulation");

  if (p.beam_test_only == 0 && p.kaskada_redo == 0)
    if (p.dyn_dis_nc or p.dyn_res_nc or p.dyn_dis_cc or p.dyn_res_cc) {
      cout << "     -> Calculating the single-pion functions..." << endl;
      singlepion(p);
    }
  if (p.kaskada_redo == 0) {
    cout << "     -> Building the target nuclei..." << endl;
    _nucleus = make_nucleus(p);
    if (p.target_type == 1)
      _mixer = new target_mixer(p);
    else
      _mixer = NULL;
    cout << "     -> Constructing the detector..." << endl;
    _detector = make_detector(p);

    cout << "     -> Creating the beam..." << endl;
    _beam = create_beam(p, _detector);
    if (_beam == NULL) {
      cerr << "No beam defined." << endl;
      exit(5);
    }
  }

  // load the input data
  cout << "     -> Loading external physical data..." << endl;
  input.initialize(p);
  input.load_data();

  cout << "     -> Configuring form factors..." << endl;
  ff_configure(p);

  cout << "     -> Extablishing the choice of dynamics..." << endl;
  // refresh_dyn(p);

  frame_bottom();
  {
    if (p.dyn_qel_cc)
      enabled_dyns.push_back(0);
    if (p.dyn_qel_nc)
      enabled_dyns.push_back(1);
    if (p.dyn_res_nc)
      enabled_dyns.push_back(3);
    if (p.dyn_res_cc)
      enabled_dyns.push_back(2);
    if (p.dyn_dis_cc)
      enabled_dyns.push_back(4);
    if (p.dyn_dis_nc)
      enabled_dyns.push_back(5);
    if (p.dyn_coh_cc)
      enabled_dyns.push_back(6);
    if (p.dyn_coh_nc)
      enabled_dyns.push_back(7);
    if (p.dyn_mec_cc)
      enabled_dyns.push_back(8);
    if (p.dyn_mec_nc)
      enabled_dyns.push_back(9);
    if (p.dyn_hyp_cc)
      enabled_dyns.push_back(10);
    if (p.dyn_lep)
      enabled_dyns.push_back(12);
    if (p.dyn_qel_el)
      enabled_dyns.push_back(20);
    if (p.dyn_res_el)
      enabled_dyns.push_back(21);
  }
}

void NuWro_metropolis::makeevent(event *e, params &p) {
  static double max_norm = 0;
  particle nu;
  int dyn = e->dyn;
  if (_detector) {
    material mat;
    do {
      // nu = _beam->shoot(1 < dyn && dyn < 6 && dismode);
      nu = _beam->shoot(dismode);
      if (nu.travelled > 0 && p.beam_weighted == 0) {
        if (nu.travelled < frandom() * max_norm)
          continue;
        if (nu.travelled > max_norm)
          max_norm = nu.travelled;
      }
      nu.travelled = 0;
      nu.r = vec(nu.r) + p.beam_offset;
      if (nu.r.x == 0 && nu.r.y == 0 && nu.r.z == 0)
        mat = _detector->getpoint();
      else
        mat = _detector->getpoint(nu.p(), nu.r);

    } while (not(mat.Z + mat.N > 0 &&
                 mat.w_density >= frandom() * _detector->max_dens()));

    /// change nucleus
    e->r = mat.r;
    p.nucleus_p = mat.Z;
    p.nucleus_n = mat.N;
    p.nucleus_E_b = 0; // Use library value for E_b
    p.nucleus_kf = 0;  // Use library value for kF
                       //		cout<<mat.Z<<' '<<mat.N<<' '<<endl;
    if (mat.Z == 0 && mat.N == 0)
      throw "Empty isotope 00";
  } else {
    nu = _beam->shoot(dismode);
    nu.r = vec(nu.r) + p.beam_offset;
  }

  if (_detector or _mixer) // _nucleus not reusable
  {
    delete _nucleus;
    _nucleus = make_nucleus(p);
    // cout<<"make_nucleus "<<_nucleus->p<<" "<<_nucleus->n<<endl;
  } else
    _nucleus->reset();
  e->in.push_back(nu); // insert neutrino
  if (dyn < 6 || (dyn >= 10 && dyn < 12) || dyn == 20) {
    // insert target nucleon
    e->in.push_back(_nucleus->get_nucleon());
    e->in[0].r = e->in[1].r;
    assert(e->in[1] * e->in[1] > 0);
  } else if (dyn >= 12 && dyn < 14) {
    // insert target electron
    e->in.push_back(particle(PDG::pdg_e, PDG::mass_e));
  }

  e->weight = 0;
  if (nu.travelled > 0)
    e->norm = nu.travelled;
  // else e->norm remains 1;

  e->flag.cc = false;
  e->flag.nc = false;

  e->flag.qel = false;
  e->flag.res = false;
  e->flag.dis = false;
  e->flag.coh = false;
  e->flag.mec = false;
  e->flag.hyp = false;
  e->flag.lep = false;

  e->flag.anty = nu.pdg < 0;

  if (p.beam_test_only) {
    e->weight = 1;
    e->out.push_back(e->in[0]);
    return;
  }
  double factor = 1.0;

  if (p.cc_smoothing and dyn == 0) // only in qel_cc
  {
    if (e->in[0].pdg > 0) {
      factor = _nucleus->frac_neutron();
      e->in[1].set_neutron();
    } else {
      factor = _nucleus->frac_proton();
      e->in[1].set_proton();
    }
  }
  e->par = p;

  if ( // (anty)-neutrino interaction
      abs(nu.pdg) == 12 or abs(nu.pdg) == 14 or abs(nu.pdg) == 16)
    switch (dyn) {
    case 0:
      e->flag.qel = e->flag.cc = true;
      if (p.dyn_qel_cc) // qel cc
      {
        if (p.sf_method > 0 and has_sf(*_nucleus, p.sf_method))
          sfevent(p, *e, *_nucleus);
        else
          qelevent1(p, *e, *_nucleus, false);
      }
      break;
    case 1:
      e->flag.qel = e->flag.nc = true;
      if (p.dyn_qel_nc) // qel nc
      {
        if (p.sf_method > 0 and has_sf(*_nucleus, p.sf_method))
          sfevent(p, *e, *_nucleus);
        else
          qelevent1(p, *e, *_nucleus, true);
      }
      break;
    case 2:
      e->flag.res = e->flag.cc = true;
      if (p.dyn_res_cc) // res cc
      {
        resevent2(p, *e, *_nucleus, true);
        if (p.pauli_blocking)
          mypauli_spp(*e, *_nucleus);
      }
      break;
    case 3:
      e->flag.res = e->flag.nc = true;
      if (p.dyn_res_nc) // res nc
      {
        resevent2(p, *e, *_nucleus, false);
        if (p.pauli_blocking)
          mypauli_spp(*e, *_nucleus);
      }
      break;
    case 4:
      e->flag.dis = e->flag.cc = true;
      if (p.dyn_dis_cc) // dis cc
      {
        disevent(p, *e, *_nucleus, true);
        if (p.pauli_blocking)
          mypauli_spp(*e, *_nucleus);
      }
      break;
    case 5:
      e->flag.dis = e->flag.nc = true;
      if (p.dyn_dis_nc) // dis nc
      {
        disevent(p, *e, *_nucleus, false);
        if (p.pauli_blocking)
          mypauli_spp(*e, *_nucleus);
      }
      break;
    case 6:
      e->flag.coh = e->flag.cc = true;
      if (p.dyn_coh_cc) // coh cc
      {
        if (p.coh_new)
          switch (p.coh_kind) {
          case 1:
            cohevent_cj(p, *e, *_nucleus, true);
            break;
          case 2:
            cohevent_bs(p, *e, *_nucleus, true);
            break;
          default:
            cohevent_bs(p, *e, *_nucleus, true);
            break;
          }
        else
          cohevent2(p, *e, *_nucleus, true);
      }
      break;
    case 7:
      e->flag.coh = e->flag.nc = true;
      if (p.dyn_coh_nc) // coh nc
      {
        if (p.coh_new)
          switch (p.coh_kind) {
          case 1:
            cohevent_cj(p, *e, *_nucleus, false);
            break;
          case 2:
            cohevent_bs(p, *e, *_nucleus, false);
            break;
          default:
            cohevent_bs(p, *e, *_nucleus, false);
            break;
          }
        else
          cohevent2(p, *e, *_nucleus, false);
      }
      break;
    case 8:
      e->flag.mec = e->flag.cc = true;
      if (p.dyn_mec_cc) // mec cc
      // if( nu.pdg>0 || !(p.mec_kind==3) )// al flavor states/antineutrinos
      // available
      {
        if (_nucleus->A() <= 1)
          break;
        switch (p.mec_kind) {
        case 1:
          mecevent_tem(p, *e, *_nucleus, true);
          break;
        case 2:
          mecevent2(p, *e, *_nucleus, true, false);
          break;
        case 3:
          mecevent_Nieves(p, *e, *_nucleus, true);
          break;
        case 4:
          mecevent2(p, *e, *_nucleus, true, true);
          break;
        case 5:
          mecevent_SuSA(p, *e, *_nucleus, true);
          break;
        default:
          mecevent_tem(p, *e, *_nucleus, true);
          break;
        }
        for (int i = 0; i < e->out.size(); i++) {
          e->out[i].r = e->in[1].r;
          e->out[i].set_momentum(e->out[i].p().fromZto(e->in[0].p()));
        }
      }
      break;
    case 9:
      e->flag.mec = e->flag.nc = true;
      if (p.dyn_mec_nc)      // mec nc
        if (p.mec_kind == 1) // only TEM for NC
        {
          if (_nucleus->A() <= 1)
            break;
          switch (p.mec_kind) {
          case 1:
            mecevent_tem(p, *e, *_nucleus, false);
            break;
          default:
            mecevent_tem(p, *e, *_nucleus, false);
            break;
          }
          for (int i = 0; i < e->out.size(); i++) {
            e->out[i].r = e->in[1].r;
            e->out[i].set_momentum(e->out[i].p().fromZto(e->in[0].p()));
          }
        }
      break;
    case 10:
      e->flag.hyp = e->flag.cc = true;
      if (p.dyn_hyp_cc) // qel hyperon
      {
        hypevent(p, *e, *_nucleus);
      }
      break;
      // case 12:
      // 	e->flag.lep=true; //->flag.cc=true;
      // 	if (p.dyn_lep) // Neutrino-lepton
      // 	{
      // 		lepevent (p, *e); //, true);
      // 	}
      // 	break;
    }
  else if (e->in[0].pdg == 11) // electron scattering
  {
    switch (dyn) {
    case 20:
      // TODO: introduce a new flag el!
      e->flag.qel = e->flag.nc = true;
      /*if(p.eel_alg=="old")
e_el_event(p,*e,*_nucleus,false);
      else
if(p.eel_alg=="fast")
e_el_event2orig(p,*e,*_nucleus,false);
else   // all remaining algorithms
e_el_event2(p,*e,*_nucleus,false); */
      if (p.sf_method > 0 and has_sf(*_nucleus, p.sf_method))
        sfevent(p, *e, *_nucleus);
      else
        qelevent1(p, *e, *_nucleus, true);
      break;
      // case 21:
      // 	e->flag.nc=true;
      //     if(p.eel_theta_lab>0)
      //                  e_spp_event(p,*e,*_nucleus,false);
      //     else // use negative theta to test new implementation
      //                  e_spp_event3(p,*e,*_nucleus,false);
      //     break;
    }
  }
  e->weight *= factor;

  if (e->weight == 0) {
    e->out.clear();

    e->out.push_back(e->in[0]);
    e->out.push_back(e->in[1]);
  }
  //      e->check();
} // end of makeevent

void NuWro_metropolis::finishevent(event *e, params &p) {
  for (int i = 0; i < 1 /* e->in.size()*/; i++) {
    e->in[i].endproc = e->dyn;
    registration(e->all, e->in[i]);
  }
  for (int i = 0; i < e->out.size(); i++) {
    e->out[i].mother = 0;
    registration(e->all, e->out[i]);
  }
  if (p.beam_test_only)
    return;

  for (int j = 0; j < e->in.size(); j++) {
    particle p = e->in[j];
    p.endproc = e->dyn;
    registration(e->all, p);
  }

  // e->pr=_nucleus->Zr(); 	// 1. po co to?
  // e->nr=_nucleus->Nr(); 	// 2. powoduje break, segmentation fault

  // copy particle from out to post if coherent interaction

  if (!e->flag.coh && !e->flag.lep) {
    kaskada k(p, *e, &input);
    k.kaskadaevent(); // runs only if p.kaskada_on is true
  } else
  //	if(e->post.size()==0)   // copy out to post if no fsi
  {
    for (int j = 0; j < e->out.size(); j++) {
      particle p = e->out[j];
      p.endproc = e->dyn;
      registration(e->all, p);
      e->post.push_back(p);
    }
  }
} // end of finishevent

void NuWro_metropolis::real_events(params &p) {
  // dismode = true;
  if (p.number_of_events < 1)
    return;

  frame_bottom();

  frame_top("Run real events");

  std::unique_ptr<event> u_e(new event);
  event *e = u_e.get();

  string output = a.output;
  int l = output.length();
  if (l < 5 || string(".root") != output.c_str() + l - 5)
    output = output + ".root";
  TFile *ff = new TFile(output.c_str(), "recreate");
  TTree *tf = new TTree("treeout", "Tree of events");
  tf->Branch("e", "event", &e);
  std::unordered_map<int, double> channel_weight_sum{},
      channel_weight_sum_fraction{};
  std::unordered_map<int, double> channel_count_final{};
  // std::random_device rd;
  // std::mt19937 gen = std::mt19937(rd());
  std::uniform_int_distribution<> dis(0, enabled_dyns.size() - 1);
  int accepted_count{};
  for (int i{}; i < p.number_of_events; i++) {
    bool accepted = false;
    for (int j{}; j < 100; j++) {
      auto e = new event();
      e->dyn = enabled_dyns[dis(gen)];
      makeevent(e, p);
      auto thisbias = 1 / e->in[0].t;
      auto biased_weight = e->weight * thisbias;
      if (isnan(thisbias) || isnan(biased_weight)) { // ignore NaNs
        i--;
        continue;
      }
      channel_weight_sum[e->dyn] += biased_weight;
      channel_weight_sum_fraction[e->dyn] += thisbias;
      accepted |= sampler.update_state(e);
    }
    accepted_count += accepted;
    *e = sampler.get_state();
    channel_count_final[e->dyn]++;
    finishevent(e, p);
    tf->Fill();
  }
  double overall_xsec{};
  for (auto &i : enabled_dyns) {
    if (channel_weight_sum[i])
      std::cout << "Channel " << i << " weight avg: "
                << channel_weight_sum[i] / channel_weight_sum_fraction[i]
                << '\n'
                << " count: " << channel_count_final[i] << '\n'
                << '\n';
    overall_xsec += channel_weight_sum[i] / channel_weight_sum_fraction[i];
  }
  std::cout << "Overall acceptance: "
            << (double)accepted_count / p.number_of_events << '\n';
  // TVectorD xsecs(1);
  // xsecs[0] = overall_xsec;
  tf->GetUserInfo()->Add(new TParameter<double>("xsec", overall_xsec));
  ff->Write();
  ff->Close();
  delete ff;
  frame_top("Finalize the simulation");
  cout << "        "
       << "-> Generated the output file: \"" << output << "\"" << endl;
  frame_bottom();
}

event NuWro_metropolis::get_event() {
  std::uniform_int_distribution<> dis(0, enabled_dyns.size() - 1);
  for (int j{}; j < 100; j++) {
    auto e = new event();
    e->dyn = enabled_dyns[dis(gen)];
    makeevent(e, p);
    sampler.update_state(e);
  }
  event e = sampler.get_state();
  finishevent(&e, p);
  return e;
}

void NuWro_metropolis::kaskada_redo(string input, string output) {
  event *e = new event;

  TFile *fi = new TFile(input.c_str());
  TTree *ti = (TTree *)fi->Get("treeout");
  if (ti == NULL) {
    cerr << "tree \"treeout\" not found in file \"" << input << "\"" << endl;
    exit(7);
  }
  ti->SetBranchAddress("e", &e);

  TFile *ff = new TFile(output.c_str(), "recreate");
  TTree *tf = new TTree("treeout", "Tree of events");
  tf->Branch("e", "event", &e);

  int nn = ti->GetEntries();
  for (int i = 0; i < nn; i++) {
    //		e = new event();
    ti->GetEntry(i);
    e->clear_fsi();
    finishevent(e, p);
    tf->Fill();
    //		delete e;
    // if(i%1000==0)
    // cout<<i/1000<<"/"<<nn/1000<<"\r"<<endl;
    // raport(i + 1, nn, " % events processed...", 100, e->dyn,
    // bool(a.progress));
  }
  cout << endl;
  fi->Close();
  delete fi;

  ff->Write();
  ff->Close();
  delete e;
  delete ff;
  cout << "Output: \"" << output << "\"" << endl;
}

// #include "UserAction.h"

void NuWro_metropolis::main(int argc, char **argv) {
  shhpythiaitokay_();
  try {
    init(argc, argv);
    if (p.kaskada_redo == 1)
      kaskada_redo(a.output, string(a.output) + ".fsi.root");
    else {
      real_events(p);
    }
    genrand_write_state();
  } catch (string s) {
    cout << s << endl;
  } catch (char const *s) {
    cout << s << endl;
  } catch (...) {
    cout << "NuWro failed" << endl;
  }
}
