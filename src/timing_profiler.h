#ifndef _timing_profiler_h_
#define _timing_profiler_h_

#include <chrono>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

using namespace std;

class TimingProfiler {
private:
    using Clock = std::chrono::high_resolution_clock;
    using TimePoint = std::chrono::time_point<Clock>;
    using Duration = std::chrono::duration<double>;

    // Timestamps
    TimePoint process_start_, process_end_;
    TimePoint test_run_start_, test_run_end_;
    TimePoint real_run_start_, real_run_end_;
    TimePoint copy_events_start_, copy_events_end_;

    // Per-channel timing for real run
    std::vector<TimePoint> channel_start_times_;
    std::vector<double> channel_elapsed_;  // in seconds

    // MH timing
    TimePoint mh_start_, mh_end_;

    // Preparation timing (from process start to test run or MH start)
    TimePoint preparation_end_;
    bool preparation_recorded_ = false;

    int num_channels_ = 0;
    bool test_run_active_ = false;
    bool real_run_active_ = false;
    bool copy_events_active_ = false;
    bool mh_active_ = false;

public:
    void init_process() {
        process_start_ = Clock::now();
    }

    void end_process() {
        process_end_ = Clock::now();
    }

    void init_test_run() {
        test_run_start_ = Clock::now();
        test_run_active_ = true;
    }

    void end_test_run() {
        test_run_end_ = Clock::now();
        test_run_active_ = false;
    }

    void init_real_run(int num_channels) {
        num_channels_ = num_channels;
        channel_start_times_.resize(num_channels);
        channel_elapsed_.assign(num_channels, 0.0);
        real_run_start_ = Clock::now();
        real_run_active_ = true;
    }

    void start_channel(int channel_idx) {
        if (channel_idx >= 0 && channel_idx < num_channels_) {
            channel_start_times_[channel_idx] = Clock::now();
        }
    }

    void end_channel(int channel_idx) {
        if (channel_idx >= 0 && channel_idx < num_channels_) {
            auto end = Clock::now();
            channel_elapsed_[channel_idx] += std::chrono::duration<double>(end - channel_start_times_[channel_idx]).count();
        }
    }

    void end_real_run() {
        real_run_end_ = Clock::now();
        real_run_active_ = false;
    }

    void init_copy_events() {
        copy_events_start_ = Clock::now();
        copy_events_active_ = true;
    }

    void end_copy_events() {
        copy_events_end_ = Clock::now();
        copy_events_active_ = false;
    }

    void init_mh() {
        mh_start_ = Clock::now();
        mh_active_ = true;
    }

    void end_mh() {
        mh_end_ = Clock::now();
        mh_active_ = false;
    }

    void end_preparation() {
        if (!preparation_recorded_) {
            preparation_end_ = Clock::now();
            preparation_recorded_ = true;
        }
    }

    double get_process_time() const {
        if (process_start_.time_since_epoch().count() == 0) return 0.0;
        auto end = (process_end_.time_since_epoch().count() == 0) ? Clock::now() : process_end_;
        return std::chrono::duration<double>(end - process_start_).count();
    }

    double get_test_run_time() const {
        if (!test_run_active_ && test_run_start_.time_since_epoch().count() == 0) return 0.0;
        auto end = test_run_active_ ? Clock::now() : test_run_end_;
        return std::chrono::duration<double>(end - test_run_start_).count();
    }

    double get_real_run_time() const {
        if (!real_run_active_ && real_run_start_.time_since_epoch().count() == 0) return 0.0;
        auto end = real_run_active_ ? Clock::now() : real_run_end_;
        return std::chrono::duration<double>(end - real_run_start_).count();
    }

    double get_copy_events_time() const {
        if (!copy_events_active_ && copy_events_start_.time_since_epoch().count() == 0) return 0.0;
        auto end = copy_events_active_ ? Clock::now() : copy_events_end_;
        return std::chrono::duration<double>(end - copy_events_start_).count();
    }

    double get_mh_time() const {
        if (!mh_active_ && mh_start_.time_since_epoch().count() == 0) return 0.0;
        auto end = mh_active_ ? Clock::now() : mh_end_;
        return std::chrono::duration<double>(end - mh_start_).count();
    }

    double get_preparation_time() const {
        if (!preparation_recorded_) return 0.0;
        return std::chrono::duration<double>(preparation_end_ - process_start_).count();
    }

    double get_channel_time(int channel_idx) const {
        if (channel_idx < 0 || channel_idx >= num_channels_) return 0.0;
        return channel_elapsed_[channel_idx];
    }

    void print_report(ostream& out, bool use_mh, const class chooser* proc = nullptr) {
        double total_time = get_process_time();

        // Helper lambda to format table lines
        auto line = [&](char c, int width) {
            out << string(width, c);
        };

        out << endl;
        line('_', 78);
        out << endl;
        out << " |" << endl;

        if (use_mh) {
            out << " |  Timing Profile (Metropolis-Hastings)" << endl;
        } else {
            out << " |  Timing Profile" << endl;
        }

        line('_', 78);
        out << endl;
        out << " |" << endl;

        // Header row
        out << " |  " << left << setw(32) << "Section"
            << " | " << right << setw(10) << "Time (s)"
            << " | " << right << setw(12) << "Percentage"
            << " |" << endl;

        line('-', 78);
        out << endl;
        out << " |" << endl;

        // Preparation time (shown for both modes)
        double prep_time = get_preparation_time();
        if (prep_time > 0) {
            double pct = total_time > 0 ? (prep_time / total_time * 100.0) : 0.0;
            out << " |  " << left << setw(32) << "Preparation"
                << " | " << right << setw(10) << fixed << setprecision(2) << prep_time
                << " | " << right << setw(11) << setprecision(1) << pct << "%"
                << " |" << endl;
        }

        if (use_mh) {
            double mh_time = get_mh_time();
            double pct = total_time > 0 ? (mh_time / total_time * 100.0) : 0.0;
            out << " |  " << left << setw(32) << "MH Generation"
                << " | " << right << setw(10) << fixed << setprecision(2) << mh_time
                << " | " << right << setw(11) << setprecision(1) << pct << "%"
                << " |" << endl;
        } else {
            // Test run time
            double test_time = get_test_run_time();
            if (test_time > 0) {
                double pct = total_time > 0 ? (test_time / total_time * 100.0) : 0.0;
                out << " |  " << left << setw(32) << "Test Run"
                    << " | " << right << setw(10) << fixed << setprecision(2) << test_time
                    << " | " << right << setw(11) << setprecision(1) << pct << "%"
                    << " |" << endl;
            }

            // Per-channel real run times
            double real_run_total = 0.0;
            for (int i = 0; i < num_channels_; i++) {
                double ch_time = channel_elapsed_[i];
                real_run_total += ch_time;
                if (ch_time > 0 && proc != nullptr) {
                    double pct = total_time > 0 ? (ch_time / total_time * 100.0) : 0.0;
                    string label = "  Real Run - " + proc->label(i);
                    out << " |  " << left << setw(32) << label
                        << " | " << right << setw(10) << fixed << setprecision(2) << ch_time
                        << " | " << right << setw(11) << setprecision(1) << pct << "%"
                        << " |" << endl;
                }
            }

            // Real run total
            if (real_run_total > 0) {
                double pct = total_time > 0 ? (real_run_total / total_time * 100.0) : 0.0;
                out << " |  " << left << setw(32) << "Real Run Total"
                    << " | " << right << setw(10) << fixed << setprecision(2) << real_run_total
                    << " | " << right << setw(11) << setprecision(1) << pct << "%"
                    << " |" << endl;
            }

            // Copy events time
            double copy_time = get_copy_events_time();
            if (copy_time > 0) {
                double pct = total_time > 0 ? (copy_time / total_time * 100.0) : 0.0;
                out << " |  " << left << setw(32) << "Event Copying"
                    << " | " << right << setw(10) << fixed << setprecision(2) << copy_time
                    << " | " << right << setw(11) << setprecision(1) << pct << "%"
                    << " |" << endl;
            }
        }

        // Separator before total
        line('-', 78);
        out << endl;
        out << " |" << endl;

        // Total process time
        out << " |  " << left << setw(32) << "Total Process Time"
            << " | " << right << setw(10) << fixed << setprecision(2) << total_time
            << " | " << right << setw(11) << setprecision(1) << 100.0 << "%"
            << " |" << endl;

        out << " |" << endl;
        line('_', 78);
        out << endl << endl;

        // Reset formatting
        out << defaultfloat << setprecision(6);
    }
};

#endif
