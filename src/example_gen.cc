#include "libnuwro.h"
#include "jednostki.h"  // GeV, cm2 unit constants
#include <iostream>
#include <iomanip>

// Channel index → human-readable label
static const char *dyn_label(int dyn) {
    switch (dyn) {
    case 0:  return "QEL-CC";
    case 1:  return "QEL-NC";
    case 2:  return "RES-CC";
    case 3:  return "RES-NC";
    case 4:  return "DIS-CC";
    case 5:  return "DIS-NC";
    case 6:  return "COH-CC";
    case 7:  return "COH-NC";
    case 8:  return "MEC-CC";
    case 9:  return "MEC-NC";
    case 10: return "HYP-CC";
    case 12: return "LEP";
    case 20: return "QEL-el";
    case 21: return "RES-el";
    default: return "unknown";
    }
}

int main() {
    libnuwro::Generator gen;
    gen.from_file("params.txt")
       .set("use_mh", "1")
       .set("number_of_events", "10")
       .initialize();

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "\n"
              << std::setw(4)  << "#"
              << std::setw(10) << "channel"
              << std::setw(12) << "Enu [GeV]"
              << std::setw(12) << "Q2 [GeV2]"
              << std::setw(8)  << "#post"
              << std::setw(16) << "weight [cm2]"
              << "\n"
              << std::string(62, '-') << "\n";

    for (int i = 0; i < 10; i++) {
        event e = gen.gen_event();

        double Enu  = e.in[0].E() / GeV;
        double Q2   = -e.q2()    / (GeV * GeV);
        int    npost = static_cast<int>(e.post.size());
        double wgt  = e.weight   / cm2;

        std::cout << std::setw(4)  << i
                  << std::setw(10) << dyn_label(e.dyn)
                  << std::setw(12) << Enu
                  << std::setw(12) << Q2
                  << std::setw(8)  << npost
                  << std::setw(16) << wgt
                  << "\n";
    }
    std::cout << "\n";
    return 0;
}
