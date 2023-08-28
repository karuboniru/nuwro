#include <nuwro_metropolis.h>

int main(int argc, char **argv) {

  cout << R"(
  ____________________________________________________________________________
 |                                                                            |
 |                                                                            |
 |                                      `.``      `.       .-.        `       |
 |                                  `-/+os+s  ./ohmN:.    sNNNy:`   .---.+`   |
 |    |\ |     |  |  _  _           :oooyysy: +oodMMd-`  .MMMMM/-   -----s/   |
 |    | \| |_| |/\| |  (_)           `.`oyy+d`   `mMMo-   yMMMyo`   .----d.   |
 |             __        __   __        .yyyoo    :MMN:.  :MMho.    `---h:    |
 |              _) /|   /  \ (__\        :yyoh-    sMMh-`.mMho-     ---h:     |
 |             /__  | . \__/  __/         oyy+h    `mMM+-mMho-     ---h:      |
 |                                        .yyys+    -MMNmMh+-     ---h:       |
 |                                         :yyod.    sMMMh+-     ---h:        |
 |   Wrocław Neutrino Event Generator       oyy+h    `mMh-+.    ---h:         |
 |   https://github.com/NuWro/nuwro         .yyss/  .s/y--oh   .--y-          |
 |                                           :yy+d`.sy+----d/ .--y-           |
 |   J. T. Sobczyk et al.                     oyy++syoy`.--:s.--y-            |
 |   Institute of Theoretical Physics         .yyssyoy.  ------y-             |
 |   University of Wrocław                     :yyyoh.   `----y-              |
 |   Poland                                     osoh.     .-:y-               |
 |                                              `-:.       .:-                |
 |                                                                            |
 |____________________________________________________________________________|
             )"
       << endl;
  NuWro_metropolis nuwro_metropolis;
  nuwro_metropolis.main(argc, argv);
  return 0;
}