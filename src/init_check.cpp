#include "init_check.hpp"

namespace lvx {

  // ------------------------------------------------------------
  void init_check1() {
    using namespace boost::filesystem;

    path p("lavax.conf");

    if (!(exists(p) && is_regular_file(p))) {
      std::cout << "Couldn't find file 'lavax.conf'. Creating default" << "\n";
      copy_file(path(DATADIR "/lavax.conf"), path("./lavax.conf"), copy_option::overwrite_if_exists);
    }

    p = path("predictor.in");

    if (!(exists(p) && is_regular_file(p))) {
      std::cout << "Couldn't find file 'predictor.in'. Creating default" << "\n";
      copy_file(path(DATADIR "/predictor.in"), path("./predictor.in"), copy_option::overwrite_if_exists); // current_path()
    }

    std::cout << "Remember to add a LAMMPS potential file as specified in 'predictor.in' "
              << "and a properly concatenated POTCAR file for VASP" << "\n";
  }

  // ------------------------------------------------------------
  bool init_check2(const std::string& init_poscar) {
    using namespace boost::filesystem;
    std::vector<std::string> files = {"INCAR", "KPOINTS", "POTCAR"};
    files.push_back(init_poscar);

    bool success = true;
    for (std::string s : files) {
      path p(s);
      if (!exists(p)) {
        std::cout << "Please provide file '" << s << "'" << "\n";
        success = false;
      }
    }
    return success;
  }
  

}

