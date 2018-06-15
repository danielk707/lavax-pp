#include "parse_config_file.hpp"
#include "lavax.hpp"
#include "init_check.hpp"
#include <iostream>
#include <fstream>

using namespace lvx;

quantity<atomic_mass_unit> operator "" _u(long double mass) {
  return static_cast<double>(mass) * u;
}

quantity<angstrom_unit> operator "" _Angstrom(long double length) {
  return static_cast<double>(length) * angstrom;
}

int main(int argc, char *argv[])
{

  std::cout << "-------------------- parse_config_file --------------------" << "\n";
  std::fstream f("lavax.conf");
  
  auto dict = parse_config_file(f);

  for (auto v : dict) {
    std::cout << v.first << "\t\t\t" << v.second << "\n";
  }

  std::cout << "-------------------- create_crystal --------------------" << "\n";

  auto P = create_crystal(crystal_structure::BCC, 3.18_Angstrom, 2, 2, 2);

  for (auto p : P) {
    std::cout << p.getPos() << "\n";
  }

  std::cout << "-------------------- init_check1 --------------------" << "\n";
  init_check1();

  std::cout << "-------------------- init_check2 --------------------" << "\n";
  init_check2("POSCAR_init");
  
  return 0;
}
