#include <iostream>
// #include "small_scale_units.hpp"
#include "lavax.hpp"
#include <boost/units/cmath.hpp>

int main(int argc, char *argv[])
{
  using namespace lvx;
  
  if (argc != 6) {
    std::cout << "usage: PKAvel <mass (u)> <x-dir> <y-dir> <z-dir> <energy (eV)>" << "\n";
    return -1;
  }
  
  quantity<atomic_mass_unit> m = std::stod(argv[1]) * u;

  vec3_dimless dir;
  dir = std::stod(argv[2]),
        std::stod(argv[3]),
        std::stod(argv[4]);

  quantity<electronvolt_unit> Ek = std::stod(argv[5]) * eV;

  quantity<velocity_unit> v = sqrt(2.0 * static_cast<quantity<energy_unit>>(Ek)/m);

  v = v / norm(dir);

  vec3_velocity vel(v*dir[0], v*dir[1], v*dir[2]);

  std::cout << vel << "\n";
  return 0;
}

