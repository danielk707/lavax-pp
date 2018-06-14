#include <blitz/array.h>
#include <cmath>
#include <boost/units/conversion.hpp>
#include <boost/units/static_constant.hpp>
#include <boost/units/systems/si/length.hpp>
#include <boost/units/systems/si/velocity.hpp>
#include <boost/units/systems/si/io.hpp>
#include <boost/units/systems/si/mass.hpp>
#include <boost/units/base_units/metric/angstrom.hpp>
#include <boost/units/systems/si/codata/physico-chemical_constants.hpp>
#include <boost/units/base_units/us/pound.hpp>

#include <boost/units/base_units/metric/angstrom.hpp>
typedef boost::units::metric::angstrom_base_unit::unit_type angstrom_unit;
BOOST_UNITS_STATIC_CONSTANT(angstrom, angstrom_unit);
BOOST_UNITS_STATIC_CONSTANT(angstroms, angstrom_unit);

BZ_USING_NAMESPACE(blitz)

using namespace boost::units;
using namespace boost::units::si;

// template<typename T>
// class TD;

namespace lvx {

template<typename U>
quantity<U>
abs(blitz::TinyVector<quantity<U>,3> v) {
  return sqrt((v[0]*v[0] + v[1]*v[1] + v[2]*v[2]).value()) * (typename U::unit());
}

}

int main(int argc, char *argv[])
{
  blitz::TinyVector<double,3> a;
  a = 1.0, 2.0, 3.0;

  blitz::firstIndex i;
  std::cout << sqrt(blitz::sum(a(i)*a(i))) << "\n";

  std::cout << a << "\n";
  std::cout << fabs(a)[0] << "\n";

  // ------------------------------

  blitz::TinyVector<quantity<angstrom_unit>,3> b;
  b = 1.0 * angstrom,
      2.0 * angstrom,
      4.0 * angstrom;

  std::cout << b << "\n";
  std::cout << lvx::abs(b) << "\n";

  // blitz::sum(b(i)*b(i));
  
  
  return 0;
}
