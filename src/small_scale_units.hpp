#ifndef SMALL_SCALE_UNITS
#define SMALL_SCALE_UNITS

#include <boost/units/conversion.hpp>
#include <boost/units/static_constant.hpp>
#include <boost/units/systems/si/length.hpp>
#include <boost/units/systems/si/velocity.hpp>
#include <boost/units/systems/si/io.hpp>
#include <boost/units/systems/si/mass.hpp>
#include <boost/units/base_units/metric/angstrom.hpp>
#include <boost/units/systems/si/codata/physico-chemical_constants.hpp>
#include <boost/units/base_units/us/pound.hpp>

#include <blitz/array.h>

// typedef boost::units::us::pound_base_unit::unit_type pound_unit;
// BOOST_UNITS_DEFINE_CONVERSION_FACTOR(pound_unit, boost::units::si::mass, float, 2.2f);

struct atomic_mass_base_unit :
  boost::units::base_unit<atomic_mass_base_unit, boost::units::mass_dimension, 2> {};

template<>
struct boost::units::base_unit_info<atomic_mass_base_unit>
{
  static std::string name()               { return "atomic mass"; }
  static std::string symbol()             { return "u"; }
};

using atomic_system = boost::units::make_system<atomic_mass_base_unit>::type;

using atomic_mass_unit = boost::units::unit<boost::units::mass_dimension, atomic_system>;

// using atomic_mass_dim  = boost::units::derived_dimension<boost::units::mass_base_dimension,1>::type;

// using my_system = make_system<metric::angstrom_base_unit,
//                               femtosecond_base_unit>::type;

// using atomic_mass_unit = boost::units::unit<atomic_mass_dim, boost::units::si::system>;
BOOST_UNITS_STATIC_CONSTANT(u, atomic_mass_unit);

BOOST_UNITS_DEFINE_CONVERSION_FACTOR(atomic_mass_base_unit, si::mass, double,
                                     boost::units::si::constants::codata::m_u.value().value());

#include <boost/units/base_units/metric/angstrom.hpp>
typedef boost::units::metric::angstrom_base_unit::unit_type angstrom_unit;
BOOST_UNITS_STATIC_CONSTANT(angstrom, angstrom_unit);
BOOST_UNITS_STATIC_CONSTANT(angstroms, angstrom_unit);

typedef boost::units::si::dimensionless::unit_type dimless_unit;
BOOST_UNITS_STATIC_CONSTANT(dimensionless, dimless_unit);

// using boost::units::derived_dimension<boost::units::length_base_dimension

namespace lvx {
  using namespace boost::units;

  using femtosecond_base_unit = scaled_base_unit<si::second_base_unit,
                                                 scale<10, static_rational<-15>>>;
  
  using my_system = make_system<metric::angstrom_base_unit,
                                femtosecond_base_unit>::type;

  using velocity = unit<velocity_dimension, my_system>;

  using velocity_unit = velocity::unit_type;
  BOOST_UNITS_STATIC_CONSTANT(angstrom_per_fs, velocity_unit);

  // BOOST_UNITS_DEFINE_CONVERSION_FACTOR(atomic_mass_unit::unit_type, si::mass, double, si::constants::codata::m_u);
  
  // using vec3_angstrom = Eigen::Matrix<units::length::angstrom_t, 3, 1>;
  // using vec3_length   = Eigen::Matrix<quantity<length>, 3, 1>;

  using vec3_length   = quantity<si::length,   blitz::Array<double,1>>;
  // using vec3_velocity = quantity<si::velocity, blitz::Array<double,1>>;
  

  using vec3_angstrom = blitz::TinyVector<quantity<angstrom_unit>,     3>;
  using vec3_velocity = blitz::TinyVector<quantity<velocity>,          3>;
  using vec3_dimless  = blitz::TinyVector<quantity<si::dimensionless>, 3>;

  // class Angstrom {
  //   double raw_length;
  // public:
  //   class DoubleIsAngstrom{};
  //   explicit constexpr Angstrom(DoubleIsAngstrom, double lng) : raw_length(lng) {}
  //   operator double () const { return raw_length; }
  // };

  // constexpr Angstrom operator "" _A(long double length) {
  //   return Angstrom{Angstrom::DoubleIsAngstrom{}, static_cast<double>(length)};
  // }

}
#endif
