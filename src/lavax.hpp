#ifndef LAVAX_INC // Inclusion gaurd
#define LAVAX_INC

#include <cmath>
#include <sstream>
#include <fstream>
#include <set>
#include <iomanip>
#include <algorithm>
#include <regex>
#include <numeric>
#include <tuple>
// #include <armadillo>
// #include <Eigen/Core>
#include "small_scale_units.hpp"

#include <blitz/array.h>
// #include <units.h>

namespace lvx {
  // using namespace arma;
  // using namespace units;
  using namespace boost::units;
  // using namespace boost::units::si;
  // boost::units::si::constants::codata::m_u

  using vec3_length   = quantity<si::length,   blitz::Array<double,1>>;
  // using vec3_velocity = quantity<si::velocity, blitz::Array<double,1>>;

  using vec3_angstrom = blitz::TinyVector<quantity<angstrom_unit>,     3>;
  using vec3_velocity = blitz::TinyVector<quantity<velocity_unit>,     3>;
  using vec3_dimless  = blitz::TinyVector<quantity<si::dimensionless>, 3>;

  template <typename T>
  std::ostream& operator<<(std::ostream &strm,
                           const blitz::TinyVector<quantity<T>,3>& a) {
    return strm << std::fixed << std::setprecision(8)
                << a(0).value() << " " << a(1).value() << " " << a(2).value();
  }

  // ------------------------------------------------------------
  class Particle_v2 {
  private:
    unsigned int type = 0;
  public:
    vec3_angstrom pos;
    vec3_velocity vel;
    
    Particle_v2(vec3_angstrom pos, vec3_velocity vel) : pos(pos), vel(vel) {
      // this->vel = vel;
    }
    Particle_v2(unsigned int type,
                vec3_angstrom pos,
                vec3_velocity vel) : type(type), pos(pos), vel(vel) {}

    vec3_angstrom getPos() const { return pos; }
    vec3_velocity getVel() const { return vel; }
  };

  class atomic_element_info {
  public:
    std::string symbol;
    int         atom_type;
    std::string vasp_symbol_good;
    std::string vasp_symbol_bad;
    std::string vasp_potential_file_good;
    std::string vasp_potential_file_bad;
    quantity<atomic_mass_unit> mass;
  };

  class atomic_particle: public Particle_v2 {
  public:
    atomic_particle(vec3_angstrom pos, vec3_velocity vel) : Particle_v2(pos,vel) {}
    std::shared_ptr<atomic_element_info> element_info;
    bool high_prec = false;
  };

  class simulation_cell_v3 {
  public:
    std::vector<std::shared_ptr<lvx::atomic_element_info> >
                                 elements_info;
    std::vector<atomic_particle> particles;
    std::vector<std::tuple<std::string,int,bool> > vasp_symbol_count_helper;
    // int index_by_vasp_good(std::string) {
      
    // }
  };

  template<typename T, typename U>
  blitz::TinyVector<quantity<typename U::unit_type>,3>
  read_vec3(T&& strm, U unit) {
    // std::stringstream strm(str);
    double x, y, z;
    strm >> x;
    strm >> y;
    strm >> z;

    blitz::TinyVector<quantity<typename U::unit_type>,3> rtn;
    rtn = x * unit,
          y * unit,
          z * unit;
    
    return rtn;
  }

  template<typename U>
  blitz::TinyVector<quantity<typename U::unit_type>,3>
  read_vec3(std::string& str, U unit) {
    return read_vec3(std::stringstream(str), unit);
  }

  // Ugly, but works for our purposes:
  template<typename U>
  quantity<U>
  norm(blitz::TinyVector<quantity<U>,3> v) {
    return sqrt((v[0]*v[0] + v[1]*v[1] + v[2]*v[2]).value()) * (typename U::unit());
  }

  enum class crystal_structure { BCC, FCC };
  
  std::vector<Particle_v2> create_crystal(crystal_structure cstruct,
                                          quantity<angstrom_unit> lattice_constant,
                                          int I, int J, int K);

  void parse_init_poscar(std::ifstream& poscar,
                         quantity<angstrom_unit>& lattice_constant,
                         vec3_dimless& a1,
                         vec3_dimless& a2,
                         vec3_dimless& a3,
                         simulation_cell_v3& sim_cell);

  std::set<int> parse_lammps_neighbor(std::ifstream& file, quantity<angstrom_unit> cutoff);
  
  std::set<int> parse_lammps_neighbor(std::ifstream& file,
                                      quantity<angstrom_unit> cutoff,
                                      int max_potential_switch,
                                      int max_vasp_nsw,
                                      int count_good_prev,
                                      int& NSW);
  
  std::string make_lammps_data(quantity<angstrom_unit> xhi,
                               quantity<angstrom_unit> yhi,
                               quantity<angstrom_unit> zhi,
                               simulation_cell_v3& sim_cell,
                               bool include_velocity = true);

  std::string make_poscar(const quantity<angstrom_unit>& lattice_constant,
                          const vec3_dimless& a1,
                          const vec3_dimless& a2,
                          const vec3_dimless& a3,
                          const simulation_cell_v3& sim_cell);
  
  // vec3_velocity PKAvel(const quantity<atomic_mass_unit>& mass,
  //                      const vec3_dimless& dir,
  //                      const quantity<electron_volt_unit>& energy);

}

#endif
