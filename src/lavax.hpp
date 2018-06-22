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

  class atom_species {
  public:
    std::string symbol;
    std::string vasp_good_pot;
    std::string vasp_bad_pot;
    quantity<atomic_mass_unit> mass;
    std::vector<Particle_v2> good_particles;
    std::vector<Particle_v2> bad_particles;
  // public:

  //   void set_good_potential(std::string buffer) {
  //     vasp_good_pot = buffer;
  //   }

  //   void set_bad_potential(std::string buffer) {
  //     vasp_bad_pot = buffer;
  //   }

  //   void set_mass(quantity<atomic_mass_unit> mass) {
  //     m_mass = mass;
  //   }

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

    // int index_by_vasp_good(std::string) {
      
    // }
  };

  class atomic_element {
  public:
    std::string symbol;
    std::string vasp_symbol_good;
    std::string vasp_symbol_bad;
    std::string vasp_potential_file_good;
    std::string vasp_potential_file_bad;
    // std::pair<int,int> vasp_indices;
    int vasp_num_good;
    int vasp_num_bad;
    quantity<atomic_mass_unit> mass;
    
    int number_atoms() {
      return vasp_num_good + vasp_num_bad;
    }
    // int number_atoms() {
    //   return std::accumulate(vasp_indices.begin(),
    //                          vasp_indices.end(), 0);
    // }
  };

  class simulation_cell_v2 {
  public:
    std::vector<atomic_element> elements;
    std::vector<Particle_v2>    particles;

    // int index_by_vasp_good(std::string) {
      
    // }
  };

  class simulation_cell {
  public:
    std::vector<std::string> vasp_potential_symbols;
    std::vector<std::string> vasp_potential_file;
    std::vector<int>         vasp_indices;
    std::vector<int>         lammps_indices;
    std::vector<Particle_v2> particles;
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

  void parse_poscar(std::ifstream& poscar,
                    quantity<angstrom_unit>& lattice_constant,
                    vec3_dimless& a1,
                    vec3_dimless& a2,
                    vec3_dimless& a3,
                    std::vector<Particle_v2>& P);

  void parse_init_poscar(std::ifstream& poscar,
                    quantity<angstrom_unit>& lattice_constant,
                    vec3_dimless& a1,
                    vec3_dimless& a2,
                    vec3_dimless& a3,
                    std::vector<atom_species>& P);

  void parse_init_poscar(std::ifstream& poscar,
                    quantity<angstrom_unit>& lattice_constant,
                    vec3_dimless& a1,
                    vec3_dimless& a2,
                    vec3_dimless& a3,
                    simulation_cell_v2& sim_cell);

  void parse_init_poscar(std::ifstream& poscar,
                    quantity<angstrom_unit>& lattice_constant,
                    vec3_dimless& a1,
                    vec3_dimless& a2,
                    vec3_dimless& a3,
                    simulation_cell_v3& sim_cell);

  std::set<int> parse_lammps_neighbor(std::ifstream& file, quantity<angstrom_unit> cutoff);
                    
  
  std::vector<std::vector<int>> parse_lammps_neighbor(std::ifstream& file);


  std::string make_lammps_data(quantity<angstrom_unit> xhi,
                               quantity<angstrom_unit> yhi,
                               quantity<angstrom_unit> zhi,
                               std::vector<atom_species>& P, bool incl_vel = true);

  std::string make_lammps_data(quantity<angstrom_unit> xhi,
                               quantity<angstrom_unit> yhi,
                               quantity<angstrom_unit> zhi,
                               simulation_cell_v2& sim_cell,
                               bool include_velocity = true);

  std::string make_lammps_data(quantity<angstrom_unit> xhi,
                               quantity<angstrom_unit> yhi,
                               quantity<angstrom_unit> zhi,
                               simulation_cell_v3& sim_cell,
                               bool include_velocity = true);

  std::string make_poscar(const quantity<angstrom_unit>& lattice_constant,
                          const vec3_dimless& a1,
                          const vec3_dimless& a2,
                          const vec3_dimless& a3,
                          const simulation_cell& sim_cell);
  
  // vec3_velocity PKAvel(const quantity<atomic_mass_unit>& mass,
  //                      const vec3_dimless& dir,
  //                      const quantity<electron_volt_unit>& energy) {

  // }

}

#endif
