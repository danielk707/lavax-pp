#ifndef LAVAX_INC // Inclusion guard
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
#include <array>
#include "small_scale_units.hpp"

template <typename T>
std::ostream& operator<<(std::ostream &strm,
                         const std::array<quantity<T>,3>& a) {
  return strm << std::fixed << std::setprecision(8)
              << a[0].value() << " " << a[1].value() << " " << a[2].value();
}

namespace lvx {
  using namespace boost::units;

  using vec3_angstrom = std::array<quantity<angstrom_unit>,3>;
  using vec3_velocity = std::array<quantity<velocity_unit>,3>;
  using vec3_dimless  = std::array<quantity<si::dimensionless>,3>; 

  // ------------------------------------------------------------
  class Particle {
  private:
    unsigned int type = 0;
  public:
    vec3_angstrom pos;
    vec3_velocity vel;
    
    Particle(vec3_angstrom pos, vec3_velocity vel) : pos(pos), vel(vel) {}
    
    Particle(unsigned int type,
             vec3_angstrom pos,
             vec3_velocity vel) : type(type), pos(pos), vel(vel) {}

    vec3_angstrom getPos() const { return pos; }
    vec3_velocity getVel() const { return vel; }
  };

  class atomic_element_info {
  public:
    std::string symbol;
    int         atom_type;
    std::string vasp_symbol_hard;
    std::string vasp_symbol_soft;
    std::string vasp_potential_file_hard;
    std::string vasp_potential_file_soft;
    quantity<atomic_mass_unit> mass;
  };

  class atomic_particle: public Particle {
  public:
    atomic_particle(vec3_angstrom pos, vec3_velocity vel) : Particle(pos,vel) {}
    std::shared_ptr<atomic_element_info> element_info;
    bool is_hard = false;
    int idx = -1;
  };

  class simulation_cell {
  public:
    std::vector<std::shared_ptr<lvx::atomic_element_info> > elements_info;
    std::vector<atomic_particle> particles;
    std::vector<std::tuple<std::string,int,bool> > vasp_symbol_count_helper;
    std::vector<int> poscar_part_order;
    std::set<int>    lammps_selected_idx;
  };

  template<typename T, typename U>
  std::array<quantity<typename U::unit_type>,3>
  read_vec3(T&& strm, U unit) {
    // std::stringstream strm(str);
    double x, y, z;
    strm >> x;
    strm >> y;
    strm >> z;

    std::array<quantity<typename U::unit_type>,3> rtn;
    rtn[0] = x * unit,
    rtn[1] = y * unit,
    rtn[2] = z * unit;
    
    return rtn;
  }

  template<typename U>
  std::array<quantity<typename U::unit_type>,3>
  read_vec3(std::string& str, U unit) {
    return read_vec3(std::stringstream(str), unit);
  }
  
  // Ugly, but works for our purposes:
  template<typename U>
  quantity<U>
  norm(std::array<quantity<U>,3> v) {
    return sqrt((v[0]*v[0] + v[1]*v[1] + v[2]*v[2]).value()) * (typename U::unit());
  }

  enum class crystal_structure { BCC, FCC };
  
  std::vector<Particle> create_crystal(crystal_structure cstruct,
                                       quantity<angstrom_unit> lattice_constant,
                                       int I, int J, int K);

  void parse_poscar(std::ifstream& poscar,
                    quantity<angstrom_unit>& lattice_constant,
                    vec3_dimless& a1,
                    vec3_dimless& a2,
                    vec3_dimless& a3,
                    simulation_cell& sim_cell);

  std::set<int> parse_lammps_neighbor(std::ifstream& file, quantity<angstrom_unit> cutoff);
  
  std::set<int> parse_lammps_neighbor(std::ifstream& file,
                                      quantity<angstrom_unit> cutoff,
                                      int max_potential_switch,
                                      int max_vasp_nsw,
                                      int& NSW);
  
  std::string make_lammps_data(quantity<angstrom_unit> xhi,
                               quantity<angstrom_unit> yhi,
                               quantity<angstrom_unit> zhi,
                               simulation_cell& sim_cell,
                               bool include_velocity = true);

  std::string make_poscar(const quantity<angstrom_unit>& lattice_constant,
                          const vec3_dimless& a1,
                          const vec3_dimless& a2,
                          const vec3_dimless& a3,
                          simulation_cell& sim_cell);

  quantity<atomic_mass_unit> get_mass_from_POTCAR(std::istream& file);

  bool parse_INCAR(int& NSW, quantity<femtosecond_unit>& POTIM);

  std::vector<std::string> parse_POTCAR();

  void concat_POTCAR(simulation_cell& sim_cell);

  std::vector<std::shared_ptr<lvx::atomic_element_info> >
  create_atomic_catalog(std::map<std::string, std::string>& conf_data);

  std::stringstream replace_all_in_file(std::istream& file,
                                        std::regex rgx,
                                        std::string replacement);

  bool backup_files(const std::vector<std::string>& file_names,
                    int unique_idx, int lavax_iterations);

}

#endif
