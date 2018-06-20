#ifndef LAVAX_INC // Inclusion gaurd
#define LAVAX_INC

#include <cmath>
#include <sstream>
#include <fstream>
#include <set>
#include <iomanip>
#include <algorithm>
#include <regex>
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

  // ------------------------------------------------------------
  class Vector3 {
  private:
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
  public:
    Vector3(double x, double y, double z) : x(x), y(y), z(z) {}
    Vector3(){};

    // explicit operator std::string() const {
    //   return std::to_string(42);
    // }

    Vector3 operator-(const Vector3& b) {
      return Vector3(x-b.x, y-b.y, z-b.z);
    }
    
    Vector3 operator+(const Vector3& b) {
      return Vector3(x+b.x, y+b.y, z+b.z);
    }

    Vector3 operator*(double b) {
      return Vector3(x*b, y*b, z*b);
    }

    double operator*(const Vector3& b) {
      return x*b.x + y*b.y + z*b.z;
    }
    
    friend double abs(const Vector3&);
    friend std::ostream& operator<<(std::ostream&, const Vector3&);
  };

  std::ostream& operator<<(std::ostream &strm, const Vector3& a);
  double abs(const Vector3& v);
  
  // ------------------------------------------------------------
  class Particle {
  private:
    Vector3 pos;
    Vector3 vel;

    std::vector<int> neigh_idx;
    const unsigned int idx;
    
  public:
    Particle(Vector3 pos, Vector3 vel) : idx(0), pos(pos), vel(vel) {}
    Particle(int idx, Vector3 pos, Vector3 vel) : idx(idx), pos(pos), vel(vel) {}

    void update(Vector3 pos, Vector3 vel, const std::vector<int>& neighbors) {
      update(pos, vel);
      neigh_idx = neighbors;
    }

    void update(Vector3 pos, Vector3 vel) {
      this->pos = pos;
      this->vel = vel;
    }
    
    Vector3 getPosition() const { return pos; }
    Vector3 getVelocity() const { return vel; }
  
  };


  template <typename T>
  std::ostream& operator<<(std::ostream &strm,
                           const blitz::TinyVector<quantity<T>,3>& a) {
    return strm << std::fixed << std::setprecision(8)
                << a(0).value() << " " << a(1).value() << " " << a(2).value();
  }

  // ------------------------------------------------------------
  class Particle_v2 {
    
  private:
    vec3_angstrom pos;
    vec3_velocity vel;
    unsigned int type = 0;
  public:
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

  // ------------------------------------------------------------
  class Vector3;
  Vector3 read_vector3(std::istream&);
  Vector3 read_vector3(std::string&);
  vec3_angstrom read_vec3_angstrom(std::istream& strm);
  vec3_angstrom read_vec3_angstrom(std::string& line);

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
  
  
  std::vector<Particle> make_bcc(double lattice_constant, int I, int J, int K);

  enum class crystal_structure { BCC, FCC };
  
  std::vector<Particle_v2> create_crystal(crystal_structure cstruct,
                                          quantity<angstrom_unit> lattice_constant,
                                          int I, int J, int K);
  
  void parse_poscar(std::ifstream& poscar, double& lattice_constant,
                    Vector3& a1, Vector3& a2, Vector3& a3, std::vector<Particle>& v);

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


  std::set<int> parse_lammps_neighbor(std::ifstream& file, quantity<angstrom_unit> cutoff);
                    
  
  std::vector<std::vector<int>> parse_lammps_neighbor(std::ifstream& file);
  
  std::string make_poscar(double lattice_constant, Vector3 a1, Vector3 a2, Vector3 a3,
                          int I, int J, int K);

  // ------------------------------------------------------------
  template<class InputIt>
  std::string make_lammps_data(double mass, double xhi, double yhi, double zhi,
                               InputIt first_atom, InputIt last_atom, bool incl_vel = true) {
    std::stringstream ss;
    ss.precision(8);
    ss << std::fixed;

    ss << "# W Crystal bcc #\n\n";
    
    ss << std::distance(first_atom, last_atom) << " atoms\n";
    ss << 1 << " atom types\n\n";

    ss << 0.0 << " " << xhi << " xlo xhi\n";
    ss << 0.0 << " " << yhi << " ylo yhi\n";
    ss << 0.0 << " " << zhi << " zlo zhi\n";

    ss << "\nMasses\n\n";
    ss << 1 << " " << mass << "\n";

    ss << "\nAtoms\n\n";

    int i = 1;
    std::for_each(first_atom, last_atom, [&ss, &i] (const Particle& p) {
        ss << i << " " << 1 << " " << p.getPosition() << "\n";
        ++i;
      });

    if (incl_vel) {
      ss << "\nVelocities\n\n";
  
      i = 1;
      std::for_each(first_atom, last_atom, [&ss, &i] (const Particle& p) {
          ss << i << " " << p.getVelocity() << "\n";
          ++i;
        });
    }

    return ss.str();
  }
  
  template<class InputIt>
  std::string make_poscar(double lattice_constant, Vector3 a1, Vector3 a2, Vector3 a3,
                          InputIt first_atom, InputIt last_atom) {
    std::stringstream ss;

    ss << "BCC Xx \n";
    ss << lattice_constant << "\n";
    ss << a1 << "\n";
    ss << a2 << "\n";
    ss << a3 << "\n";
    ss << std::distance(first_atom, last_atom) << "\n";
    ss << "Direct\n";

    std::for_each(first_atom, last_atom, [&ss] (const Particle& p) {
        ss << p.getPosition() << "\n";
      });
  
    ss << "\n";

    std::for_each(first_atom, last_atom, [&ss] (const Particle& p) {
        ss << p.getVelocity() << "\n";
      });
  
    return ss.str();
  }

  template<class InputIt>
  std::string make_lammps_data(quantity<atomic_mass_unit> mass,
                               quantity<angstrom_unit> xhi,
                               quantity<angstrom_unit> yhi,
                               quantity<angstrom_unit> zhi,
                               InputIt first_atom, InputIt last_atom, bool incl_vel = true) {
    std::stringstream ss;
    ss.precision(8);
    ss << std::fixed;

    ss << "# W Crystal bcc #\n\n";
    
    ss << std::distance(first_atom, last_atom) << " atoms\n";
    ss << 1 << " atom types\n\n";

    ss << 0.0 << " " << xhi.value() << " xlo xhi\n";
    ss << 0.0 << " " << yhi.value() << " ylo yhi\n";
    ss << 0.0 << " " << zhi.value() << " zlo zhi\n";

    ss << "\nMasses\n\n";
    ss << 1 << " " << mass.value() << "\n";

    ss << "\nAtoms\n\n";

    int i = 1;
    std::for_each(first_atom, last_atom, [&ss, &i] (const Particle_v2& p) {
        ss << i << " " << 1 << " " << p.getPos() << "\n";
        ++i;
      });

    if (incl_vel) {
      ss << "\nVelocities\n\n";
  
      i = 1;
      std::for_each(first_atom, last_atom, [&ss, &i] (const Particle_v2& p) {
          ss << i << " " << p.getVel() << "\n";
          ++i;
        });
    }

    return ss.str();
  }

  std::string make_lammps_data(quantity<angstrom_unit> xhi,
                               quantity<angstrom_unit> yhi,
                               quantity<angstrom_unit> zhi,
                               std::vector<atom_species>& P, bool incl_vel = true);

  std::string make_lammps_data(quantity<angstrom_unit> xhi,
                               quantity<angstrom_unit> yhi,
                               quantity<angstrom_unit> zhi,
                               simulation_cell_v2& sim_cell,
                               bool include_velocity = true);
  
  template<class InputIt>
  std::string make_poscar(const quantity<angstrom_unit>& lattice_constant,
                          const vec3_dimless& a1,
                          const vec3_dimless& a2,
                          const vec3_dimless& a3,
                          InputIt first_atom, InputIt last_atom) {
    std::stringstream ss;
    // ss << std::fixed;
    
    ss << "BCC Xx \n";
    ss << lattice_constant.value() << "\n";
    ss << a1 << "\n";
    ss << a2 << "\n";
    ss << a3 << "\n";
    ss << std::distance(first_atom, last_atom) << "\n";
    ss << ((typeid(first_atom->getPos()) == typeid(vec3_angstrom)) ? "Cartesian\n" : "Direct\n");

    std::for_each(first_atom, last_atom, [&ss] (const Particle_v2& p) {
        ss << p.getPos() << "\n";
      });
  
    ss << "\n";

    std::for_each(first_atom, last_atom, [&ss] (const Particle_v2& p) {
        ss << p.getVel() << "\n";
      });
  
    return ss.str();
  }

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
