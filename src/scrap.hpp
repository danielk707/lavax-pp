// namespace Eigen {
// template<> struct NumTraits<units::length::angstrom_t>
//  : NumTraits<double> // permits to get the epsilon, dummy_precision, lowest, highest functions
// {
//   typedef units::length::angstrom_t Real;
//   typedef units::length::angstrom_t NonInteger;
//   typedef units::length::angstrom_t Nested;
//   enum {
//     IsComplex = 0,
//     IsInteger = 0,
//     IsSigned = 1,
//     RequireInitialization = 1,
//     ReadCost = 1,
//     AddCost = 3,
//     MulCost = 3
//   };
// };

//   template<> struct NumTraits<boost::units::quantity<boost::units::si::length>>
//  : NumTraits<double> // permits to get the epsilon, dummy_precision, lowest, highest functions
// {
//   typedef units::length::angstrom_t Real;
//   typedef units::length::angstrom_t NonInteger;
//   typedef units::length::angstrom_t Nested;
//   enum {
//     IsComplex = 0,
//     IsInteger = 0,
//     IsSigned = 1,
//     RequireInitialization = 1,
//     ReadCost = 1,
//     AddCost = 3,
//     MulCost = 3
//   };
// };
// }

namespace lvx {

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

  std::vector<Particle> make_bcc(double lattice_constant, int I, int J, int K);
  vec3_angstrom read_vec3_angstrom(std::istream& strm);
  vec3_angstrom read_vec3_angstrom(std::string& line);
  Vector3 read_vector3(std::istream&);
  Vector3 read_vector3(std::string&);

  void parse_poscar(std::ifstream& poscar, double& lattice_constant,
                    Vector3& a1, Vector3& a2, Vector3& a3, std::vector<Particle>& v);

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

}
