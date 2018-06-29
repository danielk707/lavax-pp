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

  class simulation_cell {
  public:
    std::vector<std::string> vasp_potential_symbols;
    std::vector<std::string> vasp_potential_file;
    std::vector<int>         vasp_indices;
    std::vector<int>         lammps_indices;
    std::vector<Particle_v2> particles;
  };
  
  class simulation_cell_v2 {
  public:
    std::vector<atomic_element> elements;
    std::vector<Particle_v2>    particles;

    // int index_by_vasp_good(std::string) {
      
    // }
  };

  std::vector<Particle> make_bcc(double lattice_constant, int I, int J, int K);
  vec3_angstrom read_vec3_angstrom(std::istream& strm);
  vec3_angstrom read_vec3_angstrom(std::string& line);
  Vector3 read_vector3(std::istream&);
  Vector3 read_vector3(std::string&);
  
  std::vector<std::vector<int>> parse_lammps_neighbor(std::ifstream& file);
  
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

  std::string make_poscar(const quantity<angstrom_unit>& lattice_constant,
                          const vec3_dimless& a1,
                          const vec3_dimless& a2,
                          const vec3_dimless& a3,
                          const simulation_cell& sim_cell);
  
  // ------------------------------------------------------------

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

// #ifdef ACTIVE

//   std::vector<lvx::atom_species> w;

//   for (int i = 0; i < v.size()/2; i++) {

//     lvx::atom_species a;
//     a.vasp_good_pot = v[i];
//     a.vasp_bad_pot  = v[i + v.size()/2];
//     std::stringstream ss(a.vasp_good_pot);

//     std::string line;
//     while (std::getline(ss, line)) {
//       std::regex rgx(".*POMASS\\s*=\\s*(\\d+\\.\\d+).*");
//       std::smatch matches;

//       if (std::regex_search(line, matches, rgx)) {
//         a.mass = std::stod(matches[1]) * u;
//         break;
//       }
//     }

//     ss = std::stringstream(a.vasp_good_pot);

//     while (std::getline(ss, line)) {
//       std::regex rgx(".*VRHFIN\\s*=\\s*(\\w+)\\s*:.*");
//       std::smatch matches;

//       if (std::regex_search(line, matches, rgx)) {
//         a.symbol = matches[1];
//         break;
//       }
//     }

//     std::cout << a.symbol << " " << a.mass << "\n";

//     w.push_back(std::move(a));
//   }
//   using namespace lvx;

//   quantity<angstrom_unit> latt_const;
//   vec3_dimless a1, a2, a3;
//   // std::cout << v[0] << "\n";
//   std::ifstream file(conf_data["INIT_POSCAR"]);
//   parse_init_poscar(file,
//                     latt_const, a1, a2, a3, w);

//   std::cout << latt_const << "\n";
//   std::cout << a1 << "\n";
//   std::cout << a2 << "\n";
//   std::cout << a3 << "\n";
//   std::cout << "\n";
//   // std::cout << w[1].bad_particles.size() << "\n";
//   for (auto a : w[0].bad_particles) {
//     std::cout << a.getPos() << "\n";
//   }

//   std::cout << lvx::make_lammps_data(1.0 * angstrom, 1.0 * angstrom, 1.0 * angstrom, w, true) << "\n";

//   bool USE_ADAPTIVE_TIMESTEP;
//   std::istringstream is(conf_data["USE_ADAPTIVE_TIMESTEP"]);
//   is >> std::boolalpha >> USE_ADAPTIVE_TIMESTEP;

//   std::cout << w[0].bad_particles[6].getVel() << "\n";
//   // std::cout << abs(w[0].bad_particles[0].getVel())) << "\n";
  
//   for (int i = 0; i < std::stoi(conf_data["LVX_ITERATIONS"]); i++) {
//     if (USE_ADAPTIVE_TIMESTEP) {
//       std::vector<Particle_v2> temp;
//       for (auto a : w) {
//         temp.insert(temp.begin(), a.bad_particles.begin(), a.bad_particles.end());
//         temp.insert(temp.begin(), a.good_particles.begin(), a.good_particles.end());
//       }
//       auto max_vel = std::max_element(temp.begin(), temp.end(),
//                                       [] (const Particle_v2& p1, const Particle_v2& p2) {
//                                         return norm(p1.getVel()) < norm(p2.getVel());
//                                       });
//       std::cout << max_vel->getVel() << "\n";
                                                                            
//     }
//   }
// #endif
  // vec3_angstrom v1(units::length::angstrom_t(3.0),
  //                  units::length::angstrom_t(3.0),
  //                  units::length::angstrom_t(3.0));

  // vec3_angstrom v2(units::length::angstrom_t(3.0),
  //                  units::length::angstrom_t(3.0),
  //                  units::length::angstrom_t(3.0));

  // std::cout << v1.x()*3.0 << "\n";

  // using namespace boost::units;
  // using namespace boost::units::si;
  // using namespace boost::units::metric;
  // quantity<length> dx(meter * 2.0);

  // using vec3_length = Eigen::Matrix<quantity<length>, 3, 1>;

  // vec3_length v3(2.0*meter, 3.0*meter, 2.0*meter);
  // vec3_length v4(2.0*meter, 3.0*meter, 2.0*meter);

  // Eigen::Matrix<quantity<length>,3,3> A;
  // A << 1*meter, 2*meter, 3*meter,
  //   4*meter, 5*meter, 6*meter,
  //   7*meter, 8*meter, 9*meter;
  
  // std::cout << 3.0*dx << "\n";
  // // std::cout << A*v3 << "\n";

  // auto l = A*v3;
  // // std::cout << l.x() << "\n";

  // // std::thread t( []() {system("ls");} );
  // // t.join();
