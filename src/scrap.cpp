#include "scrap.cpp"

namespace lvx {
  // // ------------------------------------------------------------
  // std::vector<Particle> make_bcc(double lattice_constant, int I, int J, int K) {
  //   double& L = lattice_constant;
  //   std::vector<Particle> list;

  //   int idx = 1;
  //   for (int i = 0; i < I; ++i) {
  //     for (int j = 0; j < J; ++j) {
  //       for (int k = 0; k < K; ++k) {
  //         list.emplace_back(idx, Vector3(i*L, j*L, k*L), Vector3(0,0,0));
  //         ++idx;
  //         list.emplace_back(idx, Vector3((i+0.5)*L, (j+0.5)*L, (k+0.5)*L), Vector3(0,0,0));
  //         ++idx;
  //       }
  //     }
  //   }
  //   return list;
  // }

  // std::string make_poscar(double lattice_constant, Vector3 a1, Vector3 a2, Vector3 a3,
  //                         int I, int J, int K) {
  //   std::stringstream ss;

  //   ss << "BCC Xx " << 2*I*J*K <<"\n";
  //   ss << lattice_constant << "\n";
  //   ss << a1 << "\n";
  //   ss << a2 << "\n";
  //   ss << a3 << "\n";
  //   ss << 2*I*J*K  << "\n";
  //   ss << "Direct\n";

  //   ss << std::fixed << std::setprecision(8);
    
  //   auto P = make_bcc(lattice_constant, I, J, K);
    
  //   for (const auto& p : P) {
  //     ss << p.getPosition() << "\n";
  //   }
      
  //   return ss.str();
  // }

  // void parse_poscar(std::ifstream& poscar, double& lattice_constant,
  //                   Vector3& a1, Vector3& a2, Vector3& a3, std::vector<Particle>& v) {
  //   poscar.ignore(1024, '\n'); // Ignore first line

  //   poscar >> lattice_constant;
  //   a1 = read_vector3(poscar);
  //   a2 = read_vector3(poscar);
  //   a3 = read_vector3(poscar);

  //   for (int i = 0; i < 3; i++)
  //     poscar.ignore(1024, '\n');

  //   std::string line;
  //   // Regex for floating point number:
  //   std::regex reg("[+-]?(?=[.]?[0-9])[0-9]*(?:[.][0-9]*)?(?:[Ee][+-]?[0-9]+)?");

  //   // Read all the position vectors:
  //   std::vector<Vector3> pos_vec;
  //   while (getline(poscar, line)) {
  //     if (std::regex_search(line, reg)) { // If line contains floating point, assume it has an entire vector:
  //       pos_vec.push_back(read_vector3(line));
  //     } else
  //       break; // If we hit a blank line
  //   }

  //   // Read all the velocity vectors:
  //   std::vector<Vector3> vel_vec;
  //   while (getline(poscar, line)) {
  //     if (std::regex_search(line, reg)) {
  //       vel_vec.push_back(read_vector3(line));
  //     } else
  //       break;
  //   }

  //   // Construct all the particles (from positions and velocities) and add them to v:
  //   for (int i = 0; i < pos_vec.size(); i++) {
  //     v.emplace_back(pos_vec[i], vel_vec[i]);
  //   }
  // }

  // std::vector<std::vector<int>> parse_lammps_neighbor(std::ifstream& file) {
  //   std::regex reg("ITEM\\: ENTRIES c_neigh\\[1\\] c_neigh\\[2\\]");
  //   std::string line;
  //   while (getline(file, line)) {
  //     if (std::regex_search(line, reg)) {
  //       std::cout << "FOUND" << "\n";
  //       break;
  //     }
  //   }

  //   std::vector<std::vector<int>> indicies;

  //   int i = 0, i_prev = 0, j = 0;
    
  //   std::vector<int> l;
  //   getline(file, line);
  //   std::stringstream ss(line);
  //   ss >> i; i_prev = i;
  //   ss >> j;
  //   l.push_back(i);
  //   l.push_back(j);
    
  //   while (file.peek() != EOF) {

  //     while (getline(file, line)) {
  //       ss = std::stringstream(line);
  //       ss >> i;
  //       ss >> j;
        
  //       if (i != i_prev) {
  //         i_prev = i;
  //         break;
  //       }
        
  //       l.push_back(j);
  //       i_prev = i;
  //     }
  //     indicies.push_back(l);
  //     l.clear();
  //     l.push_back(i); l.push_back(j);
  //   }
  //   return indicies;
  // }

  void parse_poscar(std::ifstream& poscar,
                    quantity<angstrom_unit>& lattice_constant,
                    vec3_dimless& a1,
                    vec3_dimless& a2,
                    vec3_dimless& a3,
                    std::vector<Particle_v2>& P) {
    
    poscar.ignore(1024, '\n'); // Ignore first line
    double L;
    poscar >> L;
    
    lattice_constant = L * angstrom;
    a1 = read_vec3(poscar, dimensionless);
    a2 = read_vec3(poscar, dimensionless);
    a3 = read_vec3(poscar, dimensionless);

    // std::cout << transform(lattice_constant, a1, a2, a3, a1) << "\n";
    
    bool using_Direct = false;

    std::string line;

    // while (getline(poscar, line)) {
      // std::regex rgx("(\\d+)");
    //   std::smatch matches;
    // }

    // Check if we are using Direct or Cartesian coordinates:
    while (getline(poscar, line)) {
      if (std::regex_match(line, std::regex("C.*"))) {
        using_Direct = false;
        break;
      }
      else if (std::regex_match(line, std::regex("D.*"))) {
        using_Direct = true;
        break;
      }
    }

    // Regex for floating point number:
    std::regex reg("[+-]?(?=[.]?[0-9])[0-9]*(?:[.][0-9]*)?(?:[Ee][+-]?[0-9]+)?");

    // Read all the position vectors:
    std::vector<vec3_angstrom> pos_vec;
    while (getline(poscar, line)) {
      if (std::regex_search(line, reg)) {
        if (using_Direct) {
          vec3_dimless v = read_vec3(line, dimensionless);
          pos_vec.push_back(transform(lattice_constant, a1, a2, a3, v));
        } else
          pos_vec.push_back(read_vec3(line, angstrom));
      } else
        break; // If we hit a blank line
    }

    // Read all the velocity vectors:
    std::vector<vec3_velocity> vel_vec;
    while (getline(poscar, line)) {
      if (std::regex_search(line, reg)) {
        vel_vec.push_back(read_vec3(line, angstrom_per_fs));
      } else
        break;
    }

    // Construct all the particles (from positions and velocities) and add them to P:
    for (int i = 0; i < pos_vec.size(); i++) {
      P.emplace_back(pos_vec[i], vel_vec[i]);
    }
    
  }

  void parse_init_poscar(std::ifstream& poscar,
                         quantity<angstrom_unit>& lattice_constant,
                         vec3_dimless& a1,
                         vec3_dimless& a2,
                         vec3_dimless& a3,
                         std::vector<atom_species>& P) {
    _parse_poscar_header(poscar, lattice_constant, a1, a2, a3);

    std::string line;

    while (getline(poscar, line)) {
      if (std::regex_match(line, std::regex("\\s*\\d+[\\s\\d]*"))) {
        break;
      }
    }

    std::regex rgx("(\\d+)");
    std::vector<int> num_atom_species;
    
    std::for_each(std::sregex_iterator(line.begin(), line.end(), rgx),
                  std::sregex_iterator(),
                  [&num_atom_species] (const std::smatch& m) {
                    num_atom_species.push_back(std::stoi(m.str(1)));
                    // std::cout << std::stoi(m.str(1)) << "\n";
                  });

    // for (auto i : num_atom_species) {
    //   std::cout << i << "\n";
    // }

    // while (getline(poscar, line)) {
    // std::d::regex rgx("(\\d+)");
    //   std::smatch matches;
    // }

    // Check if we are using Direct or Cartesian coordinates:
    bool using_Direct = false;
    while (getline(poscar, line)) {
      if (std::regex_match(line, std::regex("C.*"))) {
        using_Direct = false;
        break;
      }
      else if (std::regex_match(line, std::regex("D.*"))) {
        using_Direct = true;
        break;
      }
    }

    // Regex for floating point number:
    std::regex reg("[+-]?(?=[.]?[0-9])[0-9]*(?:[.][0-9]*)?(?:[Ee][+-]?[0-9]+)?");

    // Read all the position vectors:
    std::vector<vec3_angstrom> pos_vec;
    while (getline(poscar, line)) {
      if (std::regex_search(line, reg)) {
        if (using_Direct) {
          vec3_dimless v = read_vec3(line, dimensionless);
          pos_vec.push_back(transform(lattice_constant, a1, a2, a3, v));
        } else
          pos_vec.push_back(read_vec3(line, angstrom));
      } else
        break; // If we hit a blank line
    }

    // Read all the velocity vectors:
    std::vector<vec3_velocity> vel_vec;
    while (getline(poscar, line)) {
      if (std::regex_search(line, reg)) {
        vel_vec.push_back(read_vec3(line, angstrom_per_fs));
      } else
        break;
    }

    // for (auto a : P) {
    int k = 0;
    for (int i = 0; i < num_atom_species.size(); ++i) {
      for (int j = k; j < k+num_atom_species[i]; ++j) {
        P[i].soft_particles.emplace_back(pos_vec[j], vel_vec[j]);
      }
      k += num_atom_species[i];
    }
        
    // }

    // // Construct all the particles (from positions and velocities) and add them to P:
    // for (int i = 0; i < pos_vec.size(); i++) {
    //   P.emplace_back(pos_vec[i], vel_vec[i]);
    // }
    
  }

  void parse_init_poscar(std::ifstream& poscar,
                         quantity<angstrom_unit>& lattice_constant,
                         vec3_dimless& a1,
                         vec3_dimless& a2,
                         vec3_dimless& a3,
                         simulation_cell_v2& sim_cell) {
    
    _parse_poscar_header(poscar, lattice_constant, a1, a2, a3);

    // std::smatch matches;
    std::map<std::string,int> dict;

    for (int i = 0; i < sim_cell.elements.size(); ++i) {
      dict[sim_cell.elements[i].vasp_symbol_hard] = i;
      dict[sim_cell.elements[i].vasp_symbol_soft] = i;
    }

    std::vector<std::string> vasp_symbols;
    
    std::string line;
    while (getline(poscar, line)) {
      if (std::regex_match(line, std::regex("\\s*\\d+[\\s\\d]*"))) {
        break;
      }
      // This can be simplified if you think about it:
      if (std::regex_match(line, std::regex("\\w+.*"))) {
        std::regex rgx("\\w+");
        for( std::sregex_iterator it(line.begin(), line.end(), rgx), it_end;
             it != it_end; ++it ) {
          vasp_symbols.push_back((*it)[0]);
          std::cout << "DEBUGGG " << (*it)[0] << "\n";
        }
        // continue;
      }
    }
    std::regex rgx("(\\d+)");
    std::vector<int> num_atom_species;
    
    std::for_each(std::sregex_iterator(line.begin(), line.end(), rgx),
                  std::sregex_iterator(),
                  [&num_atom_species] (const std::smatch& m) {
                    num_atom_species.push_back(std::stoi(m.str(1)));
                    // std::cout << std::stoi(m.str(1)) << "\n";
                  });

    // for (int i = 0; i < num_atom_species.size(); ++i) {
    //   atomic_element ae;
    //   ae.vasp_indices.push_back(num_atom_species[i]);
    //   sim_cell.elements.push_back(ae);
    // }
    
    std::cout << vasp_symbols.size() << "\n";
    std::cout << num_atom_species.size() << "\n";
    
    int j = 0;
    for (auto s : vasp_symbols) {
      if (s == sim_cell.elements[dict[s]].vasp_symbol_soft) {
        sim_cell.elements[dict[s]].vasp_num_soft = num_atom_species[j];
      }
      else if (s == sim_cell.elements[dict[s]].vasp_symbol_hard) {
        sim_cell.elements[dict[s]].vasp_num_hard = num_atom_species[j];
      }
      ++j;
    }

    // Check if we are using Direct or Cartesian coordinates:
    bool using_Direct = false;
    while (getline(poscar, line)) {
      if (std::regex_match(line, std::regex("C.*"))) {
        using_Direct = false;
        break;
      }
      else if (std::regex_match(line, std::regex("D.*"))) {
        using_Direct = true;
        break;
      }
    }

    // Regex for floating point number:
    std::regex reg("[+-]?(?=[.]?[0-9])[0-9]*(?:[.][0-9]*)?(?:[Ee][+-]?[0-9]+)?");

    // Read all the position vectors:
    std::vector<vec3_angstrom> pos_vec;
    while (getline(poscar, line)) {
      if (std::regex_search(line, reg)) {
        if (using_Direct) {
          vec3_dimless v = read_vec3(line, dimensionless);
          pos_vec.push_back(transform(lattice_constant, a1, a2, a3, v));
        } else
          pos_vec.push_back(read_vec3(line, angstrom));
      } else
        break; // If we hit a blank line
    }

    // Read all the velocity vectors:
    std::vector<vec3_velocity> vel_vec;
    while (getline(poscar, line)) {
      if (std::regex_search(line, reg)) {
        vel_vec.push_back(read_vec3(line, angstrom_per_fs));
      } else
        break;
    }

    sim_cell.particles.clear();
    for (int i = 0; i < pos_vec.size(); i++) {
      sim_cell.particles.emplace_back(pos_vec[i], vel_vec[i]);
    }
  }

  std::string make_poscar(const quantity<angstrom_unit>& lattice_constant,
                          const vec3_dimless& a1,
                          const vec3_dimless& a2,
                          const vec3_dimless& a3,
                          const simulation_cell& sim_cell) {
    std::stringstream ss;
    // ss << std::fixed;
    
    ss << "BCC Xx \n";
    ss << lattice_constant.value() << "\n";
    ss << a1 << "\n";
    ss << a2 << "\n";
    ss << a3 << "\n";

    for (int i = 0; i < sim_cell.vasp_indices.size(); ++i) {
      if (sim_cell.vasp_indices[i] > 0) {
        ss << sim_cell.vasp_potential_symbols[i] << " ";
      }
    }
    ss << "\n";
    for (int i = 0; i < sim_cell.vasp_indices.size(); ++i) {
      if (sim_cell.vasp_indices[i] > 0) {
        ss << sim_cell.vasp_indices[i] << " ";
      }
    }
    ss << "\n";
    
    // ss << std::distance(first_atom, last_atom) << "\n";
    // ss << ((typeid(first_atom->getPos()) == typeid(vec3_angstrom)) ? "Cartesian\n" : "Direct\n");
    ss << ((typeid(sim_cell.particles[0].getPos()) == typeid(vec3_angstrom)) ? "Cartesian\n" : "Direct\n");

    for (const auto& p : sim_cell.particles) {
      ss << p.getPos() << "\n";
    }
    ss << "\n";
    for (const auto& p : sim_cell.particles) {
      ss << p.getVel() << "\n";
    }
    
    return ss.str();
  }

  std::string make_lammps_data(quantity<angstrom_unit> xhi,
                               quantity<angstrom_unit> yhi,
                               quantity<angstrom_unit> zhi,
                               std::vector<atom_species>& P, bool incl_vel) {
    std::stringstream ss;
    ss.precision(8);
    ss << std::fixed;

    ss << "# W Crystal bcc #\n\n";

    int num_part = 0;
    std::for_each(P.begin(), P.end(),
                  [&num_part] (const atom_species& a) {
                    num_part += a.soft_particles.size() + a.hard_particles.size();
                  });

    ss << num_part << " atoms\n";
    ss << P.size() << " atom types\n\n";

    ss << 0.0 << " " << xhi.value() << " xlo xhi\n";
    ss << 0.0 << " " << yhi.value() << " ylo yhi\n";
    ss << 0.0 << " " << zhi.value() << " zlo zhi\n";

    ss << "\nMasses\n\n";
    for (int i = 0; i < P.size(); i++) {
      ss << i+1 << " " << P[i].mass.value() << "\n";
    }
    
    ss << "\nAtoms\n\n";
    int k = 1;
    for (int i = 0; i < P.size(); i++) {
      for (int j = 0; j < P[i].soft_particles.size(); ++j) {
        ss << k << " " << i+1 << " " << P[i].soft_particles[j].getPos() << "\n";
        k++;
      }
      for (int j = 0; j < P[i].hard_particles.size(); ++j) {
        ss << k << " " << i+1 << " " << P[i].hard_particles[j].getPos() << "\n";
        k++;
      }
    }

    if (incl_vel) {
      ss << "\nVelocities\n\n";

      k = 1;
      for (int i = 0; i < P.size(); i++) {
        for (int j = 0; j < P[i].soft_particles.size(); ++j) {
          ss << k << " " << P[i].soft_particles[j].getVel() << "\n";
          k++;
        }
        for (int j = 0; j < P[i].hard_particles.size(); ++j) {
          ss << k << " " << P[i].hard_particles[j].getVel() << "\n";
          k++;
        }
      }
    }
    
    return ss.str();
  }

  std::string make_lammps_data(quantity<angstrom_unit> xhi,
                               quantity<angstrom_unit> yhi,
                               quantity<angstrom_unit> zhi,
                               simulation_cell_v2& sim_cell,
                               bool include_velocity) {
    std::stringstream ss;
    ss.precision(8);
    ss << std::fixed;

    ss << "# W Crystal bcc #\n\n";

    ss << sim_cell.particles.size() << " atoms\n";
    ss << sim_cell.elements.size() << " atom types\n\n";

    ss << 0.0 << " " << xhi.value() << " xlo xhi\n";
    ss << 0.0 << " " << yhi.value() << " ylo yhi\n";
    ss << 0.0 << " " << zhi.value() << " zlo zhi\n";

    ss << "\nMasses\n\n";
    for (int i = 0; i < sim_cell.elements.size(); ++i) {
      ss << i+1 << " " << sim_cell.elements[i].mass.value() << "\n";
    }

    ss << "\nAtoms\n\n";
    // int k = 1;
    // for (int i = 0; i < sim_cell.elements.size(); ++i) {
    //   sim_cell.elements[i].number_atoms()
    // }

    // int j_upper = sim_cell.elements[0].number
    std::vector<int> num_atoms;
    for (int i = 0; i < sim_cell.elements.size(); ++i) {
      num_atoms.push_back(sim_cell.elements[i].number_atoms());
    }
    std::vector<int> indices;
    std::partial_sum(num_atoms.begin(), num_atoms.end(),
                     std::back_inserter(indices));

    std::cout << "INDICES: ";
    for (auto i : num_atoms)
      std::cout << i << " ";
    std::cout << "\n";
    
    for (int i = 0; i < sim_cell.particles.size(); ++i) {
      int j = 0;
      if (i >= indices[j]) {
        j++;
      }
      
      ss << i+1 << " " << j+1 << " " << sim_cell.particles[i].getPos() << "\n";
    }

    if (include_velocity) {
      ss << "\nVelocities\n\n";
      for (int i = 0; i < sim_cell.particles.size(); ++i) {
        // int j = 0;
        // if (i >= indices[j])
        //   ++j;
      
        ss << i+1 << " " << sim_cell.particles[i].getVel() << "\n";
      }
    }
    
    return ss.str();
  }

  std::ostream& operator<<(std::ostream &strm, const Vector3& a) {
    return strm << std::setprecision(8) << a.x << " " << a.y << " " << a.z;
  }

  double abs(const Vector3& v) {
    return sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
  }

  vec3_angstrom read_vec3_angstrom(std::istream& strm) {
    double x, y, z;
    strm >> x;
    strm >> y;
    strm >> z;

    vec3_angstrom rtn;
    rtn = x * angstrom,
          y * angstrom,
          z * angstrom;

    return rtn;
  }

  vec3_angstrom read_vec3_angstrom(std::string& line) {
    std::stringstream ss(line);
    return read_vec3_angstrom(ss);
  }
}
