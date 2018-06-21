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
