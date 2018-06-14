#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <regex>
#include <thread>
// #include <cstdlib>
#include <complex>
#include <unordered_map>
#include <boost/units/systems/si/length.hpp>
#include <boost/units/systems/si/io.hpp>
#include <boost/units/base_units/metric/angstrom.hpp>
#include "lavax.hpp"
#include "parse_config_file.hpp"
#include "init_check.hpp"

namespace lvx {
  double get_mass_from_POTCAR() {
    std::fstream file("POTCAR");
    std::string line;

    while (std::getline(file, line)) {
      std::regex rgx(".*POMASS\\s*=\\s*(\\d+\\.\\d+).*");
      std::smatch matches;

      if (std::regex_search(line, matches, rgx)) {
        return std::stod(matches[1]);
      }
    }
    return 0.0;
  }

  bool parse_INCAR(int& NSW, double& POTIM) {
    std::fstream file("INCAR");
    std::regex rgx1(".*NSW\\s*=\\s*(\\d+).*");
    std::regex rgx2(".*POTIM\\s*=\\s*(\\d+\\.\\d+).*");

    std::string line;
    bool b1 = false;
    bool b2 = false;

    while (std::getline(file, line)) {
      std::smatch matches;

      if (std::regex_search(line, matches, rgx1)) {
        NSW = std::stoi(matches[1]);
        b1 = true;
        // std::cout << std::stoi(matches[1]) << "\n";
      }

      if (std::regex_search(line, matches, rgx2)) {
        POTIM = std::stod(matches[1]);
        b2 = true;
        // std::cout << std::stod(matches[1]) << "\n";
      }
    }
    return b1 && b2;
  }

  std::vector<std::string> parse_POTCAR() {
    std::fstream file("POTCAR");
    std::vector<std::string> buffers;
    std::string line;
    
    while (std::getline(file, line)) {
      std::stringstream ss;
      ss << line << '\n';
      
      while (std::getline(file, line)) {
        ss << line << '\n';
        if (std::regex_match(line, std::regex(".*End of Dataset.*")))
          break;
      }
      buffers.push_back(ss.str());
    }
    return buffers;
  }

  void replace_all_in_file(std::fstream& file, std::regex rgx, std::string replacement) {
    std::string line;
    std::stringstream ss;
    
    while (std::getline(file, line)) {
      if (std::regex_match(line, rgx)) {
        ss << std::regex_replace(line, rgx, replacement) << "\n";
      } else
        ss << line << "\n";
    }
    // file.close();
    
    file << ss.rdbuf();
    file.close();
  }
  
  
}

// #define DATADIR 

int main(int argc, char *argv[]) {

  lvx::init_check1();

  std::fstream f("lavax.conf");
  
  auto conf_data = lvx::parse_config_file(f);
  
  if (!lvx::init_check2(conf_data["INIT_POSCAR"])) {
    return false;
  }

  std::cout << lvx::get_mass_from_POTCAR()  << "\n";

  int NSW;
  double POTIM;

  lvx::parse_INCAR(NSW, POTIM);
  std::cout << NSW << " " << POTIM << "\n";

  auto v = lvx::parse_POTCAR();

  std::vector<lvx::atom_species> w;

  for (int i = 0; i < v.size()/2; i++) {

    lvx::atom_species a;
    a.vasp_good_pot = v[i];
    a.vasp_bad_pot  = v[i + v.size()/2];
    std::stringstream ss(a.vasp_good_pot);

    std::string line;
    while (std::getline(ss, line)) {
      std::regex rgx(".*POMASS\\s*=\\s*(\\d+\\.\\d+).*");
      std::smatch matches;

      if (std::regex_search(line, matches, rgx)) {
        a.mass = std::stod(matches[1]) * u;
        break;
      }
    }

    ss = std::stringstream(a.vasp_good_pot);

    while (std::getline(ss, line)) {
      std::regex rgx(".*VRHFIN\\s*=\\s*(\\w+)\\s*:.*");
      std::smatch matches;

      if (std::regex_search(line, matches, rgx)) {
        a.symbol = matches[1];
        break;
      }
    }

    std::cout << a.symbol << " " << a.mass << "\n";

    w.push_back(std::move(a));
  }
  using namespace lvx;

  quantity<angstrom_unit> latt_const;
  vec3_dimless a1, a2, a3;
  // std::cout << v[0] << "\n";
  std::ifstream file(conf_data["INIT_POSCAR"]);
  parse_init_poscar(file,
                    latt_const, a1, a2, a3, w);

  std::cout << latt_const << "\n";
  std::cout << a1 << "\n";
  std::cout << a2 << "\n";
  std::cout << a3 << "\n";
  std::cout << "\n";
  // std::cout << w[1].bad_particles.size() << "\n";
  for (auto a : w[0].bad_particles) {
    std::cout << a.getPos() << "\n";
  }

  std::cout << lvx::make_lammps_data(1.0 * angstrom, 1.0 * angstrom, 1.0 * angstrom, w, true) << "\n";

  bool USE_ADAPTIVE_TIMESTEP;
  std::istringstream is(conf_data["USE_ADAPTIVE_TIMESTEP"]);
  is >> std::boolalpha >> USE_ADAPTIVE_TIMESTEP;

  std::cout << w[0].bad_particles[6].getVel() << "\n";
  // std::cout << abs(w[0].bad_particles[0].getVel())) << "\n";
  
  for (int i = 0; i < std::stoi(conf_data["LVX_ITERATIONS"]); i++) {
    if (USE_ADAPTIVE_TIMESTEP) {
      std::vector<Particle_v2> temp;
      for (auto a : w) {
        temp.insert(temp.begin(), a.bad_particles.begin(), a.bad_particles.end());
        temp.insert(temp.begin(), a.good_particles.begin(), a.good_particles.end());
      }
      auto max_vel = std::max_element(temp.begin(), temp.end(),
                                      [] (const Particle_v2& p1, const Particle_v2& p2) {
                                        return norm(p1.getVel()) < norm(p2.getVel());
                                      });
      std::cout << max_vel->getVel() << "\n";
                                                                            
    }
  }
  
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
  
  // std::cin.get();
  
  // std::ifstream poscar("POSCAR");

  // double latt_con;
  // Vector3 a1, a2, a3;
  // std::vector<Particle> particles;
  // parse_poscar(poscar, latt_con, a1, a2, a3, particles);

  // auto P = make_bcc(3.18, 4, 4, 4);

  // std::ofstream lmp_data_file("W_crystal.dat");
  // lmp_data_file << make_lammps_data(183.85, 3.18*4, 3.18*4, 3.18*4, P.begin(), P.end());
  // lmp_data_file.close();

  // // std::thread t1( [] () { system("} );
  // // system("mpirun -np 4 lmp_mpi -in calc_neigh.in");

  // // std::ifstream neigh_file("neigh.0");
  // // auto neigh_list = qc::parse_lammps_neighbor(neigh_file);
  

  // std::vector<Vector3> latt_vec;
  // latt_vec.push_back(a1);
  // latt_vec.push_back(a2);
  // latt_vec.push_back(a3);
  
  // std::cout << "-------------------- make_poscar --------------------\n";
  // std::cout << make_poscar(latt_con, latt_vec[0], latt_vec[1], latt_vec[2],
  //                          particles.begin(), particles.end()) << "\n";

  // std::cout << "-------------------- make_lammps_data --------------------\n";
  // std::cout << make_lammps_data(183.85, 3, 3, 3, particles.begin(), particles.end()) << "\n";

  // // std::thread t( []() {system("ls");} );

  // // t.join();

  // std::ifstream neigh("neigh.0");

  // auto vec = parse_lammps_neighbor(neigh);

  // std::for_each(vec.begin(), vec.end(),
  //               [] (const std::vector<int>& v) {
  //                 std::for_each(v.begin(), v.end(), [] (int i) { std::cout << i << " "; });
  //                 std::cout << "\n";
  //               });
  // // std::cout << static_cast<std::string>(v) << "\n";
  return 0;
}
