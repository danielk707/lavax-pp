#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <regex>
#include <thread>
#include <set>
// #include <cstdlib>
#include <complex>
#include <unordered_map>
#include <boost/units/systems/si/length.hpp>
#include <boost/units/systems/si/io.hpp>
#include <boost/units/base_units/metric/angstrom.hpp>
#include "lavax.hpp"
#include "parse_config_file.hpp"
#include "init_check.hpp"
#include <boost/format.hpp>

namespace lvx {
  quantity<atomic_mass_unit> get_mass_from_POTCAR() {
    std::fstream file("POTCAR");
    std::string line;

    while (std::getline(file, line)) {
      std::regex rgx(".*POMASS\\s*=\\s*(\\d+\\.\\d+).*");
      std::smatch matches;

      if (std::regex_search(line, matches, rgx)) {
        return std::stod(matches[1]) * u;
      }
    }
    return 0.0 * u;
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

  std::stringstream replace_all_in_file(std::istream& file,
                                        std::regex rgx,
                                        std::string replacement) {
    std::string line;
    std::stringstream ss;
    
    while (std::getline(file, line)) {
      if (std::regex_match(line, rgx)) {
        ss << std::regex_replace(line, rgx, replacement) << "\n";
      } else
        ss << line << "\n";
    }
    // file.close();
    
    // file << ss.rdbuf();
    // file.close();
    return ss;
  }

  void update_LAMMPS_script(std::string lammps_potential_file,
                            std::string lammps_atomic_symbol,
                            int vasp_nsw, double vasp_potim) {
    std::fstream f("predictor.in");
    std::stringstream ss;
    ss = replace_all_in_file(f, std::regex("pair_coeff.*"),
      std::string("pair_coeff * * ") + lammps_potential_file + " " + lammps_atomic_symbol);
    // ss.seekg(0);
    ss = replace_all_in_file(ss, std::regex("run.*"),
      std::string("run ") + std::to_string((vasp_nsw+3)*5));
    ss = replace_all_in_file(ss, std::regex("timestep.*"),
      std::string("timestep ") + std::to_string(vasp_potim/5));
    f.close();
    f.open("predictor.in", std::fstream::in | std::fstream::out | std::fstream::trunc);
    f << ss.str();
    f.close();
  }
  
  bool backup_files(const std::vector<std::string>& file_names,
                    int unique_idx, int lvx_iterations) {
    std::string fmt_str = std::string("%1$0") +
      std::to_string(std::to_string(lvx_iterations).length()) + "d";
    boost::format fmt(fmt_str);

    fmt % unique_idx;
    std::string folder_name = std::string("./RUN") + fmt.str();
    std::cout << folder_name << "\n";

    using namespace boost::filesystem;
    
    path p(folder_name);

    if (!exists(p)) 
      create_directory(p);
    
    for (const auto& file_name : file_names)
      copy_file(path(std::string("./") + file_name),
                path(folder_name + "/" + file_name),
                copy_option::overwrite_if_exists);
    return true;
  }

  std::set<int> predict_neighbors(std::string lammps_command, double cut_off_dist) {
    std::string cmd = lammps_command + " -in predictor.in > /dev/null";
    system(cmd.c_str());
    std::ifstream f("neigh.dump");
    return parse_lammps_neighbor(f, cut_off_dist * angstrom);
  }

  std::set<int> predict_neighbors(std::string lammps_command,
                                  double cut_off_dist,
                                  int max_potential_switch,
                                  int max_vasp_nsw,
                                  int count_good_prev,
                                  int& NSW) {
    std::string cmd = lammps_command + " -in predictor.in > /dev/null";
    system(cmd.c_str());
    std::ifstream f("neigh.dump");
    return parse_lammps_neighbor(f, cut_off_dist * angstrom, max_potential_switch,
                                 max_vasp_nsw, count_good_prev, NSW);
  }

  void concat_POTCAR(simulation_cell_v3& sim_cell) {
    std::ofstream writer("POTCAR");
    for (const auto& e : sim_cell.vasp_symbol_count_helper) {
      if (std::get<1>(e) != 0) {
        auto itr = std::find_if(sim_cell.elements_info.begin(),
                                sim_cell.elements_info.end(),
                                [&e] (std::shared_ptr<lvx::atomic_element_info> p) {
                                  return (std::get<2>(e) ? p->vasp_symbol_good :
                                                           p->vasp_symbol_bad) == std::get<0>(e);
                                });
        std::get<2>(e) ? (writer << (*itr)->vasp_potential_file_good) :
                         (writer << (*itr)->vasp_potential_file_bad);
      }
    }
    writer.close();
  }

  std::vector<std::shared_ptr<lvx::atomic_element_info> >
  create_atomic_catalog(std::map<std::string, std::string>& conf_data) {
    std::vector<std::shared_ptr<lvx::atomic_element_info> > atomic_catalog;

    for (int i = 0; conf_data.find(std::string("ATOMIC_SYMBOL_") + std::to_string(i))
           != conf_data.end(); ++i) {
    
      std::shared_ptr<lvx::atomic_element_info> aep(new lvx::atomic_element_info);

      aep->symbol           = conf_data[std::string("ATOMIC_SYMBOL_") + std::to_string(i)];
      aep->mass   = std::stod(conf_data[std::string("ATOMIC_MASS_")   + std::to_string(i)]) * u;
      aep->atom_type        = i+1;
      aep->vasp_symbol_good = conf_data[std::string("VASP_POTENTIAL_SYMBOL_GOOD_") + std::to_string(i)];
      aep->vasp_symbol_bad  = conf_data[std::string("VASP_POTENTIAL_SYMBOL_BAD_")  + std::to_string(i)];
    
      std::ifstream f(conf_data[std::string("VASP_POTENTIAL_FILE_GOOD_")  + std::to_string(i)]);
      aep->vasp_potential_file_good.assign((std::istreambuf_iterator<char>(f)),
                                           std::istreambuf_iterator<char>());
      f.close();
      f.open(conf_data[std::string("VASP_POTENTIAL_FILE_BAD_") + std::to_string(i)]);
      aep->vasp_potential_file_bad.assign((std::istreambuf_iterator<char>(f)),
                                          std::istreambuf_iterator<char>());
      f.close();
      atomic_catalog.push_back(std::move(aep));
    }
    return atomic_catalog;
  }
}

void PRINT_SET(const std::set<int>& s) {
  for (auto e : s)
    std::cout << e << " ";
  std::cout << "\n";
}

// #define DATADIR 

int main(int argc, char *argv[]) {

  lvx::init_check1();

  std::fstream f("lavax.conf");
  
  auto conf_data = lvx::parse_config_file(f);
  
  if (!lvx::init_check2(conf_data["INIT_POSCAR"])) {
    return false;
  }
  
  std::cout << lvx::get_mass_from_POTCAR() << "\n";

  int NSW;
  double POTIM;

  lvx::parse_INCAR(NSW, POTIM);
  std::cout << NSW << " " << POTIM << "\n";

  lvx::update_LAMMPS_script(conf_data["LAMMPS_POTENTIAL_FILE"],
                            conf_data["LAMMPS_ATOMIC_SYMBOL_0"],
                            NSW, POTIM);

  int lvx_iterations = std::stoi(conf_data["LVX_ITERATIONS"]);

  bool use_adaptive_timestep;
  std::istringstream is(conf_data["USE_ADAPTIVE_POTIM"]);
  is >> std::boolalpha >> use_adaptive_timestep;

  auto v = lvx::parse_POTCAR();

  auto atomic_catalog = lvx::create_atomic_catalog(conf_data);
  
  boost::units::quantity<angstrom_unit> latt_const;
  lvx::vec3_dimless       a1, a2, a3;
  std::ifstream           init_poscar(conf_data["INIT_POSCAR"]);
  lvx::simulation_cell_v3 sim_cell;

  sim_cell.elements_info = atomic_catalog;
  
  lvx::parse_init_poscar(init_poscar, latt_const, a1, a2, a3, sim_cell);

  std::cout << latt_const << "\n";
  std::cout << a1 << "\n";
  std::cout << a2 << "\n";
  std::cout << a3 << "\n";

  for (auto& e : sim_cell.particles)
    std::cout << e.getPos() << "\n";

  int count_good_prev = 0;
  
  for (int i = 0; i < lvx_iterations; ++i) {
    if (use_adaptive_timestep) {
      auto itr =
      std::max_element(sim_cell.particles.begin(),
                       sim_cell.particles.end(),
                       [] (const lvx::atomic_particle& a,
                           const lvx::atomic_particle& b) {
                         return lvx::norm(a.getVel()) < lvx::norm(b.getVel());
                       });
      lvx::vec3_velocity max_vel = itr->getVel();
      auto speed = lvx::norm(max_vel);

      double dt = std::min(std::stod(conf_data["MAX_DISTANCE"])/speed.value(),
                           std::stod(conf_data["MAX_POTIM"]));

      std::fstream f("INCAR");
      std::stringstream ss = lvx::replace_all_in_file(f, std::regex(".*POTIM.*"),
                                                      std::string("POTIM = ") + std::to_string(dt));
      f << ss.str(); f.close();
      f.open("predictor.in");
      ss = lvx::replace_all_in_file(f,  std::regex("^timestep.*"),
                                        std::string("timestep ") + std::to_string(dt/5));
      ss = lvx::replace_all_in_file(ss, std::regex(".*run.*"),
                                         std::string("run ") + std::to_string((NSW+3)*5));
      f << ss.str(); f.close();
      std::cout << "POTIM = " << dt << "\n";
    }
    std::ofstream writer("W_crystal.dat");
    std::string content = lvx::make_lammps_data(latt_const * a1(0),
                                                latt_const * a2(1),
                                                latt_const * a3(2),
                                                sim_cell);
    writer << content;
    writer.close();

    std::set<int> si = lvx::predict_neighbors(conf_data["LAMMPS_COMMAND"],
                                              std::stod(conf_data["POTENTIAL_DEPARTURE_DISTANCE"]),
                                              std::stoi(conf_data["MAX_POTENTIAL_SUBSTITUTIONS"]),
                                              std::stoi(conf_data["MAX_VASP_NSW"]),
                                              count_good_prev, NSW);
    count_good_prev = si.size();
    std::cout << "LAMMPS prediction DONE\n";
    std::cout << "Neighbor index: ";
    PRINT_SET(si);
    std::cout << "NSW = " << NSW << "\n";
    std::fstream f("INCAR");
    std::stringstream ss =
      lvx::replace_all_in_file(f, std::regex(".*NSW.*"),
                               std::string("NSW = ") + std::to_string(NSW));
    f << ss.str(); f.close();

    for (int j = 0; j < sim_cell.particles.size(); ++j) {
      sim_cell.particles[j].high_prec = false;
    }
    
    for (auto j : si) {
      sim_cell.particles[j-1].high_prec = true;
    }

    writer.open("POSCAR");
    content = lvx::make_poscar(latt_const, a1, a2, a3, sim_cell);
    writer << content;
    writer.close();

    lvx::concat_POTCAR(sim_cell);

    system(conf_data["VASP_COMMAND"].c_str());
    std::cout << "VASP run DONE\n";

    std::ifstream contcar("CONTCAR");
    lvx::parse_init_poscar(contcar, latt_const, a1, a2, a3, sim_cell);
    contcar.close();

    std::vector<std::string> files = {"XDATCAR", "CONTCAR", "CHG",
                                      "CHGCAR", "DOSCAR", "EIGENVAL",
                                      "OSZICAR", "PCDAT", "vasprun.xml",
                                      "OUTCAR", "INCAR", "WAVECAR",
                                      "IBZKPT", "POSCAR"};
    lvx::backup_files(files, i, lvx_iterations);
    std::cout << "--------------------\n";
  }

#ifdef ACTIVE

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
#endif
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
