#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <iterator>
#include <regex>
#include <thread>
#include <chrono>
#include <set>
// #include <cstdlib>
#include <complex>
#include <unordered_map>
#include <boost/units/systems/si/length.hpp>
#include <boost/units/systems/si/io.hpp>
#include <boost/units/base_units/metric/angstrom.hpp>
#include <boost/filesystem.hpp>
#include "lavax.hpp"
#include "parse_config_file.hpp"
#include "init_check.hpp"

namespace lvx {

  void update_LAMMPS_script(std::string lammps_potential_file,
                            std::string lammps_atomic_symbol,
                            std::string lammps_pair_style,
                            int vasp_nsw,
                            quantity<femtosecond_unit> vasp_potim) {
    
    std::fstream f("predictor.in");
    std::stringstream ss;
    ss = replace_all_in_file(f, std::regex("pair_coeff.*"),
      std::string("pair_coeff * * ") + lammps_potential_file + " " + lammps_atomic_symbol);
    // ss.seekg(0);
    ss = replace_all_in_file(ss, std::regex("pair_style.*"),
      std::string("pair_style ") + lammps_pair_style);
    // ss = replace_all_in_file(ss, std::regex("run.*"),
    //   std::string("run ") + std::to_string((vasp_nsw+3)*5));
    // ss = replace_all_in_file(ss, std::regex("timestep.*"),
    //   std::string("timestep ") + std::to_string(vasp_potim.value()/5));
    ss = replace_all_in_file(ss, std::regex("run.*"),
      std::string("run ") + std::to_string(vasp_nsw));
    ss = replace_all_in_file(ss, std::regex("timestep.*"),
      std::string("timestep ") + std::to_string(vasp_potim.value()));
    f.close();
    f.open("predictor.in", std::fstream::in | std::fstream::out | std::fstream::trunc);
    f << ss.str();
    f.close();
  }

  std::set<int> predict_neighbors(std::string lammps_command,
                                  bool lammps_hide_output,
                                  double cut_off_dist) {
    
    std::string cmd = lammps_command + " -in predictor.in" +
      (lammps_hide_output ? " > /dev/null" : "");
    system(cmd.c_str());
    std::ifstream f("neigh.dump");
    return parse_lammps_neighbor(f, cut_off_dist * angstrom);
  }

  std::set<int> predict_neighbors(std::string lammps_command,
                                  double cut_off_dist,
                                  bool lammps_hide_output,
                                  int max_potential_switch,
                                  int max_vasp_nsw,
                                  int& NSW) {
    
    std::string cmd = lammps_command + " -in predictor.in" +
      (lammps_hide_output ? " > /dev/null" : "");
    system(cmd.c_str());
    std::ifstream f("neigh.dump");
    return parse_lammps_neighbor(f, cut_off_dist * angstrom, max_potential_switch,
                                 max_vasp_nsw, NSW);
  }

  struct config_data {
    std::string VASP_COMMAND;
    std::string LAMMPS_COMMAND;
    std::string INIT_POSCAR;
    int         LAVAX_ITERATIONS;
    std::string LAMMPS_POTENTIAL_FILE;
    std::string LAMMPS_PAIR_STYLE;

    quantity<angstrom_unit> POTENTIAL_DEPARTURE_DISTANCE;
  
    bool                       USE_ADAPTIVE_POTIM;
    quantity<angstrom_unit>    MAX_DISTANCE;
    quantity<femtosecond_unit> MAX_POTIM;

    int MAX_POTENTIAL_SUBSTITUTIONS;
    int MAX_VASP_NSW;

    bool LAMMPS_HIDE_OUTPUT;
    bool VASP_HIDE_OUTPUT;
    bool REFORMAT_XDATCAR;
  };
  
}

void PRINT_SET(const std::set<int>& s) {
  for (auto e : s)
    std::cout << e << " ";
  std::cout << "\n";
}

int main(int argc, char *argv[]) {

  using namespace boost::units;

  lvx::init_check1();

  std::fstream f("lavax.conf");
  
  auto conf_data = lvx::parse_config_file(f);

  lvx::config_data conf;
  conf.VASP_COMMAND                 = conf_data["VASP_COMMAND"];
  conf.LAMMPS_COMMAND               = conf_data["LAMMPS_COMMAND"];
  conf.INIT_POSCAR                  = conf_data["INIT_POSCAR"];
  conf.LAVAX_ITERATIONS             = std::stoi(conf_data["LAVAX_ITERATIONS"]);
  conf.LAMMPS_POTENTIAL_FILE        = conf_data["LAMMPS_POTENTIAL_FILE"];
  conf.LAMMPS_PAIR_STYLE            = conf_data["LAMMPS_PAIR_STYLE"];
  conf.POTENTIAL_DEPARTURE_DISTANCE = std::stod(conf_data["POTENTIAL_DEPARTURE_DISTANCE"]) * angstrom;
  std::istringstream is(conf_data["USE_ADAPTIVE_POTIM"]);
  is >> std::boolalpha >> conf.USE_ADAPTIVE_POTIM;
  conf.MAX_DISTANCE                 = std::stod(conf_data["MAX_DISTANCE"]) * angstrom;
  conf.MAX_POTIM                    = std::stod(conf_data["MAX_POTIM"]) * femtosecond;
  conf.MAX_POTENTIAL_SUBSTITUTIONS  = std::stoi(conf_data["MAX_POTENTIAL_SUBSTITUTIONS"]);
  conf.MAX_VASP_NSW                 = std::stoi(conf_data["MAX_VASP_NSW"]);

  is.str(conf_data["LAMMPS_HIDE_OUTPUT"]);
  is >> std::boolalpha >> conf.LAMMPS_HIDE_OUTPUT;

  is.str(conf_data["VASP_HIDE_OUTPUT"]);
  is >> std::boolalpha >> conf.VASP_HIDE_OUTPUT;

  is.str(conf_data["REFORMAT_XDATCAR"]);
  is >> std::boolalpha >> conf.REFORMAT_XDATCAR;
  
  if (!lvx::init_check2(conf.INIT_POSCAR)) {
    return false;
  }
    
  int NSW;
  quantity<femtosecond_unit> POTIM;

  if (!lvx::parse_INCAR(NSW, POTIM)) {
    std::cerr << "Failed to parse INCAR!" << "\n";
    std::cerr << "... must contain fields NSW and POTIM" << std::endl;
  }
  std::cout << "NSW = " << NSW << ", POTIM = " << POTIM << std::endl;

  std::string lammps_atomic_symbols;
  for (int i = 0; conf_data.find(std::string("LAMMPS_ATOMIC_SYMBOL_") + std::to_string(i))
         != conf_data.end(); ++i) {
    lammps_atomic_symbols += conf_data[std::string("LAMMPS_ATOMIC_SYMBOL_") + std::to_string(i)] + " ";
  }

  lvx::update_LAMMPS_script(conf.LAMMPS_POTENTIAL_FILE,
                            lammps_atomic_symbols,
                            conf_data["LAMMPS_PAIR_STYLE"],
                            conf.MAX_VASP_NSW, POTIM);

  auto atomic_catalog = lvx::create_atomic_catalog(conf_data);
  
  quantity<angstrom_unit> latt_const;
  lvx::vec3_dimless       a1, a2, a3;
  std::ifstream           init_poscar(conf_data["INIT_POSCAR"]);
  lvx::simulation_cell    sim_cell;

  sim_cell.elements_info = atomic_catalog;
  
  lvx::parse_poscar(init_poscar, latt_const, a1, a2, a3, sim_cell);

  // using namespace lvx;
  std::cout << latt_const << "\n";
  std::cout << a1 << "\n";
  std::cout << a2 << "\n";
  std::cout << a3 << std::endl;

  bool deleted_prev = false;
  
  for (int i = 0; i < conf.LAVAX_ITERATIONS; ++i) {
    if (conf.USE_ADAPTIVE_POTIM) {
      auto itr =
      std::max_element(sim_cell.particles.begin(),
                       sim_cell.particles.end(),
                       [] (const lvx::atomic_particle& a,
                           const lvx::atomic_particle& b) {
                         return lvx::norm(a.getVel()) < lvx::norm(b.getVel());
                       });
      lvx::vec3_velocity max_vel = itr->getVel();
      auto speed = lvx::norm(max_vel);

      auto dt = std::min(conf.MAX_DISTANCE / speed, conf.MAX_POTIM);

      std::fstream f("INCAR");
      std::stringstream ss =
        lvx::replace_all_in_file(f, std::regex(".*POTIM.*"),
                                     std::string("POTIM = ") + std::to_string(dt.value()));
      f.close();
      f.open("INCAR", std::fstream::in | std::fstream::out | std::fstream::trunc);
      f << ss.str(); f.close();
      f.open("predictor.in");
      // ss = lvx::replace_all_in_file(f,  std::regex("^timestep.*"),
      //                                   std::string("timestep ") + std::to_string(dt.value()/5));
      // ss = lvx::replace_all_in_file(ss, std::regex(".*run.*"),
      //                                    std::string("run ") + std::to_string((NSW+3)*5));
      ss = lvx::replace_all_in_file(f,  std::regex("^timestep.*"),
                                        std::string("timestep ") + std::to_string(dt.value()));
      ss = lvx::replace_all_in_file(ss, std::regex(".*run.*"),
                                         std::string("run ") + std::to_string(conf.MAX_VASP_NSW));
      f.close();
      f.open("predictor.in", std::fstream::in | std::fstream::out | std::fstream::trunc);
      f << ss.str(); f.close();
      std::cout << "POTIM = " << dt << "\n";
    }
    std::ofstream writer("W_crystal.dat");
    std::string content = lvx::make_lammps_data(latt_const*a1[0],
                                                latt_const*a2[1],
                                                latt_const*a3[2],
                                                sim_cell);
    writer << content;
    writer.close();

    std::set<int> si = lvx::predict_neighbors(conf_data["LAMMPS_COMMAND"],
                                              conf.POTENTIAL_DEPARTURE_DISTANCE.value(),
                                              conf.LAMMPS_HIDE_OUTPUT,
                                              conf.MAX_POTENTIAL_SUBSTITUTIONS,
                                              conf.MAX_VASP_NSW,
                                              NSW);
    std::cout << "LAMMPS prediction DONE" << std::endl;
    std::cout << "Neighbor index: ";
    PRINT_SET(si);
    std::cout << "NSW = " << NSW << std::endl;

    std::set<int> diff;
    std::set_difference(si.begin(), si.end(),
                        sim_cell.lammps_selected_idx.begin(), sim_cell.lammps_selected_idx.end(),
                        std::inserter(diff, diff.end()));
    sim_cell.lammps_selected_idx = si;

    std::fstream f("INCAR");
    std::stringstream ss;
    if (diff.size() < conf.MAX_POTENTIAL_SUBSTITUTIONS) {

      ss = lvx::replace_all_in_file(f, std::regex(".*NSW.*"),
                                    std::string("NSW = ") + std::to_string(NSW));
      
      for (int j = 0; j < sim_cell.particles.size(); ++j)
        sim_cell.particles[j].is_hard = false;
    
      for (auto j : si)
        sim_cell.particles[j].is_hard = true;

      deleted_prev = false;
    } else {
      ss = lvx::replace_all_in_file(f, std::regex(".*NSW.*"),
                                    std::string("NSW = ") + std::to_string(conf.MAX_VASP_NSW));
      
      for (int j = 0; j < sim_cell.particles.size(); ++j)
        sim_cell.particles[j].is_hard = true;

      if (!deleted_prev) {
        namespace bf =  boost::filesystem;

        bf::path p("WAVECAR");

        if (bf::remove(p))
          std::cout << "REMOVED WAVECAR"  << "\n";
        else
          std::cout << "FAILED TO REMOVE WAVECAR" << "\n";

        deleted_prev = true;
      }
    }
    
    f.close();
    f.open("INCAR", std::fstream::in | std::fstream::out | std::fstream::trunc);
    f << ss.str(); f.close();

    writer.open("POSCAR");
    content = lvx::make_poscar(latt_const, a1, a2, a3, sim_cell);
    writer << content;
    writer.close();

    lvx::concat_POTCAR(sim_cell);

    std::this_thread::sleep_for(std::chrono::seconds(2));

    std::string cmd = conf_data["VASP_COMMAND"] + (conf.VASP_HIDE_OUTPUT ? " > /dev/null" : "");
    system(cmd.c_str());
    std::cout << "VASP run DONE" << std::endl;

    // ----------------------------------------

    if (conf.REFORMAT_XDATCAR) {
      std::string str1 = "\\s*";
      std::string str2;
      std::string str3;
    
      for (const auto& e : sim_cell.vasp_symbol_count_helper) {
        if (std::get<1>(e) != 0) {
          str1 += std::get<0>(e) + "\\s*";
        }
        str2 += std::get<0>(e) + " ";
        str3 += std::to_string(std::get<1>(e)) + " ";
      }
    
      f.open("XDATCAR");
      ss.str("");
    
      std::regex  rgx(str1);
      std::string line;
      while (std::getline(f, line)) {
        if (std::regex_match(line, rgx)) {
          ss << std::regex_replace(line, rgx, str2) << "\n";
          std::getline(f, line);
          ss << str3 << "\n";
        } else
          ss << line << "\n";
      }

      f.close();
      f.open("XDATCAR", std::fstream::in | std::fstream::out | std::fstream::trunc);
      f << ss.str();
      f.close();
    }

    std::ifstream contcar("CONTCAR");
    lvx::parse_poscar(contcar, latt_const, a1, a2, a3, sim_cell);
    contcar.close();

    std::vector<std::string> files = {"XDATCAR", "CONTCAR", "CHG",
                                      "CHGCAR", "DOSCAR", "EIGENVAL",
                                      "OSZICAR", "PCDAT", "vasprun.xml",
                                      "OUTCAR", "INCAR", "IBZKPT", "POSCAR"};
    lvx::backup_files(files, i, conf.LAVAX_ITERATIONS);
    std::cout << "--------------------" << std::endl;
    std::this_thread::sleep_for(std::chrono::seconds(2));
  }

  return 0;
}
