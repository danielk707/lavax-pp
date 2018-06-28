#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <regex>
#include <thread>
#include <cstdlib>
#include <complex>
#include <map>
#include "lavax.hpp"

namespace lvx {

  std::string make_poscar(const quantity<angstrom_unit>& lattice_constant,
                          const vec3_dimless& a1,
                          const vec3_dimless& a2,
                          const vec3_dimless& a3,
                          const simulation_cell_v3& sim_cell) {
    std::stringstream ss;
    // ss << std::fixed;
    
    ss << "BCC Xx \n";
    ss << lattice_constant.value() << "\n";
    ss << a1 << "\n";
    ss << a2 << "\n";
    ss << a3 << "\n";

    std::vector<std::tuple<std::string,int,bool> > symbols_count;
    for (auto e : sim_cell.elements_info) {
      int num_bad = 
      std::count_if(sim_cell.particles.cbegin(),
                    sim_cell.particles.cend(),
                    [&e] (const atomic_particle& a) {
                      return (a.element_info->vasp_symbol_bad == e->vasp_symbol_bad) &&
                        !a.high_prec;
                    });

      int num_good = 
      std::count_if(sim_cell.particles.cbegin(),
                    sim_cell.particles.cend(),
                    [&e] (const atomic_particle& a) {
                      return (a.element_info->vasp_symbol_good == e->vasp_symbol_good) &&
                        a.high_prec;
                    });
      symbols_count.push_back(std::make_tuple(e->vasp_symbol_bad,  num_bad,  false));
      symbols_count.push_back(std::make_tuple(e->vasp_symbol_good, num_good, true));      
    }

    for (const auto& sc : symbols_count) {
      if (std::get<1>(sc) != 0)
        ss << std::get<0>(sc) << " ";
    }

    ss << "\n";

    for (const auto& sc : symbols_count) {
      if (std::get<1>(sc) != 0)
        ss << std::get<1>(sc) << " ";
    }

    ss << "\n";
    ss << "Cartesian\n";

    for (const auto& sc : symbols_count) {
      for (const auto& p : sim_cell.particles) {
        if        ( std::get<2>(sc) &&  p.high_prec && p.element_info->vasp_symbol_good == std::get<0>(sc)) {
          ss << p.getPos() << "\n";
        } else if (!std::get<2>(sc) && !p.high_prec && p.element_info->vasp_symbol_bad  == std::get<0>(sc)) {
          ss << p.getPos() << "\n";
        }
      }
    }

    ss << "\n";
    for (const auto& sc : symbols_count) {
      for (const auto& p : sim_cell.particles) {
        if        ( std::get<2>(sc) &&  p.high_prec && p.element_info->vasp_symbol_good == std::get<0>(sc)) {
          ss << p.getVel() << "\n";
        } else if (!std::get<2>(sc) && !p.high_prec && p.element_info->vasp_symbol_bad  == std::get<0>(sc)) {
          ss << p.getVel() << "\n";
        }
      }
    }
    sim_cell.vasp_symbol_count_helper.assign(symbols_count.begin(),
                                             symbols_count.end());
    return ss.str();
  }

  std::vector<Particle_v2> create_crystal(crystal_structure cstruct,
                                          quantity<angstrom_unit> lattice_constant,
                                          int I, int J, int K) {
    auto& L = lattice_constant;
    
    vec3_velocity zero_vel;
    zero_vel = 0.0 * angstrom_per_fs,
               0.0 * angstrom_per_fs,
               0.0 * angstrom_per_fs;

    std::vector<Particle_v2> P;
    
    if (cstruct == crystal_structure::BCC) {
      P.reserve(2*I*J*K);
      for (int i = 0; i < I; ++i) {
        for (int j = 0; j < J; ++j) {
          for (int k = 0; k < K; ++k) {
            P.emplace_back(vec3_angstrom(L*(i+0.0), L*(j+0.0), L*(k+0.0)), zero_vel);
            P.emplace_back(vec3_angstrom(L*(i+0.5), L*(j+0.5), L*(k+0.5)), zero_vel);
          }
        }
      }
      return P;
    }
    return P;
  }

  vec3_angstrom transform(const quantity<angstrom_unit>& lattice_constant,
                          const vec3_dimless& a1,
                          const vec3_dimless& a2,
                          const vec3_dimless& a3,
                          const vec3_dimless& v) {
    
    blitz::TinyMatrix<quantity<si::dimensionless>, 3, 3> A;
    blitz::TinyVector<quantity<angstrom_unit>,1> L;
    L = lattice_constant;
    
    A = a1(0), a2(0), a3(0),
        a1(1), a2(1), a3(1),
        a1(2), a2(2), a3(2);

    // std::cout << A << "\n";
    
    // using namespace blitz::tensor;
    blitz::firstIndex i;
    blitz::secondIndex j;
    auto temp = blitz::sum(A(i,j)*v(j), j);

    vec3_angstrom rtn;
    rtn[0] = temp(0) * lattice_constant;
    rtn[1] = temp(1) * lattice_constant;
    rtn[2] = temp(2) * lattice_constant;

    return rtn;
  }

  std::set<int> parse_lammps_neighbor(std::ifstream& file,
                                      quantity<angstrom_unit> cutoff) {
    
    std::regex reg(".*ITEM\\: ENTRIES c_distance\\[1\\] c_neigh\\[1\\] c_neigh\\[2\\].*");
    std::set<int> indicies;
    std::string line;

    while (std::getline(file, line)) {
      if (std::regex_match(line, reg)) {
        double r;
        file >> r;

        while (file.good()) {
          if (r <= cutoff.value()) {
            int i;
            file >> i;
            if (file.good()) indicies.insert(i);
            file >> i;
            if (file.good()) indicies.insert(i);
          } else
            break;
          
          file >> r;
        }
        file.clear();
      }
    }
    return indicies;
  }

  std::set<int> parse_lammps_neighbor(std::ifstream& file,
                                      quantity<angstrom_unit> cutoff,
                                      int max_potential_switch,
                                      int max_vasp_nsw,
                                      int count_good_prev,
                                      int& NSW) {
    std::set<int> indices;
    // int iteration;
    std::string line;

    while (std::getline(file, line)) {
      if (std::regex_match(line, std::regex(".*ITEM\\: ENTRIES c_distance\\[1\\] "
                                            "c_neigh\\[1\\] c_neigh\\[2\\].*"))) {
        // std::cout << line << "\n";
        while (std::getline(file, line)) {
          std::smatch matches;

          // Regex for floating point number:
          std::regex reg3("([+-]?(?=[.]?[0-9])[0-9]*(?:[.][0-9]*)?(?:[Ee][+-]?[0-9]+)?)\\s+"
                          "(\\d+)\\s+(\\d+)");

          if (std::regex_search(line, matches, reg3)) {
            if (std::stod(matches[1]) <= cutoff.value()) {
              indices.insert(std::stoi(matches[2]));
              indices.insert(std::stoi(matches[3]));
              if (abs(indices.size()-count_good_prev) >= max_potential_switch)
                break;
            }
          } else
            break;
        }
      }
      if (std::regex_match(line, std::regex("ITEM\\: TIMESTEP.*"))) {
        std::getline(file, line);
        std::smatch matches;

        if (std::regex_search(line, matches, std::regex("(\\d+)"))) {
          NSW = std::stoi(matches[1]);
          if (NSW >= max_vasp_nsw)
            break;
          // std::cout << NSW << "\n";
        }
      }
      if (abs(indices.size()-count_good_prev) >= max_potential_switch)
        break;
    }
    return indices;
  }

  void _parse_poscar_header(std::istream& poscar,
                            quantity<angstrom_unit>& lattice_constant,
                            vec3_dimless& a1,
                            vec3_dimless& a2,
                            vec3_dimless& a3) {
    poscar.ignore(1024, '\n'); // Ignore first line
    double L;
    poscar >> L;
    
    lattice_constant = L * angstrom;
    a1 = read_vec3(poscar, dimensionless);
    a2 = read_vec3(poscar, dimensionless);
    a3 = read_vec3(poscar, dimensionless);
  }

  void parse_init_poscar(std::ifstream& poscar,
                         quantity<angstrom_unit>& lattice_constant,
                         vec3_dimless& a1,
                         vec3_dimless& a2,
                         vec3_dimless& a3,
                         simulation_cell_v3& sim_cell) {
    
    _parse_poscar_header(poscar, lattice_constant, a1, a2, a3);

    // std::smatch matches;
    std::map<std::string,std::pair<int, bool>> dict;

    for (int i = 0; i < sim_cell.elements_info.size(); ++i) {
      dict[sim_cell.elements_info[i]->vasp_symbol_good] = std::make_pair(i,false);
      dict[sim_cell.elements_info[i]->vasp_symbol_bad]  = std::make_pair(i,true);
    }
    // std::cout << "GOT HERE" << "\n";
    std::vector<std::string> vasp_symbols;
    
    std::string line;
    while (getline(poscar, line)) {
      if (std::regex_match(line, std::regex("\\s*\\d+[\\s\\d]*"))) {
        break;
      }
      // This can be simplified if you think about it:
      if (std::regex_match(line, std::regex("\\s*\\w+.*"))) {
        std::regex rgx("\\w+");
        for( std::sregex_iterator it(line.begin(), line.end(), rgx), it_end;
             it != it_end; ++it ) {
          vasp_symbols.push_back((*it)[0]);
          std::cout << "DEBUGGG " << (*it)[0] << "\n";
        }
        // continue;
      }
    }

    // Parse line containing number of atoms of each potential:
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
    
    // std::cout << vasp_symbols.size() << "\n";
    // std::cout << num_atom_species.size() << "\n";
    
    // int j = 0;
    // for (auto s : vasp_symbols) {
    //   if (s == sim_cell.elements[dict[s]].vasp_symbol_bad) {
    //     sim_cell.elements[dict[s]].vasp_num_bad = num_atom_species[j];
    //   }
    //   else if (s == sim_cell.elements[dict[s]].vasp_symbol_good) {
    //     sim_cell.elements[dict[s]].vasp_num_good = num_atom_species[j];
    //   }
    //   ++j;
    // }

    std::vector<int> indices;
    std::partial_sum(num_atom_species.begin(), num_atom_species.end(),
                     std::back_inserter(indices));

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
    int k = 0;
    for (int i = 0; i < pos_vec.size(); i++) {
      atomic_particle p(pos_vec[i], vel_vec[i]);
      if (i >= indices[k]) {
        ++k;
      }
        
      p.element_info = sim_cell.elements_info[dict[vasp_symbols[k]].first];
      p.high_prec    = dict[vasp_symbols[k]].second;
      sim_cell.particles.push_back(p);
      // sim_cell.particles.emplace_back(pos_vec[i], vel_vec[i]);
    }
  }
  
  std::string make_lammps_data(quantity<angstrom_unit> xhi,
                               quantity<angstrom_unit> yhi,
                               quantity<angstrom_unit> zhi,
                               simulation_cell_v3& sim_cell,
                               bool include_velocity) {
    std::stringstream ss;
    ss.precision(8);
    ss << std::fixed;

    ss << "# W Crystal bcc #\n\n";

    ss << sim_cell.particles.size()     << " atoms\n";
    ss << sim_cell.elements_info.size() << " atom types\n\n";

    ss << 0.0 << " " << xhi.value() << " xlo xhi\n";
    ss << 0.0 << " " << yhi.value() << " ylo yhi\n";
    ss << 0.0 << " " << zhi.value() << " zlo zhi\n";

    ss << "\nMasses\n\n";
    for (int i = 0; i < sim_cell.elements_info.size(); ++i) {
      ss << sim_cell.elements_info[i]->atom_type << " "
         << sim_cell.elements_info[i]->mass.value() << "\n";
    }

    ss << "\nAtoms\n\n";
    for (int i = 0; i < sim_cell.particles.size(); ++i) {
      ss << i+1 << " " << sim_cell.particles[i].element_info->atom_type
         << " " << sim_cell.particles[i].getPos() << "\n";
    }

    if (include_velocity) {
      ss << "\nVelocities\n\n";
      for (int i = 0; i < sim_cell.particles.size(); ++i) {
        ss << i+1 << " " << sim_cell.particles[i].getVel() << "\n";
      }
    }
    
    return ss.str();
  }

    // vec3_velocity PKAvel(const quantity<atomic_mass_unit>& mass,
    //                      const vec3_dimless& dir,
    //                      const quantity<electron_volt_unit>& energy) {

    // }
}
