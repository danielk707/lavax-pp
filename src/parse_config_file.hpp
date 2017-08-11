#pragma once
#include <istream>
#include <map>
#include <string>
#include <sstream>
#include <regex>
#include <iostream>

namespace lvx {
  // ------------------------------------------------------------
  std::map<std::string, std::string> parse_config_file(std::istream& conf_file) {
    std::map<std::string, std::string> dict;
    std::string line;
    
    while (std::getline(conf_file, line)) {
      // std::istringstream strm_line(line);
      // std::string key;

      std::regex rgx("\\s*(\\S*)\\s*=\\s*(\\S*)\\s*#?.*");
      std::smatch matches;
      if (std::regex_search(line, matches, rgx)) {
        if (matches.size() == 3) {
          dict[matches[1]] = matches[2];
        }
        // for (auto v : matches) {
        //   std::cout << "REG   " << v << "\n";
        // }
      }
      
      // if (std::getline(strm_line, key, '=')) {
      //   std::string value;
        
      //   if (std::getline(strm_line, value))
      //     dict[key] = value;
      // }
    }
    return dict;
  }
}
