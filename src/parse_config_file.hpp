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
      
      std::regex rgx("^\\s*([^\\s#]+)\\s*=\\s*([^\\s#][^#]*\\b)\\s*#?.*");
      std::smatch matches;
      
      if (std::regex_search(line, matches, rgx)) {
        if (matches.size() == 3) {
          dict[matches[1]] = matches[2];
        }
      }
    }
    return dict;
  }
}
