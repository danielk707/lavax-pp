#include "parse_config_file.hpp"
#include <iostream>
#include <fstream>

int main(int argc, char *argv[])
{
  using namespace lvx;
  std::cout << "parse_config_file" << "\n";
  // std::string config = "url = http://example.com\n"
  //   "file = main.exe   # comment\n"
  //   "true = 0";
  std::fstream f("lavax.conf");
  
  // std::istringstream is_file(config);

  auto dict = parse_config_file(f);

  for (auto v : dict) {
    std::cout << "KEY: " << v.first << "\t\t\t" << v.second << "\n";
  }
  
  return 0;
}
