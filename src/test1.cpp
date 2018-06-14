#include <iostream>
#include <regex>
#include <vector>
#include <string>
#include <algorithm>

int main(int argc, char *argv[])
{
  std::string str(" 12 3 1\n");

  std::regex rgx("(\\d+)");

  std::regex_iterator<std::string::iterator> rit (str.begin(), str.end(), rgx);
  std::regex_iterator<std::string::iterator> rend;

  while (rit!=rend) {
    std::cout << rit->str() << std::endl;
    ++rit;
  }

  // std::vector<int> w;
  // std::for_each(std::sregex_token_iterator(str.begin(), str.end(), rgx, -1),
  //               std::sregex_token_iterator(), [] ( )

  std::for_each(std::sregex_iterator(str.begin(), str.end(), rgx),
                std::sregex_iterator(),
                [] (const std::smatch& m) {
                  std::cout << std::stoi(m.str(1)) << "\n";
                });

    
    
  // std::smatch matches;

  // std::regex_search(str, matches, rgx);

  // for (auto v : matches) {
  //   std::cout << v << "\n";
  // }
  
  return 0;
}
