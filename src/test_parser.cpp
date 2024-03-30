#include <string_view>
#include <map>
#include <vector>
#include <toml.hpp>

using namespace std;

struct Material {
  double F;
  double D; // diffusion coefficient
  vector<double> xs_a; // absorbsion XS
  vector<double> xs_f; // fission XS
  vector<double> chi;
};

struct Config {
  map<string, Material> materials;
};

Config parse(string_view toml) {
  Config config;

  return config;
}

int parse_input_1() {
  static constexpr std::string_view some_toml = R"(
    [geometry]
    mesh = "cube.msh"
    boundaries = { S1 = "dirichlet" }

    [cross-sections]
    MAT1 = { F = 0, D = 1.5, xs_absorption = 0.2, xs_fission = 0.3, chi = 1.0 }

    [params]
    calculation = "k_eff"
    energy-groups = 1
  )";

  return 0;
}

int main() {
  parse_input_1();
  return 0;
}
