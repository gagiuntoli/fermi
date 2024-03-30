
#include <string>
#include <string_view>
#include <map>
#include <vector>
#include <assert.h>
#include <iostream>
#include <optional>
#include <toml.hpp>

using namespace std;

struct Material {
  double F;
  double D; // diffusion coefficient
  vector<double> xs_a; // absorbsion XS
  vector<double> xs_f; // fission XS
  vector<double> chi;
};

enum BoundaryCondition {
  Dirichlet,
  Neumann,
};

struct Config {
  string mesh_file;
  map<string, BoundaryCondition> boundaries;
  map<string, Material> materials;

  static Config parse(string_view toml);
};

Config Config::parse(string_view toml_string) {
  Config config;

  toml::table tbl = toml::parse(toml_string);

  config.mesh_file = tbl["geometry"]["mesh"].value_or<std::string>("");

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

  Config config = Config::parse(some_toml);

  assert(config.mesh_file == "cube.msh");

  return 0;
}

int main() {
  parse_input_1();
  return 0;
}
