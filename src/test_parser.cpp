
#include <string>
#include <string_view>
#include <map>
#include <vector>
#include <assert.h>
#include <iostream>
#include <optional>
#include <toml.hpp>

using namespace std;

template<typename K, typename V>
bool mapsAreEqual(const std::map<K, V>& map1, const std::map<K, V>& map2) {
	return map1.size() == map2.size() && std::equal(map1.begin(), map1.end(), map2.begin());
}

template<class K, class V> void printMap(map<K, V> m) {
    typename map<K, V>::iterator it = m.begin();
 
    while (it != m.end()) {
        cout << "Key: " << it->first
             << ", Value: " << it->second << endl;
        ++it;
    }
}

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

  auto m = tbl["geometry"]["boundaries"];
  
  if (m.is_table()) {
    for (auto [k, v] : *m.as_table()) {
      string key = static_cast<string>(k);
      string condition = *m[k].value<string>();
      if (condition == "dirichlet") {
        config.boundaries[key] = Dirichlet;
      }else if (condition == "neumann") {
        config.boundaries[key] = Neumann;
      }
    }
  }

  return config;
}

int parse_input_mesh_file() {
  static constexpr std::string_view some_toml = R"(
    [geometry]
    mesh = "cube.msh"
  )";

  Config config = Config::parse(some_toml);

  assert(config.mesh_file == "cube.msh");

  static constexpr std::string_view some_toml_wo_mesh_file = R"(
    [geometry]
  )";

  Config config_wo_mesh_file = Config::parse(some_toml_wo_mesh_file);

  assert(config_wo_mesh_file.mesh_file == "");
  return 0;
}

int parse_input_boundaries() {
  static constexpr std::string_view some_toml = R"(
    [geometry]
    boundaries = { S1 = "dirichlet", WALL8 = "neumann" }
  )";

  Config config = Config::parse(some_toml);

  map<string, BoundaryCondition> expected = {{"S1", Dirichlet}, {"WALL8", Neumann}};

  assert(mapsAreEqual(config.boundaries, expected));
  return 0;
}

  // static constexpr std::string_view some_toml = R"(
  //   [geometry]
  //   mesh = "cube.msh"
  //   boundaries = { S1 = "dirichlet" }
  //
  //   [cross-sections]
  //   MAT1 = { F = 0, D = 1.5, xs_absorption = 0.2, xs_fission = 0.3, chi = 1.0 }
  //
  //   [params]
  //   calculation = "k_eff"
  //   energy-groups = 1
  // )";

int main() {
  parse_input_mesh_file();
  parse_input_boundaries();
  return 0;
}
