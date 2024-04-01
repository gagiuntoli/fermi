
#include <string>
#include <string_view>
#include <map>
#include <vector>
#include <array>
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
        cout << "Key: " << it->first << ", Value: " << it->second << endl;
        ++it;
    }
}

struct Material {
  vector<double> D;    // diffusion coefficient
  vector<double> xs_a; // absorbsion XS
  vector<double> xs_f; // fission XS
  vector<double> xs_s; // scattering XS
  vector<double> chi;  // factor determing with % of all fissions ends in the groups
  //
  bool operator==(const Material& other) const {
        return D == other.D && xs_a == other.xs_a && xs_f == other.xs_f && xs_s == other.xs_s && chi == other.chi;
  }
};

enum Calculation {
  Keff
};

enum BoundaryCondition {
  Dirichlet,
  Neumann,
};

struct Config {
  string mesh_file;
  map<string, BoundaryCondition> boundaries;
  map<string, Material> materials;
  Calculation calculation;
  uint groups;

  static optional<Config> parse(string_view toml);
};

optional<Config> Config::parse(string_view toml_string) {
  Config config;

  toml::table tbl = toml::parse(toml_string);

  auto mesh_file = tbl["geometry"]["mesh"];
  if (!mesh_file.is_value()) {
    cerr << "Input error: No mesh file path found" << endl;
    return {};
  }
  config.mesh_file = mesh_file.value<string>().value();

  auto boundaries = tbl["geometry"]["boundaries"];
  if (!boundaries.is_table()) {
    cerr << "Input error: No boundaries table" << endl;
    return {};
  }
  
  for (auto [k, v] : *boundaries.as_table()) {
    string key = static_cast<string>(k);
    if (key.empty()) {
        cerr << "Input error: Empty key in boundaries" << endl;
        return {};
    }
    string condition = *boundaries[k].value<string>();
    if (condition == "dirichlet") {
      config.boundaries[key] = Dirichlet;
    }else if (condition == "neumann") {
      config.boundaries[key] = Neumann;
    } else {
      return {};
    }
  }

  auto materials = tbl["materials"];
  for (auto [k, v] : *materials.as_table()) {
    string key = static_cast<string>(k);
    if (key.empty()) {
        cerr << "Input error: Empty key in boundaries" << endl;
        return {};
    }
    auto material = *materials[k].as_table();

    for (auto&& elem : *material["D"].as_array()) {
      elem.visit([&config, key](auto&& el) noexcept {
        if constexpr (toml::is_number<decltype(el)>)
          config.materials[key].D.push_back(*el);
      });
    }

    for (auto&& elem : *material["chi"].as_array()) {
      elem.visit([&config, key](auto&& el) noexcept {
        if constexpr (toml::is_number<decltype(el)>)
          config.materials[key].chi.push_back(*el);
      });
    }

    for (auto&& elem : *material["xs_a"].as_array()) {
      elem.visit([&config, key](auto&& el) noexcept {
        if constexpr (toml::is_number<decltype(el)>)
          config.materials[key].xs_a.push_back(*el);
      });
    }

    for (auto&& elem : *material["xs_f"].as_array()) {
      elem.visit([&config, key](auto&& el) noexcept {
        if constexpr (toml::is_number<decltype(el)>)
          config.materials[key].xs_f.push_back(*el);
      });
    }

    for (auto&& elem : *material["xs_s"].as_array()) {
      elem.visit([&config, key](auto&& el) noexcept {
        if constexpr (toml::is_number<decltype(el)>)
          config.materials[key].xs_s.push_back(*el);
      });
    }
  }

  auto simulation_parameters = tbl["simulation-parameters"];
  if (!simulation_parameters.is_table()) {
    cerr << "Input error: No simulation-parameters table" << endl;
    return {};
  }

  auto calculation = simulation_parameters["calculation"];
  if (!calculation.is_value()) {
    cerr << "Input error: No calculation type specified: calculation = \"k_eff\"" << endl;
    return {};
  }

  if (calculation.value<string>().value() == "k_eff") {
    config.calculation = Keff;
  } else {
    cerr << "Input error: calculation type not recognized, values available: \"k_eff\"" << endl;
    return {};
  }

  auto groups = simulation_parameters["groups"];
  if (!groups.is_value()) {
    cerr << "Input error: No energy groups specified: groups = 1|2|..." << endl;
    return {};
  }
  config.groups = groups.value<uint>().value();

  return config;
}

int parse_input_mesh_file() {
  static constexpr std::string_view some_toml = R"(
    [geometry]
    mesh = "cube.msh"
    boundaries = { S1 = "dirichlet", WALL8 = "neumann" }

    [materials]
    MAT1 = { D = [1.5], xs_a = [0.2], xs_f = [0.3], xs_s = [1.0], chi = [1.0] }

    [simulation-parameters]
    calculation = "k_eff"
    groups = 1
  )";

  optional<Config> result = Config::parse(some_toml);
  if (!result) assert(false);

  auto config = result.value();

  assert(config.mesh_file == "cube.msh");

  static constexpr std::string_view some_toml_wo_mesh_file = R"(
    [geometry]
  )";

  result = Config::parse(some_toml_wo_mesh_file);
  assert(!result.has_value());

  return 0;
}

int parse_input_boundaries() {
  static constexpr std::string_view some_toml = R"(
    [geometry]
    mesh = "cube.msh"
    boundaries = { S1 = "dirichlet", WALL8 = "neumann" }

    [materials]
    MAT1 = { D = [1.5], xs_a = [0.2], xs_f = [0.3], xs_s = [1.0], chi = [1.0] }

    [simulation-parameters]
    calculation = "k_eff"
    groups = 1
  )";

  optional<Config> result = Config::parse(some_toml);
  if (!result) assert(false);

  Config config = result.value();

  map<string, BoundaryCondition> expected = {{"S1", Dirichlet}, {"WALL8", Neumann}};

  assert(mapsAreEqual(config.boundaries, expected));
  return 0;
}

int parse_cross_sections() {
  static constexpr std::string_view some_toml = R"(
    [geometry]
    mesh = "cube.msh"
    boundaries = { S1 = "dirichlet", WALL8 = "neumann" }

    [materials]
    MAT1 = { D = [1.5], xs_a = [0.2], xs_f = [0.3], xs_s = [1.0], chi = [1.0] }
    MAT2 = { D = [3.1], xs_a = [0.3], xs_f = [1.3], xs_s = [1.2], chi = [1.3] }

    [simulation-parameters]
    calculation = "k_eff"
    groups = 1
  )";

  optional<Config> result = Config::parse(some_toml);
  if (!result) assert(false);

  Config config = result.value();

  auto mat1 = Material {.D = {1.5}, .xs_a = {0.2}, .xs_f = {0.3}, .xs_s = {1.0}, .chi = {1.0}};
  auto mat2 = Material {.D = {3.1}, .xs_a = {0.3}, .xs_f = {1.3}, .xs_s = {1.2}, .chi = {1.3}};

  map<string, Material> expected = {{"MAT1", mat1}, {"MAT2", mat2}};
  assert(mapsAreEqual(config.materials, expected));

  return 0;
}

int parse_simulation_parameters() {
  static constexpr std::string_view some_toml = R"(
    [geometry]
    mesh = "cube.msh"
    boundaries = { S1 = "dirichlet", WALL8 = "neumann" }

    [materials]
    MAT1 = { D = [1.5], xs_a = [0.2], xs_f = [0.3], xs_s = [1.0], chi = [1.0] }
    MAT2 = { D = [3.1], xs_a = [0.3], xs_f = [1.3], xs_s = [1.2], chi = [1.3] }

    [simulation-parameters]
    calculation = "k_eff"
    groups = 1
  )";

  optional<Config> result = Config::parse(some_toml);
  if (!result) assert(false);

  Config config = result.value();
  assert(config.calculation == Keff);
  assert(config.groups == 1);

  return 0;
}

int main() {
  parse_input_mesh_file();
  parse_input_boundaries();
  parse_cross_sections();
  parse_simulation_parameters();
  return 0;
}
