#include <string_view>
#include <assert.h>
#include <iostream>
#include <toml.hpp>
#include <fermi.hpp>
#include <utils.hpp>
#include <parser.hpp>

using namespace std;

int test_parse_input_mesh_file() {
  static constexpr std::string_view some_toml = R"(
    [geometry]
    mesh = "cube.msh"
    boundaries = { S1 = "dirichlet", WALL8 = "neumann" }

    [materials]
    MAT1 = { D = [1.5], xs_a = [0.2], xs_f = [0.3], xs_s = [1.0], chi = [1.0] }
  )";

  optional<Config> result = Config::parse(some_toml);
  if (!result) assert(false);

  auto config = result.value();

  assert(config.mesh_file == "cube.msh");
  return 0;
}

int test_parse_input_boundaries() {
  static constexpr std::string_view some_toml = R"(
    [geometry]
    mesh = "cube.msh"
    boundaries = { S1 = "dirichlet", WALL8 = "neumann" }

    [materials]
    MAT1 = { D = [1.5], xs_a = [0.2], xs_f = [0.3], xs_s = [1.0], chi = [1.0] }
  )";

  optional<Config> result = Config::parse(some_toml);
  if (!result) assert(false);

  Config config = result.value();

  map<string, BoundaryCondition> expected = {{"S1", Dirichlet}, {"WALL8", Neumann}};

  assert(mapsAreEqual(config.boundaries, expected));
  return 0;
}

int test_parse_cross_sections() {
  static constexpr std::string_view some_toml = R"(
    [geometry]
    mesh = "cube.msh"
    boundaries = { S1 = "dirichlet", WALL8 = "neumann" }

    [materials]
    MAT1 = { D = [1.5], xs_a = [0.2], xs_f = [0.3], xs_s = [1.0], chi = [1.0] }
    MAT2 = { D = [3.1], xs_a = [0.3], xs_f = [1.3], xs_s = [1.2], chi = [1.3] }
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

int test_parse_bad_inputs() {
  { // miss geometry
    static constexpr std::string_view some_toml = R"(
      [materials]
      MAT1 = { D = [1.5], xs_a = [0.2], xs_f = [0.3], xs_s = [1.0], chi = [1.0] }
    )";

    optional<Config> result = Config::parse(some_toml);
    assert(!result.has_value());
  }

  { // miss materials
    static constexpr std::string_view some_toml = R"(
      [geometry]
      mesh = "cube.msh"
      boundaries = { S1 = "dirichlet", WALL8 = "neumann" }
    )";

    optional<Config> result = Config::parse(some_toml);
    assert(!result.has_value());
  }

  { // empty mesh file path
    static constexpr std::string_view some_toml = R"(
      [geometry]
      mesh = ""
      boundaries = { S1 = "dirichlet", WALL8 = "neumann" }

      [materials]
      MAT1 = { D = [1.5], xs_a = [0.2], xs_f = [0.3], xs_s = [1.0], chi = [1.0] }
    )";

    optional<Config> result = Config::parse(some_toml);
    assert(!result.has_value());
  }

  { // BC conditions are not dirichlet | neumann
    static constexpr std::string_view some_toml = R"(
      [geometry]
      mesh = "cube.msh"
      boundaries = { S1 = "dirichlet", WALL8 = "neumann_1" }

      [materials]
      MAT1 = { D = [1.5], xs_a = [0.2], xs_f = [0.3], xs_s = [1.0], chi = [1.0] }
    )";

    optional<Config> result = Config::parse(some_toml);
    assert(!result.has_value());
  }

  { // The materials are missing xs_f
    static constexpr std::string_view some_toml = R"(
      [geometry]
      mesh = "cube.msh"
      boundaries = { S1 = "dirichlet", WALL8 = "neumann" }

      [materials]
      MAT1 = { D = [1.5], xs_a = [0.2], xs_s = [1.0], chi = [1.0] }
    )";

    optional<Config> result = Config::parse(some_toml);
    assert(!result.has_value());
  }

  { // The materials are having a string instead of a double for xs_f
    static constexpr std::string_view some_toml = R"(
      [geometry]
      mesh = "cube.msh"
      boundaries = { S1 = "dirichlet", WALL8 = "neumann" }

      [materials]
      MAT1 = { D = [1.5], xs_a = [0.2], xs_f = ["a"], xs_s = [1.0], chi = [1.0] }
    )";

    optional<Config> result = Config::parse(some_toml);
    assert(!result.has_value());
  }

  return 0;
}

int main() {
  test_parse_input_mesh_file();
  test_parse_input_boundaries();
  test_parse_cross_sections();
  test_parse_bad_inputs();
  return 0;
}
