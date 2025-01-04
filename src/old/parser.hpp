#include <fermi.hpp>
#include <iostream>
#include <optional>
#include <string>
#include <vector>

using namespace std;

class Config {
   private:
    static vector<double> parse_array_from_table(toml::v3::table table,
                                                 string key);

   public:
    string mesh_file;
    map<string, BoundaryCondition> boundaries;
    map<string, Material> materials;
    uint groups = 0;

    static optional<Config> parse(string_view toml);
};

vector<double> Config::parse_array_from_table(toml::v3::table table,
                                              string key) {
    auto arr = table[key];
    if (!arr.is_array()) {
        cerr << "Input error: " << key << " not found or it is not an array"
             << endl;
        return {};
    }

    vector<double> result;
    for (auto &&elem : *arr.as_array()) {
        elem.visit([&result, key](auto &&el) noexcept {
            if constexpr (toml::is_number<decltype(el)>)
                result.push_back(*el);
        });
    }
    return result;
}

optional<Config> Config::parse(string_view toml_string) {
    Config config;
    toml::table tbl = toml::parse(toml_string);

    auto mesh_file = tbl["geometry"]["mesh"];
    if (!mesh_file.is_value()) {
        cerr << "Input error: No mesh file path found" << endl;
        return {};
    }

    string mesh_file_str = mesh_file.value<string>().value();
    if (mesh_file_str.empty()) {
        cerr << "Input error: mesh file path is empty" << endl;
        return {};
    }

    config.mesh_file = mesh_file_str;

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
        } else if (condition == "neumann") {
            config.boundaries[key] = Neumann;
        } else {
            return {};
        }
    }

    auto materials = tbl["materials"];
    if (!materials.is_table()) {
        cerr << "Input error: No materials table" << endl;
        return {};
    }

    for (auto [k, v] : *materials.as_table()) {
        string key = static_cast<string>(k);
        if (key.empty()) {
            cerr << "Input error: Empty key in boundaries" << endl;
            return {};
        }
        auto material = *materials[k].as_table();

        vector<double> vec_double;
        vec_double = Config::parse_array_from_table(material, "D");

        // the energy groups correspond to the first D length
        if (0 < vec_double.size()) {
            if (config.groups == 0) {
                config.groups = vec_double.size();
            }
        } else {
            cerr << "Input error: D needs at least 1 element" << endl;
            return {};
        }

        copy(vec_double.begin(), vec_double.end(),
             back_inserter(config.materials[key].D));

        vec_double = Config::parse_array_from_table(material, "chi");
        copy(vec_double.begin(), vec_double.end(),
             back_inserter(config.materials[key].chi));

        vec_double = Config::parse_array_from_table(material, "xs_a");
        copy(vec_double.begin(), vec_double.end(),
             back_inserter(config.materials[key].xs_a));

        vec_double = Config::parse_array_from_table(material, "xs_f");
        copy(vec_double.begin(), vec_double.end(),
             back_inserter(config.materials[key].xs_f));

        vec_double = Config::parse_array_from_table(material, "xs_s");
        copy(vec_double.begin(), vec_double.end(),
             back_inserter(config.materials[key].xs_s));

        if (config.materials[key].D.size() != config.groups ||
            config.materials[key].chi.size() != config.groups ||
            config.materials[key].xs_a.size() != config.groups ||
            config.materials[key].xs_f.size() != config.groups ||
            config.materials[key].xs_s.size() != config.groups * config.groups) {
            cerr << "Input error: lenght in material parameters is wrong" << endl;
            return {};
        }
    }

    return config;
}
