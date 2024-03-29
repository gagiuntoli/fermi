#include "assert.h"
#include "toml.h"
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

bool test_toml_equal_tables_names_are_equal() {
  int quantity = 3;
  char *name_0 = "servers";
  char *name_1 = "ports";
  char *name_2 = "databases";

  TomlTableNames *table_names_1 = (TomlTableNames *) malloc(sizeof(TomlTableNames));
  table_names_1->quantity = quantity;
  table_names_1->names = (char**) malloc(quantity * sizeof(char *));
  table_names_1->names[0] = strdup(name_0);
  table_names_1->names[1] = strdup(name_1);
  table_names_1->names[2] = strdup(name_2);

  TomlTableNames *table_names_2 = (TomlTableNames*) malloc(sizeof(TomlTableNames));
  table_names_2->quantity = quantity;
  table_names_2->names = (char**) malloc(quantity * sizeof(char *));
  table_names_2->names[0] = strdup(name_0);
  table_names_2->names[1] = strdup(name_1);
  table_names_2->names[2] = strdup(name_2);

  assert(toml_equal_table_names(table_names_1, table_names_2), __func__,
         __FILE__, __LINE__);

  return true;
}

bool test_toml_equal_tables_names_are_different() {
  int quantity = 3;
  char *name_0 = "servers";
  char *name_1 = "ports";
  char *name_2 = "databases";

  TomlTableNames *table_names_1 = (TomlTableNames*) malloc(sizeof(TomlTableNames));
  table_names_1->quantity = quantity;
  table_names_1->names = (char**) malloc(quantity * sizeof(char *));
  table_names_1->names[0] = strdup(name_0);
  table_names_1->names[1] = strdup(name_1);
  table_names_1->names[2] = strdup(name_2);

  char *name_2_diff = "database";

  TomlTableNames *table_names_2 = (TomlTableNames*) malloc(sizeof(TomlTableNames));
  table_names_2->quantity = quantity;
  table_names_2->names = (char**) malloc(quantity * sizeof(char *));
  table_names_2->names[0] = strdup(name_0);
  table_names_2->names[1] = strdup(name_1);
  table_names_2->names[2] = strdup(name_2_diff);

  assert(!toml_equal_table_names(table_names_1, table_names_2), __func__,
         __FILE__, __LINE__);

  return true;
}

bool test_toml_equal_tables_names_are_different_in_number_of_fields() {
  int quantity1 = 3;
  char *name_0 = "servers";
  char *name_1 = "ports";
  char *name_2 = "databases";

  TomlTableNames *table_names_1 = (TomlTableNames*) malloc(sizeof(TomlTableNames));
  table_names_1->quantity = quantity1;
  table_names_1->names = (char**) malloc(quantity1 * sizeof(char *));
  table_names_1->names[0] = strdup(name_0);
  table_names_1->names[1] = strdup(name_1);
  table_names_1->names[2] = strdup(name_2);

  int quantity2 = 2;
  TomlTableNames *table_names_2 = (TomlTableNames*) malloc(sizeof(TomlTableNames));
  table_names_2->quantity = quantity2;
  table_names_2->names = (char**) malloc(quantity2 * sizeof(char *));
  table_names_2->names[0] = strdup(name_0);
  table_names_2->names[1] = strdup(name_1);

  assert(!toml_equal_table_names(table_names_1, table_names_2), __func__,
         __FILE__, __LINE__);

  return true;
}

// int test_toml_get_table_names() {
//
//   char *toml_stream =
//     "[database]\n"
//     "server = \"192.168.1.1\"\n"
//     "port = 5432\n"
//     "[external-server]\n"
//     "server = \"192.168.1.2\"\n"
//     "port = 5433\n";
//
//   TomlTableNames table_names = get_table_names(toml_stream);
//   TomlTableNames expected_table_names = { };
//
//   toml_assert_table_name_equals(tables_names, )
// }

int main() {
  test_toml_equal_tables_names_are_equal();
  test_toml_equal_tables_names_are_different();
  test_toml_equal_tables_names_are_different_in_number_of_fields();
  // assert(test_toml_get_table_names());
}
