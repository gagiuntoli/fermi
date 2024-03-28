
#include "toml.h"
#include <stdbool.h>
#include <string.h>

bool toml_equal_table_names(TomlTableNames *tn1, TomlTableNames *tn2) {
  int quantity1 = tn1->quantity;
  int quantity2 = tn2->quantity;

  if (quantity1 != quantity2) {
    return false;
  }

  for (int i = 0; i < quantity1; i++) {
    if (strcmp(tn1->names[i], tn2->names[i]) != 0) {
      return false;
    }
  }

  return true;
}
