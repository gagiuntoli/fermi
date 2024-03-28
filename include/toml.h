#ifndef _TOML_H_
#define _TOML_H_

#include <stdbool.h>

typedef struct {
  int quantity;
  char **names;
} TomlTableNames;

bool toml_equal_table_names(TomlTableNames *tn1, TomlTableNames *tn2);
TomlTableNames *toml_parse_table_names(char *toml_stream);

#endif
