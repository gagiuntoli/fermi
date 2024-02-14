

#include "fermi.h"
#include <stdio.h>
#include <stdlib.h>

int main() {
  printf("Running test success\n");
  return 1;
}

typedef struct {
} Config;

Config *parse_input_1(FILE *fd) {

  Config *config = malloc(sizeof(Config));
  if (fd == NULL)
    return NULL;

  return config;
}
