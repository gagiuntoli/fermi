#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

int assert(bool target, const char *message, const char *file_name,
           int line_number) {
  if (!target) {
    printf("%s: %d - %s\n", file_name, line_number, message);
    exit(1);
  }
  return 0;
}
