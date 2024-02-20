#ifndef _ASSERT_H_
#define _ASSERT_H_

#include <stdbool.h>

int assert(bool target, char *message, char *file_name, int line_number);

#endif
