#ifndef _ASSERT_H_
#define _ASSERT_H_

#include <stdbool.h>

int assert(bool target, const char *function_name, const char *file_name, int line_number);

#endif
