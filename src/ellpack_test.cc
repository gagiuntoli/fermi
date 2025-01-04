#include "ellpack.h"

#include <gtest/gtest.h>

TEST(EllpackTest, DotProduct) {
  std::vector<double> x = {1.0, 2.0, 3.0};
  std::vector<double> y = {4.0, 5.0, 6.0};
  double result = dot(y, x, 3);
  EXPECT_DOUBLE_EQ(result, 32.0);  // 1*4 + 2*5 + 3*6 = 32
}
