#include "gtest/gtest.h"

// Helper function for implementing ASSERT_NEAR.
template <typename T>
::testing::AssertionResult NearPredFormat(const char* expr1,
                                     const char* expr2,
                                     const char* abs_error_expr,
                                     T val1,
                                     T val2,
                                     double abs_error) {
  const double diff = std::abs(val1 - val2);
  if (diff <= abs_error) return ::testing::AssertionSuccess();

  return ::testing::AssertionFailure()
      << "The difference between " << expr1 << " and " << expr2
      << " is " << diff << ", which exceeds " << abs_error_expr << ", where\n"
      << expr1 << " evaluates to " << val1 << ",\n"
      << expr2 << " evaluates to " << val2 << ".";
}

#undef EXPECT_NEAR
#define EXPECT_NEAR(val1, val2, abs_error) \
  EXPECT_PRED_FORMAT3(NearPredFormat, val1, val2, abs_error)


