#include <algebra.h>
#include <gtest/gtest.h>
#include <types.h>

#include <array>
#include <cstdint>

template <std::size_t N>
constexpr auto assert_equal(const std::array<value_t, N>& expected,
                            const std::array<value_t, N>& actual) -> void
{
  for (std::size_t i = 0; i < N; ++i)
  {
    ASSERT_TRUE(almost_equal(expected.at(i), actual.at(i)));
  }
}

TEST(Algebra_cardano, example_1)
{
  const std::array<value_t, 3> expected{1.L, 2.L, 3.L};
  const auto actual = cardano({1.L, -6.L, 11.L, -6.L});

  assert_equal(expected, actual);
}

TEST(Algebra_cardano, example_2)
{
  const std::array<value_t, 3> expected{-3.L, -2.L, -1.L};
  const auto actual = cardano({1.L, 6.L, 11.L, 6.L});

  assert_equal(expected, actual);
}
