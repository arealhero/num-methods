#include <algebra/polynomial.h>
#include <core/compare.h>
#include <gtest/gtest.h>

#include <array>
#include <cstdint>

template <std::size_t N>
constexpr auto assert_equal(const std::array<double, N>& expected,
                            const std::array<double, N>& actual) -> void
{
  for (std::size_t i = 0; i < N; ++i)
  {
    ASSERT_TRUE(almost_equal(expected.at(i), actual.at(i)));
  }
}

TEST(Algebra_cardano, example_1)
{
  const std::array<double, 3> expected{1., 2., 3.};
  const auto actual = cardano({1., -6., 11., -6.});

  assert_equal(expected, actual);
}

TEST(Algebra_cardano, example_2)
{
  const std::array<double, 3> expected{-3., -2., -1.};
  const auto actual = cardano({1., 6., 11., 6.});

  assert_equal(expected, actual);
}
