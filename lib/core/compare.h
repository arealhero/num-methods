#pragma once

#include "types.h"

#include <cmath>
#include <concepts>

// FIXME: add `almost_zero`

template <std::floating_point Type>
auto almost_equal(const Type x, const Type y, const u32 ulp = 4) -> bool
{
  constexpr auto eps = std::numeric_limits<Type>::epsilon();

  // See https://en.cppreference.com/w/cpp/types/numeric_limits/epsilon
  return std::fabs(x - y) <= eps * std::fabs(x + y) * ulp ||
         std::fabs(x - y) < std::numeric_limits<Type>::min();
}

template <std::integral Type>
auto almost_equal(const Type x,
                  const Type y,
                  [[maybe_unused]] const u32 ulp = 4) -> bool
{
  return x == y;
}
