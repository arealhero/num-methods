#pragma once

#include <cmath>
#include <concepts>
#include <limits>
#include <memory>
#include <utility>

template <typename Type>
using Ptr = std::shared_ptr<Type>;

template <typename Type, typename... Args>
auto make(Args&&... args) -> Ptr<Type>
{
  return std::make_shared<Type>(std::forward<Args>(args)...);
}

template <std::floating_point Type>
auto almost_equal(const Type x, const Type y, const std::size_t ulp = 4) -> bool
{
  constexpr auto eps = std::numeric_limits<Type>::epsilon();

  // See https://en.cppreference.com/w/cpp/types/numeric_limits/epsilon
  return std::fabs(x - y) <= eps * std::fabs(x + y) * ulp
         || std::fabs(x - y) < std::numeric_limits<Type>::min();
}

template <std::integral Type>
auto almost_equal(const Type x,
                  const Type y,
                  [[maybe_unused]] const std::size_t ulp = 4) -> bool
{
  return x == y;
}
