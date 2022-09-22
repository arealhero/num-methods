#pragma once

#include <cassert>
#include <concepts>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <type_traits>

namespace uni::math
{
template <typename T>
concept decimal = std::integral<T> || std::floating_point<T>;

constexpr auto abs(decimal auto value)
{
  assert(value != std::numeric_limits<decltype(value)>::min());
  return (value > 0 ? value : -value);
}

constexpr auto power(const auto& base, int32_t pow)
{
  std::remove_cv_t<decltype(base)> result{1};

  if (pow < 0)
  {
    while (pow < 0)
    {
      result /= base;
      ++pow;
    }
  }
  else if (pow > 0)
  {
    while (pow > 0)
    {
      result *= base;
      --pow;
    }
  }
  else
  {
    assert(base != 0 && "0 in the 0th power is undefined");
  }

  return result;
}

constexpr auto exp(const std::floating_point auto x, const double eps = 1e-6)
{
  std::remove_cv_t<decltype(x)> prev = 1, next = 0;
  long k = 0;

  do
  {
    prev += next;
    next = /* TODO */ 0;
    throw std::runtime_error{"Not implemented yet"};
  } while (abs(next) > eps);

  return prev;
}

constexpr auto sqrt(const std::floating_point auto x, const double eps = 1e-6)
{
  assert(x >= 0);

  std::remove_cv_t<decltype(x)> prev, next = x;

  do
  {
    prev = next;
    next = (prev + x / prev) / 2;
  } while (abs(prev - next) > eps);

  return next;
}

}  // namespace uni::math
