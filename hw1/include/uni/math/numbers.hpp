#pragma once

#include <concepts>
#include <limits>
#include <cstdint>
#include <cassert>

namespace uni::math
{
  template<typename T>
  concept decimal = std::integral<T> || std::floating_point<T>;

  constexpr auto abs(decimal auto value)
  {
    assert(value != std::numeric_limits<decltype(value)>::min());
    return (value > 0 ? value : -value);
  }

  template<typename T>
  constexpr T power(const T base, int32_t pow) {
    T result{1};

    if (pow < 0) {
      while (pow < 0) {
	result /= base;
	++pow;
      }
    } else if (pow > 0) {
      while (pow > 0) {
	result *= base;
	--pow;
      }
    } else {
      assert(base != 0 && "0 in the 0th power is undefined");
    }

    return result;
  }

  constexpr auto sqrt(std::floating_point auto x, double eps = 1e-6)
  {
    assert(x >= 0);

    decltype(x) prev, next = x;

    do
      {
	prev = next;
	next = (prev + x / prev) / 2;
      }
    while (abs(prev - next) > eps);

    return next;
  }

}
