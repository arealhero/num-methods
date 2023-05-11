#pragma once

#include <types.h>

#include <cmath>
#include <cstdint>

[[maybe_unused]] constexpr value_t ALPHA = 2.L / 5.L;
[[maybe_unused]] constexpr value_t BETA = 0.L;

constexpr value_t ORIG_A = 2.1L;
constexpr value_t ORIG_B = 3.3L;
constexpr value_t EXACT_VALUE = 2.95073L;
constexpr value_t EXACT_VALUE_WITH_P = 4.46151L;
constexpr std::size_t MAX_PARTITIONS = 1'000;

constexpr auto f(value_t x) -> value_t
{
  return 4.5L * std::cos(7.L * x) * std::exp(-2.L / 3.L * x) +
         1.4L * std::sin(1.5L * x) * std::exp(-1.L / 3.L * x) + 3.L;
}
