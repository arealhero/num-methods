#pragma once

#include <cmath>
#include <cstdint>

#include <core/types.h>

[[maybe_unused]] constexpr double ALPHA = 2. / 5.;
[[maybe_unused]] constexpr double BETA = 0.;

constexpr double ORIG_A = 2.1;
constexpr double ORIG_B = 3.3;
constexpr double EXACT_VALUE = 2.95073;
constexpr double EXACT_VALUE_WITH_P = 4.46151;
constexpr u32 MAX_PARTITIONS = 1'000;

using FunctionType = double (*)(double x);

constexpr auto f(double x) -> double
{
  return 4.5 * std::cos(7. * x) * std::exp(-2. / 3. * x) +
         1.4 * std::sin(1.5 * x) * std::exp(-1. / 3. * x) + 3.;
}
