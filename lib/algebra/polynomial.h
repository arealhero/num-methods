#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <numbers>

/**
 * Solves ax^3 + bx^2 + cx + d = 0
 * coeffs: { a, b, c, d }
 **/
constexpr auto cardano(const std::array<double, 4>& coeffs)
    -> std::array<double, 3>
{
  constexpr auto pi = std::numbers::pi;
  [[maybe_unused]] constexpr auto eps = 1e-10;

  const auto a = coeffs.at(0);
  const auto b = coeffs.at(1);
  const auto c = coeffs.at(2);
  const auto d = coeffs.at(3);

  const auto p = (3 * a * c - b * b) / (9 * a * a);
  const auto q =
      (2 * b * b * b - 9 * a * b * c + 27 * a * a * d) / (54 * a * a * a);

  // FIXME: change to almost_zero (see core/utils.h)
  assert(q * q + p * p * p < eps);  // Three real distinct roots

  const auto r = ((q < 0) ? -1 : 1) * std::sqrt(std::abs(p));
  const auto phi = std::acos(q / (r * r * r));

  std::array<double, 3> roots = {
      -2 * r * std::cos(phi / 3) - b / (3 * a),
      2 * r * std::cos(pi / 3 - phi / 3) - b / (3 * a),
      2 * r * std::cos(pi / 3 + phi / 3) - b / (3 * a),
  };
  std::sort(roots.begin(), roots.end());

  return roots;
}
