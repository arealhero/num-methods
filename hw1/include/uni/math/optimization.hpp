#pragma once

#include "matrix.hpp"
#include "numbers.hpp"

namespace uni::math::optimization
{
  template<std::size_t Size>
  constexpr vector<Size> steepest_descent(const matrix<Size, Size>& A,
					  const vector<Size>& b,
					  const double eps)
  {
    auto delta_container = A.get_diagonal_dominance_factor();
    if (!delta_container.has_value())
      assert(0 && "A should be diagonally dominant");

    auto delta = delta_container.value();

    auto x = vector<Size>::zero();
    auto q = A * x + b;

    while (q.norm_2() / delta > eps)
      {
	auto q_t = q.T();
	auto mu = -((q_t * q) / (q_t * (A * q)));
	x += mu * q;
	q = A * x + b;
      }

    return x;
  }

  template<std::size_t Size>
  constexpr vector<Size> coordinate_descent(const matrix<Size, Size>& A,
					    const vector<Size>& b,
					    const double eps)
  {
    auto delta_container = A.get_diagonal_dominance_factor();
    if (!delta_container.has_value())
      assert(0 && "A should be diagonally dominant");

    auto delta = delta_container.value();

    auto x = vector<Size>::zero();
    auto gradient = A * x + b;

    std::size_t i = 0;
    auto q = vector<Size>::e(i);

    while (gradient.norm_2() / delta > eps)
      {
	auto q_t = q.T();
	auto mu = -((q_t * gradient) / (q_t * (A * q)));
	x += mu * q;
	gradient = A * x + b;

	i = (i + 1) % Size;
	q = vector<Size>::e(i);
      }

    return x;
  }
}
