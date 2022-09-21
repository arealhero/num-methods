#pragma once

#include "matrix.hpp"
#include "numbers.hpp"

#include <iostream>

namespace uni::math::sle
{
template <std::size_t Size>
constexpr vector<Size> seidel(const matrix<Size>& A,
                              const vector<Size>& b,
                              const double eps)
{
  auto c = matrix<Size>::zero();
  auto d = vector<Size>::zero();

  for (std::size_t i = 0; i < Size; ++i)
  {
    auto dividor = A.at(i, i);
    d.at(i) = b.at(i) / dividor;

    for (std::size_t j = 0; j < Size; ++j)
    {
      if (i == j)
        continue;
      c.at(i, j) = -A.at(i, j) / dividor;
    }
  }

  auto prev = vector<Size>{1};
  auto x = d;
  auto diff = vector<Size>{1};

  double norm = c.norm_1();
  double factor = uni::math::abs(norm / (1 - norm));

  do
  {
    prev = x;
    x = c * x + d;
    diff = x - prev;
  } while (factor * diff.norm_1() > eps);

  return x;
}

template <std::size_t Size>
constexpr void gauss_down(matrix<Size>& a, vector<Size>& b)
{
  for (std::size_t index = 0; index < Size; ++index)
  {
    // 1. Find max element in column
    std::size_t row_to_swap = index;
    auto max_element = a.at(index, index);
    for (std::size_t row = index + 1; row < Size; ++row)
    {
      auto element = a.at(row, index);
      if (abs(element) > abs(max_element))
      {
        max_element = element;
        row_to_swap = row;
      }
    }

    // 2. Swap current row with that row
    if (row_to_swap != index)
    {
      a.swap_rows(index, row_to_swap);
      b.swap_rows(index, row_to_swap);
    }

    // 3. Normalize row
    auto divider = a.at(index, index);

    for (std::size_t col = index; col < Size; ++col)
      a.at(index, col) /= divider;

    b.at(index) /= divider;

    // 4. Subtract that row from the rows below
    for (std::size_t row = index + 1; row < Size; ++row)
    {
      auto factor = a.at(row, index);

      for (std::size_t col = index; col < Size; ++col)
        a.at(row, col) -= factor * a.at(index, col);

      b.at(row) -= factor * b.at(index);
    }
  }
}

template <std::size_t Size>
constexpr void gauss_up(matrix<Size>& a, vector<Size>& b)
{
  for (std::size_t index = Size - 1; index != 0; --index)
  {
    for (std::size_t row = index; row != 0; --row)
    {
      auto factor = a.at(row - 1, index);

      for (std::size_t col = index; col < Size; ++col)
        a.at(row - 1, col) -= factor * a.at(index, col);

      b.at(row - 1) -= factor * b.at(index);
    }
  }
}

template <std::size_t Size>
constexpr vector<Size> gauss(matrix<Size> a, vector<Size> b)
{
  gauss_down(a, b);
  gauss_up(a, b);

  return b;
}

}  // namespace uni::math::sle
