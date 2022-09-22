#pragma once

#include "numbers.hpp"

#include <array>
#include <cassert>
#include <cstdint>
#include <initializer_list>
#include <optional>
#include <ostream>
#include <stdexcept>
#include <type_traits>

namespace uni::math
{
template <std::size_t Rows, std::size_t Cols = Rows>
class matrix
{
 public:
  using value_type = double;

  constexpr explicit matrix(const value_type value = 0)
  {
    apply_all([value](auto& element) { element = value; });
  }
  constexpr matrix(const std::initializer_list<value_type>& numbers)
  {
    std::size_t i = 0;
    for (const auto& number : numbers)
    {
      m_data[i++] = number;
    }
  }

  constexpr value_type& at(std::size_t i, std::size_t j)
  {
    assert(i < Rows && j < Cols);
    return m_data[i * Cols + j];
  }
  constexpr const value_type& at(std::size_t i, std::size_t j) const
  {
    assert(i < Rows && j < Cols);
    return m_data[i * Cols + j];
  }

  constexpr value_type& at(std::size_t index)
    requires(Rows == 1 || Cols == 1)
  {
    assert(index < Rows * Cols);
    return m_data[index];
  }
  constexpr const value_type& at(std::size_t index) const
    requires(Rows == 1 || Cols == 1)
  {
    assert(index < Rows * Cols);
    return m_data[index];
  }

  constexpr auto& operator+=(const value_type number)
  {
    apply_all([number](auto& element) { element += number; });
    return *this;
  }
  constexpr auto& operator-=(const value_type number)
  {
    apply_all([number](auto& element) { element -= number; });
    return *this;
  }
  constexpr auto& operator*=(const value_type number)
  {
    apply_all([number](auto& element) { element *= number; });
    return *this;
  }
  constexpr auto& operator/=(const value_type number)
  {
    apply_all([number](auto& element) { element /= number; });
    return *this;
  }

  constexpr auto& operator+=(const matrix<Rows, Cols>& other)
  {
    for (std::size_t row = 0; row < Rows; ++row)
      for (std::size_t col = 0; col < Cols; ++col)
        at(row, col) += other.at(row, col);
    return *this;
  }
  constexpr auto& operator-=(const matrix<Rows, Cols>& other)
  {
    return *this += -other;
  }

  constexpr value_type norm_1() const
  {
    value_type result{0};

    for (std::size_t j = 0; j < Cols; ++j)
    {
      value_type sum{0};

      for (std::size_t i = 0; i < Rows; ++i)
      {
        sum += uni::math::abs(at(i, j));
      }

      if (sum > result)
      {
        result = sum;
      }
    }

    return result;
  }

  constexpr value_type norm_2() const
  {
    value_type result{0};

    for (std::size_t j = 0; j < Cols; ++j)
    {
      value_type sum{0};

      for (std::size_t i = 0; i < Rows; ++i)
      {
        auto element = at(i, j);
        sum += element * element;
      }

      if (sum > result)
        result = sum;
    }

    return sqrt(result);
  }

  constexpr matrix<Rows, Cols> operator-() const
  {
    auto result = *this;
    result.apply_all([](auto& element) { element = -element; });
    return result;
  }

  constexpr operator value_type() const requires(Rows == 1 && Cols == 1)
  {
    return at(0, 0);
  }

  constexpr auto T() const
  {
    auto result = matrix<Cols, Rows>();
    for (std::size_t row = 0; row < Rows; ++row)
      for (std::size_t col = 0; col < Cols; ++col)
        result.at(col, row) = at(row, col);
    return result;
  }

  constexpr void swap_rows(std::size_t i, std::size_t j)
  {
    assert(i < Rows && j < Rows);

    if (i == j)
      return;

    for (std::size_t col = 0; col < Cols; ++col)
    {
      auto temp = at(i, col);
      at(i, col) = at(j, col);
      at(j, col) = temp;
    }
  }

  constexpr std::optional<value_type> get_diagonal_dominance_factor() const
  {
    value_type min_factor = abs(at(0, 0));

    for (std::size_t row = 0; row < Rows; ++row)
    {
      value_type row_factor = 0;
      for (std::size_t col = 0; col < Cols; ++col)
      {
        auto current_value = abs(at(row, col));
        if (row == col)
          row_factor += current_value;
        else
          row_factor -= current_value;
      }

      if (row_factor < 0)
        return std::nullopt;

      if (row_factor < min_factor)
        min_factor = row_factor;
    }

    return min_factor;
  }

  constexpr bool is_diagonally_dominant() const
  {
    return get_diagonal_dominance_factor().has_value();
  }

  constexpr matrix<Rows, Cols> inv() const
  {
    throw new std::runtime_error{"Not implemented yet"};
  }

  constexpr static matrix<Rows, Cols> zero() { return matrix<Rows, Cols>{}; }

  constexpr static matrix<Rows, Cols> identity()
    requires(Rows == Cols)
  {
    auto result = matrix<Rows, Cols>(0);
    for (std::size_t i = 0; i < Rows; ++i)
      result.at(i, i) = 1;
    return result;
  }

  constexpr static matrix<Rows, Cols> e(std::size_t i)
    requires(Rows == 1 || Cols == 1)
  {
    assert(i < Rows * Cols);
    auto result = matrix<Rows, Cols>::zero();
    result.at(i) = 1;
    return result;
  }

 private:
  value_type m_data[Rows * Cols];

  constexpr void apply_all(auto&& func)
  {
    for (std::size_t i = 0; i < Rows * Cols; ++i)
      func(m_data[i]);
  }
};

template <std::size_t Size, bool is_row = false>
using vector = matrix<is_row ? 1 : Size, is_row ? Size : 1>;

template <std::size_t Rows, std::size_t Cols>
constexpr matrix<Rows, Cols> operator+(const matrix<Rows, Cols>& lhs,
                                       const matrix<Rows, Cols>& rhs)
{
  matrix<Rows, Cols> result = lhs;
  return result += rhs;
}
template <std::size_t Rows, std::size_t Cols>
constexpr matrix<Rows, Cols> operator-(const matrix<Rows, Cols>& lhs,
                                       const matrix<Rows, Cols>& rhs)
{
  return lhs + -rhs;
}

template <std::size_t Rows, std::size_t K, std::size_t Cols>
constexpr matrix<Rows, Cols> operator*(const matrix<Rows, K>& lhs,
                                       const matrix<K, Cols>& rhs)
{
  matrix<Rows, Cols> result;
  for (std::size_t row = 0; row < Rows; ++row)
    for (std::size_t col = 0; col < Cols; ++col)
      for (std::size_t k = 0; k < K; ++k)
        result.at(row, col) += lhs.at(row, k) * rhs.at(k, col);
  return result;
}

template <std::size_t Rows, std::size_t Cols>
constexpr matrix<Rows, Cols> operator*(
    const typename matrix<Rows, Cols>::value_type number,
    const matrix<Rows, Cols>& mat)
{
  matrix<Rows, Cols> result = mat;
  return result *= number;
}

template <std::size_t Rows, std::size_t Cols>
std::ostream& operator<<(std::ostream& out,
                         const uni::math::matrix<Rows, Cols>& matrix)
{
  for (std::size_t i = 0; i < Rows; ++i)
  {
    for (std::size_t j = 0; j < Cols; ++j)
      out << matrix.at(i, j) << ' ';
    out << '\n';
  }
  return out;
}
}  // namespace uni::math
