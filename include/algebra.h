#pragma once

#include "types.h"
#include "utils.h"

#include <Eigen/Dense>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <functional>
#include <initializer_list>
#include <numbers>
#include <optional>
#include <ostream>
#include <vector>

class Matrix
{
 public:
  constexpr explicit Matrix(const std::size_t rows,
                            const std::size_t cols,
                            const value_t value = 0)
      : m_data(rows * cols, value), m_rows(rows), m_cols(cols)
  {
  }

  // FIXME: remove NOLINT
  // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
  constexpr Matrix(const std::size_t rows,
                   const std::size_t cols,
                   const std::initializer_list<value_t>& numbers)
      : m_data(numbers), m_rows(rows), m_cols(cols)
  {
    assert(m_data.size() == m_rows * m_cols);
  }

  [[nodiscard]] constexpr auto rows() const -> std::size_t { return m_rows; }
  [[nodiscard]] constexpr auto cols() const -> std::size_t { return m_cols; }
  [[nodiscard]] constexpr auto is_square() const -> bool
  {
    return rows() == cols();
  }

  [[nodiscard]] constexpr auto size() const -> std::size_t
  {
    // FIXME: add support for vectors (rows and cols)
    assert(is_square());
    return m_rows;
  }

  [[nodiscard]] constexpr auto at(std::size_t i, std::size_t j) -> value_t&
  {
    assert(i < rows() && j < cols());
    return m_data[i * cols() + j];
  }
  [[nodiscard]] constexpr auto at(std::size_t i, std::size_t j) const
      -> const value_t&
  {
    assert(i < rows() && j < cols());
    return m_data[i * cols() + j];
  }

  [[nodiscard]] constexpr auto at(std::size_t index) -> value_t&
  {
    assert((rows() == 1 || cols() == 1) && index < rows() * cols());
    return m_data[index];
  }
  [[nodiscard]] constexpr auto at(std::size_t index) const -> const value_t&
  {
    assert((rows() == 1 || cols() == 1) && index < rows() * cols());
    return m_data[index];
  }

  constexpr auto operator+=(const value_t number) -> Matrix&
  {
    apply_all([number](auto& element) { element += number; });
    return *this;
  }
  constexpr auto operator-=(const value_t number) -> Matrix&
  {
    apply_all([number](auto& element) { element -= number; });
    return *this;
  }
  constexpr auto operator*=(const value_t number) -> Matrix&
  {
    apply_all([number](auto& element) { element *= number; });
    return *this;
  }
  constexpr auto operator/=(const value_t number) -> Matrix&
  {
    apply_all([number](auto& element) { element /= number; });
    return *this;
  }

  constexpr auto operator+=(const Matrix& other) -> Matrix&
  {
    assert(rows() == other.rows() && cols() == other.cols());

    for (std::size_t row = 0; row < rows(); ++row)
    {
      for (std::size_t col = 0; col < cols(); ++col)
      {
        at(row, col) += other.at(row, col);
      }
    }

    return *this;
  }
  constexpr auto operator-=(const Matrix& other) -> Matrix&
  {
    return *this += -other;
  }

  [[nodiscard]] constexpr auto norm_1() const -> value_t
  {
    value_t result{0};

    for (std::size_t j = 0; j < cols(); ++j)
    {
      value_t sum{0};

      for (std::size_t i = 0; i < rows(); ++i)
      {
        sum += std::abs(at(i, j));
      }

      if (sum > result)
      {
        result = sum;
      }
    }

    return result;
  }

  [[nodiscard]] constexpr auto norm_2() const -> value_t
  {
    value_t result{0};

    for (std::size_t j = 0; j < cols(); ++j)
    {
      value_t sum{0};

      for (std::size_t i = 0; i < rows(); ++i)
      {
        auto element = at(i, j);
        sum += element * element;
      }

      if (sum > result)
      {
        result = sum;
      }
    }

    return std::sqrt(result);
  }

  [[nodiscard]] constexpr auto operator-() const -> Matrix
  {
    auto result = *this;
    result.apply_all([](auto& element) { element = -element; });
    return result;
  }

  constexpr operator value_t() const
  {
    assert(rows() == 1 && cols() == 1);
    if (rows() != 1 || cols() != 1)
    {
      throw std::runtime_error{"cannot convert Matrix to value_t"};
    }

    return at(0, 0);
  }

  [[nodiscard]] constexpr auto T() const -> Matrix
  {
    auto result = Matrix(cols(), rows());
    for (std::size_t row = 0; row < rows(); ++row)
    {
      for (std::size_t col = 0; col < cols(); ++col)
      {
        result.at(col, row) = at(row, col);
      }
    }
    return result;
  }

  constexpr void swap_rows(std::size_t i, std::size_t j)
  {
    assert(i < rows() && j < rows());

    if (i == j)
    {
      return;
    }

    for (std::size_t col = 0; col < cols(); ++col)
    {
      auto temp = at(i, col);
      at(i, col) = at(j, col);
      at(j, col) = temp;
    }
  }

  [[nodiscard]] constexpr auto get_data() const -> std::vector<value_t>
  {
    return m_data;
  }

  [[nodiscard]] constexpr auto get_diagonal_dominance_factor() const
      -> std::optional<value_t>
  {
    value_t min_factor = std::abs(at(0, 0));

    for (std::size_t row = 0; row < rows(); ++row)
    {
      value_t row_factor = 0;
      for (std::size_t col = 0; col < cols(); ++col)
      {
        auto current_value = std::abs(at(row, col));
        if (row == col)
        {
          row_factor += current_value;
        }
        else
        {
          row_factor -= current_value;
        }
      }

      if (row_factor < 0)
      {
        return std::nullopt;
      }

      if (row_factor < min_factor)
      {
        min_factor = row_factor;
      }
    }

    return min_factor;
  }

  [[nodiscard]] constexpr auto is_diagonally_dominant() const -> bool
  {
    return get_diagonal_dominance_factor().has_value();
  }

  [[nodiscard]] constexpr auto inv() const -> Matrix;

  constexpr static auto zero(const std::size_t size) -> Matrix
  {
    return zero(size, size);
  }
  constexpr static auto zero(const std::size_t rows, const std::size_t cols)
      -> Matrix
  {
    return Matrix(rows, cols);
  }

  constexpr static auto identity(const std::size_t size) -> Matrix
  {
    auto result = zero(size);
    for (std::size_t i = 0; i < size; ++i)
    {
      result.at(i, i) = 1;
    }
    return result;
  }

  // FIXME: remove NOLINT
  // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
  constexpr static auto e(const std::size_t size, const std::size_t i) -> Matrix
  {
    assert(i < size);
    auto result = zero(size);
    result.at(i) = 1;
    return result;
  }
  constexpr static auto e_T(const std::size_t size, const std::size_t i)
      -> Matrix
  {
    return e(size, i).T();
  }

 private:
  std::vector<value_t> m_data;
  std::size_t m_rows;
  std::size_t m_cols;

  constexpr void apply_all(auto&& func)
  {
    std::for_each(m_data.begin(), m_data.end(), func);
  }
};

constexpr auto operator+(const Matrix& lhs, const Matrix& rhs) -> Matrix
{
  assert(lhs.rows() == rhs.rows() && lhs.cols() == rhs.cols());
  auto result = lhs;
  return result += rhs;
}
constexpr auto operator-(const Matrix& lhs, const Matrix& rhs) -> Matrix
{
  return lhs + -rhs;
}

constexpr auto operator==(const Matrix& lhs, const Matrix& rhs) -> bool
{
  if (lhs.rows() != rhs.rows() || lhs.cols() != rhs.cols())
  {
    return false;
  }

  for (std::size_t row = 0; row < lhs.rows(); ++row)
  {
    for (std::size_t col = 0; col < lhs.cols(); ++col)
    {
      if (almost_equal(lhs.at(row, col), rhs.at(row, col)))
      {
        return false;
      }
    }
  }

  return true;
}

constexpr auto operator*(const Matrix& lhs, const Matrix& rhs) -> Matrix
{
  assert(lhs.cols() == rhs.rows());
  const auto K = lhs.cols();

  auto result = Matrix::zero(lhs.rows(), rhs.cols());
  for (std::size_t row = 0; row < lhs.rows(); ++row)
  {
    for (std::size_t col = 0; col < rhs.cols(); ++col)
    {
      for (std::size_t k = 0; k < K; ++k)
      {
        result.at(row, col) += lhs.at(row, k) * rhs.at(k, col);
      }
    }
  }

  return result;
}

constexpr auto operator*(const value_t number, const Matrix& mat) -> Matrix
{
  auto result = mat;
  return result *= number;
}

constexpr auto gauss_down(Matrix& A, Matrix& b) -> void;
constexpr auto gauss_up(Matrix& A, Matrix& b) -> void;

[[nodiscard]] constexpr auto gauss(
    // NOLINTNEXTLINE(bugprone-easily-swappable-parameters,-warnings-as-errors)
    const Matrix& A,
    const Matrix& b) -> Matrix
{
  auto A_copy = A;
  auto b_copy = b;

  gauss_down(A_copy, b_copy);
  gauss_up(A_copy, b_copy);

  return b_copy;
}

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters,-warnings-as-errors)
constexpr auto gauss_down(Matrix& A, Matrix& b) -> void
{
  for (std::size_t index = 0; index < A.size(); ++index)
  {
    // 1. Find max element in column
    std::size_t row_to_swap = index;
    auto max_element = A.at(index, index);
    for (std::size_t row = index + 1; row < A.size(); ++row)
    {
      auto element = A.at(row, index);
      if (std::abs(element) > std::abs(max_element))
      {
        max_element = element;
        row_to_swap = row;
      }
    }

    // 2. Swap current row with that row
    if (row_to_swap != index)
    {
      A.swap_rows(index, row_to_swap);
      b.swap_rows(index, row_to_swap);
    }

    // 3. Normalize row
    auto divider = A.at(index, index);

    for (std::size_t col = index; col < A.size(); ++col)
    {
      A.at(index, col) /= divider;
    }

    b.at(index) /= divider;

    // 4. Subtract that row from rows below
    for (std::size_t row = index + 1; row < A.size(); ++row)
    {
      auto factor = A.at(row, index);

      for (std::size_t col = index; col < A.size(); ++col)
      {
        A.at(row, col) -= factor * A.at(index, col);
      }

      b.at(row) -= factor * b.at(index);
    }
  }
}

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters,-warnings-as-errors)
constexpr auto gauss_up(Matrix& A, Matrix& b) -> void
{
  for (std::size_t index = A.size() - 1; index != 0; --index)
  {
    for (std::size_t row = index; row != 0; --row)
    {
      auto factor = A.at(row - 1, index);

      for (std::size_t col = index; col < A.size(); ++col)
      {
        A.at(row - 1, col) -= factor * A.at(index, col);
      }

      b.at(row - 1) -= factor * b.at(index);
    }
  }
}

/* // FIXME: remove NOLINT */
/* // NOLINTNEXTLINE(bugprone-easily-swappable-parameters) */
/* auto solve_sle(const Matrix& A, const Matrix& b) -> Matrix */
/* { */
/*   using Eigen::MatrixXd; */
/*   using Eigen::VectorXd; */
/*   using Eigen::Index; */

/*   MatrixXd eigen_A(A.rows(), A.cols()); */
/*   VectorXd eigen_b(b.rows()); */

/*   for (std::size_t row = 0; row < A.rows(); ++row) */
/*   { */
/*     for (std::size_t col = 0; col < A.cols(); ++col) */
/*     { */
/*       eigen_A(Index(row), Index(col)) = A.at(row, col); */
/*     } */

/*     eigen_b(row) = b.at(row); */
/*   } */
/* } */

// FIXME: remove NOLINT
// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
constexpr Matrix seidel(const Matrix& A,
                        const Matrix& b,
                        const value_t eps = 1e-6L)
{
  assert(A.is_square() && A.rows() == b.rows());

  const auto n = A.rows();

  auto c = Matrix(n, n);
  auto d = Matrix(n, 1);

  for (std::size_t i = 0; i < n; ++i)
  {
    auto dividor = A.at(i, i);
    d.at(i) = b.at(i) / dividor;

    for (std::size_t j = 0; j < n; ++j)
    {
      if (i == j)
      {
        continue;
      }
      c.at(i, j) = -A.at(i, j) / dividor;
    }
  }

  auto prev = Matrix(n, 1, 1.L);
  auto x = d;
  auto diff = Matrix(n, 1, 1.L);

  value_t norm = c.norm_1();
  value_t factor = std::abs(norm / (1 - norm));

  do
  {
    prev = x;
    x = c * x + d;
    diff = x - prev;
  } while (factor * diff.norm_1() > eps);

  return x;
}

/**
 * Solves ax^3 + bx^2 + cx + d = 0
 * coeffs: { a, b, c, d }
 **/
constexpr auto cardano(const std::array<value_t, 4>& coeffs)
    -> std::array<value_t, 3>
{
  constexpr auto pi = std::numbers::pi_v<value_t>;
  [[maybe_unused]] constexpr auto eps = 1e-10L;

  const auto a = coeffs.at(0);
  const auto b = coeffs.at(1);
  const auto c = coeffs.at(2);
  const auto d = coeffs.at(3);

  const auto p = (3 * a * c - b * b) / (9 * a * a);
  const auto q =
      (2 * b * b * b - 9 * a * b * c + 27 * a * a * d) / (54 * a * a * a);

  assert(q * q + p * p * p < eps);  // Three real distinct roots

  const auto r = ((q < 0) ? -1 : 1) * std::sqrt(std::abs(p));
  const auto phi = std::acos(q / (r * r * r));

  std::array<value_t, 3> roots = {
      -2 * r * std::cos(phi / 3) - b / (3 * a),
      2 * r * std::cos(pi / 3 - phi / 3) - b / (3 * a),
      2 * r * std::cos(pi / 3 + phi / 3) - b / (3 * a),
  };
  std::sort(roots.begin(), roots.end());

  return roots;
}
