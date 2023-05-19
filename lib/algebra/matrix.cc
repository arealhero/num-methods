#include "matrix.h"

#include <core/compare.h>

Matrix::Matrix(const u64 rows, const u64 cols, const double value)
    : m_matrix(Eigen::MatrixXd::Constant(static_cast<i64>(rows),
                                         static_cast<i64>(cols),
                                         value))
{
}

// FIXME: remove NOLINT
// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
Matrix::Matrix(const u64 rows,
               const u64 cols,
               const std::initializer_list<double>& numbers)
    : m_matrix(static_cast<i64>(rows), static_cast<i64>(cols))
{
  assert(numbers.size() == rows * cols && "not enough elements provided");

  for (auto number : numbers)
  {
    m_matrix << number;
  }
}

Matrix::Matrix(Eigen::MatrixXd matrix) : m_matrix(std::move(matrix)) {}

[[nodiscard]] auto Matrix::zero(const u64 size) -> Matrix
{
  return zero(size, size);
}
[[nodiscard]] auto Matrix::zero(const u64 rows, const u64 cols) -> Matrix
{
  return {rows, cols};
}

[[nodiscard]] auto Matrix::identity(const u64 size) -> Matrix
{
  auto result = zero(size);
  for (u64 i = 0; i < size; ++i)
  {
    result.at(i, i) = 1;
  }
  return result;
}

// FIXME: remove NOLINT
// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
[[nodiscard]] auto Matrix::e(const u64 size, const u64 i) -> Matrix
{
  assert(i < size);
  auto result = zero(size);
  result.at(i) = 1;
  return result;
}

[[nodiscard]] auto Matrix::rows() const -> u64
{
  return static_cast<u64>(m_matrix.rows());
}
[[nodiscard]] auto Matrix::cols() const -> u64
{
  return static_cast<u64>(m_matrix.cols());
}
[[nodiscard]] auto Matrix::is_square() const -> bool
{
  return rows() == cols();
}

[[nodiscard]] auto Matrix::at(u64 row, u64 col) -> double&
{
  return m_matrix(static_cast<i64>(row), static_cast<i64>(col));
}
[[nodiscard]] auto Matrix::at(u64 row, u64 col) const -> const double&
{
  return m_matrix(static_cast<i64>(row), static_cast<i64>(col));
}

// FIXME: remove code duplication
[[nodiscard]] auto Matrix::at(u64 index) -> double&
{
  if (rows() == 1)
  {
    return m_matrix(0, static_cast<i64>(index));
  }

  if (cols() == 1)
  {
    return m_matrix(static_cast<i64>(index), 0);
  }

  assert(0 && "at(index) is only available for vectors");
  std::abort();  // FIXME: change to std::unreachable (C++23)
}
[[nodiscard]] auto Matrix::at(u64 index) const -> const double&
{
  if (rows() == 1)
  {
    return m_matrix(1, static_cast<i64>(index));
  }

  if (cols() == 1)
  {
    return m_matrix(static_cast<i64>(index), 1);
  }

  assert(0 && "at(index) is only available for vectors");
  std::abort();  // FIXME: change to std::unreachable (C++23)
}

Matrix::operator double() const
{
  assert(rows() == 1 && cols() == 1);
  if (rows() != 1 || cols() != 1)
  {
    throw std::runtime_error{"cannot convert Matrix to double"};
  }

  return at(0, 0);
}

auto Matrix::operator+=(const double number) -> Matrix&
{
  m_matrix.array() += number;
  return *this;
}
auto Matrix::operator-=(const double number) -> Matrix&
{
  m_matrix.array() -= number;
  return *this;
}
auto Matrix::operator*=(const double number) -> Matrix&
{
  m_matrix.array() *= number;
  return *this;
}
auto Matrix::operator/=(const double number) -> Matrix&
{
  m_matrix.array() /= number;
  return *this;
}

auto Matrix::operator+=(const Matrix& other) -> Matrix&
{
  assert(rows() == other.rows() && cols() == other.cols());

  for (u64 row = 0; row < rows(); ++row)
  {
    for (u64 col = 0; col < cols(); ++col)
    {
      at(row, col) += other.at(row, col);
    }
  }

  return *this;
}
auto Matrix::operator-=(const Matrix& other) -> Matrix&
{
  return *this += -other;
}

[[nodiscard]] auto Matrix::operator-() const -> Matrix
{
  auto result = *this;
  result.m_matrix.array() *= -1;
  return result;
}

[[nodiscard]] auto Matrix::norm1() const -> double
{
  double result{0};

  for (u64 j = 0; j < cols(); ++j)
  {
    double sum{0};

    for (u64 i = 0; i < rows(); ++i)
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

[[nodiscard]] auto Matrix::norm2() const -> double
{
  return m_matrix.norm();
}

[[nodiscard]] auto Matrix::T() const -> Matrix
{
  auto result = Matrix(cols(), rows());
  for (u64 row = 0; row < rows(); ++row)
  {
    for (u64 col = 0; col < cols(); ++col)
    {
      // NOLINTNEXTLINE(readability-suspicious-call-argument)
      result.at(col, row) = at(row, col);
    }
  }
  return result;
}
[[nodiscard]] auto Matrix::inv() const -> Matrix
{
  return {m_matrix.inverse()};
}

void Matrix::swap_rows(u64 i, u64 j)
{
  assert(i < rows() && j < rows());

  if (i == j)
  {
    return;
  }

  for (u64 col = 0; col < cols(); ++col)
  {
    auto temp = at(i, col);
    at(i, col) = at(j, col);
    at(j, col) = temp;
  }
}

[[nodiscard]] auto Matrix::get_diagonal_dominance_factor() const
    -> std::optional<double>
{
  double min_factor = std::abs(at(0, 0));

  for (u64 row = 0; row < rows(); ++row)
  {
    double row_factor = 0;
    for (u64 col = 0; col < cols(); ++col)
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

[[nodiscard]] auto Matrix::is_diagonally_dominant() const -> bool
{
  return get_diagonal_dominance_factor().has_value();
}

auto operator+(const Matrix& lhs, const Matrix& rhs) -> Matrix
{
  assert(lhs.rows() == rhs.rows() && lhs.cols() == rhs.cols());
  auto result = lhs;
  return result += rhs;
}
auto operator-(const Matrix& lhs, const Matrix& rhs) -> Matrix
{
  return lhs + -rhs;
}

auto operator==(const Matrix& lhs, const Matrix& rhs) -> bool
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

auto operator*(const Matrix& lhs, const Matrix& rhs) -> Matrix
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

auto operator*(const double number, const Matrix& mat) -> Matrix
{
  auto result = mat;
  return result *= number;
}

auto operator/(const Matrix& mat, double number) -> Matrix
{
  auto result = mat;
  return result /= number;
}

std::ostream& operator<<(std::ostream& out, const Matrix& matrix)
{
  for (std::size_t row = 0; row < matrix.rows(); ++row)
  {
    for (std::size_t col = 0; col < matrix.cols(); ++col)
    {
      out << matrix.at(row, col) << ' ';
    }

    out << '\n';
  }

  return out;
}
