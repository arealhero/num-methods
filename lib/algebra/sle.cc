#include "sle.h"

// FIXME: remove NOLINT
// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
auto solve_sle(const Matrix& A, const Matrix& b) -> Matrix
{
  using Eigen::Index;
  using Eigen::MatrixXd;
  using Eigen::VectorXd;

  MatrixXd eigen_A(A.rows(), A.cols());
  VectorXd eigen_b(b.rows());

  for (std::size_t row = 0; row < A.rows(); ++row)
  {
    for (std::size_t col = 0; col < A.cols(); ++col)
    {
      eigen_A(Index(row), Index(col)) = static_cast<double>(A.at(row, col));
    }

    eigen_b(Index(row)) = static_cast<double>(b.at(row));
  }

  const VectorXd eigen_result = eigen_A.colPivHouseholderQr().solve(eigen_b);
  assert(eigen_b.isApprox(eigen_A * eigen_result));

  Matrix result(b.rows(), 1);

  for (std::size_t row = 0; row < A.rows(); ++row)
  {
    result.at(row) = static_cast<double>(eigen_result(Index(row)));
  }

  return result;
}

auto gauss_down(Matrix& A, Matrix& b) -> void;
auto gauss_up(Matrix& A, Matrix& b) -> void;

[[nodiscard]] auto gauss(
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
auto gauss_down(Matrix& A, Matrix& b) -> void
{
  assert(A.is_square());
  const auto size = A.rows();

  for (std::size_t index = 0; index < size; ++index)
  {
    // 1. Find max element in column
    std::size_t row_to_swap = index;
    auto max_element = A.at(index, index);
    for (std::size_t row = index + 1; row < size; ++row)
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

    for (std::size_t col = index; col < size; ++col)
    {
      A.at(index, col) /= divider;
    }

    b.at(index) /= divider;

    // 4. Subtract that row from rows below
    for (std::size_t row = index + 1; row < size; ++row)
    {
      auto factor = A.at(row, index);

      for (std::size_t col = index; col < size; ++col)
      {
        A.at(row, col) -= factor * A.at(index, col);
      }

      b.at(row) -= factor * b.at(index);
    }
  }
}

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters,-warnings-as-errors)
auto gauss_up(Matrix& A, Matrix& b) -> void
{
  assert(A.is_square());
  const auto size = A.rows();

  for (std::size_t index = size - 1; index != 0; --index)
  {
    for (std::size_t row = index; row != 0; --row)
    {
      auto factor = A.at(row - 1, index);

      for (std::size_t col = index; col < size; ++col)
      {
        A.at(row - 1, col) -= factor * A.at(index, col);
      }

      b.at(row - 1) -= factor * b.at(index);
    }
  }
}

// FIXME: remove NOLINT
// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
Matrix seidel(const Matrix& A, const Matrix& b, const double eps)
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

  double norm = c.norm1();
  double factor = std::abs(norm / (1 - norm));

  do
  {
    prev = x;
    x = c * x + d;
    diff = x - prev;
  } while (factor * diff.norm1() > eps);

  return x;
}
