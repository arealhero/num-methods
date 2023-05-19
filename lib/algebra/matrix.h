#pragma once

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wnull-dereference"
#include <Eigen/Dense>
#pragma GCC diagnostic pop

#include <core/utils.h>

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
  Matrix(u64 rows, u64 cols, double value = 0);
  Matrix(u64 rows, u64 cols, const std::initializer_list<double>& numbers);
  Matrix(Eigen::MatrixXd matrix);

  [[nodiscard]] static auto zero(u64 size) -> Matrix;
  [[nodiscard]] static auto zero(u64 rows, u64 cols) -> Matrix;

  [[nodiscard]] static auto identity(u64 size) -> Matrix;
  [[nodiscard]] static auto e(u64 size, u64 i) -> Matrix;

  [[nodiscard]] auto rows() const -> u64;
  [[nodiscard]] auto cols() const -> u64;
  [[nodiscard]] auto is_square() const -> bool;

  [[nodiscard]] auto at(u64 row, u64 col) -> double&;
  [[nodiscard]] auto at(u64 row, u64 col) const -> const double&;

  [[nodiscard]] auto at(u64 index) -> double&;
  [[nodiscard]] auto at(u64 index) const -> const double&;

  operator double() const;

  auto operator+=(double number) -> Matrix&;
  auto operator-=(double number) -> Matrix&;
  auto operator*=(double number) -> Matrix&;
  auto operator/=(double number) -> Matrix&;

  auto operator+=(const Matrix& other) -> Matrix&;
  auto operator-=(const Matrix& other) -> Matrix&;

  [[nodiscard]] auto operator-() const -> Matrix;

  [[nodiscard]] auto norm1() const -> double;
  [[nodiscard]] auto norm2() const -> double;

  [[nodiscard]] auto T() const -> Matrix;
  [[nodiscard]] auto inv() const -> Matrix;

  void swap_rows(u64 i, u64 j);

  [[nodiscard]] auto get_diagonal_dominance_factor() const
      -> std::optional<double>;
  [[nodiscard]] auto is_diagonally_dominant() const -> bool;

 private:
  Eigen::MatrixXd m_matrix;
};

auto operator+(const Matrix& lhs, const Matrix& rhs) -> Matrix;
auto operator-(const Matrix& lhs, const Matrix& rhs) -> Matrix;

auto operator==(const Matrix& lhs, const Matrix& rhs) -> bool;

auto operator*(const Matrix& lhs, const Matrix& rhs) -> Matrix;
auto operator*(double number, const Matrix& mat) -> Matrix;
auto operator/(const Matrix& mat, double number) -> Matrix;

std::ostream& operator<<(std::ostream& out, const Matrix& matrix);
