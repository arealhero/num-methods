#pragma once

#include "config.h"
#include <algebra/algebra.h>
#include <algebra/polynomial.h>
#include <algebra/sle.h>
#include <core/types.h>

#include <Eigen/Dense>
#include <array>
#include <iostream>
#include <string>

struct Interval
{
  double left;
  double right;

  [[nodiscard]] constexpr double length() const { return right - left; }
  [[nodiscard]] constexpr double middle() const { return left + length() / 2; }
};

class IIntegrator
{
 public:
  virtual ~IIntegrator() = default;

  [[nodiscard]] virtual constexpr std::string get_name() const = 0;
  [[nodiscard]] virtual constexpr double operator()(
      FunctionType x,
      const Interval& interval) const = 0;
};

class LeftRect : public IIntegrator
{
 public:
  [[nodiscard]] constexpr std::string get_name() const override
  {
    return "Левый прямоугольник";
  }

  [[nodiscard]] constexpr double operator()(
      const FunctionType f,
      const Interval& interval) const override
  {
    return interval.length() * f(interval.left);
  }
};

class MidRect : public IIntegrator
{
 public:
  [[nodiscard]] constexpr std::string get_name() const override
  {
    return "Средний прямоугольник";
  }

  [[nodiscard]] constexpr double operator()(
      const FunctionType f,
      const Interval& interval) const override
  {
    return interval.length() * f(interval.middle());
  }
};

class Trapezoid : public IIntegrator
{
 public:
  [[nodiscard]] constexpr std::string get_name() const override
  {
    return "Трапеция";
  }

  [[nodiscard]] constexpr double operator()(
      const FunctionType f,
      const Interval& interval) const override
  {
    return interval.length() / 2. * (f(interval.left) + f(interval.right));
  }
};

class Simpson : public IIntegrator
{
 public:
  [[nodiscard]] constexpr std::string get_name() const override
  {
    return "Симпсон";
  }

  [[nodiscard]] constexpr double operator()(
      const FunctionType f,
      const Interval& interval) const override
  {
    return interval.length() / 6. *
           (f(interval.left) + 4 * f(interval.middle()) + f(interval.right));
  }
};

class NewtonCotes : public IIntegrator
{
 public:
  explicit constexpr NewtonCotes(const u32 n = 3) : m_n(n)
  {
    assert(m_n > 1 && "n should be greater than 1");
  }

  [[nodiscard]] constexpr std::string get_name() const override
  {
    return "Ньютон-Котс";
  }

  [[nodiscard]] double operator()(const FunctionType f,
                                  const Interval& interval) const override
  {
    Matrix mu(m_n, 1);
    const auto a = interval.left;
    const auto b = interval.right;

    for (u32 i = 0; i < m_n; ++i)
    {
      const auto power = i - ALPHA + 1;
      mu.at(i) =
          (std::pow(b - ORIG_A, power) - std::pow(a - ORIG_A, power)) / power;
    }

    Matrix t(m_n, 1);
    const auto step = interval.length() / (m_n - 1);
    for (u32 i = 0; i < m_n; ++i)
    {
      t.at(i) = a - ORIG_A + i * step;
    }

    Matrix ts(m_n, m_n);
    for (std::size_t row = 0; row < m_n; ++row)
    {
      for (std::size_t col = 0; col < m_n; ++col)
      {
        ts.at(row, col) = std::pow(t.at(col), row);
      }
    }

    const auto coeffs = solve_sle(ts, mu);

    Matrix values(m_n, 1);
    for (std::size_t i = 0; i < m_n; ++i)
    {
      values.at(i) = f(t.at(i) + ORIG_A);
    }

    auto result = coeffs.T() * values;
    return result;
  }

  [[nodiscard]] std::pair<double, Matrix> run_with_coeffs(
      const FunctionType f,
      const Interval& interval) const
  {
    Matrix mu(m_n, 1);
    const auto a = interval.left;
    const auto b = interval.right;

    for (u32 i = 0; i < m_n; ++i)
    {
      const auto power = i - ALPHA + 1;
      mu.at(i) =
          (std::pow(b - ORIG_A, power) - std::pow(a - ORIG_A, power)) / power;
    }

    Matrix t(m_n, 1);
    const auto step = interval.length() / (m_n - 1);
    for (u32 i = 0; i < m_n; ++i)
    {
      t.at(i) = a - ORIG_A + i * step;
    }

    Matrix ts(m_n, m_n);
    for (std::size_t row = 0; row < m_n; ++row)
    {
      for (std::size_t col = 0; col < m_n; ++col)
      {
        ts.at(row, col) = std::pow(t.at(col), row);
      }
    }

    auto coeffs = solve_sle(ts, mu);
    Matrix values(m_n, 1);
    for (std::size_t i = 0; i < m_n; ++i)
    {
      values.at(i) = f(t.at(i) + ORIG_A);
    }

    auto result = coeffs.T() * values;
    return {result, coeffs};
  }

 private:
  u32 m_n;
};

class Gauss : public IIntegrator
{
 public:
  [[nodiscard]] constexpr std::string get_name() const override
  {
    return "Гаусс";
  }

  [[nodiscard]] double operator()(const FunctionType f,
                                  const Interval& interval) const override
  {
    constexpr std::size_t n = 3;
    const auto a = interval.left;
    const auto b = interval.right;

    auto mu = Matrix(2 * n, 1);
    for (u32 i = 0; i < 2 * n; ++i)
    {
      const auto power = i - ALPHA + 1;
      mu.at(i) =
          (std::pow(b - ORIG_A, power) - std::pow(a - ORIG_A, power)) / power;
    }

    // TODO: come up with better names
    auto A = Matrix(n, n);
    auto B = Matrix(n, 1);
    for (std::size_t row = 0; row < n; ++row)
    {
      for (std::size_t col = 0; col < n; ++col)
      {
        A.at(row, col) = mu.at(row + col);
      }

      B.at(row) = -mu.at(n + row);
    }

    const auto solution = solve_sle(A, B);
    auto polynomial_coefficients = std::array<double, n + 1>{1.L};
    for (std::size_t i = 0; i < n; ++i)
    {
      polynomial_coefficients.at(i + 1) = solution.at(n - i - 1);
    }

    const auto roots = cardano(polynomial_coefficients);

    for (std::size_t row = 0; row < n; ++row)
    {
      for (std::size_t col = 0; col < n; ++col)
      {
        A.at(row, col) = std::pow(roots.at(col), row);
      }

      B.at(row) = mu.at(row);
    }

    const auto coeffs = solve_sle(A, B);
    Matrix values(n, 1);
    for (std::size_t i = 0; i < n; ++i)
    {
      values.at(i) = f(roots.at(i) + ORIG_A);
    }

    return coeffs.T() * values;
  }
};
