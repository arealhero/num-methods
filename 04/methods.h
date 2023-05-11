#pragma once

#include "config.h"
#include <algebra.h>
#include <types.h>

#include <Eigen/Dense>
#include <array>
#include <iostream>
#include <string>

struct Interval
{
  value_t left;
  value_t right;

  [[nodiscard]] constexpr value_t length() const { return right - left; }
  [[nodiscard]] constexpr value_t middle() const { return left + length() / 2; }
};

class IIntegrator
{
 public:
  virtual ~IIntegrator() = default;

  [[nodiscard]] virtual constexpr std::string get_name() const = 0;
  [[nodiscard]] virtual constexpr value_t operator()(
      func_t x,
      const Interval& interval) const = 0;
};

class LeftRect : public IIntegrator
{
 public:
  [[nodiscard]] constexpr std::string get_name() const override
  {
    return "Левый прямоугольник";
  }

  [[nodiscard]] constexpr value_t operator()(
      const func_t f,
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

  [[nodiscard]] constexpr value_t operator()(
      const func_t f,
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

  [[nodiscard]] constexpr value_t operator()(
      const func_t f,
      const Interval& interval) const override
  {
    return interval.length() / 2.L * (f(interval.left) + f(interval.right));
  }
};

class Simpson : public IIntegrator
{
 public:
  [[nodiscard]] constexpr std::string get_name() const override
  {
    return "Симпсон";
  }

  [[nodiscard]] constexpr value_t operator()(
      const func_t f,
      const Interval& interval) const override
  {
    return interval.length() / 6.L *
           (f(interval.left) + 4 * f(interval.middle()) + f(interval.right));
  }
};

class NewtonCotes : public IIntegrator
{
 public:
  explicit constexpr NewtonCotes(const std::size_t n = 3) : m_n(n) {}

  [[nodiscard]] constexpr std::string get_name() const override
  {
    return "Ньютон-Котс";
  }

  [[nodiscard]] value_t operator()(const func_t f,
                                   const Interval& interval) const override
  {
    Matrix mu(m_n, 1);
    const auto a = interval.left;
    const auto b = interval.right;

    for (std::size_t i = 0; i < m_n; ++i)
    {
      const auto power = i - ALPHA + 1;
      mu.at(i) =
          (std::pow(b - ORIG_A, power) - std::pow(a - ORIG_A, power)) / power;
    }

    Matrix t(m_n, 1);
    const auto step = interval.length() / (m_n - 1);
    for (std::size_t i = 0; i < m_n; ++i)
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

    const auto coeffs = gauss(ts, mu);
    std::cout << get_name() << ": coeffs =";
    for (std::size_t i = 0; i < m_n; ++i)
    {
      std::cout << ' ' << coeffs.at(i);
    }
    std::cout << '\n';

    Matrix values(m_n, 1);
    for (std::size_t i = 0; i < m_n; ++i)
    {
      values.at(i) = f(t.at(i) + ORIG_A);
    }

    auto result = coeffs.T() * values;
    return result;
  }

  [[nodiscard]] std::pair<value_t, Matrix> run_with_coeffs(
      const func_t f,
      const Interval& interval) const
  {
    Matrix mu(m_n, 1);
    const auto a = interval.left;
    const auto b = interval.right;

    for (std::size_t i = 0; i < m_n; ++i)
    {
      const auto power = i - ALPHA + 1;
      mu.at(i) =
          (std::pow(b - ORIG_A, power) - std::pow(a - ORIG_A, power)) / power;
    }

    Matrix t(m_n, 1);
    const auto step = interval.length() / (m_n - 1);
    for (std::size_t i = 0; i < m_n; ++i)
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

    auto coeffs = seidel(ts, mu, 1e-8L);
    Matrix values(m_n, 1);
    for (std::size_t i = 0; i < m_n; ++i)
    {
      values.at(i) = f(t.at(i) + ORIG_A);
    }

    auto result = coeffs.T() * values;
    return {result, coeffs};
  }

 private:
  std::size_t m_n;
};

class Gauss : public IIntegrator
{
 public:
  [[nodiscard]] constexpr std::string get_name() const override
  {
    return "Гаусс";
  }

  [[nodiscard]] value_t operator()(const func_t f,
                                   const Interval& interval) const override
  {
    constexpr std::size_t n = 3;
    const auto a = interval.left;
    const auto b = interval.right;

    auto mu = Matrix(2 * n, 1);
    for (std::size_t i = 0; i < 2 * n; ++i)
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

    const auto solution = seidel(A, B, 1e-8L);
    auto polynomial_coefficients = std::array<value_t, n + 1>{1.L};
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

    const auto coeffs = seidel(A, B, 1e-8L);
    Matrix values(n, 1);
    for (std::size_t i = 0; i < n; ++i)
    {
      values.at(i) = f(roots.at(i) + ORIG_A);
    }

    return coeffs.T() * values;
  }
};
