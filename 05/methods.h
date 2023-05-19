#pragma once

#include "config.h"
#include "functions.h"
#include <algebra/algebra.h>
#include <core/types.h>

#include <cmath>
#include <iostream>

class IODESolver
{
 public:
  virtual ~IODESolver() = default;

  [[nodiscard]] virtual Matrix operator()(Config& config) const = 0;
};

class IRungeKutta : public IODESolver
{
 public:
  enum class Type
  {
    FixedStep,
    AutoStep,
  };

  constexpr explicit IRungeKutta(const Type type_ = Type::FixedStep)
      : type(type_)
  {
  }

  [[nodiscard]] Matrix operator()(Config& config) const final
  {
    auto& f = config.f;

    auto x0 = config.x0;
    auto y0 = config.y0;

    auto xk = config.xk;

    auto eps = config.eps;

    auto x = x0;
    auto y = y0;
    auto y_halved = y0;

    auto f_value = (*f)(x0, y0);
    bool should_calculate_f = false;

    auto h = calculate_starting_step(f_value, f, x0, y0, xk, eps);
    auto n = std::floor((xk - x0) / h);
    h = (xk - x0) / n;

    std::cout << "h: " << h << " (approx. # of steps: " << n << ")\n";

    if (type == Type::FixedStep)
    {
      config.h = h;
    }

    while (x < xk)
    {
      if (should_calculate_f)
      {
        f_value = (*f)(x, y);
      }

      std::cout << "\nx: " << x << ", xk: " << xk << " (h: " << h << ")\n";

      auto next_y = make_step(f, x, y, h, f_value);

      switch (type)
      {
        case Type::FixedStep:
        {
          y_halved = make_step(f, x, y_halved, h / 2);
          y_halved = make_step(f, x + h / 2, y_halved, h / 2);

          const auto total_error =
              (y_halved - next_y).norm2() / (1 - 1. / (1 << s()));

          std::cout << "Оценка полной погрешности: " << total_error << '\n';
          std::cout << "y: " << y.T() << "y_halved: " << y_halved.T() << '\n';

          y = next_y;
          x += h;

          should_calculate_f = true;
        }
        break;

        case Type::AutoStep:
        {
          const auto bottom = eps / (1 << (s() + 1));
          const auto middle = eps;
          const auto top = eps * (1 << s());

          auto half_step_y = make_step(f, x, y, h / 2, f_value);
          auto full_step_y = make_step(f, x + h / 2, half_step_y, h / 2);

          auto local_error =
              (full_step_y - next_y).norm2() / (1 - 1. / (1 << s()));

          std::cout << "Оценка локальной погрешности: " << local_error
                    << "\nbottom: " << bottom << ", middle: " << middle
                    << ", top: " << top << '\n';

          while (local_error > top)
          {
            std::cout << "CASE 1: decreasing h & recalculating (before: " << h
                      << ", after: " << h / 2 << ")\n";
            h /= 2;

            next_y = half_step_y;

            half_step_y = make_step(f, x, y, h / 2, f_value);
            full_step_y = make_step(f, x + h / 2, half_step_y, h / 2);

            local_error =
                (full_step_y - next_y).norm2() / (1 - 1. / (1 << s()));

            std::cout << "\nОценка локальной погрешности: " << local_error
                      << '\n';
          }

          if (middle < local_error && local_error <= top)
          {
            std::cout << "CASE 2: decreasing h (before: " << h
                      << ", after: " << h / 2 << ")\n";
            y = full_step_y;
            x += h;
            h /= 2;
            should_calculate_f = true;
          }
          else if (bottom <= local_error && local_error <= middle)
          {
            std::cout << "CASE 3: doing nothing\n";
            y = next_y;
            x += h;
            // h = h;
            should_calculate_f = true;
          }
          else  // (local_error < bottom)
          {
            std::cout << "CASE 4: increasing h (before: " << h
                      << ", after: " << h * 2 << ")\n";
            y = next_y;
            x += h;
            h *= 2;
            should_calculate_f = true;
          }

          if (x < xk && x + h > xk)
          {
            std::cout << "x + h > xk, clamping h (before: " << h
                      << ", after: " << (xk - x) << ")\n";
            h = xk - x;
          }

          std::cout << "number of calls: " << f->get_number_of_calls() << '\n';
          std::cout << "total number of calls: "
                    << f->get_total_number_of_calls() << '\n';

          f->reset_counter();
        }
        break;
      }
    }

    std::cout << "x: " << x << ", |xk - x| = " << std::abs(xk - x) << '\n';

    return y;
  }

 protected:
  [[nodiscard]] virtual auto make_step(const FunctionPtr& f,
                                       double x,
                                       const Matrix& y,
                                       double h,
                                       const Matrix& f_value) const
      -> Matrix = 0;

  [[nodiscard]] virtual auto m() const -> u32 = 0;
  [[nodiscard]] virtual auto s() const -> u32 = 0;

 private:
  Type type;

  [[nodiscard]] double calculate_starting_step(
      const Matrix& f0,
      const FunctionPtr& f,
      const double x0,
      const Matrix& y0,
      // FIXME: remove NOLINT
      // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
      const double xk,
      const double eps) const
  {
    auto f_value = f0;

    auto calculate_h = [&]()
    {
      const auto delta =
          std::pow(1. / std::max(std::abs(x0), std::abs(xk)), s() + 1) +
          std::pow(f_value.norm2(), s() + 1);

      return std::pow(eps / delta, 1. / (s() + 1));
    };

    auto h = calculate_h();

    const auto threshold = 1e-3;
    if (f_value.norm2() <= threshold)
    {
      auto next_y = y0 + h * f_value;
      f_value = (*f)(x0 + h, next_y);
      auto next_h = calculate_h();
      h = std::min(h, next_h);
    }

    return h;
  }

  [[nodiscard]] auto make_step(const FunctionPtr& f,
                               const double x,
                               const Matrix& y,
                               const double h) const -> Matrix
  {
    auto f_value = (*f)(x, y);
    return make_step(f, x, y, h, f_value);
  }
};

class RungeKutta2 : public IRungeKutta
{
 public:
  constexpr explicit RungeKutta2(const double c2_,
                                 const Type type_ = Type::FixedStep)
      : IRungeKutta(type_), c2(c2_)
  {
  }

 protected:
  [[nodiscard]] auto make_step(const FunctionPtr& f,
                               const double x,
                               const Matrix& y,
                               const double h,
                               const Matrix& f_value) const -> Matrix override
  {
    const double a21 = c2;
    const double b2 = 1. / (2. * c2);
    const double b1 = 1 - b2;

    const auto k1 = h * f_value;
    const auto k2 = h * (*f)(x + c2 * h, y + a21 * k1);

    return y + b1 * k1 + b2 * k2;
  }

  [[nodiscard]] auto m() const -> u32 override { return 2; }
  [[nodiscard]] auto s() const -> u32 override { return 2; }

 private:
  double c2;
};

class RungeKutta4 : public IRungeKutta
{
 protected:
  [[nodiscard]] auto make_step(const FunctionPtr& f,
                               const double x,
                               const Matrix& y,
                               const double h,
                               const Matrix& f_value) const -> Matrix override
  {
    const auto k1 = h * f_value;
    const auto k2 = h * (*f)(x + 0.5 * h, y + 0.5 * k1);
    const auto k3 = h * (*f)(x + 0.5 * h, y + 0.5 * k2);
    const auto k4 = h * (*f)(x + h, y + k3);

    return y + (k1 + 2. * k2 + 2. * k3 + k4) / 6.;
  }

  [[nodiscard]] auto m() const -> u32 override { return 4; }
  [[nodiscard]] auto s() const -> u32 override { return 4; }
};
