#pragma once

#include "functions.h"
#include "loggers.h"
#include <core/utils.h>

#include <cmath>
#include <cstdint>

// Вариант 9
constexpr double xi = 1. / 11.;

constexpr double A_VALUE = 1. / 15.;
constexpr double B_VALUE = 1. / 20.;

struct Config
{
  explicit Config(const double eps_)
  {
    static constexpr double pi = std::numbers::pi_v<double>;

    Matrix S(2, 2);

    S.at(0, 1) = A_VALUE;
    S.at(1, 0) = -B_VALUE;

    f = make<LinearFunction>(S);

    x0 = 0.;
    y0.at(0) = B_VALUE * pi;
    y0.at(1) = A_VALUE * pi;

    xk = pi;

    eps = eps_;
  }

  FunctionPtr f{};
  LoggerPtr logger{};

  double x0{};
  Matrix y0 = Matrix(2, 1);

  double xk{};

  double eps{};
};
