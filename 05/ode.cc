#include "config.h"
#include "methods.h"
#include <algebra/matrix.h>
#include <algebra/sle.h>
#include <core/utils.h>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <vector>

static constexpr double pi = std::numbers::pi_v<double>;

[[nodiscard]] static auto exact_solution(double x, double x0, const Matrix& y0)
    -> Matrix
{
  const auto A = A_VALUE;
  const auto B = B_VALUE;
  const auto sqrtAB = std::sqrt(A * B);

  auto A_coeffs = Matrix(2, 2);
  A_coeffs.at(0, 0) = sqrtAB / B * std::sin(sqrtAB * x0);
  A_coeffs.at(0, 1) = -sqrtAB / B * std::cos(sqrtAB * x0);
  A_coeffs.at(1, 0) = std::cos(sqrtAB * x0);
  A_coeffs.at(1, 1) = std::sin(sqrtAB * x0);

  auto C = gauss(A_coeffs, y0);

  Matrix y(2, 1);

  const auto arg = sqrtAB * x;
  y.at(0) = C.at(0) * sqrtAB * std::sin(arg) / B -
            C.at(1) * sqrtAB * std::cos(arg) / B;
  y.at(1) = C.at(0) * std::cos(arg) + C.at(1) * std::sin(arg);

  return y;
}

void print_usage();

void first_part();
void second_part();
void third_part_fixed();
void third_part_automatic();

auto main(int argc, char* const* argv) -> int
{
  if (argc == 1)
  {
    print_usage();
    return EXIT_FAILURE;
  }

  int c;
  // NOLINTNEXTLINE(concurrency-mt-unsafe)
  while ((c = getopt(argc, argv, "123h")) != -1)
  {
    switch (c)
    {
      case '1':
        first_part();
        break;

      case '2':
        second_part();
        break;

      case '3':
        third_part_fixed();
        third_part_automatic();
        break;

      case 'h':
      case '?':
        print_usage();
        return EXIT_FAILURE;
    }
  }

  return EXIT_SUCCESS;
}

void print_usage()
{
  std::cout << "usage: ode [-1|-2|-3|-h]\n"
               "\t-1\trun first part\n"
               "\t-2\trun second part\n"
               "\t-3\trun third part\n"
               "\t-h\tprint this message\n";
}

void first_part()
{
  auto config = Config{1e-4};
  config.logger = make<DummyLogger>();

  auto method = make<RungeKutta2>(xi);
  const auto result = (*method)(config);

  const auto exact = exact_solution(config.xk, config.x0, config.y0);
  const auto error = (exact - result).norm2();

  std::cout << "\nexpected:\n"
            << exact << '\n'
            << "got:\n"
            << result << '\n'
            << "error: " << error << '\n';
}

void second_part()
{
  auto config = Config{1e-5};
  config.logger = make<DummyLogger>();

  auto method = make<RungeKutta2>(xi, IRungeKutta::Type::AutoStep);
  const auto result = (*method)(config);

  const auto exact = exact_solution(config.xk, config.x0, config.y0);
  const auto error = (exact - result).norm2();

  std::cout << "\nexpected:\n"
            << exact << '\n'
            << "got:\n"
            << result << '\n'
            << "error: " << error << '\n';

  std::cout << "total number of calls: "
            << config.f->get_total_number_of_calls() << '\n';
}

void third_part_fixed()
{
  static constexpr auto eps = 1e-4;

  auto rk2_config = Config{eps};
  rk2_config.logger = make<VectorLogger>();

  auto rk2 = make<RungeKutta2>(xi);
  auto rk2_y = (*rk2)(rk2_config);

  {
    std::ofstream fout("data/rk2-total-error.csv");

    fout << "x,RK2 error\n";
    auto points = rk2_config.logger->get_points();
    auto rk2_total_errors = rk2_config.logger->get_total_errors();
    for (auto&& [x, y] : points)
    {
      const auto error =
          (exact_solution(x, rk2_config.x0, rk2_config.y0) - y).norm2();
      fout << x << ',' << error << '\n';
    }

    fout.close();
  }

  // ------------------------------------------------------------------------ //

  auto rk4_config = Config{eps};
  rk4_config.logger = make<VectorLogger>();

  auto rk4 = make<RungeKutta4>();
  auto rk4_y = (*rk4)(rk4_config);

  {
    std::ofstream fout("data/rk4-total-error.csv");

    fout << "x,RK4 error\n";
    auto points = rk4_config.logger->get_points();
    auto rk4_total_errors = rk4_config.logger->get_total_errors();
    for (auto&& [x, y] : points)
    {
      const auto error =
          (exact_solution(x, rk4_config.x0, rk4_config.y0) - y).norm2();
      fout << x << ',' << error << '\n';
    }

    fout.close();
  }
}

void third_part_automatic()
{
  struct CallsInfo
  {
    // FIXME: remove NOLINT
    // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
    CallsInfo(double eps_, double calls_) : eps(eps_), calls(calls_) {}

    double eps;
    double calls;
  };

  std::vector<CallsInfo> rk2_calls{};
  std::vector<CallsInfo> rk4_calls{};

  for (const double eps : {1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-8, 1e-12})
  {
    auto rk2_config = Config{eps};
    rk2_config.logger = make<VectorLogger>();

    auto rk2 = make<RungeKutta2>(xi, IRungeKutta::Type::AutoStep);
    auto rk2_y = (*rk2)(rk2_config);

    rk2_calls.emplace_back(eps, rk2_config.f->get_total_number_of_calls());

    {
      std::ofstream step_fout("data/rk2-auto-step-" + std::to_string(eps) +
                              ".csv");
      step_fout << "x,RK2 h\n";

      auto rk2_step_values = rk2_config.logger->get_step_values();
      for (auto&& [x, h] : rk2_step_values)
      {
        step_fout << x << ',' << h << '\n';
      }

      step_fout.close();
    }

    {
      std::ofstream local_error_fout("data/rk2-local-error-" +
                                     std::to_string(eps) + ".csv");
      local_error_fout << "x,one,RK2 local error proportion\n";

      auto rk2_local_errors = rk2_config.logger->get_local_errors();
      auto rk2_points = rk2_config.logger->get_points();

      double x_prev = rk2_config.x0;
      Matrix y_prev = rk2_config.y0;
      for (u64 i = 0; i < rk2_points.size(); ++i)
      {
        auto&& [x, y] = rk2_points.at(i);
        auto&& [_, local_error] = rk2_local_errors.at(i);

        auto exact_y = exact_solution(x, x_prev, y_prev);
        auto true_local_error = (exact_y - y).norm2();

        local_error_fout << x << ',' << 1. << ','
                         << (true_local_error / local_error) << '\n';

        x_prev = x;
        y_prev = y;
      }

      local_error_fout.close();
    }

    // ---------------------------------------------------------------------- //

    auto rk4_config = Config{eps};
    rk4_config.logger = make<VectorLogger>();

    auto rk4 = make<RungeKutta4>(IRungeKutta::Type::AutoStep);
    auto rk4_y = (*rk4)(rk4_config);

    rk4_calls.emplace_back(eps, rk4_config.f->get_total_number_of_calls());

    {
      std::ofstream step_fout("data/rk4-auto-step-" + std::to_string(eps) +
                              ".csv");

      step_fout << "x,RK4 h\n";

      auto rk4_step_values = rk4_config.logger->get_step_values();
      for (auto&& [x, h] : rk4_step_values)
      {
        step_fout << x << ',' << h << '\n';
      }

      step_fout.close();
    }

    {
      std::ofstream local_error_fout("data/rk4-local-error-" +
                                     std::to_string(eps) + ".csv");
      local_error_fout << "x,one,RK4 local error proportion\n";

      auto rk4_local_errors = rk4_config.logger->get_local_errors();
      auto rk4_points = rk4_config.logger->get_points();

      double x_prev = rk4_config.x0;
      Matrix y_prev = rk4_config.y0;
      for (u64 i = 0; i < rk4_points.size(); ++i)
      {
        auto&& [x, y] = rk4_points.at(i);
        auto&& [_, local_error] = rk4_local_errors.at(i);

        auto exact_y = exact_solution(x, x_prev, y_prev);
        auto true_local_error = (exact_y - y).norm2();

        local_error_fout << x << ',' << 1. << ','
                         << (true_local_error / local_error) << '\n';

        x_prev = x;
        y_prev = y;
      }

      local_error_fout.close();
    }
  }

  {
    std::ofstream fout("data/function-calls.csv");

    fout << "eps,RK2 calls,RK4 calls\n";
    for (std::size_t i = 0; i < rk2_calls.size(); ++i)
    {
      auto&& [eps, calls2] = rk2_calls.at(i);
      auto&& [_, calls4] = rk4_calls.at(i);
      fout << eps << ',' << calls2 << ',' << calls4 << '\n';
    }

    fout.close();
  }
}
