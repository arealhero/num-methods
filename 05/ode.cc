#include "config.h"
#include "methods.h"
#include <algebra/matrix.h>
#include <core/utils.h>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <vector>

static constexpr double pi = std::numbers::pi_v<double>;

[[nodiscard]] static auto exact_solution(double x) -> Matrix
{
  Matrix result(2, 1);

  const auto A = A_VALUE;
  const auto B = B_VALUE;
  const auto sqrtAB = std::sqrt(A * B);
  const auto arg = sqrtAB * x;

  result.at(0) = A * pi * sqrtAB * std::sin(arg) / B + B * pi * std::cos(arg);
  result.at(1) = A * pi * std::cos(arg) - B * B * pi * std::sin(arg) / sqrtAB;

  return result;
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

  const auto exact = exact_solution(config.xk);
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

  const auto exact = exact_solution(config.xk);
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
    auto rk2_total_errors = rk2_config.logger->get_total_errors();
    for (auto&& [x, error] : rk2_total_errors)
    {
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
    auto rk4_total_errors = rk4_config.logger->get_total_errors();
    for (auto&& [x, error] : rk4_total_errors)
    {
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

  for (const double eps : {1e-1, 1e-2, 1e-3, 1e-4, 1e-5})
  {
    auto rk2_config = Config{eps};
    rk2_config.logger = make<VectorLogger>();

    auto rk2 = make<RungeKutta2>(xi, IRungeKutta::Type::AutoStep);
    auto rk2_y = (*rk2)(rk2_config);

    rk2_calls.emplace_back(eps, rk2_config.f->get_total_number_of_calls());

    {
      std::ofstream step_fout("data/rk2-auto-step-" + std::to_string(eps) +
                              ".csv");
      std::ofstream error_fout("data/rk2-local-error-" + std::to_string(eps) +
                               ".csv");

      step_fout << "x,RK2 h\n";
      error_fout << "x,RK2 error\n";

      auto rk2_step_values = rk2_config.logger->get_step_values();
      auto rk2_local_errors = rk2_config.logger->get_local_errors();

      for (u64 i = 0; i < rk2_step_values.size(); ++i)
      {
        auto&& [x, h] = rk2_step_values.at(i);
        step_fout << x << ',' << h << '\n';
        std::cout << x << ',' << h << '\n';

        auto&& [_, error] = rk2_local_errors.at(i);

        error_fout << x << ',' << (error * h * eps) << '\n';
      }

      step_fout.close();
      error_fout.close();
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
      std::ofstream error_fout("data/rk4-local-error-" + std::to_string(eps) +
                               ".csv");

      step_fout << "x,RK4 h\n";
      error_fout << "x,RK4 error\n";

      auto rk4_step_values = rk4_config.logger->get_step_values();
      auto rk4_local_errors = rk4_config.logger->get_local_errors();

      for (u64 i = 0; i < rk4_step_values.size(); ++i)
      {
        auto&& [x, h] = rk4_step_values.at(i);
        step_fout << x << ',' << h << '\n';

        auto&& [_, error] = rk4_local_errors.at(i);

        error_fout << x << ',' << (error * h * eps) << '\n';
      }

      step_fout.close();
      error_fout.close();
    }
  }

  {
    std::ofstream fout("data/rk2-function-calls.csv");

    fout << "eps,calls\n";
    for (auto&& [eps, calls] : rk2_calls)
    {
      fout << eps << ',' << calls << '\n';
    }

    fout.close();
  }

  {
    std::ofstream fout("data/rk4-function-calls.csv");

    fout << "eps,calls\n";
    for (auto&& [eps, calls] : rk4_calls)
    {
      fout << eps << ',' << calls << '\n';
    }

    fout.close();
  }
}
