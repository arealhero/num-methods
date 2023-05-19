#include "config.h"
#include "methods.h"
#include <algebra/matrix.h>
#include <core/utils.h>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>

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
void third_part();

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
        third_part();
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

void third_part()
{
  std::ofstream fout("part-3.csv");

  auto config = Config{1e-5};

  auto runge_kutta_2_fixed = make<RungeKutta2>(xi);
  auto runge_kutta_4_fixed = make<RungeKutta4>();
}
