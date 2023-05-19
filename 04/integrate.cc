#include "config.h"
#include "methods.h"
#include <core/types.h>
#include <core/utils.h>
#include <unistd.h>

#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <string>
#include <utility>
#include <vector>

using IntegratorPair = std::pair<Ptr<IIntegrator>, Ptr<IIntegrator>>;

struct Run
{
  // FIXME: remove NOLINT
  // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
  constexpr Run(const std::size_t power_, const double h_, const double S_)
      : power(power_), h(h_), S(S_)
  {
  }

  std::size_t power;
  double h;
  double S;
  double error = 0;
};

auto print_usage() -> void;

auto first_part() -> void;
auto second_part() -> void;

auto penalty_part_1() -> void;
auto penalty_part_2() -> void;
auto penalty_part_3() -> void;
auto penalty_part_4() -> void;

constexpr auto run_method(FunctionType f, const Ptr<IIntegrator>& method, u32 n)
    -> double;
constexpr auto run_methods_interleaved(FunctionType f,
                                       const IntegratorPair& methods,
                                       u32 n) -> double;

auto calculate_optimum_step(FunctionType f,
                            const Ptr<IIntegrator>& method,
                            double eps) -> double;

auto run_method_with_eps(FunctionType f,
                         const Ptr<IIntegrator>& method,
                         double eps,
                         bool choose_optimum_step = false) -> Run;

auto run_method_with_step(FunctionType f,
                          const Ptr<IIntegrator>& method,
                          double step) -> double;
constexpr auto run_methods_with_step_interleaved(FunctionType f,
                                                 const IntegratorPair& methods,
                                                 double step) -> double;

constexpr auto aitken(double S_1, double S_2, double S_3, double L) -> double;
constexpr auto aitken(const std::vector<Run>& runs, double L) -> double;

auto main(int argc, char* const* argv) -> int
{
  if (argc == 1)
  {
    print_usage();
    return EXIT_FAILURE;
  }

  int c;
  while ((c = getopt(argc, argv, "fsp:h")) != -1)
  {
    switch (c)
    {
      case 'f':
        first_part();
        break;

      case 's':
        second_part();
        break;

      case 'p':
      {
        const auto opt = std::string{optarg};

        if (opt == std::string{"1"})
        {
          penalty_part_1();
        }
        else if (opt == std::string{"2"})
        {
          penalty_part_2();
        }
        else if (opt == std::string{"3"})
        {
          penalty_part_3();
        }
        else if (opt == std::string{"4"})
        {
          penalty_part_4();
        }
        else if (opt == std::string{"all"})
        {
          penalty_part_1();
          std::cout << "\n\n";
          penalty_part_2();
          std::cout << "\n\n";
          penalty_part_3();
          std::cout << "\n\n";
          penalty_part_4();
        }
        else
        {
          print_usage();
          return EXIT_FAILURE;
        }
      }
      break;

      case 'h':
      case '?':
        print_usage();
        return EXIT_FAILURE;
    }
  }
}

auto print_usage() -> void
{
  std::cout << "usage: integrate [-f|-s|-p <num>|-h]\n"
               "\t-f\trun first part\n"
               "\t-s\trun second part\n"
               "\t-p <num>\trun penalty part <num>\n"
               "\t-h\tprint this message\n";
}

auto first_part() -> void
{
  const std::vector<Ptr<IIntegrator>> methods = {
      make<LeftRect>(),
      make<MidRect>(),
      make<Trapezoid>(),
      make<Simpson>(),
  };

  const std::vector<Ptr<IIntegrator>> weighted_methods = {
      make<NewtonCotes>(),
      make<Gauss>(),
  };

  std::ofstream fout("first_part.csv");

  fout << "N";
  for (const auto& method : methods)
  {
    fout << ',' << method->get_name();
  }
  for (const auto& method : weighted_methods)
  {
    fout << ',' << method->get_name();
  }
  fout << '\n';

  for (u32 n = 2; n < MAX_PARTITIONS; ++n)
  {
    fout << std::fixed << n << std::scientific << std::setprecision(8);
    std::cout << n << '\r' << std::flush;
    for (const auto& method : methods)
    {
      const auto sum = run_method(f, method, n);
      fout << ',' << std::abs(EXACT_VALUE - sum);
    }

    for (const auto& method : weighted_methods)
    {
      const auto sum = run_method(f, method, n);
      fout << ',' << std::abs(EXACT_VALUE_WITH_P - sum);
    }
    fout << '\n';
  }

  fout.close();
  std::cout << "Файл 'first_part.csv' создан.\n";
}

auto second_part() -> void
{
  constexpr double eps = 1e-6;

  const std::vector<Ptr<IIntegrator>> weighted_methods = {
      make<NewtonCotes>(),
      make<Gauss>(),
  };

  for (const auto& method : weighted_methods)
  {
    std::cout << " --- " << method->get_name() << " ---\n";
    std::cout << "Запуск без вычисления оптимального шага\n";
    const auto run = run_method_with_eps(f, method, eps);
    std::cout << "\nИтог:\n - h = " << run.h << "\n - S = " << run.S
              << "\n - error = " << run.error
              << "\n - abs_error = " << std::abs(EXACT_VALUE_WITH_P - run.S)
              << "\n - abs_error (S+error) = "
              << std::abs(EXACT_VALUE_WITH_P - run.S - run.error) << "\n\n";

    std::cout << "Запуск с вычислением оптимального шага\n";
    const auto optimum_run = run_method_with_eps(f, method, eps, true);
    std::cout << "\nИтог:\n - h = " << optimum_run.h
              << "\n - S = " << optimum_run.S
              << "\n - error = " << optimum_run.error << "\n - abs_error = "
              << std::abs(EXACT_VALUE_WITH_P - optimum_run.S)
              << "\n - abs_error (S+error) = "
              << std::abs(EXACT_VALUE_WITH_P - optimum_run.S -
                          optimum_run.error)
              << "\n\n";
  }
}

auto penalty_part_1() -> void
{
  std::cout << "--- ЧАСТЬ 1 ---\n";
  std::size_t n = 2;
  bool has_negative_coeffs = false;
  while (!has_negative_coeffs)
  {
    ++n;

    const auto method = make<NewtonCotes>(n);
    const auto [value, coeffs] = method->run_with_coeffs(f, {ORIG_A, ORIG_B});

    std::cout << "coeffs:";
    for (std::size_t i = 0; i < coeffs.rows(); ++i)
    {
      std::cout << ' ' << coeffs.at(i);
    }
    std::cout << '\n';

    for (std::size_t i = 0; i < coeffs.rows(); ++i)
    {
      has_negative_coeffs |= coeffs.at(i) < 0;
    }
  }

  std::cout << "n: " << n << '\n';
}

auto penalty_part_2() -> void
{
  std::cout << "--- ЧАСТЬ 2 ---\n";
  std::ofstream fout("newton_cotes_comparison.csv");

  constexpr std::size_t max_partitions = 40;
  fout << "N,Составная,Малая\n";
  for (u32 i = 2; i < max_partitions; ++i)
  {
    const std::vector<Ptr<IIntegrator>> weighted_methods = {
        make<NewtonCotes>(),
        make<NewtonCotes>(2 * i + 1),
    };

    fout << std::fixed << i << std::scientific << std::setprecision(8);
    std::cout << i << '\r' << std::flush;

    for (const auto& method : weighted_methods)
    {
      const auto sum = run_method(f, method, i);
      fout << ',' << std::abs(EXACT_VALUE_WITH_P - sum);
    }
    fout << '\n';
  }

  fout.close();
  std::cout << "Файл 'newton_cotes_comparison.csv' создан.\n";
}

auto penalty_part_3() -> void
{
  std::cout << "--- ЧАСТЬ 3 ---\n";

  const IntegratorPair same_degree_methods = {
      make<Gauss>(),
      make<NewtonCotes>(6),
  };

  const IntegratorPair different_degree_methods = {
      make<Gauss>(),
      make<NewtonCotes>(8),
  };

  std::vector<Run> same_degree_runs;
  std::vector<Run> different_degree_runs;
  std::vector<std::optional<double>> same_degree_convergence_rate;
  std::vector<std::optional<double>> different_degree_convergence_rate;

  constexpr double L = 2;
  auto h = (ORIG_B - ORIG_A);
  for (std::size_t i = 0; i < 9; ++i)
  {
    std::cout << (i + 1) << ". h: " << h << std::endl;
    same_degree_runs.emplace_back(
        L, h, run_methods_with_step_interleaved(f, same_degree_methods, h));
    different_degree_runs.emplace_back(
        L,
        h,
        run_methods_with_step_interleaved(f, different_degree_methods, h));
    h /= L;

    if (i >= 2)
    {
      // Процесс Эйткена
      same_degree_convergence_rate.emplace_back(aitken(same_degree_runs, L));
      different_degree_convergence_rate.emplace_back(
          aitken(different_degree_runs, L));
    }
    else
    {
      same_degree_convergence_rate.emplace_back(std::nullopt);
      different_degree_convergence_rate.emplace_back(std::nullopt);
    }
  }

  std::ofstream fout("penalty_3.csv");

  fout << "L,h,Одинаковые АСТ,m,Разные АСТ,m\n";

  for (std::size_t i = 0; i < same_degree_runs.size(); ++i)
  {
    const auto same_degree_run = same_degree_runs.at(i);
    const auto same_degree_convergence = same_degree_convergence_rate.at(i);

    const auto different_degree_run = different_degree_runs.at(i);
    const auto different_degree_convergence =
        different_degree_convergence_rate.at(i);

    std::cout << std::fixed << i + 1 << std::scientific << std::setprecision(8)
              << '\t' << same_degree_run.h << '\t' << same_degree_run.S << '\t'
              << (same_degree_convergence
                      ? std::to_string(same_degree_convergence.value())
                      : "-\t")
              << '\t' << different_degree_run.S << '\t'
              << (different_degree_convergence
                      ? std::to_string(different_degree_convergence.value())
                      : "-\t")
              << '\n';

    fout << std::fixed << i + 1 << std::scientific << std::setprecision(8)
         << ',' << same_degree_run.h << ',' << same_degree_run.S << ','
         << (same_degree_convergence
                 ? std::to_string(same_degree_convergence.value())
                 : "-")
         << ',' << different_degree_run.S << ','
         << (different_degree_convergence
                 ? std::to_string(different_degree_convergence.value())
                 : "-")
         << '\n';
  }

  fout.close();
  std::cout << "Файл 'penalty_3.csv' создан.\n";
}

auto penalty_part_4() -> void
{
  std::cout << "--- ЧАСТЬ 4 ---\n";
  constexpr double eps = 1e-4;

  const std::vector<Ptr<IIntegrator>> methods = {
      make<LeftRect>(),
      make<MidRect>(),
      make<Trapezoid>(),
      make<Simpson>(),
  };

  std::cout << "Точное значение: " << EXACT_VALUE << '\n';

  for (const auto& method : methods)
  {
    std::cout << "--- " << method->get_name() << " ---\n";
    const auto run = run_method_with_eps(f, method, eps);
    std::cout << "\nИтог:\n - h = " << run.h << "\n - S = " << run.S
              << "\n - error = " << run.error
              << "\n - abs_error = " << std::abs(EXACT_VALUE - run.S) << "\n\n";
  }
}

// ----------------------------------------------------------------------------
// |                              Implementation                              |
// ----------------------------------------------------------------------------

constexpr auto run_method(const FunctionType f,
                          const Ptr<IIntegrator>& method,
                          u32 n) -> double
{
  const auto step = (ORIG_B - ORIG_A) / n;

  double sum = 0.L;
  for (u32 i = 0; i < n; ++i)
  {
    const auto left = ORIG_A + i * step;
    const auto right = left + step;
    sum += (*method)(f, {left, right});
  }

  return sum;
}

constexpr auto run_methods_interleaved(const FunctionType f,
                                       const IntegratorPair& methods,
                                       const u32 n) -> double
{
  const auto step = (ORIG_B - ORIG_A) / n;

  double sum = 0.L;
  for (u32 i = 0; i < n; ++i)
  {
    const auto left = ORIG_A + i * step;
    const auto right = left + step;
    if (i % 2 == 0)
    {
      sum += (*methods.first)(f, {left, right});
    }
    else
    {
      sum += (*methods.second)(f, {left, right});
    }
  }

  return sum;
}

auto calculate_optimum_step(FunctionType f,
                            const Ptr<IIntegrator>& method,
                            const double eps) -> double
{
  std::cout << "Вычисление оптимального шага\n";
  std::vector<Run> runs;

  constexpr double L = 2;

  u32 partitions = 1;
  for (u32 i = 0; i < 3; ++i)
  {
    const double h = (ORIG_B - ORIG_A) / partitions;
    runs.emplace_back(i + 1, h, run_method(f, method, partitions));
    partitions *= 2;
  }

  const double m = aitken(runs, L);
  std::cout << "   m = " << m << '\n';

  auto& run = runs.back();
  run.error = std::abs((run.S - runs.at(1).S) / (std::pow(L, m) - 1));

  return run.h * std::pow(eps / run.error, 1. / m);
}

auto run_method_with_eps(const FunctionType f,
                         const Ptr<IIntegrator>& method,
                         const double eps,
                         const bool choose_optimum_step) -> Run
{
  constexpr double L = 2;

  std::vector<Run> runs;

  double h = 0;
  if (choose_optimum_step)
  {
    h = calculate_optimum_step(f, method, eps);
    std::cout << "Оптимальный шаг: " << h << '\n';
  }
  else
  {
    h = ORIG_B - ORIG_A;
  }

  for (std::size_t i = 0; i < 2; ++i)
  {
    runs.emplace_back(i, h, run_method_with_step(f, method, h));
    h /= L;
  }

  std::size_t r = 2;
  do
  {
    runs.emplace_back(r, h, run_method_with_step(f, method, h));
    h /= L;

    // Процесс Эйткена
    const double m = aitken(runs, L);
    std::cout << r - 1 << ". m = " << m << '\n';

    // Метод Ричардсона
    Matrix A(r, r);
    Matrix b(r, 1);
    for (u32 row = 0; row < r; ++row)
    {
      for (u32 col = 0; col < r; ++col)
      {
        A.at(row, col) = std::pow(runs.at(row + 1).h, m + col) -
                         std::pow(runs.at(row).h, m + col);
      }

      b.at(row) = runs.at(row + 1).S - runs.at(row).S;
    }

    const auto C = solve_sle(A, b);

    for (auto& run : runs)
    {
      run.error = 0;

      for (u32 row = 0; row < C.rows(); ++row)
      {
        run.error += C.at(row) * std::pow(run.h, m + row);
      }

      if (std::abs(run.error) < eps)
      {
        return run;
      }
    }

    r += 1;
  } while (true);
}

auto run_method_with_step(const FunctionType f,
                          const Ptr<IIntegrator>& method,
                          const double step) -> double
{
  double sum = 0.L;

  auto left = ORIG_A;
  while (left < ORIG_B)
  {
    const auto right = std::min(ORIG_B, left + step);
    sum += (*method)(f, {left, right});
    left = right;
  }

  return sum;
}

constexpr auto run_methods_with_step_interleaved(const FunctionType f,
                                                 const IntegratorPair& methods,
                                                 const double step) -> double
{
  double sum = 0.L;

  std::size_t i = 0;
  auto left = ORIG_A;
  while (left < ORIG_B)
  {
    const auto right = std::min(ORIG_B, left + step);
    if (i % 2 == 0)
    {
      sum += (*methods.first)(f, {left, right});
    }
    else
    {
      sum += (*methods.second)(f, {left, right});
    }
    left = right;
    ++i;
  }

  return sum;
}

constexpr auto aitken(double S_1, double S_2, double S_3, double L) -> double
{
  return -(std::log(std::abs((S_3 - S_2) / (S_2 - S_1)))) / (std::log(L));
}
constexpr auto aitken(const std::vector<Run>& runs, double L) -> double
{
  const auto size = runs.size();
  assert(size >= 3);
  return aitken(
      runs.at(size - 3).S, runs.at(size - 2).S, runs.at(size - 1).S, L);
}
