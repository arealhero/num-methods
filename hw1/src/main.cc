#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <uni/math/math.hpp>

constexpr std::size_t N = 30;
constexpr static double eps = 1e-6;

// clang-format off
constexpr static uni::math::matrix<3> A = {
  4.,  1.,                  1.,
  1.,  2. * (3 + 0.1 * N), -1.,
  1., -1.,                  2. * (4 + 0.1 * N),
};
constexpr static uni::math::vector<3> b = {
   1.,
  -2.,
   3.,
};
// clang-format on

constexpr static auto f(const uni::math::vector<3>& x)
{
  auto x_t = x.T();
  return 0.5 * (x_t * (A * x)) + x_t * b + N;
}

template <std::size_t Rows, std::size_t Cols>
void print_matrix(const uni::math::matrix<Rows, Cols>& mat, const std::size_t precision = 7)
{
  std::cout << std::setprecision(8);
  for (std::size_t row = 0; row < Rows; ++row)
  {
    for (std::size_t col = 0; col < Cols; ++col)
    {
      std::cout << std::right << std::setw(precision + 8) << mat.at(row, col) << ' ';
    }
    std::cout << '\n';
  }
}

int main()
{
  constexpr static std::size_t precision = 7;
  std::cout << std::scientific << std::setprecision(precision);
  std::cout << "eps: " << eps << ", N: " << N << '\n';
  std::cout << "\nx_0:\n";
  print_matrix(uni::math::vector<3>::zero(), precision);
  std::cout << "\nA:\n";
  print_matrix(A, precision);
  std::cout << "\nb:\n";
  print_matrix(b, precision);
  std::cout << '\n';

  constexpr static auto gauss_x = uni::math::sle::gauss(A, -b);
  std::cout << "Метод Гаусса\n\nx:\n";
  print_matrix(gauss_x, precision);
  std::cout << "f(x) = \t\t" << f(gauss_x) << '\n';

  constexpr static auto steepest_x = uni::math::optimization::steepest_descent(A, b, eps);

  std::cout << "\nМНГС\n\nx*:\n";
  print_matrix(steepest_x, precision);
  std::cout << "f(x*) = \t" << f(steepest_x) << '\n';
  std::cout << "||x - x*||_2: \t" << (gauss_x - steepest_x).norm_2() << "\n\n";

  constexpr static auto coordinate_x = uni::math::optimization::coordinate_descent(A, b, eps);
  std::cout << "МНПС\n\nx*:\n";
  print_matrix(coordinate_x, precision);
  std::cout << "f(x*) = \t" << f(coordinate_x) << '\n';
  std::cout << "||x - x*||_2: \t" << (gauss_x - coordinate_x).norm_2() << "\n\n";

  return EXIT_SUCCESS;
}
