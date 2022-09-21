#include <cstdio>
#include <cstdlib>
#include <uni/math/math.hpp>

constexpr std::size_t N = 3;
constexpr static double eps = 1e-6;

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

int main()
{

  std::cout << "eps: " << eps << '\n'
	    << "\nA:\n" << A << '\n'
	    << "b:\n" << b << '\n';

  constexpr static auto gauss_x = uni::math::sle::gauss(A, -b);
  std::cout << "gauss_x\n" << gauss_x << '\n';

  constexpr static auto gauss_test = A * gauss_x;
  std::cout << "gauss_test\n" << gauss_test << '\n';

  auto steepest_x = uni::math::optimization::steepest_descent(A, b, eps);
  std::cout << "steepest_x\n" << steepest_x
	    << "(diff: " << (gauss_x - steepest_x).norm_2() << ")\n\n";

  auto coordinate_x = uni::math::optimization::coordinate_descent(A, b, eps);
  std::cout << "coordinate_x\n" << coordinate_x
	    << "(diff: " << (gauss_x - coordinate_x).norm_2() << ")\n\n";

  return EXIT_SUCCESS;
}
