#pragma once

#include "matrix.h"

[[nodiscard]] auto solve_sle(const Matrix& A, const Matrix& b) -> Matrix;

[[nodiscard]] auto gauss(const Matrix& A, const Matrix& b) -> Matrix;

[[nodiscard]] Matrix seidel(const Matrix& A,
                            const Matrix& b,
                            double eps = 1e-6);
