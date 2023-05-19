#pragma once

#include <algebra/algebra.h>
#include <core/types.h>

#include <algorithm>
#include <cstdint>
#include <vector>

class Counter
{
 public:
  Counter(u32 count_ = 0) : count(count_) {}

  constexpr Counter& operator++()
  {
    ++count;
    return *this;
  }

  [[nodiscard]] constexpr auto get() const -> u32 { return count; }
  constexpr auto reset() -> void { count = 0; }

 private:
  u32 count{0};
};

class IFunction
{
 public:
  virtual ~IFunction() = default;

  auto operator()(double x, const Matrix& y) -> Matrix
  {
    ++total_counter;
    ++counter;

    return call(x, y);
  }

  auto get_number_of_calls() -> u32 { return counter.get(); }
  auto reset_counter() -> void { counter.reset(); }

  auto get_total_number_of_calls() -> u32 { return total_counter.get(); }

 protected:
  virtual auto call(double x, const Matrix& y) -> Matrix = 0;

 private:
  Counter total_counter{0};
  Counter counter{0};
};

class LinearFunction : public IFunction
{
 public:
  explicit LinearFunction(Matrix A_) : A(std::move(A_)) {}

  auto call([[maybe_unused]] const double x, const Matrix& y) -> Matrix override
  {
    return A * y;
  }

 private:
  Matrix A;
};

using FunctionPtr = std::shared_ptr<IFunction>;
