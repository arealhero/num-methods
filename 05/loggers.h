#pragma once

#include <algebra/matrix.h>

#include <memory>
#include <vector>

class ILogger
{
 public:
  template <typename Type = double>
  struct Point
  {
    double x;
    Type value;
  };

  virtual ~ILogger() = default;

  virtual void insert_point(const Point<Matrix>& point) = 0;

  virtual void insert_total_error(const Point<double>& point) = 0;
  virtual void insert_local_error(const Point<double>& point) = 0;

  virtual void insert_step_value(const Point<double>& point) = 0;

  [[nodiscard]] virtual std::vector<Point<Matrix>> get_points() const = 0;

  [[nodiscard]] virtual std::vector<Point<double>> get_total_errors() const = 0;
  [[nodiscard]] virtual std::vector<Point<double>> get_local_errors() const = 0;

  [[nodiscard]] virtual std::vector<Point<double>> get_step_values() const = 0;
};

class DummyLogger : public ILogger
{
 public:
  void insert_point([[maybe_unused]] const Point<Matrix>& point) override {}

  void insert_total_error([[maybe_unused]] const Point<double>& point) override
  {
  }
  void insert_local_error([[maybe_unused]] const Point<double>& point) override
  {
  }

  void insert_step_value([[maybe_unused]] const Point<double>& point) override
  {
  }

  [[nodiscard]] std::vector<Point<Matrix>> get_points() const override
  {
    return {};
  }

  [[nodiscard]] std::vector<Point<double>> get_total_errors() const override
  {
    return {};
  }
  [[nodiscard]] std::vector<Point<double>> get_local_errors() const override
  {
    return {};
  }

  [[nodiscard]] std::vector<Point<double>> get_step_values() const override
  {
    return {};
  }
};

class VectorLogger : public ILogger
{
 public:
  void insert_point(const Point<Matrix>& point) override
  {
    points.push_back(point);
  }

  void insert_total_error(const Point<double>& point) override
  {
    total_errors.push_back(point);
  }
  void insert_local_error(const Point<double>& point) override
  {
    local_errors.push_back(point);
  }

  void insert_step_value(const Point<double>& point) override
  {
    step_values.push_back(point);
  }

  [[nodiscard]] std::vector<Point<Matrix>> get_points() const override
  {
    return points;
  }

  [[nodiscard]] std::vector<Point<double>> get_total_errors() const override
  {
    return total_errors;
  }
  [[nodiscard]] std::vector<Point<double>> get_local_errors() const override
  {
    return local_errors;
  }

  [[nodiscard]] std::vector<Point<double>> get_step_values() const override
  {
    return step_values;
  }

 private:
  std::vector<Point<Matrix>> points{};
  std::vector<Point<double>> total_errors{};
  std::vector<Point<double>> local_errors{};
  std::vector<Point<double>> step_values{};
};

using LoggerPtr = std::shared_ptr<ILogger>;
