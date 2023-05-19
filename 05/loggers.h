#pragma once

#include <memory>
#include <vector>

class ILogger
{
 public:
  struct Point
  {
    double x;
    double value;
  };

  virtual ~ILogger() = default;

  virtual void insert_total_error(const Point& point) = 0;
  virtual void insert_local_error(const Point& point) = 0;

  virtual void insert_step_value(const Point& point) = 0;

  [[nodiscard]] virtual std::vector<Point> get_total_errors() const = 0;
  [[nodiscard]] virtual std::vector<Point> get_local_errors() const = 0;

  [[nodiscard]] virtual std::vector<Point> get_step_values() const = 0;
};

class DummyLogger : public ILogger
{
 public:
  void insert_total_error([[maybe_unused]] const Point& point) override {}
  void insert_local_error([[maybe_unused]] const Point& point) override {}

  void insert_step_value([[maybe_unused]] const Point& point) override {}

  [[nodiscard]] std::vector<Point> get_total_errors() const override
  {
    return {};
  }
  [[nodiscard]] std::vector<Point> get_local_errors() const override
  {
    return {};
  }

  [[nodiscard]] std::vector<Point> get_step_values() const override
  {
    return {};
  }
};

class VectorLogger : public ILogger
{
 public:
  void insert_total_error([[maybe_unused]] const Point& point) override
  {
    total_errors.push_back(point);
  }
  void insert_local_error([[maybe_unused]] const Point& point) override
  {
    local_errors.push_back(point);
  }

  void insert_step_value([[maybe_unused]] const Point& point) override
  {
    step_values.push_back(point);
  }

  [[nodiscard]] std::vector<Point> get_total_errors() const override
  {
    return total_errors;
  }
  [[nodiscard]] std::vector<Point> get_local_errors() const override
  {
    return local_errors;
  }

  [[nodiscard]] std::vector<Point> get_step_values() const override
  {
    return step_values;
  }

 private:
  std::vector<Point> total_errors{};
  std::vector<Point> local_errors{};
  std::vector<Point> step_values{};
};

using LoggerPtr = std::shared_ptr<ILogger>;
