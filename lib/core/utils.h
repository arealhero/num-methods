#pragma once

#include "types.h"

#include <cmath>
#include <concepts>
#include <limits>
#include <memory>
#include <utility>

template <typename Type>
using Ptr = std::shared_ptr<Type>;

template <typename Type, typename... Args>
auto make(Args&&... args) -> Ptr<Type>
{
  return std::make_shared<Type>(std::forward<Args>(args)...);
}
