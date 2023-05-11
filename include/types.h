#pragma once

using value_t = long double;
using func_t = value_t (*)(value_t x);
using method_t = value_t (*)(func_t, value_t a, value_t b);
