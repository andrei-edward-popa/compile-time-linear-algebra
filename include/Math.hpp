#pragma once

#include <limits>

#include "Definition.hpp"

namespace cte::math {

template<typename Type>
constexpr Type abs(const Type value) noexcept {
    return value == Type{} ? Type{} : value < Type{} ? - value : value;
}

template<typename Type, Integer IntegerType>
constexpr double pow(const Type base, const IntegerType exponent) noexcept {
    if (base == Type{0}) return base;
    double result = 1.0;
    double cbase = exponent < 0 ? 1 / static_cast<double>(base) : static_cast<double>(base);
    for (IntegerType index = 0; index < abs(exponent); index++) {
        result *= cbase;
    }
    return result;
}

template<FloatingPoint FloatingPointType>
constexpr FloatingPointType sqrt(FloatingPointType value) noexcept {
    constexpr auto sqrtNewtonRaphson = [](const auto sqrtNewtonRaphson, FloatingPointType value, FloatingPointType current, FloatingPointType previous) -> FloatingPointType {
        return current == previous ? current : sqrtNewtonRaphson(sqrtNewtonRaphson, value, 0.5 * (current + value / current), current);
    };
    return value >= 0 && value < std::numeric_limits<FloatingPointType>::infinity() ? sqrtNewtonRaphson(sqrtNewtonRaphson, value, value, 0) : std::numeric_limits<FloatingPointType>::quiet_NaN();
}

}

