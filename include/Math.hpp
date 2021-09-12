#pragma once

#include <limits>

#include "Definition.hpp"

namespace cte::math {

template<FloatingPoint FloatingPointType, Integer IntegerType>
constexpr FloatingPointType pow(const FloatingPointType base, const IntegerType exponent) noexcept {
    if (base == FloatingPointType{0}) return base;
    FloatingPointType result = FloatingPointType{1};
    FloatingPointType cbase = exponent < 0 ? 1 / static_cast<FloatingPointType>(base) : static_cast<FloatingPointType>(base);
    for (IntegerType index = 0; index < abs(exponent); index++) {
        result *= cbase;
    }
    return result;
}

template<Complex ComplexType, Integer IntegerType>
constexpr ComplexType pow(const ComplexType base, const IntegerType exponent) noexcept {
    using ComplexFloatingPointType = ComplexType::value_type;
    if (base == ComplexType{0, 0}) return base;
    ComplexType result{1, 0};
    ComplexType cbase = exponent < 0 ? std::complex{1 / static_cast<ComplexFloatingPointType>(base.real()), 1 / static_cast<ComplexFloatingPointType>(base.imag())} 
                                                   : std::complex{static_cast<ComplexFloatingPointType>(base.real()), static_cast<ComplexFloatingPointType>(base.imag())};
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

template<Integer IntegerType>
constexpr IntegerType abs(const IntegerType value) noexcept {
    return value < IntegerType{} ? - value : value;
}

template<FloatingPoint FloatingPointType>
constexpr FloatingPointType abs(const FloatingPointType value) noexcept {
    return value < FloatingPointType{} ? - value : value;
}

template<Complex ComplexType>
constexpr auto abs(const ComplexType value) noexcept -> ComplexType::value_type {
    return cte::math::sqrt(static_cast<ComplexType::value_type>(value.real() * value.real() + value.imag() * value.imag()));
}

}

