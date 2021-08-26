#pragma once

#include <limits>
#include <numeric>
#include <complex>
#include <type_traits>

template<typename... Types> void ignore(const Types&...) { }

template<typename Type>
concept Arithmetic = std::is_arithmetic<Type>::value && !std::is_same<Type, char>::value;

template<typename Type>
concept Integer = std::is_integral<Type>::value && !std::is_same<Type, char>::value;

template <typename charType, typename Traits>
std::basic_ostream<charType, Traits>& operator<<(std::basic_ostream<charType, Traits> &os, const std::complex<double>& value) {
    os << value.real();
    if (value.imag() != 0) {
        os << std::showpos << value.imag() << "*i";
    }
    os << std::noshowpos;
    return os;
}

namespace cte::math {

template<typename Type>
constexpr Type abs(const Type value) noexcept {
    return value == Type{} ? Type{} : value < Type{} ? - value : value;
}

template<typename Type, Integer IntegerType>
constexpr double pow(const Type base, const IntegerType exponent) {
    double result = 1.0;
    double cbase = exponent < 0 ? 1 / static_cast<double>(base) : static_cast<double>(base);
    for (IntegerType index = 0; index < abs(exponent); index++) {
        result *= cbase;
    }
    return result;
}

constexpr double sqrt(double value) {
    constexpr auto sqrtNewtonRaphson = [](const auto sqrtNewtonRaphson, double value, double current, double previous) -> double {
        return current == previous ? current : sqrtNewtonRaphson(sqrtNewtonRaphson, value, 0.5 * (current + value / current), current);
    };
    return value >= 0 && value < std::numeric_limits<double>::infinity() ? sqrtNewtonRaphson(sqrtNewtonRaphson, value, value, 0) : std::numeric_limits<double>::quiet_NaN();
}

}

namespace cte::utils {

constexpr uint64_t factorial(uint64_t n) {
    return (n <= 1) ? 1 : n * factorial(n - 1);
}

template<std::size_t Elem>
constexpr std::array<std::array<uint8_t, Elem>, factorial(Elem)> computePermutations() {
    std::array<uint8_t, Elem> permutation;
    std::array<std::array<uint8_t, Elem>, factorial(Elem)> result;
    std::iota(permutation.begin(), permutation.end(), 0);
    uint64_t index = 0;
    do {
        result[index++] = permutation;
    } while(std::next_permutation(permutation.begin(), permutation.end()));
    return result;
}

template<std::size_t Elem>
constexpr int8_t permutationSign(const std::array<uint8_t, Elem>& permutation) {
    uint64_t counter = 0;
    for (std::size_t i = 0; i < permutation.size(); i++) {
        for (std::size_t j = i + 1; j < permutation.size(); j++) {
            if (permutation[i] > permutation[j]) {
                counter++;
            }
        }
    }
    if (counter % 2 == 0) return 1;
    else return -1;
}

}

