#pragma once

#include <type_traits>

template<typename... Types> void ignore(const Types&...) { }

template<typename Type>
concept Arithmetic = std::is_arithmetic<Type>::value && !std::is_same<Type, char>::value;

template<typename Type>
concept Integer = std::is_integral<Type>::value && !std::is_same<Type, char>::value;

template<typename Type>
concept FloatingPoint = std::is_floating_point<Type>::value;

namespace cte::mat {

template<Arithmetic, std::size_t, std::size_t> struct Matrix;

template<typename... Types, size_t... Sizes>
Matrix(Types (&&...arr)[Sizes]) -> Matrix<std::common_type_t<Types...>, sizeof...(Sizes), std::get<0>(std::forward_as_tuple(Sizes...))>;

template<std::size_t> constexpr auto reductionBareissAlgorithmMapping();
template<Arithmetic, std::size_t> constexpr auto matrixPolynomialComposition(auto);
template<auto> constexpr auto reductionBareissAlgorithm();
template<auto> constexpr auto getCharacteristicPolynomial();
template<auto> constexpr auto calculateEigenvalues();
template<auto> constexpr auto calculateEigenvectors();
template<Arithmetic, std::size_t, std::size_t, Arithmetic, std::size_t, std::size_t> struct QRType;
template<std::size_t N> Matrix<uint16_t, N, N - 1> linearEquationMapping();
template<Arithmetic ArithmeticType, std::size_t Size> constexpr std::array<long double, Size> solveDeterminateSystem(const Matrix<ArithmeticType, Size, Size>& coeffs, const std::array<ArithmeticType, Size> freeTerms);

template<typename M1, typename M2>
QRType(M1, M2) -> QRType<typename M1::type, M1::getRows(), M1::getCols(), typename M2::type, M2::getRows(), M2::getCols()>;

}

namespace cte::poly {

template<Arithmetic, std::size_t> struct Polynomial;

template<typename... Types>
Polynomial(Types...) -> Polynomial<std::common_type_t<Types...>, sizeof...(Types) - 1>;

constexpr std::size_t countLeadingZeros(const auto&);
constexpr auto reduceDegree(const auto&);

template<auto, auto> constexpr auto add();
template<auto, auto> constexpr auto sub();
template<auto, auto> constexpr auto mult();
template<auto, auto> constexpr auto div();
template<auto> constexpr std::size_t countLeadingZeros();
template<auto> constexpr auto reduceDegree();

template<Arithmetic, std::size_t, Arithmetic, std::size_t> struct DivType;

template<typename P1, typename P2>
DivType(P1, P2) -> DivType<typename P1::type, P1::getDegree(), typename P2::type, P2::getDegree()>;

}

