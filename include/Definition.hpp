#pragma once

#include <complex>
#include <type_traits>

template<typename... Types> void ignore(const Types&...) { }

template<typename Type>
concept Integer = std::is_integral<Type>::value && !std::is_same<Type, char>::value;

template<typename Type>
concept FloatingPoint = std::is_floating_point<Type>::value;

template <FloatingPoint FloatingPointType>
struct is_complex : std::false_type {};

template <FloatingPoint FloatingPointType>
struct is_complex<std::complex<FloatingPointType>> : std::true_type {};

template<typename Type>
concept Complex = is_complex<Type>::value;

namespace cte::mat {

template<FloatingPoint, std::size_t, std::size_t> struct Matrix;

template<typename... Types, size_t... Sizes>
Matrix(Types (&&...arr)[Sizes]) -> Matrix<std::common_type_t<Types...>, sizeof...(Sizes), std::get<0>(std::forward_as_tuple(Sizes...))>;

template<std::size_t> auto linearEquationMapping();
template<auto> constexpr auto calculateEigenvalues();
template<auto> constexpr auto calculateEigenvectors();
template<auto> constexpr auto reductionBareissAlgorithm();
template<auto> constexpr auto getCharacteristicPolynomial();
template<std::size_t> constexpr auto reductionBareissAlgorithmMapping();
template<FloatingPoint, std::size_t> constexpr auto determinant(const auto&);
template<FloatingPoint, std::size_t> constexpr auto gaussDecomposition(const auto&);
template<FloatingPoint, std::size_t> constexpr auto matrixPolynomialComposition(const auto&);
template<FloatingPoint, std::size_t> constexpr auto solveDeterminateSystem(const auto&, const auto&);

template<FloatingPoint, std::size_t, std::size_t, FloatingPoint, std::size_t, std::size_t> struct QRType;

template<typename M1, typename M2>
QRType(M1, M2) -> QRType<typename M1::type, M1::getRows(), M1::getCols(), typename M2::type, M2::getRows(), M2::getCols()>;

}

namespace cte::poly {

template<FloatingPoint, std::size_t> struct Polynomial;

template<typename... Types>
Polynomial(Types...) -> Polynomial<std::common_type_t<Types...>, sizeof...(Types) - 1>;

constexpr auto reduceDegree(const auto&);
template<auto, auto> constexpr auto add();
template<auto, auto> constexpr auto sub();
template<auto, auto> constexpr auto mult();
template<auto, auto> constexpr auto div();
template<auto> constexpr auto reduceDegree();
constexpr std::size_t countLeadingZeros(const auto&);
template<auto> constexpr std::size_t countLeadingZeros();

template<FloatingPoint, std::size_t, FloatingPoint, std::size_t> struct DivType;

template<typename P1, typename P2>
DivType(P1, P2) -> DivType<typename P1::type, P1::getDegree(), typename P2::type, P2::getDegree()>;

}

