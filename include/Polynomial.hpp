#pragma once

#include <array>
#include <complex>
#include <cassert>

#include "Math.hpp"
#include "Exception.hpp"
#include "Definition.hpp"
#include "CustomFormatter.hpp"

namespace cte::poly {

template<Arithmetic ArithmeticType, std::size_t Degree>
struct Polynomial {

    typedef ArithmeticType type;

    constexpr Polynomial() {
        mCoeffs.fill(ArithmeticType{0});
    }

    constexpr Polynomial(std::initializer_list<ArithmeticType>&& param) {
        std::copy(param.begin(), param.end(), mCoeffs.begin());
    }

    constexpr Polynomial(const Polynomial& rhs) {
        mCoeffs = rhs.mCoeffs;
    }

    constexpr Polynomial(Polynomial&& rhs) {
        mCoeffs = std::move(rhs.mCoeffs);
    }

    constexpr ArithmeticType& operator[](const std::size_t index) {
        if (index >= Degree + 1) throw cte::OutOfRangeException(index);
        return mCoeffs[index];
    }

    constexpr const ArithmeticType& operator[](const std::size_t index) const {
        if (index >= Degree + 1) throw cte::OutOfRangeException(index);
        return mCoeffs[index];
    }

    constexpr Polynomial& operator=(const Polynomial& rhs) {
        mCoeffs = rhs.mCoeffs;
        return *this;
    }

    constexpr Polynomial& operator=(Polynomial&& rhs) {
        mCoeffs = std::move(rhs.mCoeffs);
        return *this;
    }

    template<Arithmetic OtherArithmeticType, std::size_t OtherDegree>
    constexpr operator Polynomial<OtherArithmeticType, OtherDegree>() const & {
        if (Degree != OtherDegree) throw cte::DimensionMismatchException("Invalid conversion, polynomials have different degrees", Degree, OtherDegree);
        Polynomial<OtherArithmeticType, OtherDegree> result;
        for (std::size_t i = 0; i <= Degree; i++) {
            result[i] = static_cast<OtherArithmeticType>(mCoeffs[i]);
        }
        return result;
    }

    template<Arithmetic OtherArithmeticType, std::size_t OtherDegree>
    constexpr operator Polynomial<OtherArithmeticType, OtherDegree>() && {
        if (Degree != OtherDegree) throw cte::DimensionMismatchException("Invalid conversion, polynomials have different degrees", Degree, OtherDegree);
        Polynomial<OtherArithmeticType, OtherDegree> result;
        for (std::size_t i = 0; i <= Degree; i++) {
            result[i] = static_cast<OtherArithmeticType>(std::move(mCoeffs[i]));
        }
        return result;
    }

    static constexpr std::size_t getDegree() {
        return Degree;
    }

    constexpr auto operator+(const auto& rhs) const noexcept -> Polynomial<ArithmeticType, Degree>;
    constexpr auto operator-(const auto& rhs) const noexcept -> Polynomial<ArithmeticType, Degree>;
    constexpr auto operator*(const auto& rhs) const noexcept -> Polynomial<ArithmeticType, Degree>;
    constexpr auto operator/(const auto& rhs) const noexcept -> DivType<ArithmeticType, Degree, ArithmeticType, Degree>;
    constexpr auto operator*(const ArithmeticType& constant) const noexcept -> Polynomial<ArithmeticType, Degree>;
    constexpr auto operator/(const ArithmeticType& constant) const noexcept -> Polynomial<ArithmeticType, Degree>;
    constexpr Polynomial& operator*=(const ArithmeticType constant) const noexcept;
    constexpr Polynomial& operator/=(const ArithmeticType constant) const noexcept;

    constexpr auto derivate() const noexcept -> Polynomial<ArithmeticType, static_cast<int16_t>(Degree) - 1 < 0 ? 0 : Degree - 1>;
    constexpr std::array<std::complex<long double>, Degree> roots() const noexcept;
    constexpr auto calculate(const auto& value) const noexcept -> std::common_type_t<ArithmeticType, decltype(value)>;

public:
    std::array<ArithmeticType, Degree + 1> mCoeffs;
    constexpr static int8_t sPrecision = 4;
};

template<Arithmetic AT1, std::size_t D1, Arithmetic AT2, std::size_t D2>
struct DivType {

    constexpr DivType(const auto& q, const auto& r) : mQuotient(q), mRemainder(r) {

    }

    constexpr const auto& getRemainder() const {
        return mRemainder;
    }

    constexpr const auto& getQuotient() const {
        return mQuotient;
    }

    constexpr auto& getR() {
        return mRemainder;
    }

    constexpr auto& getQ() {
        return mQuotient;
    }

private:
    Polynomial<AT1, D1> mQuotient;
    Polynomial<AT2, D2> mRemainder;
};

template<Arithmetic ArithmeticType, std::size_t Degree>
constexpr auto Polynomial<ArithmeticType, Degree>::operator+(const auto& rhs) const noexcept -> Polynomial<ArithmeticType, Degree> {
    using removed_cvref_rhs = typename std::remove_cvref_t<decltype(rhs)>::type;
    static_assert(Degree == rhs.getDegree());
    static_assert(std::is_same_v<ArithmeticType, removed_cvref_rhs>);
    Polynomial<ArithmeticType, Degree> result;
    for (std::size_t index = 0; index <= Degree; index++) {
        result[index] = mCoeffs[index] + rhs[index];
    }
    return result;
}

template<Arithmetic ArithmeticType, std::size_t Degree>
constexpr auto Polynomial<ArithmeticType, Degree>::operator-(const auto& rhs) const noexcept -> Polynomial<ArithmeticType, Degree> {
    using removed_cvref_rhs = typename std::remove_cvref_t<decltype(rhs)>::type;
    static_assert(Degree == rhs.getDegree());
    static_assert(std::is_same_v<ArithmeticType, removed_cvref_rhs>);
    auto indirect_sub = [&]() {
        Polynomial<ArithmeticType, Degree> result;
        for (std::size_t index = 0; index <= Degree; index++) {
            result[index] = mCoeffs[index] - rhs[index];
        }
        return result;
    };
    Polynomial<ArithmeticType, Degree> result = indirect_sub();
    return result;
}

template<Arithmetic ArithmeticType, std::size_t Degree>
constexpr auto Polynomial<ArithmeticType, Degree>::operator*(const auto& rhs) const noexcept -> Polynomial<ArithmeticType, Degree> {
    using removed_cvref_rhs = typename std::remove_cvref_t<decltype(rhs)>::type;
    static_assert(Degree == rhs.getDegree());
    static_assert(std::is_same_v<ArithmeticType, removed_cvref_rhs>);
    constexpr std::size_t degree = Degree + rhs.getDegree();
    Polynomial<ArithmeticType, degree> result;
    for (std::size_t index1 = 0; index1 <= Degree; index1++) {
        for (std::size_t index2 = 0; index2 <= rhs.getDegree(); index2++) {
            result[degree - index1 - index2] += mCoeffs[Degree - index1] * rhs[rhs.getDegree() - index2];
        }
    }
    Polynomial<ArithmeticType, Degree> final_result;
    for (std::size_t index = 0; index <= Degree; index++) {
        final_result[Degree - index] = result[degree - index];
    }
    return final_result;
}

template<Arithmetic ArithmeticType, std::size_t Degree>
constexpr auto Polynomial<ArithmeticType, Degree>::operator/(const auto& rhs) const noexcept -> DivType<ArithmeticType, Degree, ArithmeticType, Degree> {
    using removed_cvref_rhs = typename std::remove_cvref_t<decltype(rhs)>::type;
    static_assert(Degree == rhs.getDegree());
    static_assert(std::is_same_v<ArithmeticType, removed_cvref_rhs>);
    auto indirect_div = [&]() {
        Polynomial<ArithmeticType, Degree> q;
        Polynomial<ArithmeticType, Degree> r;
        Polynomial<ArithmeticType, Degree> copy_poly1 = *this;

        std::size_t nz1 = 0, nz2 = 0;
        for (std::size_t index = 0; index <= Degree; index++) {
            if (mCoeffs[index] != 0) {
                nz1 = index;
                break;
            }
        }
        for (std::size_t index = 0; index <= rhs.getDegree(); index++) {
            if (rhs[index] != 0) {
                nz2 = index;
                break;
            }
        }
        std::size_t p1_degree = Degree - nz1;
        std::size_t p2_degree = rhs.getDegree() - nz2;

        std::size_t degree = p1_degree > p2_degree ? p1_degree - p2_degree : 0;

        for (std::size_t q_index = nz1 + p2_degree; q_index <= degree + nz1 + p2_degree; q_index++) {
            q[q_index] = copy_poly1[q_index - p2_degree] / static_cast<long double>(rhs[nz2]);
            for (std::size_t p2_index = q_index; p2_index <= p2_degree + q_index; p2_index++) {
                copy_poly1[p2_index - p2_degree] -= rhs[p2_index - q_index + nz2] * q[q_index];
            }
        }
        return DivType(q, copy_poly1);
    };
    auto result = indirect_div();
    return DivType(result.getQuotient(), result.getRemainder());
}

template<Arithmetic ArithmeticType, std::size_t Degree>
constexpr auto Polynomial<ArithmeticType, Degree>::operator*(const ArithmeticType& constant) const noexcept -> Polynomial<ArithmeticType, Degree> {
    Polynomial<ArithmeticType, Degree> result;
    for (std::size_t index = 0; index <= Degree; index++) {
        result[index] = mCoeffs[index] * constant;
    }
    return result;
}

template<Arithmetic ArithmeticType, std::size_t Degree>
constexpr auto Polynomial<ArithmeticType, Degree>::operator/(const ArithmeticType& constant) const noexcept -> Polynomial<ArithmeticType, Degree> {
    Polynomial<ArithmeticType, Degree> result;
    for (std::size_t index = 0; index <= Degree; index++) {
        result[index] = mCoeffs[index] / constant;
    }
    return result;
}

template<Arithmetic ArithmeticType, std::size_t Degree>
constexpr Polynomial<ArithmeticType, Degree>& Polynomial<ArithmeticType, Degree>::operator*=(const ArithmeticType constant) const noexcept {
    for (std::size_t index = 0; index <= Degree; index++) {
        mCoeffs[index] = mCoeffs[index] * constant;
    }
    return *this;
}

template<Arithmetic ArithmeticType, std::size_t Degree>
constexpr Polynomial<ArithmeticType, Degree>& Polynomial<ArithmeticType, Degree>::operator/=(const ArithmeticType constant) const noexcept {
    for (std::size_t index = 0; index <= Degree; index++) {
        mCoeffs[index] = mCoeffs[index] / constant;
    }
    return *this;
}

template<Arithmetic ArithmeticType, std::size_t Degree>
constexpr auto Polynomial<ArithmeticType, Degree>::derivate() const noexcept -> Polynomial<ArithmeticType, static_cast<int16_t>(Degree) - 1 < 0 ? 0 : Degree - 1> {
    constexpr std::size_t degree = static_cast<int16_t>(Degree) - 1 < 0 ? 0 : Degree - 1;
    Polynomial<ArithmeticType, degree> result;
    for (std::size_t index = 0; index < Degree; index++) {
        result[index] = mCoeffs[index] * (Degree - index); 
    }
    return result;
}

template<Arithmetic ArithmeticType, std::size_t Degree>
constexpr std::array<std::complex<long double>, Degree> Polynomial<ArithmeticType, Degree>::roots() const noexcept {
    std::array<std::complex<long double>, Degree> z, w;
    std::complex<long double> ratio, inverse_sum;
    long double error = cte::math::pow(10, -sPrecision);
    auto derivative = derivate();
    bool check = true;   

    for (std::size_t index = 0; index < Degree; index++) {
        z[index] = std::complex<long double>(0.5L * index + 0.1L, 0.5L * index + 0.1L);
    }
    
    while (check == true) {
        check = false;
        for (std::size_t index = 0; index < Degree; index++) {
            inverse_sum = std::complex<long double>(0.0L, 0.0L);
            for (std::size_t sum_index = 0; sum_index < Degree; sum_index++) {
                inverse_sum += index != sum_index ? 1.0L / (z[index] - z[sum_index]) : 0.0L;
            }
            ratio = calculate(z[index]) / derivative.calculate(z[index]);
            w[index] = ratio / (1.0L - ratio * inverse_sum);
        }
        for (std::size_t index = 0; index < Degree; index++) {
            z[index] = z[index] - w[index];
            z[index] = cte::math::abs(z[index].imag()) < 1e-20 ? std::complex<long double>(z[index].real(), 0.0) : z[index];
            check = cte::math::abs(w[index].real()) > error || cte::math::abs(w[index].imag()) > error;
        }
    }
    return z;
}

template<Arithmetic ArithmeticType, std::size_t Degree>
constexpr auto Polynomial<ArithmeticType, Degree>::calculate(const auto& value) const noexcept -> std::common_type_t<ArithmeticType, decltype(value)> {
    std::common_type_t<ArithmeticType, decltype(value)> result{}, temp{1};
    for (std::size_t index = 0; index < Degree + 1; index++) {
        result += static_cast<decltype(value)>(mCoeffs[Degree - index]) * temp;
        temp *= value;
    }
    return result;
}

template<auto polynomial>
constexpr std::size_t countLeadingZeros() {
    std::size_t count = 0;
    for (std::size_t i = 0; i < polynomial.getDegree(); i++) {
        if (polynomial[i] == 0) { count++; }
    }
    return count;
}

constexpr std::size_t countLeadingZeros(const auto& polynomial) {
    std::size_t count = 0;
    for (std::size_t i = 0; i < polynomial.getDegree(); i++) {
        if (polynomial[i] == 0) { count++; }
    }
    return count;
}

template<auto polynomial>
constexpr auto reduceDegree() {
    using removed_cvref_type = typename std::remove_cvref_t<decltype(polynomial)>::type;
    constexpr std::size_t zeros = countLeadingZeros<polynomial>();
    constexpr std::size_t degree = polynomial.getDegree() - zeros;
    Polynomial<removed_cvref_type, degree> result;
    for (std::size_t index = 0; index <= degree; index++) {
        result[index] = polynomial[zeros + index];
    }
    return result;
}

constexpr auto reduceDegree(const auto& polynomial) {
    using removed_cvref_type = typename std::remove_cvref_t<decltype(polynomial)>::type;
    constexpr std::size_t zeros = countLeadingZeros(polynomial);
    constexpr std::size_t degree = polynomial.getDegree() - zeros;
    Polynomial<removed_cvref_type, degree> result;
    for (std::size_t index = 0; index <= degree; index++) {
        result[index] = polynomial[zeros + index];
    }
    return result;
}

template<auto polynomial1, auto polynomial2>
constexpr auto add() {
    using removed_cvref_type1 = typename std::remove_cvref_t<decltype(polynomial1)>::type;
    using removed_cvref_type2 = typename std::remove_cvref_t<decltype(polynomial2)>::type;
    constexpr std::size_t degree = polynomial1.getDegree() > polynomial2.getDegree() ? polynomial1.getDegree() : polynomial2.getDegree();
    Polynomial<std::common_type_t<removed_cvref_type1, removed_cvref_type2>, degree> result;
    for (std::size_t index = 0; index <= degree; index++) {
        result[index] = (static_cast<int>(polynomial1.getDegree() - degree + index) >= 0 ? polynomial1[polynomial1.getDegree() - degree + index] : 0) + (static_cast<int>(polynomial2.getDegree() - degree + index) >= 0 ? polynomial2[polynomial2.getDegree() - degree + index] : 0);
    }
    return result;
}

template<auto polynomial1, auto polynomial2>
constexpr auto sub() {
    using removed_cvref_type1 = typename std::remove_cvref_t<decltype(polynomial1)>::type;
    using removed_cvref_type2 = typename std::remove_cvref_t<decltype(polynomial2)>::type;
    constexpr std::size_t degree = polynomial1.getDegree() > polynomial2.getDegree() ? polynomial1.getDegree() : polynomial2.getDegree();
    constexpr auto indirect_sub = []() {
        Polynomial<std::common_type_t<removed_cvref_type1, removed_cvref_type2>, degree> result;
        for (std::size_t index = 0; index <= degree; index++) {
            result[index] = (static_cast<int>(polynomial1.getDegree() - degree + index) >= 0 ? polynomial1[polynomial1.getDegree() - degree + index] : 0) - (static_cast<int>(polynomial2.getDegree() - degree + index) >= 0 ? polynomial2[polynomial2.getDegree() - degree + index] : 0);
        }
        return result;
    };
    constexpr Polynomial<std::common_type_t<removed_cvref_type1, removed_cvref_type2>, degree> result = indirect_sub();
    constexpr std::size_t zeros = countLeadingZeros<result>();
    Polynomial<std::common_type_t<removed_cvref_type1, removed_cvref_type2>, degree - zeros> final_result;
    for (std::size_t index = 0; index <= degree - zeros; index++) {
        final_result[index] = result[index + zeros];
    }
    return final_result;
}

template<auto polynomial1, auto polynomial2>
constexpr auto mult() {
    using removed_cvref_type1 = typename std::remove_cvref_t<decltype(polynomial1)>::type;
    using removed_cvref_type2 = typename std::remove_cvref_t<decltype(polynomial2)>::type;
    constexpr std::size_t degree = polynomial1.getDegree() + polynomial2.getDegree();
    Polynomial<std::common_type_t<removed_cvref_type1, removed_cvref_type2>, degree> result;
    for (std::size_t index1 = 0; index1 <= polynomial1.getDegree(); index1++) {
        for (std::size_t index2 = 0; index2 <= polynomial2.getDegree(); index2++) {
            result[degree - index1 - index2] += polynomial1[polynomial1.getDegree() - index1] * polynomial2[polynomial2.getDegree() - index2];
        }
    }
    return result;
}

template<auto polynomial1, auto polynomial2>
constexpr auto div() {
    constexpr auto indirect_div = []() {
        constexpr std::size_t degree = polynomial1.getDegree() > polynomial2.getDegree() ? polynomial1.getDegree() - polynomial2.getDegree() : 0;
        Polynomial<long double, degree> q;
        Polynomial<long double, polynomial2.getDegree()> r;
        Polynomial<long double, polynomial1.getDegree()> copy_poly1 = polynomial1;
        for (std::size_t q_index = 0; q_index <= degree; q_index++) {
            q[q_index] = copy_poly1[q_index] / static_cast<long double>(polynomial2[0]);
            for (std::size_t p2_index = q_index; p2_index <= polynomial2.getDegree() + q_index; p2_index++) {
                copy_poly1[p2_index] -= polynomial2[p2_index - q_index] * q[q_index];
            }
        }
        return DivType(q, copy_poly1);
    };
    constexpr auto result = indirect_div();
    return DivType(result.getQuotient(), reduceDegree<result.getRemainder()>());
}

}

