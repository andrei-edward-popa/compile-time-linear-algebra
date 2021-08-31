#pragma once

#include <array>
#include <complex>
#include <cassert>

#include "Math.hpp"
#include "Definition.hpp"
#include "CustomFormatter.hpp"

namespace cte::poly {

template<Arithmetic ArithmeticType, std::size_t Degree>
struct Polynomial {

    typedef ArithmeticType type;

    constexpr Polynomial() {
        mCoeffs.fill(ArithmeticType{0});
    }

    constexpr Polynomial(const auto& rhs) {
        assert(Degree == rhs.getDegree());
        for (std::size_t index = 0; index <= Degree; index++) {
            mCoeffs[index] = static_cast<ArithmeticType>(rhs[index]);
        }
    }

    constexpr Polynomial(const Polynomial& rhs) {
        mCoeffs = rhs.mCoeffs;
    }

    constexpr Polynomial(Polynomial&& rhs) {
        mCoeffs = std::move(rhs.mCoeffs);
    }

    constexpr Polynomial(std::initializer_list<ArithmeticType>&& param) {
        assert(param.size() == Degree + 1);
        std::copy(param.begin(), param.end(), mCoeffs.begin());
    }

    constexpr auto derivate() const {
        constexpr std::size_t degree = static_cast<int16_t>(Degree) - 1 < 0 ? 0 : Degree - 1;
        Polynomial<ArithmeticType, degree> result;
        for (std::size_t index = 0; index < Degree; index++) {
            result[index] = mCoeffs[index] * (Degree - index); 
        }
        return result;
    }

    template<typename ArgType>
    constexpr ArgType calculate(const ArgType& value) const {
        ArgType result{}, temp{1};
        for (std::size_t index = 0; index < Degree + 1; index++) {
            result += static_cast<ArgType>(mCoeffs[Degree - index]) * temp;
            temp *= value;
        }
        return result;
    }

    template<typename ArgType>
    constexpr std::complex<ArgType> calculate(const std::complex<ArgType>& value) const {
        std::complex<ArgType> result{}, temp{1, 0};
        for (std::size_t index = 0; index < Degree + 1; index++) {
            result += static_cast<ArgType>(mCoeffs[Degree - index]) * temp;
            temp *= value;
            if (cte::math::abs(temp.imag()) < 1e-20) {
                temp = std::complex<ArgType>(temp.real(), ArgType{});
            }
        }
        return result;
    }

    constexpr std::array<std::complex<double>, Degree> roots() const {
        std::array<std::complex<double>, Degree> z, w;
        auto derivative = derivate();
        for (std::size_t index = 0; index < Degree; index++) {
            z[index] = std::complex<double>(0.5 * index + 0.1, 0.5 * index + 0.1);
        }
        std::complex<double> ratio, inverse_sum;
        double error = cte::math::pow(10, -sPrecision);
        bool check;
        while (true) {
            check = false;
            for (std::size_t index = 0; index < Degree; index++) {
                inverse_sum = std::complex<double>(0, 0);
                for (std::size_t sum_index = 0; sum_index < Degree; sum_index++) {
                    if (index != sum_index) {
                        inverse_sum += 1.0 / (z[index] - z[sum_index]);
                    }
                }
                ratio = calculate(z[index]) / derivative.calculate(z[index]);
                w[index] = ratio / (1.0 - ratio * inverse_sum);
            }
            for (std::size_t index = 0; index < Degree; index++) {
                z[index] = z[index] - w[index];
                if (cte::math::abs(w[index].real()) > error || cte::math::abs(w[index].imag()) > error) {
                    check = true;
                }
                if (cte::math::abs(z[index].imag()) < 1e-20) {
                    z[index] = std::complex<double>(z[index].real(), 0.0);
                }
            }
            if (check == false) {
                break;
            }
        }
        return z;
    }

    constexpr bool isNull() const {
        return Degree == 0 && mCoeffs[0] == 0;
    }

    static constexpr std::size_t getDegree() {
        return Degree;
    }

    constexpr ArithmeticType& operator[](const std::size_t index) {
        assert(index < Degree + 1);
        return mCoeffs[index];
    }

    constexpr const ArithmeticType& operator[](const std::size_t index) const {
        assert(index < Degree + 1);
        return mCoeffs[index];
    }

    constexpr Polynomial& operator=(const Polynomial& rhs) {
        mCoeffs = rhs.mCoeffs;
        return *this;
    }

    constexpr auto& operator=(const auto& rhs) {
        assert(Degree == rhs.getDegree());
        for (std::size_t index = 0; index <= Degree; index++) {
            mCoeffs[index] = static_cast<ArithmeticType>(rhs[index]);
        }
        return *this;
    }

    constexpr Polynomial& operator=(Polynomial&& rhs) {
        mCoeffs = std::move(rhs.mCoeffs);
        return *this;
    }


    template<Arithmetic ArgArithmeticType>
    constexpr Polynomial operator*(const ArgArithmeticType constant) const {
        Polynomial<std::common_type_t<ArithmeticType, ArgArithmeticType>, Degree> result;
        for (std::size_t index = 0; index <= Degree; index++) {
            result[index] = mCoeffs[index] * constant;
        }
        return result;
    }

    constexpr Polynomial& operator*=(const ArithmeticType constant) const {
        for (std::size_t index = 0; index <= Degree; index++) {
            mCoeffs[index] = mCoeffs[index] * constant;
        }
        return *this;
    }

    template<Arithmetic ArgArithmeticType>
    constexpr Polynomial operator/(const ArithmeticType constant) const {
        Polynomial<std::common_type_t<ArithmeticType, ArgArithmeticType>, Degree> result;
        for (std::size_t index = 0; index <= Degree; index++) {
            result[index] = mCoeffs[index] / constant;
        }
        return result;
    }

    constexpr Polynomial& operator/=(const ArithmeticType constant) const {
        for (std::size_t index = 0; index <= Degree; index++) {
            mCoeffs[index] = mCoeffs[index] / constant;
        }
        return *this;
    }

public:
    std::array<ArithmeticType, Degree + 1> mCoeffs;
    constexpr static int8_t sPrecision = 4;
};

template<Arithmetic AT1, Arithmetic AT2, std::size_t D1, std::size_t D2>
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

template<auto polynomial>
constexpr std::size_t countLeadingZeros() {
    std::size_t count = 0;
    for (std::size_t i = 0; i < polynomial.getDegree(); i++) {
        if (polynomial[i] == 0) { count++; }
    }
    return count;
}

template<auto polynomial>
constexpr auto reduceDegree() {
    constexpr std::size_t zeros = countLeadingZeros<polynomial>();
    constexpr std::size_t degree = polynomial.getDegree() - zeros;
    Polynomial<typename decltype(polynomial)::type, degree> result;
    for (std::size_t index = 0; index <= degree; index++) {
        result[index] = polynomial[zeros + index];
    }
    return result;
}

template<auto polynomial1, auto polynomial2>
constexpr auto add() {
    constexpr std::size_t degree = polynomial1.getDegree() > polynomial2.getDegree() ? polynomial1.getDegree() : polynomial2.getDegree();
    Polynomial<std::common_type_t<typename decltype(polynomial1)::type, typename decltype(polynomial2)::type>, degree> result;
    for (std::size_t index = 0; index <= degree; index++) {
        result[index] = (static_cast<int>(polynomial1.getDegree() - degree + index) >= 0 ? polynomial1[polynomial1.getDegree() - degree + index] : 0) + (static_cast<int>(polynomial2.getDegree() - degree + index) >= 0 ? polynomial2[polynomial2.getDegree() - degree + index] : 0);
    }
    return result;
}

template<auto polynomial1, auto polynomial2>
constexpr auto sub() {
    constexpr auto indirect_sub = []() {
        constexpr std::size_t degree = polynomial1.getDegree() > polynomial2.getDegree() ? polynomial1.getDegree() : polynomial2.getDegree();
        Polynomial<std::common_type_t<typename decltype(polynomial1)::type, typename decltype(polynomial2)::type>, degree> result;
        for (std::size_t index = 0; index <= degree; index++) {
            result[index] = (static_cast<int>(polynomial1.getDegree() - degree + index) >= 0 ? polynomial1[polynomial1.getDegree() - degree + index] : 0) - (static_cast<int>(polynomial2.getDegree() - degree + index) >= 0 ? polynomial2[polynomial2.getDegree() - degree + index] : 0);
        }
        return result;
    };
    constexpr std::size_t degree = polynomial1.getDegree() > polynomial2.getDegree() ? polynomial1.getDegree() : polynomial2.getDegree();
    constexpr Polynomial<std::common_type_t<typename decltype(polynomial1)::type, typename decltype(polynomial2)::type>, degree> result = indirect_sub();
    constexpr std::size_t zeros = countLeadingZeros<result>();
    Polynomial<std::common_type_t<typename decltype(polynomial1)::type, typename decltype(polynomial2)::type>, degree - zeros> final_result;
    for (std::size_t index = 0; index <= degree - zeros; index++) {
        final_result[index] = result[index + zeros];
    }
    return final_result;
}

template<auto polynomial1, auto polynomial2>
constexpr auto mult() {
    constexpr std::size_t degree = polynomial1.getDegree() + polynomial2.getDegree();
    Polynomial<std::common_type_t<typename decltype(polynomial1)::type, typename decltype(polynomial2)::type>, degree> result;
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
        Polynomial<double, degree> q;
        Polynomial<double, polynomial2.getDegree()> r;
        Polynomial<double, polynomial1.getDegree()> copy_poly1 = polynomial1;
        for (std::size_t q_index = 0; q_index <= degree; q_index++) {
            q[q_index] = copy_poly1[q_index] / static_cast<double>(polynomial2[0]);
            for (std::size_t p2_index = q_index; p2_index <= polynomial2.getDegree() + q_index; p2_index++) {
                copy_poly1[p2_index] -= polynomial2[p2_index - q_index] * q[q_index];
            }
        }
        return DivType(q, copy_poly1);
    };
    constexpr auto result = indirect_div();
    return DivType(result.getQuotient(), reduceDegree<result.getRemainder()>());
}

constexpr auto add_preserve_dim(auto polynomial1, auto polynomial2) {
    static_assert(polynomial1.getDegree() == polynomial2.getDegree());
    static_assert(std::is_same_v<typename decltype(polynomial1)::type, typename decltype(polynomial2)::type>);
    constexpr std::size_t degree = polynomial1.getDegree();
    Polynomial<typename decltype(polynomial1)::type, degree> result;
    for (std::size_t index = 0; index <= degree; index++) {
        result[index] = (static_cast<int>(polynomial1.getDegree() - degree + index) >= 0 ? polynomial1[polynomial1.getDegree() - degree + index] : 0) + (static_cast<int>(polynomial2.getDegree() - degree + index) >= 0 ? polynomial2[polynomial2.getDegree() - degree + index] : 0);
    }
    return result;
}

constexpr auto sub_preserve_dim(auto polynomial1, auto polynomial2) {
    static_assert(polynomial1.getDegree() == polynomial2.getDegree());
    static_assert(std::is_same_v<typename decltype(polynomial1)::type, typename decltype(polynomial2)::type>);
    constexpr std::size_t degree = polynomial1.getDegree();
    auto indirect_sub = [polynomial1, polynomial2]() {
        Polynomial<typename decltype(polynomial1)::type, degree> result;
        for (std::size_t index = 0; index <= degree; index++) {
            result[index] = (static_cast<int>(polynomial1.getDegree() - degree + index) >= 0 ? polynomial1[polynomial1.getDegree() - degree + index] : 0) - (static_cast<int>(polynomial2.getDegree() - degree + index) >= 0 ? polynomial2[polynomial2.getDegree() - degree + index] : 0);
        }
        return result;
    };
    Polynomial<typename decltype(polynomial1)::type, degree> result = indirect_sub();
    return result;
}

constexpr auto mult_preserve_dim(auto polynomial1, auto polynomial2) {
    static_assert(polynomial1.getDegree() == polynomial2.getDegree());
    static_assert(std::is_same_v<typename decltype(polynomial1)::type, typename decltype(polynomial2)::type>);
    constexpr std::size_t degree = polynomial1.getDegree() + polynomial2.getDegree();
    Polynomial<typename decltype(polynomial1)::type, degree> result;
    for (std::size_t index1 = 0; index1 <= polynomial1.getDegree(); index1++) {
        for (std::size_t index2 = 0; index2 <= polynomial2.getDegree(); index2++) {
            result[degree - index1 - index2] += polynomial1[polynomial1.getDegree() - index1] * polynomial2[polynomial2.getDegree() - index2];
        }
    }
    Polynomial<typename decltype(polynomial1)::type, polynomial1.getDegree()> final_result;
    for (std::size_t index = 0; index <= polynomial1.getDegree(); index++) {
        final_result[polynomial1.getDegree() - index] = result[degree - index];
    }
    return final_result;
}

constexpr auto div_preserve_dim(auto polynomial1, auto polynomial2) {
    static_assert(polynomial1.getDegree() == polynomial2.getDegree());
    static_assert(std::is_same_v<typename decltype(polynomial1)::type, typename decltype(polynomial2)::type>);
    auto indirect_div = [&]() {
        Polynomial<typename decltype(polynomial1)::type, polynomial1.getDegree()> q;
        Polynomial<typename decltype(polynomial1)::type, polynomial2.getDegree()> r;
        Polynomial<typename decltype(polynomial1)::type, polynomial1.getDegree()> copy_poly1 = polynomial1;

        std::size_t nz1 = 0, nz2 = 0;
        for (std::size_t index = 0; index <= polynomial1.getDegree(); index++) {
            if (polynomial1[index] != 0) {
                nz1 = index;
                break;
            }
        }
        for (std::size_t index = 0; index <= polynomial2.getDegree(); index++) {
            if (polynomial2[index] != 0) {
                nz2 = index;
                break;
            }
        }
        std::size_t p1_degree = polynomial1.getDegree() - nz1;
        std::size_t p2_degree = polynomial2.getDegree() - nz2;

        std::size_t degree = p1_degree > p2_degree ? p1_degree - p2_degree : 0;

        for (std::size_t q_index = nz1 + p2_degree; q_index <= degree + nz1 + p2_degree; q_index++) {
            q[q_index] = copy_poly1[q_index - p2_degree] / static_cast<double>(polynomial2[nz2]);
            for (std::size_t p2_index = q_index; p2_index <= p2_degree + q_index; p2_index++) {
                copy_poly1[p2_index - p2_degree] -= polynomial2[p2_index - q_index + nz2] * q[q_index];
            }
        }
        return DivType(q, copy_poly1);
    };
    auto result = indirect_div();
    return DivType(result.getQuotient(), result.getRemainder());
}

}

