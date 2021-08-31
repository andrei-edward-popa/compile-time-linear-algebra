#pragma once

#include <array>
#include <tuple>
#include <cassert>

#include "Math.hpp"
#include "Definition.hpp"
#include "CustomFormatter.hpp"

namespace cte::mat {

template<Arithmetic ArithmeticType, std::size_t Rows, std::size_t Cols>
struct Matrix {

    typedef ArithmeticType type;

    static constexpr Matrix sNull = [] -> Matrix {
        Matrix null;
        auto&& zeros = [](const ArithmeticType){ return ArithmeticType{0}; };
        null.apply(zeros);
        return null;
    }();

    static constexpr Matrix sIdentity = [] -> Matrix {
        Matrix identity = sNull;
        auto&& ones = [](const ArithmeticType){ return ArithmeticType{1}; };
        for (std::size_t index = 0; index < Rows; index++) {
            identity.apply(index, index, ones);
        }
        return identity;
    }();

    static constexpr std::size_t getRows()  {
        return Rows;
    }

    static constexpr std::size_t getCols() {
        return Cols;
    }

    constexpr std::array<ArithmeticType, Cols>& operator[](const std::size_t row) {
        assert(row < Rows);
        return mData[row];
    }

    constexpr const std::array<ArithmeticType, Cols>& operator[](const std::size_t row) const {
        assert(row < Rows);
        return mData[row];
    }

    template<typename Callable>
    constexpr void apply(Callable&& callable) {
        for (std::size_t row = 0; row < Rows; row++) {
            for (std::size_t col = 0; col < Cols; col++) {
                mData[row][col] = std::forward<Callable>(callable)(mData[row][col]);
            }
        }
    }

    template<typename Callable>
    constexpr void apply(const std::size_t row, const std::size_t col, Callable&& callable) {
        assert(col < Cols && row < Rows);
        mData[row][col] = std::forward<Callable>(callable)(mData[row][col]);
    }

    constexpr Matrix operator+(const Matrix& rhs) const {
        assert(Cols == rhs.getCols() && Rows == rhs.getRows());
        Matrix result;
        for (std::size_t row = 0; row < Rows; row++) {
            for (std::size_t col = 0; col < Cols; col++) {
                result[row][col] = mData[row][col] + rhs[row][col];
            }
        }
        return result;
    }

    constexpr Matrix operator-(const Matrix& rhs) const {
        assert(Cols == rhs.getCols() && Rows == rhs.getRows());
        Matrix result;
        for (std::size_t row = 0; row < Rows; row++) {
            for (std::size_t col = 0; col < Cols; col++) {
                result[row][col] = mData[row][col] - rhs[row][col];
            }
        }
        return result;
    }

    template<Arithmetic ArgArithmeticType, std::size_t ArgRows, std::size_t ArgCols>
    constexpr auto operator*(const Matrix<ArgArithmeticType, ArgRows, ArgCols>& rhs) const {
        assert(Cols == rhs.getRows());
        Matrix<ArithmeticType, Rows, rhs.getCols()> result;
        for(std::size_t row = 0; row < Rows; row++) {
            for(std::size_t col = 0; col < rhs.getCols(); col++) {
                result[row][col] = ArithmeticType{};
                for(std::size_t index = 0; index < Cols; index++) {
                    result[row][col] += mData[row][index] * rhs[index][col];
                }
            }
        }
        return result;
    }

    constexpr Matrix operator*(const ArithmeticType constant) const {
        Matrix<ArithmeticType, Rows, Cols> result;
        for(std::size_t row = 0; row < Rows; row++) {
            for(std::size_t col = 0; col < Cols; col++) {
                result[row][col] = mData[row][col] * constant;
            }
        }
        return result;
    }

    constexpr Matrix operator/(const ArithmeticType constant) const {
        Matrix<ArithmeticType, Rows, Cols> result;
        for(std::size_t row = 0; row < Rows; row++) {
            for(std::size_t col = 0; col < Cols; col++) {
                result[row][col] = mData[row][col] / constant;
            }
        }
        return result;
    }

    constexpr Matrix operator^(const Matrix& rhs) const {
        assert(Cols == rhs.getCols() && Rows == rhs.getRows());
        Matrix result;
        for(std::size_t row = 0; row < Rows; row++) {
            for(std::size_t col = 0; col < Cols; col++) {
                result[row][col] = mData[row][col] * rhs[row][col];
            }
        }
        return result;
    }

    constexpr Matrix& operator=(const Matrix& rhs) {
        mData = rhs.mData;
        return *this;
    }

    constexpr Matrix& operator=(const auto& rhs) {
        assert(Cols == rhs.getCols() && Rows == rhs.getRows());
        for (std::size_t row = 0; row < Rows; row++) {
            for (std::size_t col = 0; col < Cols; col++) {
                mData[row][col] = static_cast<ArithmeticType>(rhs[row][col]);
            }
        }
        return *this;
    }

    constexpr Matrix& operator=(Matrix&& rhs) {
        mData = std::move(rhs.mData);
        return *this;
    }

    constexpr Matrix() {
        for (std::size_t i = 0; i < Rows; i++) {
            mData[i] = std::array<ArithmeticType, Cols>();
        }
    }

    constexpr Matrix(std::initializer_list<std::initializer_list<ArithmeticType>>&& param) {
        assert(param.size() == Rows);
        for (std::size_t i = 0; i < Rows; i++) {
            assert(std::next(param.begin(), i)->size() == Cols);
            std::copy(std::next(param.begin(), i)->begin(), std::next(param.begin(), i)->end(), std::next(mData.begin(), i)->begin());
        }
    }

    constexpr Matrix(const Matrix& rhs) {
        mData = rhs.mData;
    }

    constexpr Matrix(const auto& rhs) {
        assert(Cols == rhs.getCols() && Rows == rhs.getRows());
        for (std::size_t row = 0; row < Rows; row++) {
            for (std::size_t col = 0; col < Cols; col++) {
                mData[row][col] = static_cast<ArithmeticType>(rhs[row][col]);
            }
        }
    }

    constexpr Matrix(Matrix&& rhs) {
        mData = std::move(rhs.mData);
    }

    constexpr ArithmeticType determinant() const {
        assert(Cols == Rows);
        double pivot, result = 1.0;
        Matrix<double, Rows, Cols> matrix(*this);
        for (std::size_t current = 0; current < Cols; current++) {
            for (std::size_t row = current + 1; row < Rows; row++) {
                pivot = matrix[row][current] / matrix[current][current];
                for (std::size_t col = 0; col < Cols; col++) {
                    matrix[row][col] -= pivot * matrix[current][col];
                }
            }
            result *= matrix[current][current];
        }
        return static_cast<ArithmeticType>(result);
    }

    constexpr Matrix<ArithmeticType, Cols, Rows> transpose() const {
        Matrix<ArithmeticType, Cols, Rows> result;
        for (std::size_t row = 0; row < Rows; row++) {
            for (std::size_t col = 0; col < Cols; col++) {
                result[col][row] = mData[row][col];
            }
        }
        return result;
    }

    constexpr auto inverse() const {
        assert(Rows == Cols);
        assert(this->determinant() != ArithmeticType{0});
        Matrix<double, Rows, Cols> result;
        double pivot;
        auto gaussian = this->concatenateCols(sIdentity);
        Matrix<double, gaussian.getRows(), gaussian.getCols()> matrix(gaussian);
        for (std::size_t current = 0; current < Cols; current++) {
            for (std::size_t row = current + 1; row < Rows; row++) {
                pivot = matrix[row][current] / matrix[current][current];
                for (std::size_t col = 0; col < 2 * Cols; col++) {
                    matrix[row][col] -= pivot * matrix[current][col];
                }
            }
        }
        for (int current = Cols - 1; current >= 0; current--) {
            for (int row = current - 1; row >= 0; row--) {
                pivot = matrix[row][current] / matrix[current][current];
                for (std::size_t col = 0; col < 2 * Cols; col++) {
                    matrix[row][col] -= pivot * matrix[current][col];
                }
            }
        }
        for (std::size_t row = 0; row < Rows; row++) {
            pivot = matrix[row][row];
            for (std::size_t col = Cols; col < 2 * Cols; col++) {
                matrix[row][col] /= pivot;
                result[row][col - Cols] = matrix[row][col];
            }
        }
        return result;
    }

    constexpr auto concatenateRows(const auto& rhs) const {
        assert(Cols == rhs.getCols());
        Matrix<ArithmeticType, Rows + rhs.getRows(), Cols> result;
        for (std::size_t row = 0; row < Rows; row++) {
            for (std::size_t col = 0; col < Cols; col++) {
                result[row][col] = mData[row][col];
            }
        }
        for (std::size_t row = 0; row < rhs.getRows(); row++) {
            for (std::size_t col = 0; col < rhs.getCols(); col++) {
                result[Rows + row][col] = rhs[row][col];
            }
        }
        return result;
    }

    constexpr auto concatenateCols(const auto& rhs) const {
        assert(Rows == rhs.getRows());
        Matrix<ArithmeticType, Rows, Cols + rhs.getCols()> result;
        for (std::size_t row = 0; row < Rows; row++) {
            for (std::size_t col = 0; col < Cols; col++) {
                result[row][col] = mData[row][col];
            }
        }
        for (std::size_t row = 0; row < rhs.getRows(); row++) {
            for (std::size_t col = 0; col < rhs.getCols(); col++) {
                result[row][Cols + col] = rhs[row][col];
            }
        }
        return result;
    }

    constexpr double norm() const {
        assert(Rows == 1 || Cols == 1);
        double result = 0;
        for (std::size_t row = 0; row < Rows; row++) {
            for (std::size_t col = 0; col < Cols; col++) {
                result += cte::math::pow(mData[row][col], 2);
            }
        }
        return cte::math::sqrt(result);
    }

    constexpr std::array<double, Rows> eigenvalues() const {
        assert(Rows == Cols);
        auto qr = QRDecomposition();
        Matrix<double, Rows, Cols> temp;
        std::array<double, Rows> result;
        for (std::size_t index = 0; index < 10; index++) {
            temp = qr[1] * qr[0];
            qr = temp.QRDecomposition();
        }
        temp = qr[0] * qr[1];
        for (std::size_t index = 0; index < Rows; index++) {
            result[index] = temp[index][index];
        }
        return result;
    }

    constexpr std::array<Matrix<double, Rows, 1>, Cols> eigenvectors() const {
        assert(Rows == Cols);
        auto qr = QRDecomposition();
        Matrix<double, Rows, Cols> temp;
        Matrix<double, Rows, Cols> result = qr[0];
        for (std::size_t index = 0; index < 10; index++) {
            temp = qr[1] * qr[0];
            qr = temp.QRDecomposition();
            result = result * qr[0];
        }
        return result.splitCols();
    }

    constexpr std::array<Matrix<ArithmeticType, 1, Cols>, Rows> splitRows() const {
        std::array<Matrix<ArithmeticType, 1, Cols>, Rows> result;
        for (std::size_t row = 0; row < Rows; row++) {
            for (std::size_t col = 0; col < Cols; col++) {
                result[row][0][col] = mData[row][col]; 
            }
        }
        return result;
    }

    constexpr std::array<Matrix<ArithmeticType, Rows, 1>, Cols> splitCols() const {
        std::array<Matrix<ArithmeticType, Rows, 1>, Cols> result;
        for (std::size_t row = 0; row < Rows; row++) {
            for (std::size_t col = 0; col < Cols; col++) {
                result[col][row][0] = mData[row][col]; 
            }
        }
        return result;
    }

    constexpr auto QRDecomposition() const;

public:
    std::array<std::array<ArithmeticType, Cols>, Rows> mData;
};

template<Arithmetic AT1, std::size_t R1, std::size_t C1, Arithmetic AT2, std::size_t R2, std::size_t C2>
struct QRType {
    constexpr QRType(const auto& q, const auto& r) : mOrthogonal(q), mUpperTriangular(r) {

    }
    constexpr const auto& getOrthogonal() const {
        return mOrthogonal;
    }
    constexpr const auto& getUpperTriangular() const {
        return mUpperTriangular;
    }
    constexpr auto& getOrthogonal() {
        return mOrthogonal;
    }
    constexpr auto& getUpperTriangular() {
        return mUpperTriangular;
    }
private:
    Matrix<AT1, R1, C1> mOrthogonal;
    Matrix<AT2, R2, C2> mUpperTriangular;
};

template<Arithmetic ArithmeticType, std::size_t Rows, std::size_t Cols>
constexpr auto Matrix<ArithmeticType, Rows, Cols>::QRDecomposition() const {
    assert(Rows == Cols);
    std::array<Matrix<ArithmeticType, Rows, 1>, Cols> columns = splitCols();
    std::array<Matrix<double, Rows, 1>, Cols> new_columns;
    Matrix<double, Rows, Cols> Q;
    double coeff;
    for (std::size_t col = 0; col < Cols; col++) {
        new_columns[col] = columns[col];
        for (std::size_t index = 0; index < col; index++) {
            coeff = (new_columns[index].transpose() * columns[col])[0][0] / (new_columns[index].transpose() * new_columns[index])[0][0];
            new_columns[col] = new_columns[col] - new_columns[index] * coeff;
        }
    }
    for (std::size_t col = 0; col < Cols; col++) {
        new_columns[col] = new_columns[col] / new_columns[col].norm();
    }
    for (std::size_t row = 0; row < Rows; row++) {
        for (std::size_t col = 0; col < Cols; col++) {
            Q[row][col] = new_columns[col][row][0];
        }
    }
    Matrix<double, Rows, Cols> R = Q.transpose() * (*this);
    return QRType(Q, R);
}

template<std::size_t Dim>
constexpr auto reductionBareissAlgorithmMapping() {
    constexpr std::size_t size = (Dim - 1) * Dim * (2 * Dim - 1) / 6;
    std::size_t v_index = 0;
    std::array<std::tuple<int, int, int, int>, size> values;
    int v1, v2, v3, v4, t1, t2 = Dim + 1;
    for (int dim = Dim - 1; dim >= 1; dim--) {
        t1 = Dim - 1 - dim;
        for (int index = 0; index < dim * dim; index++) {
            v1 = t1 * t2;
            v2 = index % dim + 1 + v1;
            v3 = (index / dim + 1) * Dim + v1;
            v4 = v2 + v3 - v1;
            values[v_index++] = std::tuple{ v1, v2, v3, v4 };
        }
    }
    return values;
}

template<Arithmetic ArithmeticType, std::size_t Dim>
constexpr auto matrixPolynomialComposition(cte::mat::Matrix<ArithmeticType, Dim, Dim> mat) {

    std::array<cte::poly::Polynomial<ArithmeticType, 2 * (Dim - 1)>, Dim * Dim> polynomials{};

    for (std::size_t row = 0; row < Dim; row++) {
        for (std::size_t col = 0; col < Dim; col++) {
            polynomials[row * Dim + col] = cte::poly::Polynomial<ArithmeticType, 2 * (Dim - 1)>{};
            polynomials[row * Dim + col][2 * (Dim - 1)] = mat[row][col];
            polynomials[row * Dim + col][2 * (Dim - 1) - 1] = (row * Dim + col) % (Dim + 1) == 0 ? ArithmeticType{-1} : ArithmeticType{0};
        }
    }

    return polynomials;

}

template<auto matrix>
constexpr auto reductionBareissAlgorithm() {
    constexpr std::size_t Dim = matrix.getRows();
    constexpr auto mapping = reductionBareissAlgorithmMapping<Dim>();
    auto polynomials = matrixPolynomialComposition(matrix);
    auto pivot = cte::poly::Polynomial<typename decltype(matrix)::type, 2 * (Dim - 1)>{};
    pivot[2 * (Dim - 1)] = 1.0l;
    std::size_t index = 0;
    for (int power = Dim - 1; power >= 1; power--) {
        for (std::size_t i = 0; i < power * power; i++) {
            polynomials[std::get<3>(mapping[index])] = cte::poly::div_preserve_dim(cte::poly::sub_preserve_dim(cte::poly::mult_preserve_dim(polynomials[std::get<0>(mapping[index])], polynomials[std::get<3>(mapping[index])]), cte::poly::mult_preserve_dim(polynomials[std::get<1>(mapping[index])], polynomials[std::get<2>(mapping[index])])), pivot).getQuotient();
            index++;
        }
        pivot = polynomials[(Dim - 1 - power) * (Dim + 1)];
    }
    return polynomials[Dim * Dim - 1];
}

template<auto matrix>
constexpr auto getCharacteristicPolynomial() {
    return cte::poly::reduceDegree<reductionBareissAlgorithm<matrix>()>();
}

template<auto matrix>
constexpr auto calculateEigenvalues() {
    constexpr auto characteristic = getCharacteristicPolynomial<matrix>();
    constexpr auto eigenvalues = characteristic.roots();
    return eigenvalues;
}

}

