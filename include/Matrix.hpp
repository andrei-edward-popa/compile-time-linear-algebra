#pragma once

#include <array>
#include <tuple>
#include <cassert>

#include "Math.hpp"
#include "Exception.hpp"
#include "Polynomial.hpp"
#include "Definition.hpp"
#include "CustomFormatter.hpp"

namespace cte::mat {

template<Arithmetic ArithmeticType, std::size_t Rows, std::size_t Cols>
struct Matrix {

    typedef ArithmeticType type;

    constexpr Matrix() {
        std::array<ArithmeticType, Cols> arr;
        arr.fill(ArithmeticType{0});
        for (std::size_t i = 0; i < Rows; i++) {
            mData[i] = arr;
        }
    }

    constexpr Matrix(std::initializer_list<std::initializer_list<ArithmeticType>>&& param) {
        for (std::size_t i = 0; i < Rows; i++) {
            std::copy(std::next(param.begin(), i)->begin(), std::next(param.begin(), i)->end(), std::next(mData.begin(), i)->begin());
        }
    }

    constexpr Matrix(const Matrix& rhs) {
        mData = rhs.mData;
    }

    constexpr Matrix(Matrix&& rhs) {
        mData = std::move(rhs.mData);
    }

    constexpr Matrix& operator=(const Matrix& rhs) {
        mData = rhs.mData;
        return *this;
    }

    constexpr Matrix& operator=(Matrix&& rhs) {
        mData = std::move(rhs.mData);
        return *this;
    }

    constexpr std::array<ArithmeticType, Cols>& operator[](const std::size_t row) {
        if (row >= Rows) throw cte::OutOfRangeException(row);
        return mData[row];
    }

    constexpr const std::array<ArithmeticType, Cols>& operator[](const std::size_t row) const {
        if (row >= Rows) throw cte::OutOfRangeException(row);
        return mData[row];
    }

    template<Arithmetic OtherArithmeticType, std::size_t OtherRows, std::size_t OtherCols>
    constexpr operator Matrix<OtherArithmeticType, OtherRows, OtherCols>() const & {
        if (Rows != OtherRows) throw cte::DimensionMismatchException("Invalid assign, matrices have different number of rows", Rows, OtherRows);
        if (Cols != OtherCols) throw cte::DimensionMismatchException("Invalid assign, matrices have different number of columns", Cols, OtherCols);
        Matrix<std::common_type_t<ArithmeticType, OtherArithmeticType>, Rows, Cols> result;
        for (std::size_t row = 0; row < Rows; row++) {
            for (std::size_t col = 0; col < Cols; col++) {
                result[row][col] = static_cast<OtherArithmeticType>(mData[row][col]);
            }
        }
        return result;
    }

    template<Arithmetic OtherArithmeticType, std::size_t OtherRows, std::size_t OtherCols>
    constexpr operator Matrix<OtherArithmeticType, OtherRows, OtherCols>() && {
        if (Rows != OtherRows) throw cte::DimensionMismatchException("Invalid assign, matrices have different number of rows", Rows, OtherRows);
        if (Cols != OtherCols) throw cte::DimensionMismatchException("Invalid assign, matrices have different number of columns", Cols, OtherCols);
        Matrix<std::common_type_t<ArithmeticType, OtherArithmeticType>, Rows, Cols> result;
        for (std::size_t row = 0; row < Rows; row++) {
            for (std::size_t col = 0; col < Cols; col++) {
                result[row][col] = static_cast<OtherArithmeticType>(std::move(mData[row][col]));
            }
        }
        return result;
    }

    static constexpr std::size_t getRows()  {
        return Rows;
    }

    static constexpr std::size_t getCols() {
        return Cols;
    }

    constexpr auto operator+(const auto& rhs) const -> Matrix<std::common_type_t<ArithmeticType, typename std::remove_cvref_t<decltype(rhs)>::type>, Rows, Cols>;
    constexpr auto operator-(const auto& rhs) const -> Matrix<std::common_type_t<ArithmeticType, typename std::remove_cvref_t<decltype(rhs)>::type>, Rows, Cols>;
    constexpr auto operator*(const auto& rhs) const -> Matrix<std::common_type_t<ArithmeticType, typename std::remove_cvref_t<decltype(rhs)>::type>, Rows, std::remove_cvref_t<decltype(rhs)>::getCols()>;
    constexpr auto operator*(const ArithmeticType& constant) const noexcept -> Matrix<ArithmeticType, Rows, Cols>;
    constexpr auto operator/(const ArithmeticType& constant) const noexcept -> Matrix<ArithmeticType, Rows, Cols>;
    constexpr auto operator^(const auto& rhs) const -> Matrix<std::common_type_t<ArithmeticType, typename std::remove_cvref_t<decltype(rhs)>::type>, Rows, Cols>;

    constexpr void apply(auto&& callable);
    constexpr void apply(const std::size_t row, const std::size_t col, auto&& callable);
    constexpr long double norm() const;
    constexpr auto gaussDecomposition() const noexcept -> Matrix<long double, Rows, Cols>;
    constexpr ArithmeticType determinant() const;
    constexpr Matrix<long double, Rows, Cols> inverse() const;
    constexpr Matrix<ArithmeticType, Cols, Rows> transpose() const noexcept;
    constexpr auto concatenateRows(const auto& rhs) const -> decltype(auto);
    constexpr auto concatenateCols(const auto& rhs) const -> decltype(auto);
    constexpr std::array<Matrix<ArithmeticType, 1, Cols>, Rows> splitRows() const noexcept;
    constexpr std::array<Matrix<ArithmeticType, Rows, 1>, Cols> splitCols() const noexcept;
    constexpr auto QRDecomposition() const -> QRType<long double, Rows, Cols, long double, Cols, Cols>;

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
constexpr auto Matrix<ArithmeticType, Rows, Cols>::operator+(const auto& rhs) const -> Matrix<std::common_type_t<ArithmeticType, typename std::remove_cvref_t<decltype(rhs)>::type>, Rows, Cols> {
    using remove_cvref_rhs = typename std::remove_cvref_t<decltype(rhs)>::type;
    if (Rows != rhs.getRows()) throw cte::DimensionMismatchException("Invalid addition, matrices have different number of rows", Rows, rhs.getRows());
    if (Cols != rhs.getCols()) throw cte::DimensionMismatchException("Invalid addition, matrices have different number of columns", Cols, rhs.getCols());
    Matrix<std::common_type_t<ArithmeticType, remove_cvref_rhs>, Rows, Cols> result;
    for (std::size_t row = 0; row < Rows; row++) {
        for (std::size_t col = 0; col < Cols; col++) {
            result[row][col] = mData[row][col] + rhs[row][col];
        }
    }
    return result;
}

template<Arithmetic ArithmeticType, std::size_t Rows, std::size_t Cols>
constexpr auto Matrix<ArithmeticType, Rows, Cols>::operator-(const auto& rhs) const -> Matrix<std::common_type_t<ArithmeticType, typename std::remove_cvref_t<decltype(rhs)>::type>, Rows, Cols> {
    using remove_cvref_rhs = typename std::remove_cvref_t<decltype(rhs)>::type;
    if (Rows != rhs.getRows()) throw cte::DimensionMismatchException("Invalid subtraction, matrices have different number of rows", Rows, rhs.getRows());
    if (Cols != rhs.getCols()) throw cte::DimensionMismatchException("Invalid subtraction, matrices have different number of columns", Cols, rhs.getCols());
    Matrix<std::common_type_t<ArithmeticType, remove_cvref_rhs>, Rows, Cols> result;
    for (std::size_t row = 0; row < Rows; row++) {
        for (std::size_t col = 0; col < Cols; col++) {
            result[row][col] = mData[row][col] - rhs[row][col];
        }
    }
    return result;
}

template<Arithmetic ArithmeticType, std::size_t Rows, std::size_t Cols>
constexpr auto Matrix<ArithmeticType, Rows, Cols>::operator*(const auto& rhs) const -> Matrix<std::common_type_t<ArithmeticType, typename std::remove_cvref_t<decltype(rhs)>::type>, Rows, std::remove_cvref_t<decltype(rhs)>::getCols()> {
    if (Cols != rhs.getRows()) throw cte::DimensionMismatchException("Invalid multiplication, first matrix has different number of columns than number of rows of second matrix", Cols, rhs.getRows());
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

template<Arithmetic ArithmeticType, std::size_t Rows, std::size_t Cols>
constexpr auto Matrix<ArithmeticType, Rows, Cols>::operator*(const ArithmeticType& constant) const noexcept -> Matrix<ArithmeticType, Rows, Cols> {
    Matrix<std::common_type_t<ArithmeticType, decltype(constant)>, Rows, Cols> result;
    for(std::size_t row = 0; row < Rows; row++) {
        for(std::size_t col = 0; col < Cols; col++) {
            result[row][col] = mData[row][col] * constant;
        }
    }
    return result;
}

template<Arithmetic ArithmeticType, std::size_t Rows, std::size_t Cols>
constexpr auto Matrix<ArithmeticType, Rows, Cols>::operator/(const ArithmeticType& constant) const noexcept -> Matrix<ArithmeticType, Rows, Cols> {
    Matrix<std::common_type_t<ArithmeticType, decltype(constant)>, Rows, Cols> result;
    for(std::size_t row = 0; row < Rows; row++) {
        for(std::size_t col = 0; col < Cols; col++) {
            result[row][col] = mData[row][col] / constant;
        }
    }
    return result;
}

template<Arithmetic ArithmeticType, std::size_t Rows, std::size_t Cols>
constexpr auto Matrix<ArithmeticType, Rows, Cols>::operator^(const auto& rhs) const -> Matrix<std::common_type_t<ArithmeticType, typename std::remove_cvref_t<decltype(rhs)>::type>, Rows, Cols> {
    using remove_cvref_rhs = typename std::remove_cvref_t<decltype(rhs)>::type;
    if (Rows != rhs.getRows()) throw cte::DimensionMismatchException("Invalid Hadamard multiplication, matrices have different number of rows", Rows, rhs.getRows());
    if (Cols != rhs.getCols()) throw cte::DimensionMismatchException("Invalid Hadamard multiplication, matrices have different number of columns", Cols, rhs.getCols());
    Matrix<std::common_type_t<ArithmeticType, remove_cvref_rhs>, Rows, Cols> result;
    for(std::size_t row = 0; row < Rows; row++) {
        for(std::size_t col = 0; col < Cols; col++) {
            result[row][col] = mData[row][col] * rhs[row][col];
        }
    }
    return result;
}

template<Arithmetic ArithmeticType, std::size_t Rows, std::size_t Cols>
constexpr void Matrix<ArithmeticType, Rows, Cols>::apply(auto&& callable) {
    for (std::size_t row = 0; row < Rows; row++) {
        for (std::size_t col = 0; col < Cols; col++) {
            mData[row][col] = std::forward<std::remove_cvref_t<decltype(callable)>>(callable)(mData[row][col]);
        }
    }
}

template<Arithmetic ArithmeticType, std::size_t Rows, std::size_t Cols>
constexpr void Matrix<ArithmeticType, Rows, Cols>::apply(const std::size_t row, const std::size_t col, auto&& callable) {
    if (row >= Rows) throw cte::OutOfRangeException(row);   
    if (col >= Cols) throw cte::OutOfRangeException(col);
    mData[row][col] = std::forward<std::remove_cvref_t<decltype(callable)>>(callable)(mData[row][col]);
}

template<Arithmetic ArithmeticType, std::size_t Rows, std::size_t Cols>
constexpr long double Matrix<ArithmeticType, Rows, Cols>::norm() const {
    if (Rows != 1UL && Cols != 1UL) throw cte::DimensionMismatchException("Invalid norm, row or column matrix needed", Rows, Cols);
    long double result = 0.0L;
    for (std::size_t row = 0; row < Rows; row++) {
        for (std::size_t col = 0; col < Cols; col++) {
            result += cte::math::pow(mData[row][col], 2);
        }
    }
    return cte::math::sqrt(result);
}

template<Arithmetic ArithmeticType, std::size_t Rows, std::size_t Cols>
constexpr auto Matrix<ArithmeticType, Rows, Cols>::gaussDecomposition() const noexcept -> Matrix<long double, Rows, Cols> {
    Matrix<long double, Rows, Cols> matrix(*this);
    long double pivot;
    for (std::size_t current = 0; current < Cols; current++) {
        for (std::size_t row = current + 1; row < Rows; row++) {
            pivot = matrix[row][current] / matrix[current][current];
            for (std::size_t col = 0; col < Cols; col++) {
                matrix[row][col] -= pivot * matrix[current][col];
            }
        }
    }
    return matrix;
}

template<Arithmetic ArithmeticType, std::size_t Rows, std::size_t Cols>
constexpr ArithmeticType Matrix<ArithmeticType, Rows, Cols>::determinant() const {
    if (Rows != Cols) throw cte::DimensionMismatchException("Invalid determinant, matrix has different dimensions", Rows, Cols);
    long double result = 1.0L;
    Matrix<long double, Rows, Cols> matrix = gaussDecomposition();
    for (std::size_t current = 0; current < Cols; current++) {
        result *= matrix[current][current];
    }
    return static_cast<ArithmeticType>(result);
}

template<Arithmetic ArithmeticType, std::size_t Rows, std::size_t Cols>
constexpr Matrix<long double, Rows, Cols> Matrix<ArithmeticType, Rows, Cols>::inverse() const {
    if (Rows != Cols) throw cte::DimensionMismatchException("Invalid inverse, matrix has different dimensions", Rows, Cols);
    if (this->determinant() == ArithmeticType{0}) throw cte::DimensionMismatchException("Invalid inverse, matrix is singular", Rows, Cols);
    long double pivot;
    Matrix<long double, Rows, Cols> result;
    
    auto&& zeros = [](const ArithmeticType){ return ArithmeticType{0}; };
    auto&& ones = [](const ArithmeticType){ return ArithmeticType{1}; };

    Matrix identity;
    identity.apply(zeros);
    for (std::size_t index = 0; index < Rows; index++) {
        identity.apply(index, index, ones);
    }
    auto gaussian = this->concatenateCols(identity);
    Matrix<long double, gaussian.getRows(), gaussian.getCols()> matrix(gaussian);
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

template<Arithmetic ArithmeticType, std::size_t Rows, std::size_t Cols>
constexpr Matrix<ArithmeticType, Cols, Rows> Matrix<ArithmeticType, Rows, Cols>::transpose() const noexcept {
    Matrix<ArithmeticType, Cols, Rows> result;
    for (std::size_t row = 0; row < Rows; row++) {
        for (std::size_t col = 0; col < Cols; col++) {
            result[col][row] = mData[row][col];
        }
    }
    return result;
}

template<Arithmetic ArithmeticType, std::size_t Rows, std::size_t Cols>
constexpr auto Matrix<ArithmeticType, Rows, Cols>::concatenateRows(const auto& rhs) const -> decltype(auto) {
    using remove_cvref_rhs = std::remove_cvref_t<decltype(rhs)>;
    if (Cols != remove_cvref_rhs::getCols()) throw cte::DimensionMismatchException("Invalid rows concatenation, matrices have different number of columns", Cols, remove_cvref_rhs::getCols());
    Matrix<ArithmeticType, Rows + remove_cvref_rhs::getCols(), Cols> result;
    for (std::size_t row = 0; row < Rows; row++) {
        for (std::size_t col = 0; col < Cols; col++) {
            result[row][col] = mData[row][col];
        }
    }
    for (std::size_t row = 0; row < remove_cvref_rhs::getRows(); row++) {
        for (std::size_t col = 0; col < remove_cvref_rhs::getCols(); col++) {
            result[Rows + row][col] = rhs[row][col];
        }
    }
    return result;
}

template<Arithmetic ArithmeticType, std::size_t Rows, std::size_t Cols>
constexpr auto Matrix<ArithmeticType, Rows, Cols>::concatenateCols(const auto& rhs) const -> decltype(auto) {
    using remove_cvref_rhs = std::remove_cvref_t<decltype(rhs)>;
    if (Rows != remove_cvref_rhs::getRows()) throw cte::DimensionMismatchException("Invalid columns concatenation, matrices have different number of rows", Rows, remove_cvref_rhs::getRows());
    Matrix<ArithmeticType, Rows, Cols + remove_cvref_rhs::getCols()> result;
    for (std::size_t row = 0; row < Rows; row++) {
        for (std::size_t col = 0; col < Cols; col++) {
            result[row][col] = mData[row][col];
        }
    }
    for (std::size_t row = 0; row < remove_cvref_rhs::getRows(); row++) {
        for (std::size_t col = 0; col < remove_cvref_rhs::getCols(); col++) {
            result[row][Cols + col] = rhs[row][col];
        }
    }
    return result;
}

template<Arithmetic ArithmeticType, std::size_t Rows, std::size_t Cols>
constexpr std::array<Matrix<ArithmeticType, 1, Cols>, Rows> Matrix<ArithmeticType, Rows, Cols>::splitRows() const noexcept {
    std::array<Matrix<ArithmeticType, 1, Cols>, Rows> result;
    for (std::size_t row = 0; row < Rows; row++) {
        for (std::size_t col = 0; col < Cols; col++) {
            result[row][0][col] = mData[row][col]; 
        }
    }
    return result;
}

template<Arithmetic ArithmeticType, std::size_t Rows, std::size_t Cols>
constexpr std::array<Matrix<ArithmeticType, Rows, 1>, Cols> Matrix<ArithmeticType, Rows, Cols>::splitCols() const noexcept {
    std::array<Matrix<ArithmeticType, Rows, 1>, Cols> result;
    for (std::size_t row = 0; row < Rows; row++) {
        for (std::size_t col = 0; col < Cols; col++) {
            result[col][row][0] = mData[row][col]; 
        }
    }
    return result;
}

template<Arithmetic ArithmeticType, std::size_t Rows, std::size_t Cols>
constexpr auto Matrix<ArithmeticType, Rows, Cols>::QRDecomposition() const -> QRType<long double, Rows, Cols, long double, Cols, Cols> {
    if (Rows != Cols) throw cte::DimensionMismatchException("Invalid QR decomposition, matrix has different dimensions", Rows, Cols);
    std::array<Matrix<ArithmeticType, Rows, 1>, Cols> columns = splitCols();
    std::array<Matrix<long double, Rows, 1>, Cols> new_columns;
    Matrix<long double, Rows, Cols> Q;
    long double coeff;
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
    Matrix<long double, Cols, Cols> R = Q.transpose() * (*this);
    return QRType(Q, R);
}

template<std::size_t N>
Matrix<uint16_t, N, N - 1> linearEquationMapping() {
    Matrix<uint16_t, N, N - 1> mapping;
    uint16_t up = N - 1, down = 1;
    for (uint16_t col = 0; col < N - 1; col++) {
        for (uint16_t row = 0; row < up; row++) {
            mapping[row][col] = col + 1;
        }
        for (uint16_t row = up; row < up + down; row++) {
            mapping[row][col] = col + 2;
        }
        up--;
        down++;
    }
    return mapping;
}

template<Arithmetic ArithmeticType, std::size_t Size>
constexpr std::array<long double, Size> solveDeterminateSystem(const Matrix<ArithmeticType, Size, Size>& coeffs, const std::array<ArithmeticType, Size> freeTerms) {
    ArithmeticType determ = coeffs.determinant();
    if (determ == ArithmeticType{0}) throw cte::DimensionMismatchException("System is not determinate because determinant is equal to 0", Size, Size);
    std::array<long double, Size> solutions;
    Matrix<ArithmeticType, Size, Size> new_matrix;
    for (std::size_t index = 0; index < Size; index++) {
        for (std::size_t row = 0; row < Size; row++) {
            for (std::size_t col = 0; col < Size; col++) {
                if (index == col) {
                    new_matrix[row][col] = freeTerms[row];
                } else {
                    new_matrix[row][col] = coeffs[row][col];
                }
            }
        }
        solutions[index] = new_matrix.determinant() / static_cast<long double>(determ);
    }
    return solutions;
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
            polynomials[std::get<3>(mapping[index])] = ((polynomials[std::get<0>(mapping[index])] * polynomials[std::get<3>(mapping[index])] - polynomials[std::get<1>(mapping[index])] * polynomials[std::get<2>(mapping[index])]) / pivot).getQuotient();
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

template<auto matrix>
constexpr auto calculateEigenvectors() {
    constexpr auto eigenvalues = calculateEigenvalues<matrix>();
    Matrix new_matrix = matrix;
    Matrix<long double, matrix.getRows() - 1, matrix.getCols() - 1> egv;
    std::array<long double, matrix.getRows() - 1> freeTerms;
    Matrix<long double, matrix.getRows(), matrix.getCols()> eigenvectors;
    for (std::size_t index = 0; index < eigenvalues.size(); index++) {
        for (std::size_t pr = 0; pr < matrix.getRows(); pr++) {
            new_matrix[pr][pr] = matrix[pr][pr] - eigenvalues[index].real();
        }
        auto gaussDecomp = new_matrix.gaussDecomposition();
        for (std::size_t row = 0; row < matrix.getRows() - 1; row++) {
            for (std::size_t col = 0; col < matrix.getCols() - 1; col++) {
                egv[row][col] = gaussDecomp[row][col];
            }
        }
        for (std::size_t i = 0; i < matrix.getRows() - 1; i++) {
            freeTerms[i] = -gaussDecomp[i][matrix.getCols() - 1];
        }
        auto solutions = solveDeterminateSystem(egv, freeTerms);
        for (std::size_t row = 0; row < freeTerms.size(); row++) {
            eigenvectors[row][index] = solutions[row];
        }
        eigenvectors[matrix.getRows() - 1][index] = 1.0l;
    }
    return eigenvectors;
}

}

