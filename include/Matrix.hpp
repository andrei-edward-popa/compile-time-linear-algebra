#pragma once

#include <array>
#include <tuple>
#include <cassert>
#include <complex>

#include "Math.hpp"
#include "Exception.hpp"
#include "Polynomial.hpp"
#include "Definition.hpp"
#include "CustomFormatter.hpp"

namespace cte::mat {

template<FloatingPoint FloatingPointType, std::size_t Rows, std::size_t Cols>
struct Matrix {

    typedef FloatingPointType type;

    constexpr Matrix() {
        std::array<FloatingPointType, Cols> arr;
        arr.fill(FloatingPointType{0});
        for (std::size_t i = 0; i < Rows; i++) {
            mData[i] = arr;
        }
    }

    constexpr Matrix(std::initializer_list<std::initializer_list<FloatingPointType>>&& param) {
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

    constexpr std::array<FloatingPointType, Cols>& operator[](const std::size_t row) {
        if (row >= Rows) throw cte::OutOfRangeException(row);
        return mData[row];
    }

    constexpr const std::array<FloatingPointType, Cols>& operator[](const std::size_t row) const {
        if (row >= Rows) throw cte::OutOfRangeException(row);
        return mData[row];
    }

    template<FloatingPoint OtherNumberType, std::size_t OtherRows, std::size_t OtherCols>
    constexpr operator Matrix<OtherNumberType, OtherRows, OtherCols>() const & {
        if (Rows != OtherRows) throw cte::DimensionMismatchException("Invalid assign, matrices have different FloatingPoint of rows", Rows, OtherRows);
        if (Cols != OtherCols) throw cte::DimensionMismatchException("Invalid assign, matrices have different FloatingPoint of columns", Cols, OtherCols);
        Matrix<std::common_type_t<FloatingPointType, OtherNumberType>, Rows, Cols> result;
        for (std::size_t row = 0; row < Rows; row++) {
            for (std::size_t col = 0; col < Cols; col++) {
                result[row][col] = static_cast<OtherNumberType>(mData[row][col]);
            }
        }
        return result;
    }

    template<FloatingPoint OtherNumberType, std::size_t OtherRows, std::size_t OtherCols>
    constexpr operator Matrix<OtherNumberType, OtherRows, OtherCols>() && {
        if (Rows != OtherRows) throw cte::DimensionMismatchException("Invalid assign, matrices have different FloatingPoint of rows", Rows, OtherRows);
        if (Cols != OtherCols) throw cte::DimensionMismatchException("Invalid assign, matrices have different FloatingPoint of columns", Cols, OtherCols);
        Matrix<std::common_type_t<FloatingPointType, OtherNumberType>, Rows, Cols> result;
        for (std::size_t row = 0; row < Rows; row++) {
            for (std::size_t col = 0; col < Cols; col++) {
                result[row][col] = static_cast<OtherNumberType>(std::move(mData[row][col]));
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

    constexpr auto operator+(const auto& rhs) const -> Matrix<std::common_type_t<FloatingPointType, typename std::remove_cvref_t<decltype(rhs)>::type>, Rows, Cols>;
    constexpr auto operator-(const auto& rhs) const -> Matrix<std::common_type_t<FloatingPointType, typename std::remove_cvref_t<decltype(rhs)>::type>, Rows, Cols>;
    constexpr auto operator*(const auto& rhs) const -> Matrix<std::common_type_t<FloatingPointType, typename std::remove_cvref_t<decltype(rhs)>::type>, Rows, std::remove_cvref_t<decltype(rhs)>::getCols()>;
    constexpr auto operator*(const FloatingPointType& constant) const noexcept -> Matrix<FloatingPointType, Rows, Cols>;
    constexpr auto operator/(const FloatingPointType& constant) const noexcept -> Matrix<FloatingPointType, Rows, Cols>;
    constexpr auto operator^(const auto& rhs) const -> Matrix<std::common_type_t<FloatingPointType, typename std::remove_cvref_t<decltype(rhs)>::type>, Rows, Cols>;

    constexpr void apply(auto&& callable);
    constexpr void apply(const std::size_t row, const std::size_t col, auto&& callable);
    constexpr FloatingPointType norm() const;
    constexpr auto gaussDecomposition() const noexcept -> Matrix<FloatingPointType, Rows, Cols>;
    constexpr FloatingPointType determinant() const;
    constexpr Matrix<FloatingPointType, Rows, Cols> inverse() const;
    constexpr Matrix<FloatingPointType, Cols, Rows> transpose() const noexcept;
    constexpr auto concatenateRows(const auto& rhs) const -> decltype(auto);
    constexpr auto concatenateCols(const auto& rhs) const -> decltype(auto);
    constexpr std::array<Matrix<FloatingPointType, 1, Cols>, Rows> splitRows() const noexcept;
    constexpr std::array<Matrix<FloatingPointType, Rows, 1>, Cols> splitCols() const noexcept;
    constexpr auto QRDecomposition() const -> QRType<FloatingPointType, Rows, Cols, FloatingPointType, Cols, Cols>;

public:
    std::array<std::array<FloatingPointType, Cols>, Rows> mData;
};

template<FloatingPoint NT1, std::size_t R1, std::size_t C1, FloatingPoint NT2, std::size_t R2, std::size_t C2>
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
    Matrix<NT1, R1, C1> mOrthogonal;
    Matrix<NT2, R2, C2> mUpperTriangular;
};

template<FloatingPoint FloatingPointType, std::size_t Rows, std::size_t Cols>
constexpr auto Matrix<FloatingPointType, Rows, Cols>::operator+(const auto& rhs) const -> Matrix<std::common_type_t<FloatingPointType, typename std::remove_cvref_t<decltype(rhs)>::type>, Rows, Cols> {
    using remove_cvref_rhs = typename std::remove_cvref_t<decltype(rhs)>::type;
    if (Rows != rhs.getRows()) throw cte::DimensionMismatchException("Invalid addition, matrices have different FloatingPoint of rows", Rows, rhs.getRows());
    if (Cols != rhs.getCols()) throw cte::DimensionMismatchException("Invalid addition, matrices have different FloatingPoint of columns", Cols, rhs.getCols());
    Matrix<std::common_type_t<FloatingPointType, remove_cvref_rhs>, Rows, Cols> result;
    for (std::size_t row = 0; row < Rows; row++) {
        for (std::size_t col = 0; col < Cols; col++) {
            result[row][col] = mData[row][col] + rhs[row][col];
        }
    }
    return result;
}

template<FloatingPoint FloatingPointType, std::size_t Rows, std::size_t Cols>
constexpr auto Matrix<FloatingPointType, Rows, Cols>::operator-(const auto& rhs) const -> Matrix<std::common_type_t<FloatingPointType, typename std::remove_cvref_t<decltype(rhs)>::type>, Rows, Cols> {
    using remove_cvref_rhs = typename std::remove_cvref_t<decltype(rhs)>::type;
    if (Rows != rhs.getRows()) throw cte::DimensionMismatchException("Invalid subtraction, matrices have different FloatingPoint of rows", Rows, rhs.getRows());
    if (Cols != rhs.getCols()) throw cte::DimensionMismatchException("Invalid subtraction, matrices have different FloatingPoint of columns", Cols, rhs.getCols());
    Matrix<std::common_type_t<FloatingPointType, remove_cvref_rhs>, Rows, Cols> result;
    for (std::size_t row = 0; row < Rows; row++) {
        for (std::size_t col = 0; col < Cols; col++) {
            result[row][col] = mData[row][col] - rhs[row][col];
        }
    }
    return result;
}

template<FloatingPoint FloatingPointType, std::size_t Rows, std::size_t Cols>
constexpr auto Matrix<FloatingPointType, Rows, Cols>::operator*(const auto& rhs) const -> Matrix<std::common_type_t<FloatingPointType, typename std::remove_cvref_t<decltype(rhs)>::type>, Rows, std::remove_cvref_t<decltype(rhs)>::getCols()> {
    if (Cols != rhs.getRows()) throw cte::DimensionMismatchException("Invalid multiplication, first matrix has different FloatingPoint of columns than FloatingPoint of rows of second matrix", Cols, rhs.getRows());
    Matrix<FloatingPointType, Rows, rhs.getCols()> result;
    for(std::size_t row = 0; row < Rows; row++) {
        for(std::size_t col = 0; col < rhs.getCols(); col++) {
            result[row][col] = FloatingPointType{};
            for(std::size_t index = 0; index < Cols; index++) {
                result[row][col] += mData[row][index] * rhs[index][col];
            }
        }
    }
    return result;
}

template<FloatingPoint FloatingPointType, std::size_t Rows, std::size_t Cols>
constexpr auto Matrix<FloatingPointType, Rows, Cols>::operator*(const FloatingPointType& constant) const noexcept -> Matrix<FloatingPointType, Rows, Cols> {
    Matrix<std::common_type_t<FloatingPointType, decltype(constant)>, Rows, Cols> result;
    for(std::size_t row = 0; row < Rows; row++) {
        for(std::size_t col = 0; col < Cols; col++) {
            result[row][col] = mData[row][col] * constant;
        }
    }
    return result;
}

template<FloatingPoint FloatingPointType, std::size_t Rows, std::size_t Cols>
constexpr auto Matrix<FloatingPointType, Rows, Cols>::operator/(const FloatingPointType& constant) const noexcept -> Matrix<FloatingPointType, Rows, Cols> {
    Matrix<std::common_type_t<FloatingPointType, decltype(constant)>, Rows, Cols> result;
    for(std::size_t row = 0; row < Rows; row++) {
        for(std::size_t col = 0; col < Cols; col++) {
            result[row][col] = mData[row][col] / constant;
        }
    }
    return result;
}

template<FloatingPoint FloatingPointType, std::size_t Rows, std::size_t Cols>
constexpr auto Matrix<FloatingPointType, Rows, Cols>::operator^(const auto& rhs) const -> Matrix<std::common_type_t<FloatingPointType, typename std::remove_cvref_t<decltype(rhs)>::type>, Rows, Cols> {
    using remove_cvref_rhs = typename std::remove_cvref_t<decltype(rhs)>::type;
    if (Rows != rhs.getRows()) throw cte::DimensionMismatchException("Invalid Hadamard multiplication, matrices have different FloatingPoint of rows", Rows, rhs.getRows());
    if (Cols != rhs.getCols()) throw cte::DimensionMismatchException("Invalid Hadamard multiplication, matrices have different FloatingPoint of columns", Cols, rhs.getCols());
    Matrix<std::common_type_t<FloatingPointType, remove_cvref_rhs>, Rows, Cols> result;
    for(std::size_t row = 0; row < Rows; row++) {
        for(std::size_t col = 0; col < Cols; col++) {
            result[row][col] = mData[row][col] * rhs[row][col];
        }
    }
    return result;
}

template<FloatingPoint FloatingPointType, std::size_t Rows, std::size_t Cols>
constexpr void Matrix<FloatingPointType, Rows, Cols>::apply(auto&& callable) {
    for (std::size_t row = 0; row < Rows; row++) {
        for (std::size_t col = 0; col < Cols; col++) {
            mData[row][col] = std::forward<std::remove_cvref_t<decltype(callable)>>(callable)(mData[row][col]);
        }
    }
}

template<FloatingPoint FloatingPointType, std::size_t Rows, std::size_t Cols>
constexpr void Matrix<FloatingPointType, Rows, Cols>::apply(const std::size_t row, const std::size_t col, auto&& callable) {
    if (row >= Rows) throw cte::OutOfRangeException(row);   
    if (col >= Cols) throw cte::OutOfRangeException(col);
    mData[row][col] = std::forward<std::remove_cvref_t<decltype(callable)>>(callable)(mData[row][col]);
}

template<FloatingPoint FloatingPointType, std::size_t Rows, std::size_t Cols>
constexpr FloatingPointType Matrix<FloatingPointType, Rows, Cols>::norm() const {
    if (Rows != 1UL && Cols != 1UL) throw cte::DimensionMismatchException("Invalid norm, row or column matrix needed", Rows, Cols);
    FloatingPointType result = FloatingPointType{0};
    for (std::size_t row = 0; row < Rows; row++) {
        for (std::size_t col = 0; col < Cols; col++) {
            result += cte::math::pow(mData[row][col], 2);
        }
    }
    return cte::math::sqrt(result);
}

template<FloatingPoint FloatingPointType, std::size_t Rows, std::size_t Cols>
constexpr auto Matrix<FloatingPointType, Rows, Cols>::gaussDecomposition() const noexcept -> Matrix<FloatingPointType, Rows, Cols> {
    Matrix<FloatingPointType, Rows, Cols> matrix(*this);
    FloatingPointType pivot;
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

template<FloatingPoint FloatingPointType, std::size_t Size>
constexpr auto gaussDecomposition(const std::array<std::array<std::complex<FloatingPointType>, Size>, Size>& coeffs) {
    std::array<std::array<std::complex<FloatingPointType>, Size>, Size> copy_coeffs(coeffs);
    std::complex<FloatingPointType> pivot;
    for (std::size_t current = 0; current < Size; current++) {
        for (std::size_t row = current + 1; row < Size; row++) {
            pivot = copy_coeffs[row][current] / copy_coeffs[current][current];
            for (std::size_t col = 0; col < Size; col++) {
                copy_coeffs[row][col] -= pivot * copy_coeffs[current][col];
            }
        }
    }
    return copy_coeffs;
}

template<FloatingPoint FloatingPointType, std::size_t Rows, std::size_t Cols>
constexpr FloatingPointType Matrix<FloatingPointType, Rows, Cols>::determinant() const {
    if (Rows != Cols) throw cte::DimensionMismatchException("Invalid determinant, matrix has different dimensions", Rows, Cols);
    FloatingPointType result = FloatingPointType{1};
    Matrix<FloatingPointType, Rows, Cols> matrix = gaussDecomposition();
    for (std::size_t current = 0; current < Cols; current++) {
        result *= matrix[current][current];
    }
    return result;
}

template<FloatingPoint FloatingPointType, std::size_t Size>
constexpr std::complex<FloatingPointType> determinant(const std::array<std::array<std::complex<FloatingPointType>, Size>, Size>& coeffs) {
    std::complex<FloatingPointType> pivot;
    std::array<std::array<std::complex<FloatingPointType>, Size>, Size> matrix = gaussDecomposition(coeffs);
    std::complex<FloatingPointType> result = std::complex<FloatingPointType>{1};
    for (std::size_t current = 0; current < Size; current++) {
        result *= matrix[current][current];
    }
    return result;
}

template<FloatingPoint FloatingPointType, std::size_t Rows, std::size_t Cols>
constexpr Matrix<FloatingPointType, Rows, Cols> Matrix<FloatingPointType, Rows, Cols>::inverse() const {
    if (Rows != Cols) throw cte::DimensionMismatchException("Invalid inverse, matrix has different dimensions", Rows, Cols);
    if (this->determinant() == FloatingPointType{0}) throw cte::DimensionMismatchException("Invalid inverse, matrix is singular", Rows, Cols);
    FloatingPointType pivot;
    Matrix<FloatingPointType, Rows, Cols> result;
    
    auto&& zeros = [](const FloatingPointType){ return FloatingPointType{0}; };
    auto&& ones = [](const FloatingPointType){ return FloatingPointType{1}; };

    Matrix identity;
    identity.apply(zeros);
    for (std::size_t index = 0; index < Rows; index++) {
        identity.apply(index, index, ones);
    }
    auto gaussian = this->concatenateCols(identity);
    Matrix<FloatingPointType, gaussian.getRows(), gaussian.getCols()> matrix(gaussian);
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

template<FloatingPoint FloatingPointType, std::size_t Rows, std::size_t Cols>
constexpr Matrix<FloatingPointType, Cols, Rows> Matrix<FloatingPointType, Rows, Cols>::transpose() const noexcept {
    Matrix<FloatingPointType, Cols, Rows> result;
    for (std::size_t row = 0; row < Rows; row++) {
        for (std::size_t col = 0; col < Cols; col++) {
            result[col][row] = mData[row][col];
        }
    }
    return result;
}

template<FloatingPoint FloatingPointType, std::size_t Rows, std::size_t Cols>
constexpr auto Matrix<FloatingPointType, Rows, Cols>::concatenateRows(const auto& rhs) const -> decltype(auto) {
    using remove_cvref_rhs = std::remove_cvref_t<decltype(rhs)>;
    if (Cols != remove_cvref_rhs::getCols()) throw cte::DimensionMismatchException("Invalid rows concatenation, matrices have different FloatingPoint of columns", Cols, remove_cvref_rhs::getCols());
    Matrix<FloatingPointType, Rows + remove_cvref_rhs::getCols(), Cols> result;
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

template<FloatingPoint FloatingPointType, std::size_t Rows, std::size_t Cols>
constexpr auto Matrix<FloatingPointType, Rows, Cols>::concatenateCols(const auto& rhs) const -> decltype(auto) {
    using remove_cvref_rhs = std::remove_cvref_t<decltype(rhs)>;
    if (Rows != remove_cvref_rhs::getRows()) throw cte::DimensionMismatchException("Invalid columns concatenation, matrices have different FloatingPoint of rows", Rows, remove_cvref_rhs::getRows());
    Matrix<FloatingPointType, Rows, Cols + remove_cvref_rhs::getCols()> result;
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

template<FloatingPoint FloatingPointType, std::size_t Rows, std::size_t Cols>
constexpr std::array<Matrix<FloatingPointType, 1, Cols>, Rows> Matrix<FloatingPointType, Rows, Cols>::splitRows() const noexcept {
    std::array<Matrix<FloatingPointType, 1, Cols>, Rows> result;
    for (std::size_t row = 0; row < Rows; row++) {
        for (std::size_t col = 0; col < Cols; col++) {
            result[row][0][col] = mData[row][col]; 
        }
    }
    return result;
}

template<FloatingPoint FloatingPointType, std::size_t Rows, std::size_t Cols>
constexpr std::array<Matrix<FloatingPointType, Rows, 1>, Cols> Matrix<FloatingPointType, Rows, Cols>::splitCols() const noexcept {
    std::array<Matrix<FloatingPointType, Rows, 1>, Cols> result;
    for (std::size_t row = 0; row < Rows; row++) {
        for (std::size_t col = 0; col < Cols; col++) {
            result[col][row][0] = mData[row][col]; 
        }
    }
    return result;
}

template<FloatingPoint FloatingPointType, std::size_t Rows, std::size_t Cols>
constexpr auto Matrix<FloatingPointType, Rows, Cols>::QRDecomposition() const -> QRType<FloatingPointType, Rows, Cols, FloatingPointType, Cols, Cols> {
    if (Rows != Cols) throw cte::DimensionMismatchException("Invalid QR decomposition, matrix has different dimensions", Rows, Cols);
    std::array<Matrix<FloatingPointType, Rows, 1>, Cols> columns = splitCols();
    std::array<Matrix<FloatingPointType, Rows, 1>, Cols> new_columns;
    Matrix<FloatingPointType, Rows, Cols> Q;
    FloatingPointType coeff;
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
    Matrix<FloatingPointType, Cols, Cols> R = Q.transpose() * (*this);
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

template<FloatingPoint FloatingPointType, std::size_t Size>
constexpr std::array<std::complex<FloatingPointType>, Size> solveDeterminateSystem(const std::array<std::array<std::complex<FloatingPointType>, Size>, Size>& coeffs, const std::array<std::complex<FloatingPointType>, Size>& freeTerms) {
    std::complex<FloatingPointType> determ = determinant(coeffs);
    if (determ == std::complex<FloatingPointType>{0}) throw cte::DimensionMismatchException("System is not determinate because determinant is equal to 0", Size, Size);
    std::array<std::complex<FloatingPointType>, Size> solutions;
    std::array<std::array<std::complex<FloatingPointType>, Size>, Size> new_coeffs;
    for (std::size_t index = 0; index < Size; index++) {
        for (std::size_t row = 0; row < Size; row++) {
            for (std::size_t col = 0; col < Size; col++) {
                if (index == col) {
                    new_coeffs[row][col] = freeTerms[row];
                } else {
                    new_coeffs[row][col] = coeffs[row][col];
                }
            }
        }
        solutions[index] = determinant(new_coeffs) / determ;
    }
    return solutions;
}

template<FloatingPoint FloatingPointType, std::size_t Size>
constexpr std::array<FloatingPointType, Size> solveDeterminateSystem(const Matrix<FloatingPointType, Size, Size>& coeffs, const std::array<FloatingPointType, Size> freeTerms) {
    FloatingPointType determ = coeffs.determinant();
    if (determ == FloatingPointType{0}) throw cte::DimensionMismatchException("System is not determinate because determinant is equal to 0", Size, Size);
    std::array<FloatingPointType, Size> solutions;
    Matrix<FloatingPointType, Size, Size> new_matrix;
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
        solutions[index] = new_matrix.determinant() / static_cast<FloatingPointType>(determ);
    }
    return solutions;
}

template<std::size_t Dim>
constexpr auto reductionBareissAlgorithmMapping() {
    constexpr std::size_t size = (Dim - 1) * Dim * (2 * Dim - 1) / 6;
    std::size_t v_index = 0;
    std::array<std::tuple<std::size_t, std::size_t, std::size_t, std::size_t>, size> values;
    std::size_t v1, v2, v3, v4, t1, t2 = Dim + 1;
    for (std::size_t dim = Dim - 1; dim >= 1; dim--) {
        t1 = Dim - 1 - dim;
        for (std::size_t index = 0; index < dim * dim; index++) {
            v1 = t1 * t2;
            v2 = index % dim + 1 + v1;
            v3 = (index / dim + 1) * Dim + v1;
            v4 = v2 + v3 - v1;
            values[v_index++] = std::tuple{ v1, v2, v3, v4 };
        }
    }
    return values;
}

template<FloatingPoint FloatingPointType, std::size_t Dim>
constexpr auto matrixPolynomialComposition(const cte::mat::Matrix<FloatingPointType, Dim, Dim>& mat) {
    std::array<cte::poly::Polynomial<FloatingPointType, 2 * (Dim - 1)>, Dim * Dim> polynomials{};
    for (std::size_t row = 0; row < Dim; row++) {
        for (std::size_t col = 0; col < Dim; col++) {
            polynomials[row * Dim + col] = cte::poly::Polynomial<FloatingPointType, 2 * (Dim - 1)>{};
            polynomials[row * Dim + col][2 * (Dim - 1)] = mat[row][col];
            polynomials[row * Dim + col][2 * (Dim - 1) - 1] = (row * Dim + col) % (Dim + 1) == 0 ? FloatingPointType{-1} : FloatingPointType{0};
        }
    }
    return polynomials;
}

template<auto matrix>
constexpr auto reductionBareissAlgorithm() {
    using FloatingPointType = typename decltype(matrix)::type;
    constexpr std::size_t Dim = matrix.getRows();
    constexpr auto mapping = reductionBareissAlgorithmMapping<Dim>();
    auto polynomials = matrixPolynomialComposition(matrix);
    auto pivot = cte::poly::Polynomial<FloatingPointType, 2 * (Dim - 1)>{};
    pivot[2 * (Dim - 1)] = FloatingPointType{1};
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
    using FloatingPointType = typename decltype(matrix)::type;
    constexpr std::size_t Size = matrix.getRows();
    constexpr auto eigenvalues = calculateEigenvalues<matrix>();
    std::array<std::array<std::complex<FloatingPointType>, Size>, Size> new_matrix;
    for (std::size_t row = 0; row < Size; row++) {
        for (std::size_t col = 0; col < Size; col++) {
            new_matrix[row][col] = std::complex<FloatingPointType>{matrix[row][col], 0};
        }
    }
    std::array<std::array<std::complex<FloatingPointType>, Size - 1>, Size - 1> egv;
    std::array<std::complex<FloatingPointType>, Size - 1> freeTerms;
    std::array<std::array<std::complex<FloatingPointType>, Size>, Size> eigenvectors;
    for (std::size_t index = 0; index < eigenvalues.size(); index++) {
        for (std::size_t pr = 0; pr < matrix.getRows(); pr++) {
            new_matrix[pr][pr] = matrix[pr][pr] - eigenvalues[index];
        }
        auto gaussDecomp = gaussDecomposition(new_matrix);
        for (std::size_t row = 0; row < matrix.getRows() - 1; row++) {
            for (std::size_t col = 0; col < matrix.getCols() - 1; col++) {
                egv[row][col] = gaussDecomp[row][col];
            }
        }
        for (std::size_t i = 0; i < matrix.getRows() - 1; i++) {
            freeTerms[i] = -gaussDecomp[i][matrix.getRows() - 1];
        }
        auto solutions = solveDeterminateSystem(egv, freeTerms);
        for (std::size_t row = 0; row < freeTerms.size(); row++) {
            eigenvectors[index][row] = solutions[row];
        }
        eigenvectors[index][matrix.getCols() - 1] = FloatingPointType{1};
    }
    return eigenvectors;
}

}

