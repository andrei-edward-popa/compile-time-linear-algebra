#include <gtest/gtest.h>

#include "Matrix.hpp"

TEST(Dimensions, TwoByThreeMatrix) {
	constexpr cte::mat::Matrix matrix{ { 1.0l, -2.0l, 7.0l }, { 3.0l, 4.0l, 3.0l } };
	constexpr std::size_t rows = matrix.getRows();
	constexpr std::size_t cols = matrix.getCols();
	constexpr std::size_t expected_rows = 2;
	constexpr std::size_t expected_cols = 3;
	EXPECT_EQ(rows, expected_rows);
	EXPECT_EQ(cols, expected_cols);
}

TEST(Dimensions, ColumnVector) {
	constexpr cte::mat::Matrix matrix{ { 1.0l }, { 3.0l }, { -11.1l } };
	constexpr std::size_t rows = matrix.getRows();
	constexpr std::size_t cols = matrix.getCols();
	constexpr std::size_t expected_rows = 3;
	constexpr std::size_t expected_cols = 1;
	EXPECT_EQ(rows, expected_rows);
	EXPECT_EQ(cols, expected_cols);
}

TEST(Dimensions, RowVector) {
	constexpr cte::mat::Matrix matrix{ { 1.0l, -2.0l, 7.0l } };
	constexpr std::size_t rows = matrix.getRows();
	constexpr std::size_t cols = matrix.getCols();
	constexpr std::size_t expected_rows = 1;
	constexpr std::size_t expected_cols = 3;
	EXPECT_EQ(rows, expected_rows);
	EXPECT_EQ(cols, expected_cols);
}

TEST(Dimensions, OneElement) {
	constexpr cte::mat::Matrix matrix{ { 1.0l } };
	constexpr std::size_t rows = matrix.getRows();
	constexpr std::size_t cols = matrix.getCols();
	constexpr std::size_t expected_rows = 1;
	constexpr std::size_t expected_cols = 1;
	EXPECT_EQ(rows, expected_rows);
	EXPECT_EQ(cols, expected_cols);
}

TEST(Equality, DifferentNumberOfRows) {
	constexpr cte::mat::Matrix matrix1{ { 1.0l, -2.0l }, { 3.0l, 4.0l } };
	constexpr cte::mat::Matrix matrix2{ { 3.0l, 2.1l }, { 7.2l, -4.0l }, { 1.1l, 3.3l } };
	constexpr bool eq_result = matrix1 == matrix2;
	constexpr bool not_eq_result = matrix1 != matrix2;
	EXPECT_FALSE(eq_result);
	EXPECT_TRUE(not_eq_result);
}

TEST(Equality, DifferentNumberOfColumns) {
	constexpr cte::mat::Matrix matrix1{ { 1.0l, -2.0l }, { 3.0l, 4.0l } };
	constexpr cte::mat::Matrix matrix2{ { 3.0l, 2.1l, 4.4l }, { 7.2l, -4.0l, 5.5l } };
	constexpr bool eq_result = matrix1 == matrix2;
	constexpr bool not_eq_result = matrix1 != matrix2;
	EXPECT_FALSE(eq_result);
	EXPECT_TRUE(not_eq_result);
}

TEST(Equality, NotEqualAndSameDimensions) {
	constexpr cte::mat::Matrix matrix1{ { 1.0l, -2.0l }, { 3.0l, 4.0l } };
	constexpr cte::mat::Matrix matrix2{ { 1.0l, -2.0l }, { 4.0l, 4.0l } };
	constexpr bool eq_result = matrix1 == matrix2;
	constexpr bool not_eq_result = matrix1 != matrix2;
	EXPECT_FALSE(eq_result);
	EXPECT_TRUE(not_eq_result);
}

TEST(Equality, EqualAndSameDimensions) {
	constexpr cte::mat::Matrix matrix1{ { 1.0l, -2.0l }, { 3.0l, 4.0l } };
	constexpr cte::mat::Matrix matrix2{ { 1.0l, -2.0l }, { 3.0l, 4.0l } };
	constexpr bool eq_result = matrix1 == matrix2;
	constexpr bool not_eq_result = matrix1 != matrix2;
	EXPECT_TRUE(eq_result);
	EXPECT_FALSE(not_eq_result);
}

TEST(Addition, TwoByTwo) {
	constexpr cte::mat::Matrix matrix1{ { 1.0l, -2.0l }, { 3.0l, 4.0l} };
	constexpr cte::mat::Matrix matrix2{ { 3.0l, 2.1l }, { 7.2l, -4.0l} };
	constexpr cte::mat::Matrix add_result = matrix1 + matrix2;
	constexpr cte::mat::Matrix expected_result{ { 4.0l, 0.1l }, { 10.2l, 0.0l } };
	EXPECT_EQ(add_result, expected_result);
}

TEST(Addition, ThreeByThree) {
	constexpr cte::mat::Matrix matrix1{ { 1.0l, -2.0l, 4.3l }, { 3.0l, 4.0l, -2.2l }, { 7.0l, -1.2l, -3.6l } };
	constexpr cte::mat::Matrix matrix2{ { 3.0l, -2.1l, -4.3l }, { 6.1l, 2.5l, 5.4l }, { -7.4l, -6.5l, -3.6l } };
	constexpr cte::mat::Matrix add_result = matrix1 + matrix2;
	constexpr cte::mat::Matrix expected_result{ { 4.0l, -4.1l, 0.0l }, { 9.1l, 6.5l, 3.2l }, { -0.4l, -7.7l, -7.2l } };
	EXPECT_EQ(add_result, expected_result);
}

TEST(Subtraction, TwoByTwo) {
	constexpr cte::mat::Matrix matrix1{ { 1.0l, -2.0l }, { 3.0l, 4.0l} };
	constexpr cte::mat::Matrix matrix2{ { 3.0l, 2.1l }, { 7.2l, -4.0l} };
	constexpr cte::mat::Matrix sub_result = matrix1 - matrix2;
	constexpr cte::mat::Matrix expected_result{ { -2.0l, -4.1l }, { -4.2l, 8.0l } };
	EXPECT_EQ(sub_result, expected_result);
}

TEST(Subtraction, ThreeByThree) {
	constexpr cte::mat::Matrix matrix1{ { 1.0l, -2.0l, 4.3l }, { 3.0l, 4.0l, -2.2l }, { 7.0l, -1.2l, -3.6l } };
	constexpr cte::mat::Matrix matrix2{ { 3.0l, -2.1l, -4.3l }, { 6.1l, 2.5l, 5.4l }, { -7.4l, -6.5l, -3.6l } };
	constexpr cte::mat::Matrix sub_result = matrix1 - matrix2;
	constexpr cte::mat::Matrix expected_result{ { -2.0l, 0.1l, 8.6l }, { -3.1l, 1.5l, -7.6l }, { 14.4l, 5.3l, 0.0l } };
	EXPECT_EQ(sub_result, expected_result);
}

TEST(UsualMultiplication, TwoByTwoAndTwoByTwo) {
	constexpr cte::mat::Matrix matrix1{ { 1.0l, -2.0l }, { 3.0l, 4.0l} };
	constexpr cte::mat::Matrix matrix2{ { 3.0l, 2.1l }, { 7.2l, -4.0l} };
	constexpr cte::mat::Matrix mult_result = matrix1 * matrix2;
	constexpr cte::mat::Matrix expected_result{ { -11.4l, 10.1l }, { 37.8l, -9.7l } };
	EXPECT_EQ(mult_result, expected_result);
}

TEST(UsualMultiplication, ThreeByTwoAndTwoByThree) {
	constexpr cte::mat::Matrix matrix1{ { 1.0l, -2.0l }, { 3.0l, 4.0l }, { 7.0l, -1.2l } };
	constexpr cte::mat::Matrix matrix2{ { 3.0l, -2.1l, -4.3l }, { 6.1l, 2.5l, 5.4l } };
	constexpr cte::mat::Matrix mult_result = matrix1 * matrix2;
	constexpr cte::mat::Matrix expected_result{ { -9.2l, -7.1l, -15.1l }, { 33.4l, 3.7l, 8.7l }, { 13.68l, -17.7l, -36.58l } };
	constexpr std::size_t rows = mult_result.getRows();
	constexpr std::size_t cols = mult_result.getCols();
	constexpr std::size_t expected_rows = 3;
	constexpr std::size_t expected_cols = 3;
	EXPECT_EQ(rows, expected_rows);
	EXPECT_EQ(cols, expected_cols);
	EXPECT_EQ(mult_result, expected_result);
}

TEST(Power, TwoByTwo) {
	constexpr cte::mat::Matrix matrix{ { 1.0l, -2.0l }, { 3.0l, 4.0l} };
	constexpr cte::mat::Matrix power_result = matrix ^ 3;
	constexpr cte::mat::Matrix expected_result{ { -35.0l, -30.0l }, { 45.0l, 10.0l } };
	EXPECT_EQ(power_result, expected_result);
}

TEST(Power, ThreeByThree) {
	constexpr cte::mat::Matrix matrix{ { 1.0l, -2.0l, -2.0l }, { 3.0l, 4.0l, -1.0l }, { 7.0l, -1.2l, 3.3l } };
	constexpr cte::mat::Matrix power_result = matrix ^ 4;
	constexpr cte::mat::Matrix expected_result{ { 125.3l, 209.496l, 239.086l }, { -414.85l, 367.348l, -176.357l }, { -736.195l, -412.8404l, 131.4561l } };
	EXPECT_EQ(power_result, expected_result);
}

TEST(InPlacePower, TwoByTwo) {
	cte::mat::Matrix matrix{ { 1.0l, -2.0l }, { 3.0l, 4.0l} };
	matrix ^= 3;
	constexpr cte::mat::Matrix expected_result{ { -35.0l, -30.0l }, { 45.0l, 10.0l } };
	EXPECT_EQ(matrix, expected_result);
}

TEST(InPlacePower, ThreeByThree) {
	cte::mat::Matrix matrix{ { 1.0l, -2.0l, -2.0l }, { 3.0l, 4.0l, -1.0l }, { 7.0l, -1.2l, 3.3l } };
	matrix ^= 4;
	constexpr cte::mat::Matrix expected_result{ { 125.3l, 209.496l, 239.086l }, { -414.85l, 367.348l, -176.357l }, { -736.195l, -412.8404l, 131.4561l } };
	EXPECT_EQ(matrix, expected_result);
}

TEST(HadamardMultiplication, TwoByTwo) {
	constexpr cte::mat::Matrix matrix1{ { 1.0l, -2.0l }, { 3.0l, 4.0l} };
	constexpr cte::mat::Matrix matrix2{ { 3.0l, 2.1l }, { 7.2l, -4.0l} };
	constexpr cte::mat::Matrix mult_result = matrix1 ^ matrix2;
	constexpr cte::mat::Matrix expected_result{ { 3.0l, -4.2l }, { 21.6l, -16.0l } };
	EXPECT_EQ(mult_result, expected_result);
}

TEST(HadamardMultiplication, ThreeByThree) {
	constexpr cte::mat::Matrix matrix1{ { 1.0l, -2.0l, -2.0l }, { 3.0l, 4.0l, -1.0l }, { 7.0l, -1.2l, 3.3l } };
	constexpr cte::mat::Matrix matrix2{ { 3.0l, -2.1l, -4.3l }, { 6.1l, 2.5l, 5.4l }, { 2.3l, -1.2l, -1.2l } };
	constexpr cte::mat::Matrix mult_result = matrix1 ^ matrix2;
	constexpr cte::mat::Matrix expected_result{ { 3.0l, 4.2l, 8.6l }, { 18.3l, 10.0l, -5.4l }, { 16.1l, 1.44l, -3.96l } };
	EXPECT_EQ(mult_result, expected_result);
}

TEST(ConstantMultiplication, TwoByTwo) {
	constexpr cte::mat::Matrix matrix{ { 1.0l, -2.0l }, { 3.0l, 4.0l} };
	constexpr cte::mat::Matrix mult_result = matrix * 3.0l;
	constexpr cte::mat::Matrix expected_result{ { 3.0l, -6.0l }, { 9.0l, 12.0l } };
	EXPECT_EQ(mult_result, expected_result);
}

TEST(ConstantMultiplication, ThreeByThree) {
	constexpr cte::mat::Matrix matrix{ { 1.0l, -2.0l, -2.0l }, { 3.0l, 4.0l, -1.0l }, { 7.0l, -1.2l, 3.3l } };
	constexpr cte::mat::Matrix mult_result = matrix * 4.0l;
	constexpr cte::mat::Matrix expected_result{ { 4.0l, -8.0l, -8.0l }, { 12.0l, 16.0l, -4.0l }, { 28.0l, -4.8l, 13.2l } };
	EXPECT_EQ(mult_result, expected_result);
}

TEST(InPlaceConstantMultiplication, TwoByTwo) {
	cte::mat::Matrix matrix{ { 1.0l, -2.0l }, { 3.0l, 4.0l} };
	matrix *= 0.5l;
	constexpr cte::mat::Matrix expected_result{ { 0.5l, -1.0l }, { 1.5l, 2.0l } };
	EXPECT_EQ(matrix, expected_result);
}

TEST(InPlaceConstantMultiplication, ThreeByThree) {
	cte::mat::Matrix matrix{ { 1.0l, -2.0l, -2.0l }, { 3.0l, 4.0l, -1.0l }, { 7.0l, -1.2l, 3.3l } };
	matrix *= 0.5l;
	constexpr cte::mat::Matrix expected_result{ { 0.5l, -1.0l, -1.0l }, { 1.5l, 2.0l, -0.5l }, { 3.5l, -0.6l, 1.65l } };
	EXPECT_EQ(matrix, expected_result);
}

TEST(ConstantDivision, TwoByTwo) {
	constexpr cte::mat::Matrix matrix{ { 1.0l, -2.0l }, { 3.0l, 4.0l} };
	constexpr cte::mat::Matrix div_result = matrix / 5.0l;
	constexpr cte::mat::Matrix expected_result{ { 0.2l, -0.4l }, { 0.6l, 0.8l } };
	EXPECT_EQ(div_result, expected_result);
}

TEST(ConstantDivision, ThreeByThree) {
	constexpr cte::mat::Matrix matrix{ { 1.0l, -2.0l, -2.0l }, { 3.0l, 4.0l, -1.0l }, { 7.0l, -1.2l, 3.3l } };
	constexpr cte::mat::Matrix div_result = matrix / 10.0l;
	constexpr cte::mat::Matrix expected_result{ { 0.1l, -0.2l, -0.2l }, { 0.3l, 0.4l, -0.1l }, { 0.7l, -0.12l, 0.33l } };
	EXPECT_EQ(div_result, expected_result);
}

TEST(InPlaceConstantDivision, TwoByTwo) {
	cte::mat::Matrix matrix{ { 1.0l, -2.0l }, { 3.0l, 4.0l} };
	matrix /= 8.0l;
	constexpr cte::mat::Matrix expected_result{ { 0.125l, -0.25l }, { 0.375l, 0.5l } };
	EXPECT_EQ(matrix, expected_result);
}

TEST(InPlaceConstantDivision, ThreeByThree) {
	cte::mat::Matrix matrix{ { 1.0l, -2.0l, -2.0l }, { 3.0l, 4.0l, -1.0l }, { 7.0l, -1.2l, 3.3l } };
	matrix /= 8.0l;
	constexpr cte::mat::Matrix expected_result{ { 0.125l, -0.25l, -0.25l }, { 0.375l, 0.5l, -0.125l }, { 0.875l, -0.15l, 0.4125l } };
	EXPECT_EQ(matrix, expected_result);
}

TEST(ApplyCallableToMatrix, TwoByTwo) {
	cte::mat::Matrix matrix{ { 1.0l, -2.0l }, { 3.0l, 4.0l} };
	matrix.apply([](const typename decltype(matrix)::type value){
		return value * 2.1l;
	});
	constexpr cte::mat::Matrix expected_result{ { 2.1l, -4.2l }, { 6.3l, 8.4l } };
	EXPECT_EQ(matrix, expected_result);
}

TEST(ApplyCallableToMatrix, TwoByFour) {
	cte::mat::Matrix matrix{ { 1.0l, -2.0l, 4.3l, -1.1l }, { 3.2l, 4.0l, 3.3l, -4.1l } };
	matrix.apply([](const typename decltype(matrix)::type value){
		return value * 10.0l;
	});
	constexpr cte::mat::Matrix expected_result{ { 10.0l, -20.0l, 43.0l, -11.0l }, { 32.0l, 40.0l, 33.0l, -41.0l } };
	EXPECT_EQ(matrix, expected_result);
}

TEST(ApplyCallableToElement, TwoByTwo) {
	cte::mat::Matrix matrix{ { 1.0l, -2.0l }, { 3.0l, 4.0l} };
	for (std::size_t i = 0; i < matrix.getRows(); i++) {
		for (std::size_t j = 0; j < matrix.getCols(); j++) {
			matrix.apply(i, j, [&](const typename decltype(matrix)::type value){
				return value * (i + j);
			});
		}
	}
	constexpr cte::mat::Matrix expected_result{ { 0.0l, -2.0l }, { 3.0l, 8.0l } };
	EXPECT_EQ(matrix, expected_result);
}

TEST(ApplyCallableToElement, ThreeByFour) {
	cte::mat::Matrix matrix{ { 1.0l, -2.0l, 4.3l, -1.1l }, { 3.2l, 4.0l, 3.3l, -4.1l }, { -1.2l, 2.0l, 3.7l, -6.1l } };
	for (std::size_t i = 0; i < matrix.getRows(); i++) {
		for (std::size_t j = 0; j < matrix.getCols(); j++) {
			matrix.apply(i, j, [&](const typename decltype(matrix)::type value){
				return value * i * j;
			});
		}
	}
	constexpr cte::mat::Matrix expected_result{ { 0.0l, 0.0l, 0.0l, 0.0l }, { 0.0l, 4.0l, 6.6l, -12.3l }, { 0.0l, 4.0l, 14.8l, -36.6l } };
	EXPECT_EQ(matrix, expected_result);
}

TEST(Normalization, TwoElementsColumnVector) {
	constexpr cte::mat::Matrix matrix{ { 3.0l}, { 4.0l } };
	constexpr typename decltype(matrix)::type norm_result = matrix.norm();
	constexpr typename decltype(matrix)::type expected_result = 5.0l;
	EXPECT_EQ(norm_result, expected_result);
}

TEST(Normalization, ThreeElementsColumnVector) {
	constexpr cte::mat::Matrix matrix{ { 2.0l}, { 3.0l }, { 6.0l } };
	constexpr typename decltype(matrix)::type norm_result = matrix.norm();
	constexpr typename decltype(matrix)::type expected_result = 7.0l;
	EXPECT_EQ(norm_result, expected_result);
}

TEST(Normalization, TwoElementsRowVector) {
	constexpr cte::mat::Matrix matrix{ { 6.0l, 8.0l } };
	constexpr typename decltype(matrix)::type norm_result = matrix.norm();
	constexpr typename decltype(matrix)::type expected_result = 10.0l;
	EXPECT_EQ(norm_result, expected_result);
}

TEST(Normalization, ThreeElementsRowVector) {
	constexpr cte::mat::Matrix matrix{ { 4.0l, 5.0l, 10.0l } };
	constexpr typename decltype(matrix)::type norm_result = matrix.norm();
	constexpr typename decltype(matrix)::type expected_result = 11.874342087037917234l;
	EXPECT_EQ(norm_result, expected_result);
}

TEST(Normalization, OneElement) {
	constexpr cte::mat::Matrix matrix{ { 16.0l } };
	constexpr typename decltype(matrix)::type norm_result = matrix.norm();
	constexpr typename decltype(matrix)::type expected_result = 16.0l;
	EXPECT_EQ(norm_result, expected_result);
}

TEST(GaussDecomposition, ThreeByTwo) {
	constexpr cte::mat::Matrix matrix{ { 1.0l, 2.0l }, { 4.0l, 6.0l }, { 7.0l, 8.0l } };
	constexpr cte::mat::Matrix gauss_result = matrix.gaussDecomposition();
	constexpr cte::mat::Matrix expected_result{ { 1.0l, 2.0l }, { 0.0l, -2.0l }, { 0.0l, 0.0l } };
	EXPECT_EQ(gauss_result, expected_result);
}

TEST(GaussDecomposition, ThreeByThree) {
	constexpr cte::mat::Matrix matrix{ { 1.0l, 2.0l, 3.0l }, { 4.0l, 5.0l, 6.0l }, { 7.0l, 8.0l, 8.0l } };
	constexpr cte::mat::Matrix gauss_result = matrix.gaussDecomposition();
	constexpr cte::mat::Matrix expected_result{ { 1.0l, 2.0l, 3.0l }, { 0.0l, -3.0l, -6.0l }, { 0.0l, 0.0l, -1.0l } };
	EXPECT_EQ(gauss_result, expected_result);
}

TEST(Determinant, TwoByTwo) {
	constexpr cte::mat::Matrix matrix{ { 1.0l, 3.0l }, { 4.0l, 6.0l } };
	constexpr typename decltype(matrix)::type determinant_result = matrix.determinant();
	constexpr typename decltype(matrix)::type expected_result = -6.0l;
	EXPECT_EQ(determinant_result, expected_result);
}

TEST(Determinant, ThreeByThree) {
	constexpr cte::mat::Matrix matrix{ { 1.0l, 2.0l, 3.0l }, { 4.0l, 5.0l, 6.0l }, { 7.0l, 8.0l, 8.0l } };
	constexpr typename decltype(matrix)::type determinant_result = matrix.determinant();
	constexpr typename decltype(matrix)::type expected_result = 3.0l;
	EXPECT_EQ(determinant_result, expected_result);
}

TEST(Inverse, ThreeByThree) {
	constexpr cte::mat::Matrix matrix{ { 1.0l, 2.0l, 3.0l }, { 4.0l, 5.0l, 6.0l }, { 7.0l, 8.0l, 8.0l } };
	constexpr cte::mat::Matrix inverse_result = matrix.inverse();
	constexpr cte::mat::Matrix expected_result{ { -2.666666666666666667l, 2.666666666666666667l, -1.0l }, { 3.333333333333333333l, -4.333333333333333333l, 2.0l }, { -1.0l, 2.0l, -1.0l } };
	EXPECT_EQ(inverse_result, expected_result);
}

TEST(Inverse, FourByFour) {
	constexpr cte::mat::Matrix matrix{ { 2.0l, 3.0l, 4.0l, 5.0l }, { 1.0l, 0.0l, 5.0l, 6.0l }, { 3.0l, -2.0l, 9.0l, 9.0l }, { 1.0l, 2.0l, 3.0l, 4.0l} };
	constexpr cte::mat::Matrix inverse_result = matrix.inverse();
	constexpr cte::mat::Matrix expected_result{ { 7.0l, 5.5l, -2.0l, -12.5l }, { -3.0l, -3.0l, 1.0l, 6.0l }, { -11.0l, -10.5l, 4.0l, 20.5l }, { 8.0l, 8.0l, -3.0l, -15.0l } };
	EXPECT_EQ(inverse_result, expected_result);
}

TEST(Transpose, TwoByFour) {
	constexpr cte::mat::Matrix matrix{ { 1.0l, 2.0l, 4.0l, 1.0l }, { 4.0l, 5.0l, -1.0l, -5.0l } };
	constexpr cte::mat::Matrix transpose_result = matrix.transpose();
	constexpr std::size_t rows = transpose_result.getRows();
	constexpr std::size_t cols = transpose_result.getCols();
	constexpr std::size_t expected_rows = matrix.getCols();
	constexpr std::size_t expected_cols = matrix.getRows();
	constexpr cte::mat::Matrix expected_result{ { 1.0l, 4.0l }, { 2.0l, 5.0l }, { 4.0l, -1.0l }, { 1.0l, -5.0l } };
	EXPECT_EQ(rows, expected_rows);
	EXPECT_EQ(cols, expected_cols);
	EXPECT_EQ(transpose_result, expected_result);
}

TEST(Transpose, ThreeByTwo) {
	constexpr cte::mat::Matrix matrix{ { 1.0l, 2.0l }, { 4.0l, 5.0l }, { 7.0l, 8.0l } };
	constexpr cte::mat::Matrix transpose_result = matrix.transpose();
	constexpr std::size_t rows = transpose_result.getRows();
	constexpr std::size_t cols = transpose_result.getCols();
	constexpr std::size_t expected_rows = matrix.getCols();
	constexpr std::size_t expected_cols = matrix.getRows();
	constexpr cte::mat::Matrix expected_result{ { 1.0l, 4.0l, 7.0l }, { 2.0l, 5.0l, 8.0l } };
	EXPECT_EQ(rows, expected_rows);
	EXPECT_EQ(cols, expected_cols);
	EXPECT_EQ(transpose_result, expected_result);
}

TEST(ConcatenateRows, ThreeByTwoAndTwoByTwo) {
	constexpr cte::mat::Matrix matrix1{ { 1.0l, 2.0l }, { 4.0l, 5.0l }, { 7.0l, 8.0l } };
	constexpr cte::mat::Matrix matrix2{ { -1.0l, -2.0l }, { -4.0l, -5.0l } };
	constexpr cte::mat::Matrix concat_result = matrix1.concatenateRows(matrix2);
	constexpr cte::mat::Matrix expected_result{ { 1.0l, 2.0l }, { 4.0l, 5.0l }, { 7.0l, 8.0l }, { -1.0l, -2.0l }, { -4.0l, -5.0l } };
	constexpr std::size_t rows = concat_result.getRows();
	constexpr std::size_t cols = concat_result.getCols();
	constexpr std::size_t expected_rows = matrix1.getRows() + matrix2.getRows();
	constexpr std::size_t expected_cols = matrix1.getCols();
	EXPECT_EQ(rows, expected_rows);
	EXPECT_EQ(cols, expected_cols);
	EXPECT_EQ(concat_result, expected_result);
}

TEST(ConcatenateRows, ThreeByFourAndTwoByFour) {
	constexpr cte::mat::Matrix matrix1{ { 1.0l, 2.0l, 1.0l, 5.0l }, { 4.0l, 5.0l, -1.0l, -2.2l }, { 7.0l, 8.0l, 11.0l, -44.2l } };
	constexpr cte::mat::Matrix matrix2{ { -1.0l, -2.0l, 4.0l, 5.0l }, { -4.0l, -5.0l, 8.0l, 11.0l } };
	constexpr cte::mat::Matrix concat_result = matrix1.concatenateRows(matrix2);
	constexpr cte::mat::Matrix expected_result{ { 1.0l, 2.0l, 1.0l, 5.0l }, { 4.0l, 5.0l, -1.0l, -2.2l }, { 7.0l, 8.0l, 11.0l, -44.2l }, { -1.0l, -2.0l, 4.0l, 5.0l }, { -4.0l, -5.0l, 8.0l, 11.0l } };
	constexpr std::size_t rows = concat_result.getRows();
	constexpr std::size_t cols = concat_result.getCols();
	constexpr std::size_t expected_rows = matrix1.getRows() + matrix2.getRows();
	constexpr std::size_t expected_cols = matrix1.getCols();
	EXPECT_EQ(rows, expected_rows);
	EXPECT_EQ(cols, expected_cols);
	EXPECT_EQ(concat_result, expected_result);
}


TEST(ConcatenateColumns, ThreeByTwoAndThreeByThree) {
	constexpr cte::mat::Matrix matrix1{ { 1.0l, 2.0l }, { 4.0l, 5.0l }, { 7.0l, 8.0l } };
	constexpr cte::mat::Matrix matrix2{ { 1.0l, -2.0l, 3.0l }, { -4.0l, 5.0l, 5.0l }, { 7.0l, 8.0l, -7.0l } };
	constexpr cte::mat::Matrix concat_result = matrix1.concatenateCols(matrix2);
	constexpr cte::mat::Matrix expected_result{ { 1.0l, 2.0l, 1.0l, -2.0l, 3.0l }, { 4.0l, 5.0l, -4.0l, 5.0l, 5.0l }, { 7.0l, 8.0l, 7.0l, 8.0l, -7.0l } };
	constexpr std::size_t rows = concat_result.getRows();
	constexpr std::size_t cols = concat_result.getCols();
	constexpr std::size_t expected_rows = matrix1.getRows();
	constexpr std::size_t expected_cols = matrix1.getCols() + matrix2.getCols();
	EXPECT_EQ(rows, expected_rows);
	EXPECT_EQ(cols, expected_cols);
	EXPECT_EQ(concat_result, expected_result);
}

TEST(ConcatenateColumns, TwoByTwoAndTwoByFour) {
	constexpr cte::mat::Matrix matrix1{ { 1.0l, 2.0l }, { 4.0l, 5.0l } };
	constexpr cte::mat::Matrix matrix2{ { 1.0l, -2.0l, -4.0l, 5.0l }, { 7.0l, 8.0l, 1.0l, 2.0l } };
	constexpr cte::mat::Matrix concat_result = matrix1.concatenateCols(matrix2);
	constexpr cte::mat::Matrix expected_result{ { 1.0l, 2.0l, 1.0l, -2.0l, -4.0l, 5.0l }, { 4.0l, 5.0l, 7.0l, 8.0l, 1.0l, 2.0l } };
	constexpr std::size_t rows = concat_result.getRows();
	constexpr std::size_t cols = concat_result.getCols();
	constexpr std::size_t expected_rows = matrix1.getRows();
	constexpr std::size_t expected_cols = matrix1.getCols() + matrix2.getCols();
	EXPECT_EQ(rows, expected_rows);
	EXPECT_EQ(cols, expected_cols);
	EXPECT_EQ(concat_result, expected_result);
}

TEST(SplitRows, ThreeByTwo) {
	constexpr cte::mat::Matrix matrix{ { 1.0l, 2.0l }, { 4.0l, 5.0l }, { 7.0l, 8.0l } };
	constexpr auto split_result = matrix.splitRows();
	constexpr cte::mat::Matrix row0{ { 1.0l, 2.0l } };
	constexpr cte::mat::Matrix row1{ { 4.0l, 5.0l } };
	constexpr cte::mat::Matrix row2{ { 7.0l, 8.0l } };
	EXPECT_EQ(split_result[0], row0);
	EXPECT_EQ(split_result[1], row1);
	EXPECT_EQ(split_result[2], row2);
}

TEST(SplitRows, FourByThree) {
	constexpr cte::mat::Matrix matrix{ { 1.0l, 2.0l, -1.0l }, { 4.0l, 5.0l, -3.3l }, { 7.0l, 8.0l, 2.2l }, { 7.0l, 8.0l, 5.3l } };
	constexpr auto split_result = matrix.splitRows();
	constexpr cte::mat::Matrix row0{ { 1.0l, 2.0l, -1.0l } };
	constexpr cte::mat::Matrix row1{ { 4.0l, 5.0l, -3.3l } };
	constexpr cte::mat::Matrix row2{ { 7.0l, 8.0l, 2.2l } };
	constexpr cte::mat::Matrix row3{ { 7.0l, 8.0l, 5.3l } };
	EXPECT_EQ(split_result[0], row0);
	EXPECT_EQ(split_result[1], row1);
	EXPECT_EQ(split_result[2], row2);
	EXPECT_EQ(split_result[3], row3);
}

TEST(SplitColumns, ThreeByTwo) {
	constexpr cte::mat::Matrix matrix{ { 1.0l, 2.0l }, { 4.0l, 5.0l }, { 7.0l, 8.0l } };
	constexpr auto split_result = matrix.splitCols();
	constexpr cte::mat::Matrix col0{ { 1.0l }, { 4.0l }, { 7.0l } };
	constexpr cte::mat::Matrix col1{ { 2.0l }, { 5.0l }, { 8.0l } };
	EXPECT_EQ(split_result[0], col0);
	EXPECT_EQ(split_result[1], col1);
}

TEST(SplitColumns, TwoByFour) {
	constexpr cte::mat::Matrix matrix{ { 1.0l, 2.0l, 3.0l, 4.0l }, { 4.0l, 5.0l, 6.0l, 7.0l } };
	constexpr auto split_result = matrix.splitCols();
	constexpr cte::mat::Matrix col0{ { 1.0l }, { 4.0l } };
	constexpr cte::mat::Matrix col1{ { 2.0l }, { 5.0l } };
	constexpr cte::mat::Matrix col2{ { 3.0l }, { 6.0l } };
	constexpr cte::mat::Matrix col3{ { 4.0l }, { 7.0l } };
	EXPECT_EQ(split_result[0], col0);
	EXPECT_EQ(split_result[1], col1);
	EXPECT_EQ(split_result[2], col2);
	EXPECT_EQ(split_result[3], col3);
}

TEST(QRDecomposition, ThreeByThree) {
	constexpr cte::mat::Matrix matrix{ { 12.0l, -51.0l, 4.0l }, { 6.0l, 167.0l, -68.0l }, { -4.0l, 24.0l, -41.0l } };
	constexpr auto qr_result = matrix.QRDecomposition();
	constexpr auto q_result = qr_result.getOrthogonal();
	constexpr auto r_result = qr_result.getUpperTriangular();
	constexpr cte::mat::Matrix expected_q_result{ { 0.857142857142857143l, -0.394285714285714286l, -0.331428571428571429l }, { 0.428571428571428571l, 0.902857142857142857l, 0.034285714285714286l }, { -0.285714285714285714l, 0.171428571428571429l, -0.942857142857142857l } };
	constexpr cte::mat::Matrix expected_r_result{ { 14.0l, 21.0l, -14.0l }, { 0.0l, 175.0l, -70.0l }, { 0.0l, 0.0l, 35.0l } };
	EXPECT_EQ(q_result, expected_q_result);
	EXPECT_EQ(r_result, expected_r_result);
}

TEST(QRDecomposition, FourByFour) {
	constexpr cte::mat::Matrix matrix{ { 2.0l, 2.0l, 3.0l, 4.0l }, { 2.0l, 3.0l, 4.0l, 5.0l }, { 10.0l, -10.0l, 11.0l, 12.0l }, { -1.0l, -2.0l, -3.0l, -4.0l } };
	constexpr auto qr_result = matrix.QRDecomposition();
	constexpr auto q_result = qr_result.getOrthogonal();
	constexpr auto r_result = qr_result.getUpperTriangular();
	constexpr cte::mat::Matrix expected_q_result{ { 0.19156525704423l, 0.533221260271359l, -0.606841785511149l, 0.557423436218682l }, { 0.19156525704423l, 0.680736786590086l, 0.00922076787101272l, -0.706976065448082l }, { 0.957826285221152l, -0.284204225017729l, 0.0401487601050427l, -0.0135956935663092l }, { -0.0957826285221152l, -0.414126156454406l, -0.793754434229843l, -0.435062194121895l } };
	constexpr cte::mat::Matrix expected_r_result{ { 10.4403065089105l, -8.42887130994614l, 12.1643938223086l, 13.6011332501404l }, { 0.0l, 6.77894744339908l, 2.43874292134262l, 3.78262289964074l }, { 0.0l, 0.0l, 1.0392573787956l, 1.27553955549035l }, { 0.0l, 0.0l, 0.0l, 0.271913871326186l } };
	EXPECT_EQ(q_result, expected_q_result);
	EXPECT_EQ(r_result, expected_r_result);
}


