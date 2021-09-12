#include <gtest/gtest.h>

#include "Matrix.hpp"

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

TEST(Subtraction, TwoByTwo) {
	constexpr cte::mat::Matrix matrix1{ { 1.0l, -2.0l }, { 3.0l, 4.0l} };
	constexpr cte::mat::Matrix matrix2{ { 3.0l, 2.1l }, { 7.2l, -4.0l} };
	constexpr cte::mat::Matrix sub_result = matrix1 - matrix2;
	constexpr cte::mat::Matrix expected_result{ { -2.0l, -4.1l }, { -4.2l, 8.0l } };
	EXPECT_EQ(sub_result, expected_result);
}

TEST(UsualMultiplication, TwoByTwo) {
	constexpr cte::mat::Matrix matrix1{ { 1.0l, -2.0l }, { 3.0l, 4.0l} };
	constexpr cte::mat::Matrix matrix2{ { 3.0l, 2.1l }, { 7.2l, -4.0l} };
	constexpr cte::mat::Matrix mult_result = matrix1 * matrix2;
	constexpr cte::mat::Matrix expected_result{ { -11.4l, 10.1l }, { 37.8l, -9.7l } };
	EXPECT_EQ(mult_result, expected_result);
}

TEST(ConstantMultiplication, TwoByTwo) {
	constexpr cte::mat::Matrix matrix{ { 1.0l, -2.0l }, { 3.0l, 4.0l} };
	constexpr cte::mat::Matrix mult_result = matrix * 3.0l;
	constexpr cte::mat::Matrix expected_result{ { 3.0l, -6.0l }, { 9.0l, 12.0l } };
	EXPECT_EQ(mult_result, expected_result);
}

TEST(HadamardMultiplication, TwoByTwo) {
	constexpr cte::mat::Matrix matrix1{ { 1.0l, -2.0l }, { 3.0l, 4.0l} };
	constexpr cte::mat::Matrix matrix2{ { 3.0l, 2.1l }, { 7.2l, -4.0l} };
	constexpr cte::mat::Matrix mult_result = matrix1 ^ matrix2;
	constexpr cte::mat::Matrix expected_result{ { 3.0l, -4.2l }, { 21.6l, -16.0l } };
	EXPECT_EQ(mult_result, expected_result);
}

TEST(ConstantDivision, TwoByTwo) {
	constexpr cte::mat::Matrix matrix{ { 1.0l, -2.0l }, { 3.0l, 4.0l} };
	constexpr cte::mat::Matrix div_result = matrix / 5.0l;
	constexpr cte::mat::Matrix expected_result{ { 0.2l, -0.4l }, { 0.6l, 0.8l } };
	EXPECT_EQ(div_result, expected_result);
}

