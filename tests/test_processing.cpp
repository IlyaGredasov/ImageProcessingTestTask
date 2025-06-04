#include <gtest/gtest.h>
#include <matrix.hpp>
#include <processing.hpp>

using namespace iptt;

TEST(Processing, GaussianFilter_Basic) {
	const Matrix<double> input{{10, 10, 10}, {10, 10, 10}, {10, 10, 10}};
	auto result = gaussian_filter(input);
	Matrix<double> expected{{10}};
	EXPECT_EQ(result, expected);
}

TEST(Processing, SobelOperator_Basic) {
	Matrix<int> input{{0, 0, 0}, {0, 255, 0}, {0, 0, 0}};
	auto result = sobel_operator(input);
	Matrix<int> expected{{0}};
	EXPECT_EQ(result, expected);
}

TEST(Processing, Dilatate_Basic) {
	Matrix input{{0, 0, 0}, {0, 255, 0}, {0, 0, 0}};
	auto result = dilate(input);
	Matrix expected{{255}};
	EXPECT_EQ(result, expected);
}

TEST(Processing, Erose_Basic) {
	Matrix<int> input{{255, 255, 255}, {255, 255, 255}, {255, 255, 255}};
	auto result = erose(input);
	Matrix<int> expected{{255}};
	EXPECT_EQ(result, expected);
}

TEST(Processing, ContrastStretch_RGBRange) {
	Matrix<RGB> input{{{10, 20, 30}, {40, 50, 60}}};
	auto stretched = contrast_stretch(input);
	Matrix<RGB> expected{{{0, 0, 0}, {255, 255, 255}}};
	EXPECT_EQ(stretched, expected);
}

TEST(Processing, Convolute_Identity) {
	Matrix<int> input{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
	Matrix<int> identity_kernel{{0, 0, 0}, {0, 1, 0}, {0, 0, 0}};
	auto result = convolute(input, identity_kernel);
	Matrix<int> expected{{5}};
	EXPECT_EQ(result, expected);
}
