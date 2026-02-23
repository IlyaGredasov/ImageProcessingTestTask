#include <gtest/gtest.h>
#include <matrix.hpp>
#include <processing.hpp>

using namespace iptt;

TEST(ProcessingBasics, ClampAndRGBOps) {
    EXPECT_EQ(clamp(-10), 0);
    EXPECT_EQ(clamp(0), 0);
    EXPECT_EQ(clamp(255), 255);
    EXPECT_EQ(clamp(300), 255);

    const RGB a{10, 20, 30};
    const RGB b{5, 10, 15};
    EXPECT_EQ(a + b, (RGB{15, 30, 45}));
    EXPECT_EQ(a - b, (RGB{5, 10, 15}));
    EXPECT_EQ(a * 2.0, (RGB{20, 40, 60}));
    EXPECT_EQ(a / 2.0, (RGB{5, 10, 15}));
    EXPECT_TRUE(a > b);
}

TEST(ProcessingBasics, BrightnessOrdering) {
    const RGB dark{0, 0, 0};
    const RGB bright{255, 255, 255};
    EXPECT_LT(brightness(dark), brightness(bright));
    EXPECT_TRUE(bright > dark);
}

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

TEST(Processing, DilateAndErose_Basic) {
    Matrix<int> input{{0, 0, 0}, {0, 255, 0}, {0, 0, 0}};
    auto dilated = dilate(input, 200);
    Matrix<int> expected_dilate{{255}};
    EXPECT_EQ(dilated, expected_dilate);

    auto eroded = erose(input, 200);
    Matrix<int> expected_erode{{0}};
    EXPECT_EQ(eroded, expected_erode);
}

TEST(Processing, Erose_AllOnes) {
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

TEST(Processing, ContrastStretch_EqualValues) {
    Matrix<RGB> input{{{10, 10, 10}, {10, 10, 10}}};
    auto stretched = contrast_stretch(input);
    EXPECT_EQ(stretched, input);
}

TEST(Processing, Convolute_Identity) {
    Matrix<int> input{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    Matrix<int> identity_kernel{{0, 0, 0}, {0, 1, 0}, {0, 0, 0}};
    auto result = convolute(input, identity_kernel);
    Matrix<int> expected{{5}};
    EXPECT_EQ(result, expected);
}

TEST(Processing, CropBorders_Basic) {
    Matrix<int> input{{0, 0, 0}, {0, 1, 0}, {0, 0, 0}};
    auto cropped = crop_borders(input, [](const int v) { return v != 0; });
    EXPECT_EQ(cropped.rows(), 1u);
    EXPECT_EQ(cropped.cols(), 1u);
    EXPECT_EQ(cropped.at(0, 0), 1);
}

TEST(Processing, CropBorders_AllZeroThrows) {
    Matrix<int> input{{0, 0}, {0, 0}};
    EXPECT_THROW(crop_borders(input, [](const int v) { return v != 0; }), std::invalid_argument);
}

TEST(Processing, RotateZeroKeepsImage) {
    Matrix<int> input{{1, 2}, {3, 4}};
    auto rotated = rotate(input, 0.0);
    EXPECT_EQ(rotated, input);
}

TEST(Processing, ScaleOneKeepsImage) {
    Matrix<int> input{{1, 2}, {3, 4}};
    auto scaled = scale(input, 1.0, 1.0);
    EXPECT_EQ(scaled, input);
}

TEST(Processing, FlipHorizontal) {
    Matrix<RGB> input{{{1, 0, 0}, {2, 0, 0}}};
    auto flipped = flip(input, true, false);
    EXPECT_EQ(flipped.rows(), 1u);
    EXPECT_EQ(flipped.cols(), 2u);
    EXPECT_EQ(flipped.at(0, 0).r, 2);
    EXPECT_EQ(flipped.at(0, 1).r, 1);
}
