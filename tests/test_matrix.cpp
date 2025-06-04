#include <gtest/gtest.h>
#include <matrix.hpp>
#include <processing.hpp>

TEST(MatrixBasics, InitListConstructionAndAt) {
	iptt::Matrix m{{1, 2, 3}, {4, 5, 6}};
	EXPECT_EQ(m.rows(), 2u);
	EXPECT_EQ(m.cols(), 3u);
	EXPECT_EQ(m.at(1, 2), 6);
}

TEST(MatrixBasics, ConstCorrectness) {
	const iptt::Matrix m{{7}};
	EXPECT_EQ(m.at(0, 0), 7);
}

TEST(MatrixBasics, ViewConstruction) {
	iptt::Matrix base{{1, 2}, {3, 4}};
	iptt::MatrixView v(base, 0, 0, 0, 1);
	EXPECT_EQ(v.rows(), 1u);
	EXPECT_EQ(v.at(0, 1), 2);
}


TEST(MatrixOps, PlusMinusSameType) {
	const iptt::Matrix a{{1, 1}, {1, 1}};
	const iptt::Matrix b{{2, 2}, {2, 2}};
	iptt::Matrix<int> c = a + b;
	EXPECT_EQ(c.at(0, 0), 3);
	c = c - b;
	EXPECT_EQ(c.at(1, 1), 1);
}

TEST(MatrixOps, PlusMixMatrixAndView) {
	iptt::Matrix a{{1, 2}, {3, 4}};
	const iptt::MatrixView<int> v(a);
	const iptt::Matrix b{{4, 3}, {2, 1}};
	iptt::Matrix<int> res = v + b;
	EXPECT_EQ(res.at(0, 1), 5);
}

TEST(MatrixOps, ScalarMulDiv) {
	const iptt::Matrix a{{1, 2}, {3, 4}};
	iptt::Matrix<int> b = a * 10;
	EXPECT_EQ(b.at(1, 0), 30);
	iptt::Matrix<int> c = b / 5;
	EXPECT_EQ(c.at(1, 1), 8);
}

TEST(MatrixOps, MatrixMultiplication) {
	const iptt::Matrix a{{1, 2, 3}, {4, 5, 6}};
	const iptt::Matrix b{{7, 8}, {9, 10}, {11, 12}};
	iptt::Matrix<int> c = a * b;
	EXPECT_EQ(c.rows(), 2u);
	EXPECT_EQ(c.cols(), 2u);
	EXPECT_EQ(c.at(0, 0), 58);
	EXPECT_EQ(c.at(1, 1), 154);
}


TEST(MatrixMove, MoveCtorAndAssign) {
	iptt::Matrix a{{42}};
	iptt::Matrix b(std::move(a));
	EXPECT_EQ(b.at(0, 0), 42);
	iptt::Matrix c{{0}};
	c = std::move(b);
	EXPECT_EQ(c.at(0, 0), 42);
}

TEST(MatrixMove, ViewConstructionOnLvalueOnly) {
	static_assert(!std::is_constructible_v<iptt::MatrixView<int>, iptt::Matrix<int>&&>,
				  "iptt::MatrixView must not be constructible from rvalue iptt::Matrix");
}

TEST(MatrixAlgo, TransformSquared) {
	iptt::Matrix m{{1, 2}, {3, 4}};
	auto r = iptt::transform(m, [](const int x) { return x * x; });
	EXPECT_EQ(r.at(1, 1), 16);
}

TEST(MatrixThrow, ShapeMismatchAdd) {
	const iptt::Matrix a{{1, 2}, {3, 4}};
	const iptt::Matrix b{{1}};
	EXPECT_THROW(static_cast<void>(a + b), std::invalid_argument);
}

TEST(MatrixThrow, ViewOutOfBounds) {
	iptt::Matrix m{{1, 2}, {3, 4}};
	EXPECT_THROW(iptt::MatrixView<int>(m, 0, 0, 2, 1), std::out_of_range);
}

TEST(MatrixEdge, RvaluePlusView) {
	iptt::Matrix a{{1, 1}, {1, 1}};
	const iptt::MatrixView<int> v(a);
	iptt::Matrix<int> res = std::move(a) + v;
	EXPECT_EQ(res.at(0, 0), 2);
}

TEST(MatrixEdge, RvalueScalarMul) {
	iptt::Matrix a{{2}};
	iptt::Matrix<int> r = std::move(a) * 5;
	EXPECT_EQ(r.at(0, 0), 10);
}

TEST(MatrixMeta, RowsColsPreserved) {
	iptt::Matrix a{{1, 2}, {3, 4}};
	auto b = a * 2;
	EXPECT_EQ(b.rows(), 2u);
	EXPECT_EQ(b.cols(), 2u);
	iptt::Matrix<int> c = a + b;
	EXPECT_EQ(c.rows(), 2u);
	EXPECT_EQ(c.cols(), 2u);
	EXPECT_EQ(c.at(0, 0), 3);
	EXPECT_EQ(c.at(1, 1), 12);
}

TEST(MatrixInplace, PlusEqualMinusEqual) {
	iptt::Matrix a{{1, 2}, {3, 4}};
	iptt::Matrix b{{1, 1}, {1, 1}};
	a += b;
	EXPECT_EQ(a.at(1, 1), 5);
	a -= b;
	EXPECT_EQ(a.at(0, 0), 1);
}

TEST(MatrixViewOps, MoveByInsideBounds) {
	iptt::Matrix m{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
	iptt::MatrixView<int> v(m, 0, 0, 1, 1);
	v.move_by(1, 1);
	EXPECT_EQ(v.at(0, 0), 5);
}

TEST(MatrixViewOps, MoveByThrowsOutOfRange) {
	iptt::Matrix<int> m(2, 2);
	iptt::MatrixView<int> v(m);
	EXPECT_THROW(v.move_by(1, 0), std::out_of_range);
}

TEST(MatrixAlgo, TransformNegate) {
	iptt::Matrix m{{1, -2}};
	auto r = iptt::transform(m, [](const int x) { return -x; });
	EXPECT_EQ(r.at(0, 0), -1);
	EXPECT_EQ(r.at(0, 1), 2);
}

TEST(MatrixInplace, DivideEqual) {
	iptt::Matrix m{{10, 20}, {30, 40}};
	m /= 10;
	EXPECT_EQ(m.at(1, 0), 3);
}

TEST(MatrixUnary, UnaryMinus) {
	iptt::Matrix m{{1, -2}};
	const iptt::MatrixView<int> v(m);
	const iptt::MatrixView v2(v);
	iptt::Matrix neg = -m;
	EXPECT_EQ(neg.at(0, 0), -1);
	EXPECT_EQ(v.at(0, 0), -1);
	EXPECT_EQ(v2.at(0, 0), -1);
}
