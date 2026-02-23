#include <gtest/gtest.h>
#include <matrix.hpp>
#include <processing.hpp>
#include <type_traits>

using namespace iptt;

TEST(MatrixBasics, InitListAndAt) {
    Matrix m{{1, 2, 3}, {4, 5, 6}};
    EXPECT_EQ(m.rows(), 2u);
    EXPECT_EQ(m.cols(), 3u);
    EXPECT_EQ(m.at(1, 2), 6);
}

TEST(MatrixBasics, InitListMismatchThrows) { EXPECT_THROW((Matrix<int>{{1, 2}, {3}}), std::invalid_argument); }

TEST(MatrixBasics, ConstCorrectness) {
    const Matrix m{{7}};
    EXPECT_EQ(m.at(0, 0), 7);
}

TEST(MatrixBasics, AtOutOfRangeThrows) {
    Matrix<int> m(2, 2);
    EXPECT_THROW(m.at(2, 0), std::out_of_range);
    EXPECT_THROW(m.at(0, 2), std::out_of_range);
}

TEST(MatrixViewBasics, ViewConstructionAndAt) {
    Matrix base{{1, 2}, {3, 4}};
    MatrixView v(base, 0, 0, 0, 1);
    EXPECT_EQ(v.rows(), 1u);
    EXPECT_EQ(v.at(0, 1), 2);
}

TEST(MatrixViewBasics, ViewInvalidRangeThrows) {
    Matrix<int> m(2, 2);
    EXPECT_THROW(MatrixView<int>(m, 1, 1, 0, 0), std::invalid_argument);
}

TEST(MatrixViewBasics, ViewOutOfBoundsThrows) {
    Matrix m{{1, 2}, {3, 4}};
    EXPECT_THROW(MatrixView<int>(m, 0, 0, 2, 1), std::out_of_range);
}

TEST(MatrixViewBasics, ConstructFromConst) {
    const Matrix<int> m{{1, 2}, {3, 4}};
    const MatrixView<int> v(m);
    EXPECT_EQ(v.at(1, 0), 3);
}

TEST(MatrixViewOps, MoveByInsideBounds) {
    Matrix m{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    MatrixView<int> v(m, 0, 0, 1, 1);
    v.move_by(1, 1);
    EXPECT_EQ(v.at(0, 0), 5);
}

TEST(MatrixViewOps, MoveByThrowsOutOfRange) {
    Matrix<int> m(2, 2);
    MatrixView<int> v(m);
    EXPECT_THROW(v.move_by(1, 0), std::out_of_range);
}

TEST(MatrixOps, PlusMinusSameType) {
    const Matrix a{{1, 1}, {1, 1}};
    const Matrix b{{2, 2}, {2, 2}};
    Matrix<int> c = a + b;
    EXPECT_EQ(c.at(0, 0), 3);
    c = c - b;
    EXPECT_EQ(c.at(1, 1), 1);
}

TEST(MatrixOps, PlusMixMatrixAndView) {
    Matrix a{{1, 2}, {3, 4}};
    const MatrixView<int> v(a);
    const Matrix b{{4, 3}, {2, 1}};
    Matrix<int> res = v + b;
    EXPECT_EQ(res.at(0, 1), 5);
}

TEST(MatrixOps, ScalarMulDiv) {
    const Matrix a{{1, 2}, {3, 4}};
    Matrix<int> b = a * 10;
    EXPECT_EQ(b.at(1, 0), 30);
    Matrix<int> c = b / 5;
    EXPECT_EQ(c.at(1, 1), 8);
}

TEST(MatrixOps, ShapeMismatchThrows) {
    Matrix a{{1, 2}, {3, 4}};
    const Matrix b{{1}};
    EXPECT_THROW(static_cast<void>(a + b), std::invalid_argument);
    EXPECT_THROW(static_cast<void>(a -= b), std::invalid_argument);
    EXPECT_THROW(static_cast<void>(a += b), std::invalid_argument);
}

TEST(MatrixOps, MatrixMultiplication) {
    const Matrix a{{1, 2, 3}, {4, 5, 6}};
    const Matrix b{{7, 8}, {9, 10}, {11, 12}};
    Matrix<int> c = a * b;
    EXPECT_EQ(c.rows(), 2u);
    EXPECT_EQ(c.cols(), 2u);
    EXPECT_EQ(c.at(0, 0), 58);
    EXPECT_EQ(c.at(1, 1), 154);
}

TEST(MatrixOps, MatrixMultiplicationMismatchThrows) {
    const Matrix a{{1, 2}, {3, 4}};
    const Matrix b{{1, 2}};
    EXPECT_THROW(static_cast<void>(a * b), std::invalid_argument);
}

TEST(MatrixMove, MoveCtorAndAssignResizes) {
    Matrix<int> a{{1, 2, 3}, {4, 5, 6}};
    Matrix<int> b(std::move(a));
    EXPECT_EQ(b.rows(), 2u);
    EXPECT_EQ(b.cols(), 3u);
    EXPECT_EQ(b.at(1, 2), 6);

    Matrix<int> c{{0}};
    c = std::move(b);
    EXPECT_EQ(c.rows(), 2u);
    EXPECT_EQ(c.cols(), 3u);
    EXPECT_EQ(c.at(0, 0), 1);
}

TEST(MatrixMove, RvalueOperators) {
    Matrix<int> a{{1, 1}, {1, 1}};
    Matrix<int> b{{2, 2}, {2, 2}};
    Matrix<int> c = std::move(a) + b;
    EXPECT_EQ(c.at(0, 0), 3);

    Matrix<int> d = std::move(b) - Matrix<int>{{1, 1}, {1, 1}};
    EXPECT_EQ(d.at(0, 0), 1);

    Matrix<int> e = std::move(d) * 2;
    EXPECT_EQ(e.at(1, 1), 2);

    Matrix<int> f = std::move(e) / 2;
    EXPECT_EQ(f.at(1, 1), 1);
}

TEST(MatrixOps, AssignFromView) {
    Matrix<int> base{{1, 2}, {3, 4}};
    MatrixView<int> v(base);
    Matrix<int> out(2, 2);
    out = v;
    EXPECT_EQ(out.at(1, 0), 3);
}

TEST(MatrixUnary, UnaryMinus) {
    Matrix m{{1, -2}};
    Matrix neg = -m;
    EXPECT_EQ(neg.at(0, 0), -1);
    EXPECT_EQ(neg.at(0, 1), 2);
}

TEST(MatrixMeta, RowsColsPreserved) {
    Matrix a{{1, 2}, {3, 4}};
    auto b = a * 2;
    EXPECT_EQ(b.rows(), 2u);
    EXPECT_EQ(b.cols(), 2u);
    Matrix<int> c = a + b;
    EXPECT_EQ(c.rows(), 2u);
    EXPECT_EQ(c.cols(), 2u);
    EXPECT_EQ(c.at(0, 0), 3);
    EXPECT_EQ(c.at(1, 1), 12);
}

TEST(MatrixViewRules, ViewConstructionOnLvalueOnly) {
    static_assert(!std::is_constructible_v<MatrixView<int>, Matrix<int>&&>,
        "iptt::MatrixView must not be constructible from rvalue iptt::Matrix");
}
