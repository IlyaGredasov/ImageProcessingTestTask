#pragma once

#include "matrix.hpp"
#include <array>
#include <cmath>

namespace iptt { // При работе с матрицами переводим их в математические координаты
struct Box {
    double min_x, min_y, max_x, max_y;
};
struct RGB {
    unsigned char r = 0, g = 0, b = 0;
};
struct AffineGeometry { // Смещение и размер при аффинном преобразовании
    size_t rows{}, cols{};
    double offsetX{}, offsetY{};
};
inline unsigned char clamp(const int value) noexcept {
    return static_cast<unsigned char>(std::max(0, std::min(255, value)));
}
inline RGB operator+(const RGB& a, const RGB& b) noexcept {
    return {clamp(a.r + b.r), clamp(a.g + b.g), clamp(a.b + b.b)};
}

inline RGB operator-(const RGB& a, const RGB& b) noexcept {
    return {clamp(a.r - b.r), clamp(a.g - b.g), clamp(a.b - b.b)};
}
inline RGB operator*(const RGB& a, const double scalar) noexcept {
    return {clamp(static_cast<int>(a.r * scalar)), clamp(static_cast<int>(a.g * scalar)),
        clamp(static_cast<int>(a.b * scalar))};
}

inline RGB operator/(const RGB& a, const double scalar) noexcept { return scalar != 0.f ? a * (1.f / scalar) : RGB{}; }
inline RGB operator*(const float scalar, const RGB& a) noexcept { return a * scalar; }
inline RGB& operator+=(RGB& a, const RGB& b) noexcept {
    a = a + b;
    return a;
}

inline RGB& operator-=(RGB& a, const RGB& b) noexcept {
    a = a - b;
    return a;
}

inline RGB& operator*=(RGB& a, const double scalar) noexcept {
    a = a * scalar;
    return a;
}

inline RGB& operator/=(RGB& a, const double scalar) noexcept {
    a = a / scalar;
    return a;
}
inline bool operator==(const RGB& a, const RGB& b) noexcept { return a.r == b.r && a.g == b.g && a.b == b.b; }

inline bool operator!=(const RGB& a, const RGB& b) noexcept { return !(a == b); }

inline int brightness(const RGB& c) noexcept { return static_cast<int>(0.299 * c.r + 0.587 * c.g + 0.114 * c.b); }

inline bool operator<(const RGB& a, const RGB& b) noexcept { return brightness(a) < brightness(b); }

inline bool operator>(const RGB& a, const RGB& b) noexcept { return brightness(a) > brightness(b); }

template <typename M, typename T, typename Func>
Matrix<std::invoke_result_t<Func, T>> transform(const GeneralMatrix<M, T>& m, const Func& function) {
    Matrix<std::invoke_result_t<Func, T>> result(m.rows(), m.cols());
    for (size_t i = 0; i < m.rows(); ++i) {
        for (size_t j = 0; j < m.cols(); ++j) {
            result.at(i, j) = function(m.at(i, j));
        }
    }
    return result;
}

inline Box compute_bbox(const size_t h, const size_t w,
    const Matrix<double>& A) { // Возвращает box с координатами углов матрицы размера (h, w) после умножения на A
    const double max_x = w > 0 ? static_cast<double>(w - 1) : 0.0;
    const double max_y = h > 0 ? static_cast<double>(h - 1) : 0.0;
    const std::array<std::array<double, 2>, 4> corners{{{0.0, 0.0}, {max_x, 0.0}, {0.0, max_y}, {max_x, max_y}}};

    Box b{1e30, 1e30, -1e30, -1e30};

    for (const auto [x, y] : corners) {
        const Matrix p = {{x}, {y}, {1.0}};
        Matrix<double> q = A * p;
        b.min_x = std::min(b.min_x, q.at(0, 0));
        b.max_x = std::max(b.max_x, q.at(0, 0));
        b.min_y = std::min(b.min_y, q.at(1, 0));
        b.max_y = std::max(b.max_y, q.at(1, 0));
    }
    return b;
}

template <typename M, typename T, typename Func, typename = std::enable_if_t<std::is_invocable_r_v<bool, Func, T>>>
Matrix<T> crop_borders(const GeneralMatrix<M, T>& m, Func&& crop_while_predicate) {
    // Срезаем границы, пока предикат истинен
    const size_t rows = m.rows();
    const size_t cols = m.cols();

    size_t top = 0;
    size_t bottom = rows;
    size_t left = 0;
    size_t right = cols;

    while (top < bottom) {
        bool non_zero = false;
        for (size_t j = 0; j < cols; ++j) {
            if (crop_while_predicate(m.at(top, j))) {
                non_zero = true;
                break;
            }
        }
        if (non_zero)
            break;
        ++top;
    }

    while (bottom > top) {
        bool non_zero = false;
        for (size_t j = 0; j < cols; ++j) {
            if (crop_while_predicate(m.at(bottom - 1, j))) {
                non_zero = true;
                break;
            }
        }
        if (non_zero)
            break;
        --bottom;
    }

    while (left < right) {
        bool non_zero = false;
        for (size_t i = top; i < bottom; ++i) {
            if (crop_while_predicate(m.at(i, left))) {
                non_zero = true;
                break;
            }
        }
        if (non_zero)
            break;
        ++left;
    }

    while (right > left) {
        bool non_zero = false;
        for (size_t i = top; i < bottom; ++i) {
            if (crop_while_predicate(m.at(i, right - 1))) {
                non_zero = true;
                break;
            }
        }
        if (non_zero)
            break;
        --right;
    }

    const size_t out_rows = bottom - top;
    const size_t out_cols = right - left;

    Matrix<T> result(out_rows, out_cols);
    result = MatrixView<T>(m, top, left, bottom - 1, right - 1);
    return result;
}

inline AffineGeometry make_output_geometry(const size_t rows, const size_t cols, const Matrix<double>& A) {
    const auto [min_x, min_y, max_x, max_y] = compute_bbox(rows, cols, A);

    const auto minX = static_cast<int>(std::floor(min_x));
    const auto maxX = static_cast<int>(std::ceil(max_x));
    const auto minY = static_cast<int>(std::floor(min_y));
    const auto maxY = static_cast<int>(std::ceil(max_y));

    AffineGeometry g;
    g.rows = maxY - minY + 1;
    g.cols = maxX - minX + 1;
    g.offsetX = minX;
    g.offsetY = minY;
    return g;
}

template <typename M, typename T>
Matrix<T> warp_affine(const GeneralMatrix<M, T>& src, const Matrix<double>& inv_A, const AffineGeometry& out_geom) {
    Matrix<T> dst(out_geom.rows, out_geom.cols);

    for (size_t row = 0; row < dst.rows(); ++row)
        for (size_t col = 0; col < dst.cols(); ++col) {
            // Переводим индексы (row, col) -> декартовы (x,y)
            const auto x = static_cast<double>(col) + out_geom.offsetX;
            const auto y = static_cast<double>(row) + out_geom.offsetY;
            // Идём обратно в оригинал
            const Matrix p{{x}, {y}, {1.0}};
            const Matrix<double> q = inv_A * p;

            const double srcX = q.at(0, 0);
            const double srcY = q.at(1, 0);

            const auto si = static_cast<int>(std::round(srcY));
            const auto sj = static_cast<int>(std::round(srcX));

            if (si >= 0 && si < static_cast<int>(src.rows()) && sj >= 0 && sj < static_cast<int>(src.cols()))
                dst.at(row, col) = src.at(si, sj);
            else
                dst.at(row, col) = T{};
        }
    return dst;
}

template <typename M, typename T>
static Matrix<T> affine_transform(const GeneralMatrix<M, T>& src, const Matrix<double>& A,
    const Matrix<double>& inv_A) {
    const auto geom = make_output_geometry(src.rows(), src.cols(), A);
    return warp_affine(src, inv_A, geom);
}

inline Matrix<double> shift_matrix(const double x, const double y) {
    return Matrix<double>{{1, 0, x}, {0, 1, y}, {0, 0, 1}};
}

inline Matrix<double> rotate_matrix(const double theta) {
    return Matrix<double>{{std::cos(theta), -std::sin(theta), 0}, {std::sin(theta), std::cos(theta), 0}, {0, 0, 1}};
}

template <typename M, typename T>
Matrix<T> rotate(const GeneralMatrix<M, T>& src, const double theta) {
    const double cx = (src.cols() - 1) / 2.0;
    const double cy = (src.rows() - 1) / 2.0;

    Matrix<double> A = shift_matrix(cx, cy) * rotate_matrix(theta) * shift_matrix(-cx, -cy);

    Matrix<double> inv_A = shift_matrix(cx, cy) * rotate_matrix(-theta) * shift_matrix(-cx, -cy);

    return affine_transform(src, A, inv_A);
}

inline Matrix<double> scale_matrix(const double sx, const double sy) {
    return Matrix<double>{{sx, 0, 0}, {0, sy, 0}, {0, 0, 1}};
}

template <typename M, typename T>
Matrix<T> scale(const GeneralMatrix<M, T>& src, const double sx, const double sy) {
    const double cx = (src.cols() - 1) / 2.0;
    const double cy = (src.rows() - 1) / 2.0;

    Matrix<double> A = shift_matrix(cx, cy) * scale_matrix(sx, sy) * shift_matrix(-cx, -cy);

    Matrix<double> inv_A = shift_matrix(cx, cy) * scale_matrix(1.0 / sx, 1.0 / sy) * shift_matrix(-cx, -cy);

    return affine_transform(src, A, inv_A);
}

inline Matrix<double> flip_matrix(const bool horizontal, const bool vertical) {
    const double x = horizontal ? -1 : 1;
    const double y = vertical ? -1 : 1;
    return Matrix<double>{{x, 0, 0}, {0, y, 0}, {0, 0, 1}};
}

template <typename M, typename T>
Matrix<T> flip(const GeneralMatrix<M, T>& src, const bool horizontal, const bool vertical) {
    const double cx = (src.cols() - 1) / 2.0;
    const double cy = (src.rows() - 1) / 2.0;
    const Matrix A = shift_matrix(cx, cy) * flip_matrix(horizontal, vertical) * shift_matrix(-cx, -cy);
    return affine_transform(src, A, A);
}

template <typename M1, typename M2, typename T1, typename T2>
T1 dot_product(const GeneralMatrix<M1, T1>& m1,
    const GeneralMatrix<M2, T2>& m2) { // Поэлементное скалярное произведение
    if (!is_same_shape(m1, m2)) {
        throw std::invalid_argument("Matrices have different shapes");
    }
    T1 sum = T1{};
    for (size_t i = 0; i < m1.rows(); ++i) {
        for (size_t j = 0; j < m1.cols(); ++j) {
            sum += m1.at(i, j) * m2.at(i, j);
        }
    }
    return sum;
}

template <typename M1, typename M2, typename T1, typename T2>
Matrix<T1> convolute(const GeneralMatrix<M1, T1>& m, const GeneralMatrix<M2, T2>& kernel) {
    const auto result_rows = m.rows() - kernel.rows() + 1;
    const auto result_cols = m.cols() - kernel.cols() + 1;
    Matrix<T1> result(result_rows, result_cols);
    MatrixView<T1> v(m, 0, 0, kernel.rows() - 1, kernel.cols() - 1);
    for (int i = 0; i < result_rows; ++i) {
        for (int j = 0; j < result_cols; ++j) {
            result.at(i, j) = dot_product(v, kernel);
            if (j < result_cols - 1) {
                v.move_by(0, 1);
            }
        }
        if (i < result_rows - 1) {
            v.move_by(1, -result_cols + 1);
        }
    }
    return result;
}

template <typename M, typename T>
Matrix<T> gaussian_filter(const GeneralMatrix<M, T>& m) {
    return convolute(m, Matrix<double>{{1, 2, 1}, {2, 4, 2}, {1, 2, 1}} / 16.0);
}

template <typename M>
Matrix<int> sobel_operator(const GeneralMatrix<M, int>& m) {
    const auto sobel_x = Matrix{{-1, 0, 1}, {-2, 0, 2}, {-1, 0, 1}};
    const auto sobel_y = Matrix{{-1, -2, -1}, {0, 0, 0}, {1, 2, 1}};
    const auto gx = convolute(m, sobel_x);
    const auto gy = convolute(m, sobel_y);
    Matrix<double> grad(gx.rows(), gx.cols());
    for (size_t i = 0; i < grad.rows(); ++i)
        for (size_t j = 0; j < grad.cols(); ++j)
            grad.at(i, j) = std::sqrt(gx.at(i, j) * gx.at(i, j) + gy.at(i, j) * gy.at(i, j));
    return transform(grad, [](const int x) { return static_cast<int>(x); });
}
template <typename M>
Matrix<int> dilate(const GeneralMatrix<M, int>& m, int thresh = 150) {
    const auto bin = transform(m, [thresh](const int x) { return static_cast<int>(x > thresh); });

    const auto conv = convolute(bin, Matrix{{1, 1, 1}, {1, 1, 1}, {1, 1, 1}});

    Matrix<int> out(conv.rows(), conv.cols());
    for (size_t i = 0; i < out.rows(); ++i)
        for (size_t j = 0; j < out.cols(); ++j)
            out.at(i, j) = (conv.at(i, j) > 0) ? 255 : 0;

    return out;
}

template <typename M>
Matrix<int> erose(const GeneralMatrix<M, int>& m, int thresh = 150) {
    const auto bin = transform(m, [thresh](const int x) { return static_cast<int>(x > thresh); });
    const auto conv = convolute(bin, Matrix{{1, 1, 1}, {1, 1, 1}, {1, 1, 1}});
    Matrix<int> out(conv.rows(), conv.cols());
    for (size_t i = 0; i < out.rows(); ++i)
        for (size_t j = 0; j < out.cols(); ++j) {
            out.at(i, j) = conv.at(i, j) == 9 ? 255 : 0;
        }
    return out;
}
template <typename M>
Matrix<int> grayscale(const GeneralMatrix<M, RGB>& m) {
    return iptt::transform(m, brightness);
}
template <typename M>
Matrix<RGB> contrast_stretch(const GeneralMatrix<M, RGB>& m) {
    RGB min_val = m.at(0, 0);
    RGB max_val = min_val;

    for (size_t i = 0; i < m.rows(); ++i) {
        for (size_t j = 0; j < m.cols(); ++j) {
            const RGB& px = m.at(i, j);
            min_val = std::min(min_val, px);
            max_val = std::max(max_val, px);
        }
    }

    if (min_val == max_val) {
        return static_cast<Matrix<RGB>>(m);
    }

    Matrix<RGB> result(m.rows(), m.cols());

    const auto r = (max_val.r == min_val.r);
    const auto g = (max_val.g == min_val.g);
    const auto b = (max_val.b == min_val.b);
    for (size_t i = 0; i < m.rows(); ++i) {
        for (size_t j = 0; j < m.cols(); ++j) {
            const RGB& src = m.at(i, j);
            RGB res{};
            res.r = r ? 0 : clamp((src.r - min_val.r) * 255 / (max_val.r - min_val.r));
            res.g = g ? 0 : clamp((src.g - min_val.g) * 255 / (max_val.g - min_val.g));
            res.b = b ? 0 : clamp((src.b - min_val.b) * 255 / (max_val.b - min_val.b));

            result.at(i, j) = res;
        }
    }

    return result;
}
} // namespace iptt
