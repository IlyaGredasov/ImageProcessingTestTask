#pragma once
#include <iostream>
#include <stdexcept>
#include <vector>


namespace iptt { // сделаем namespace ImageProcessingTestTask, чтобы не засорять std
	template <typename T>
	class Matrix; // Класс матрицы

	template <typename T>
	class MatrixView; // Вьюшка, не владеющая ресурсами

	template <typename Derived, typename T>
	class GeneralMatrix; // Родительский класс

	template <typename M1, typename M2, typename T1, typename T2>
	static bool is_same_shape(const GeneralMatrix<M1, T1>& m1, const GeneralMatrix<M2, T2>& m2) {
		return m1.rows() == m2.rows() && m1.cols() == m2.cols();
	}

	// Используем паттерн CRTP
	template <typename Derived, typename T>
	class GeneralMatrix {
	public:
		static_assert(!std::is_const_v<T>, "T must not be const"); // Используем стандартные assert как в std::vector
		static_assert(!std::is_volatile_v<T>, "T must not be volatile");

		GeneralMatrix(const size_t rows, const size_t cols) : rows_(rows), cols_(cols) {}
		size_t rows() const { return rows_; }
		size_t cols() const { return cols_; }

		// Если мы работаем с lvalue, безопасно возвращать T&
		T& at(const size_t i, const size_t j) & { return static_cast<Derived&>(*this).at_impl(i, j); }
		const T& at(const size_t i, const size_t j) const& { return static_cast<const Derived&>(*this).at_impl(i, j); }

		T* data() { return static_cast<Derived&>(*this).data_impl(); }
		const T* data() const { return static_cast<const Derived&>(*this).data_impl(); }

		// Перемещаем тяжелые T без копирования
		T at(const size_t i, const size_t j) && { return std::move(static_cast<Derived&>(*this).at_impl(i, j)); }

		// Делаем шаблонное присваивание для возможности обмена данными между наследниками
		template <typename OtherDerived>
		Derived& operator=(const GeneralMatrix<OtherDerived, T>& other) {
			if (this->data() == other.data()) {
				return *static_cast<Derived*>(this);
			}
			if (!is_same_shape(*this, other)) {
				throw std::invalid_argument("Matrices have different shapes");
			}
			for (size_t i = 0; i < rows(); ++i)
				for (size_t j = 0; j < cols(); ++j)
					at(i, j) = other.at(i, j);
			return static_cast<Derived&>(*this);
		}

		template <typename OtherDerived>
		Derived& operator=(GeneralMatrix<OtherDerived, T>&& other) {
			if (this->data() == other.data()) {
				return *static_cast<Derived*>(this);
			}
			if (!is_same_shape(*this, other)) {
				throw std::invalid_argument("Matrices have different shapes");
			}
			for (size_t i = 0; i < rows(); ++i)
				for (size_t j = 0; j < cols(); ++j)
					at(i, j) = std::move(other.at(i, j));
			return static_cast<Derived&>(*this);
		}

		template <typename OtherDerived>
		Derived& operator+=(const GeneralMatrix<OtherDerived, T>& other) {
			if (!is_same_shape(*this, other)) {
				throw std::invalid_argument("Matrices have different shapes");
			}
			for (size_t i = 0; i < rows(); ++i) {
				for (size_t j = 0; j < cols(); ++j) {
					at(i, j) += other.at(i, j);
				}
			}
			return static_cast<Derived&>(*this);
		}

		template <typename OtherDerived>
		Derived& operator-=(const GeneralMatrix<OtherDerived, T>& other) {
			if (!is_same_shape(*this, other)) {
				throw std::invalid_argument("Matrices have different shape");
			}
			for (size_t i = 0; i < rows(); ++i) {
				for (size_t j = 0; j < cols(); ++j) {
					at(i, j) -= other.at(i, j);
				}
			}
			return static_cast<Derived&>(*this);
		}

		template <typename OtherDerived>
		Derived& operator*=(const GeneralMatrix<OtherDerived, T>& other) {
			if (cols() != other.rows()) {
				throw std::invalid_argument("Matrices can't be multiplied on each other");
			}
			Matrix<T> result(rows(), other.cols());
			for (size_t i = 0; i < rows(); ++i) {
				for (size_t j = 0; j < other.cols(); ++j) {
					T sum = T{};
					for (size_t k = 0; k < cols(); ++k) {
						sum += at(i, k) * other.at(k, j);
					}
					result.at(i, j) = sum;
				}
			}
			return static_cast<Derived&>(*this = std::move(result));
		}

		template <typename S, typename = std::enable_if_t<std::is_arithmetic_v<S>>>
		Derived& operator*=(const S scalar) {
			for (size_t i = 0; i < rows(); ++i) {
				for (size_t j = 0; j < cols(); ++j) {
					at(i, j) *= scalar;
				}
			}
			return static_cast<Derived&>(*this);
		}

		template <typename S, typename = std::enable_if_t<std::is_arithmetic_v<S>>>
		Derived& operator/=(const S scalar) {
			for (size_t i = 0; i < rows(); ++i) {
				for (size_t j = 0; j < cols(); ++j) {
					at(i, j) /= scalar;
				}
			}
			return static_cast<Derived&>(*this);
		}

		template <typename OtherDerived>
		bool operator==(const GeneralMatrix<OtherDerived, T>& other) const {
			if (!is_same_shape(*this, other)) {
				return false;
			}
			for (size_t i = 0; i < rows(); ++i) {
				for (size_t j = 0; j < cols(); ++j) {
					if (at(i, j) != other.at(i, j)) {
						return false;
					}
				}
			}
			return true;
		}

		Derived& operator-() { return static_cast<Derived&>(*this *= T(-1)); }
	protected:
		size_t rows_, cols_;
	};

	template <typename T>
	class Matrix final : public GeneralMatrix<Matrix<T>, T> {
	public:
		using Base = GeneralMatrix<Matrix, T>;
		friend Base; // Разрешаем GeneralMatrix доступ к имплементации at_impl()
		Matrix(const size_t rows, const size_t cols) : Base{rows, cols}, data_(rows * cols) {}
		Matrix(const std::initializer_list<std::initializer_list<T>> init) :
			Base(init.size(), init.begin()->size()), data_(this->rows() * this->cols()) {
			size_t i = 0;
			for (const auto& row : init) {
				if (row.size() != this->cols()) {
					throw std::invalid_argument("Lists in initializer list have different size");
				}
				std::copy(row.begin(), row.end(), data_.begin() + i * this->cols());
				++i;
			}
		}

		Matrix(const Matrix&) = default; // Так как мы определяем свой operator= через шаблонный в базовом классе,
										 // компилятор не сгенерирует конструкторы по умолчанию
		Matrix(Matrix&&) = default;

		template <typename OtherDerived> // Чтобы все было по open/closed, сделаем шаблонный ctor от всех наследников
		explicit Matrix(const GeneralMatrix<OtherDerived, T>& other) :
			Base{other.rows(), other.cols()}, data_(other.rows() * other.cols()) {
			for (size_t i = 0; i < other.rows(); ++i)
				for (size_t j = 0; j < other.cols(); ++j)
					this->at(i, j) = other.at(i, j);
		}

		Matrix& operator=(const Matrix& other) { return Base::template operator= <Matrix>(other); }
		Matrix& operator=(Matrix&& other) { return Base::template operator= <Matrix>(std::move(other)); }

		template <typename OtherDerived>
		Matrix& operator=(const GeneralMatrix<OtherDerived, T>& other) {
			return Base::template operator= <OtherDerived>(other);
		}

		template <typename OtherDerived>
		Matrix& operator=(GeneralMatrix<OtherDerived, T>&& other) {
			return Base::template operator= <OtherDerived>(std::move(other));
		}

		void swap(Matrix& other) noexcept {
			std::swap(this->rows_, other.rows_);
			std::swap(this->cols_, other.cols_);
			std::swap(this->data_, other.data_);
		}

	private:
		T& at_impl(const size_t i, const size_t j) {
			if (i >= this->rows() || j >= this->cols())
				throw std::out_of_range("Matrix index out of range");
			return data_[i * this->cols() + j];
		}
		const T& at_impl(const size_t i, const size_t j) const {
			if (i >= this->rows() || j >= this->cols())
				throw std::out_of_range("Matrix index out of range");
			return data_[i * this->cols() + j];
		}
		T* data_impl() { return data_.data(); }
		const T* data_impl() const { return data_.data(); }
		std::vector<T> data_{}; // Единый вектор более cache-friendly, чем T**
	};

	template <typename T>
	class MatrixView final : public GeneralMatrix<MatrixView<T>, T> {
	public:
		using Base = GeneralMatrix<MatrixView, T>;
		friend Base;

		// Единый шаблонный конструктор для всех наследников GeneralMatrix
		template <typename OtherDerived>
		MatrixView(GeneralMatrix<OtherDerived, T>& src, const size_t i0, const size_t j0, const size_t i1,
				   const size_t j1) :
			Base(i1 - i0 + 1, j1 - j0 + 1), data_(src.data()), i0_{i0}, j0_{j0}, parent_rows_{src.rows()},
			parent_cols_{src.cols()} {
			if (i0 > i1 || j0 > j1)
				throw std::invalid_argument("Invalid range of view");
			if (i1 >= src.rows() || j1 >= src.cols())
				throw std::out_of_range("MatrixView bounds out of range of original GeneralMatrix");
		}

		template <typename OtherDerived>
		MatrixView(const GeneralMatrix<OtherDerived, T>& src, const size_t i0, const size_t j0, const size_t i1,
				   const size_t j1) :
			Base(i1 - i0 + 1, j1 - j0 + 1), data_(const_cast<T*>(src.data())), i0_{i0}, j0_{j0},
			parent_rows_{src.rows()}, parent_cols_{src.cols()} {
			if (i0 > i1 || j0 > j1)
				throw std::invalid_argument("Invalid range of view");
			if (i1 >= src.rows() || j1 >= src.cols())
				throw std::out_of_range("MatrixView bounds out of range of original GeneralMatrix");
		}

		template <typename OtherDerived>
		explicit MatrixView(GeneralMatrix<OtherDerived, T>&& src, size_t i0 = 0, size_t j0 = 0, size_t i1 = 0,
							size_t j1 = 0) = delete; // Нельзя создать view на &&

		template <typename OtherDerived>
		explicit MatrixView(GeneralMatrix<OtherDerived, T>& src, size_t i0 = 0, size_t j0 = 0) :
			MatrixView(src, i0, j0, src.rows() - 1, src.cols() - 1) {}

		template <typename OtherDerived>
		explicit MatrixView(const GeneralMatrix<OtherDerived, T>& src, size_t i0 = 0, size_t j0 = 0) :
			MatrixView(src, i0, j0, src.rows() - 1, src.cols() - 1) {}

		MatrixView& operator=(const Matrix<T>& other) { return Base::template operator= <Matrix<T>>(other); }

		MatrixView& operator=(const MatrixView& other) { return Base::template operator= <MatrixView>(other); }

		void move_by(const int x,
					 const int y) { // Метод, чтобы сдвинуть view на вектор (x, y) (работает в матричных координатах)
			const int new_i0 = static_cast<int>(i0_) + x;
			const int new_j0 = static_cast<int>(j0_) + y;
			const int new_i1 = new_i0 + static_cast<int>(this->rows()) - 1;
			const int new_j1 = new_j0 + static_cast<int>(this->cols()) - 1;

			// Проверка на выход за границу
			if (new_i0 < 0 || new_j0 < 0 || new_i1 >= static_cast<int>(parent_rows_) ||
				new_j1 >= static_cast<int>(parent_cols_)) {
				throw std::out_of_range("MatrixView has been moved out of range");
			}

			// Обновляем положение
			i0_ = static_cast<size_t>(new_i0);
			j0_ = static_cast<size_t>(new_j0);
		}


	private:
		T& at_impl(const size_t i, const size_t j) {
			if (i >= this->rows_ || j >= this->cols_) {
				throw std::out_of_range("MatrixView index out of range");
			}
			return data_[(i + i0_) * parent_cols_ + (j + j0_)];
		}
		const T& at_impl(const size_t i, const size_t j) const {
			if (i >= this->rows_ || j >= this->cols_) {
				throw std::out_of_range("MatrixView index out of range");
			}
			return data_[(i + i0_) * parent_cols_ + (j + j0_)];
		}
		T* data_impl() { return data_; }
		const T* data_impl() const { return data_; }

		T* data_;
		size_t i0_, j0_;
		const size_t parent_rows_, parent_cols_; // Чтобы знать, как вычислить адрес в data()
	};

	template <typename M1, typename M2, typename T1, typename T2>
	Matrix<T1> operator+(const GeneralMatrix<M1, T1>& m1, const GeneralMatrix<M2, T2>& m2) {
		if (!is_same_shape(m1, m2)) {
			throw std::invalid_argument("Matrices have different shape");
		}
		Matrix<T1> result(m1.rows(), m1.cols());
		for (size_t i = 0; i < result.rows(); ++i)
			for (size_t j = 0; j < result.cols(); ++j)
				result.at(i, j) = m1.at(i, j) + m2.at(i, j);
		return result;
	}

	template <typename M1, typename M2, typename T1, typename T2>
	Matrix<T1> operator+(GeneralMatrix<M1, T1>&& m1, const GeneralMatrix<M2, T2>& m2) { // Можно переиспользовать rvalue
		return m1 += m2;
	}

	template <typename M1, typename M2, typename T1, typename T2>
	Matrix<T1> operator+(const GeneralMatrix<M1, T1>& m1, GeneralMatrix<M2, T2>&& m2) { // Можно переиспользовать rvalue
		return m2 += m1; // Считаем, что сложение коммутативно
	}

	template <typename M1, typename M2, typename T1, typename T2>
	Matrix<T1> operator+(GeneralMatrix<M1, T1>&& m1, GeneralMatrix<M2, T2>&& m2) { // Можно переиспользовать rvalue
		return m1 += m2;
	}

	template <typename M1, typename M2, typename T1, typename T2>
	Matrix<T1> operator-(const GeneralMatrix<M1, T1>& m1, const GeneralMatrix<M2, T2>& m2) {
		if (!is_same_shape(m1, m2)) {
			throw std::invalid_argument("Matrices have different shape");
		}
		Matrix<T1> result(m1.rows(), m1.cols());
		for (size_t i = 0; i < result.rows(); ++i)
			for (size_t j = 0; j < result.cols(); ++j)
				result.at(i, j) = m1.at(i, j) - m2.at(i, j);
		return result;
	}

	template <typename M1, typename M2, typename T1, typename T2>
	Matrix<T1> operator-(GeneralMatrix<M1, T1>&& m1, const GeneralMatrix<M2, T2>& m2) { // Можно переиспользовать rvalue
		return -(m2 -= m1);
	}

	template <typename M1, typename M2, typename T1, typename T2>
	Matrix<T1> operator-(const GeneralMatrix<M1, T1>& m1, GeneralMatrix<M2, T2>&& m2) { // Можно переиспользовать rvalue
		return m1 -= m2;
	}

	template <typename M1, typename M2, typename T1, typename T2>
	Matrix<T1> operator-(const GeneralMatrix<M1, T1>&& m1,
						 GeneralMatrix<M2, T2>&& m2) { // Можно переиспользовать rvalue
		return m1 -= m2;
	}

	template <typename M1, typename M2, typename T1, typename T2>
	Matrix<T1> operator*(const GeneralMatrix<M1, T1>& m1, const GeneralMatrix<M2, T2>& m2) {
		if (m1.cols() != m2.rows()) {
			throw std::invalid_argument("Matrices can't be multiplied on each other");
		}
		Matrix<T1> result(m1.rows(), m2.cols());
		for (size_t i = 0; i < m1.rows(); ++i) {
			for (size_t j = 0; j < m2.cols(); ++j) {
				T1 sum = T1{};
				for (size_t k = 0; k < m1.cols(); ++k) {
					sum += m1.at(i, k) * m2.at(k, j);
				}
				result.at(i, j) = sum;
			}
		}

		return result;
	}

	template <typename M, typename S, typename T, typename = std::enable_if_t<std::is_arithmetic_v<S>>>
	Matrix<T> operator*(const GeneralMatrix<M, T>& m1, const S scalar) {
		Matrix<T> result(m1.rows(), m1.cols());
		for (size_t i = 0; i < result.rows(); ++i) {
			for (size_t j = 0; j < result.cols(); ++j) {
				result.at(i, j) = m1.at(i, j) * scalar;
			}
		}
		return result;
	}

	template <typename M, typename S, typename T, typename = std::enable_if_t<std::is_arithmetic_v<S>>>
	Matrix<T> operator*(GeneralMatrix<M, T>&& m1, const S scalar) { // Можно умножить in-place
		return m1 *= scalar;
	}

	template <typename M, typename S, typename T, typename = std::enable_if_t<std::is_arithmetic_v<S>>>
	Matrix<T> operator/(const GeneralMatrix<M, T>& m1, const S scalar) {
		Matrix<T> result(m1.rows(), m1.cols());
		for (size_t i = 0; i < result.rows(); ++i) {
			for (size_t j = 0; j < result.cols(); ++j) {
				result.at(i, j) = m1.at(i, j) / scalar;
			}
		}
		return result;
	}

	template <typename M, typename S, typename T, typename = std::enable_if_t<std::is_arithmetic_v<S>>>
	Matrix<T> operator/(GeneralMatrix<M, T>&& m1, const S scalar) {
		return m1 /= scalar;
	}

	template <typename M, typename T>
	Matrix<T> operator-(const GeneralMatrix<M, T>& m1) {
		return m1 * -1;
	}
} // namespace iptt
