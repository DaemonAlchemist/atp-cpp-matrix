#ifndef MATH_MATRIX
#define MATH_MATRIX

#ifndef ATP_MATRIX_LIB
#pragma comment(lib, "atp-matrix.lib")
#endif

#include <vector>
#include <exception>

#include <math.h>
#include <iostream>

namespace ATP
{
	namespace Math
	{
		class Matrix
		{
		public:
			class Exception : public std::exception {};
			class MismatchedDimensions : public Exception {};
			class SingularMatrix : public Exception {};
			class NonInvertibleMatrix : public Exception {};
			class NonSquareMatrix : public Exception {};
			class SubscriptOutOfRange : public Exception {};

			Matrix(unsigned int rows, unsigned int cols);
			Matrix(const Matrix& src);
			~Matrix();

			double& operator()(unsigned int row, unsigned int col);
			double operator()(unsigned int row, unsigned int col) const;
			unsigned int offset(unsigned int row, unsigned int col) const;

			unsigned int rows() const;
			unsigned int cols() const;
			void resize(unsigned int rows, unsigned int cols);

			void assertSameSize(const Matrix& other) const;

			//Operators
			const Matrix operator+(const Matrix& other) const;
			Matrix& operator+=(const Matrix& other);

			const Matrix operator-(const Matrix& other) const;
			Matrix& operator-=(const Matrix& other);

			const Matrix operator*(const Matrix& other) const;
			Matrix& operator*=(const Matrix& other);
			const Matrix operator*(const double d) const;
			Matrix& operator*=(const double d);

			const Matrix operator/(const Matrix& other) const;
			Matrix& operator/=(const Matrix& other);
			const Matrix operator/(const double d) const;
			Matrix& operator/=(const double d);

			Matrix& operator=(const Matrix& other);

			void populate(std::vector<double> values);

			//Functions
			Matrix transpose() const;
			Matrix inverse() const;
			double determinant() const;
			double cofactor(unsigned int row, unsigned int col) const;
			Matrix cofactorMatrix() const;

			//Simultaneous equations solving
			Matrix solve(Matrix values);

		private:
			void clear();

			enum DataDirection { ROW_COL, COL_ROW };

			std::vector<double> m_data;
			unsigned int m_rows, m_cols;
			DataDirection m_dir;
		};
	}
}

#endif
