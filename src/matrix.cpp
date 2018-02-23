
#include "matrix.hpp"
#include <cassert>

#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <algorithm>
#include <stdio.h>
#include <stdarg.h>

namespace ATP
{
	namespace Math
	{
		Matrix::Matrix(unsigned int rows, unsigned int cols)
		{
			resize(rows, cols);
		}

		Matrix::Matrix(const Matrix& src)
		{
			(*this) = src;
		}

		Matrix::~Matrix()
		{
			clear();
		}

		void Matrix::clear()
		{
			m_rows = 0;
			m_cols = 0;
			m_data.resize(0, 0);
		}

		double& Matrix::operator()(unsigned int row, unsigned int col)
		{
			return m_data[offset(row, col)];
		}

		double Matrix::operator()(unsigned int row, unsigned int col) const
		{
			return m_data[offset(row, col)];
		}

		unsigned int Matrix::offset(unsigned int row, unsigned int col) const
		{
			if (row >= m_rows || col >= m_cols)
			{
				throw SubscriptOutOfRange();
			}

			return m_dir == ROW_COL ? row * m_cols + col : col * m_rows + row;
		}

		unsigned int Matrix::rows() const
		{
			return m_rows;
		}

		unsigned int Matrix::cols() const
		{
			return m_cols;
		}

		void Matrix::resize(unsigned int rows, unsigned int cols)
		{
			clear();

			m_rows = rows;
			m_cols = cols;
			m_dir = ROW_COL;

			m_data.resize(m_rows*m_cols, 0.0f);
		}

		void Matrix::assertSameSize(const Matrix& other) const
		{
			if (other.m_rows != m_rows || other.m_cols != m_cols)
			{
				throw MismatchedDimensions();
			}
		}

		//Operators
		const Matrix Matrix::operator+(const Matrix& other) const
		{
			assertSameSize(other);

			Matrix sum(m_rows, m_cols);

			unsigned int row, col;
			for (row = 0; row<m_rows; row++)
			{
				for (col = 0; col<m_cols; col++)
				{
					sum(row, col) = (*this)(row, col) + other(row, col);
				}
			}

			return sum;
		}

		Matrix& Matrix::operator+=(const Matrix& other)
		{
			assertSameSize(other);

			unsigned int row, col;
			for (row = 0; row<m_rows; row++)
			{
				for (col = 0; col<m_cols; col++)
				{
					(*this)(row, col) += other(row, col);
				}
			}

			return *this;

		}

		const Matrix Matrix::operator-(const Matrix& other) const
		{
			assertSameSize(other);

			Matrix diff(m_rows, m_cols);

			unsigned int row, col;
			for (row = 0; row<m_rows; row++)
			{
				for (col = 0; col<m_cols; col++)
				{
					diff(row, col) = (*this)(row, col) - other(row, col);
				}
			}

			return diff;
		}

		Matrix& Matrix::operator-=(const Matrix& other)
		{
			assertSameSize(other);

			unsigned int row, col;
			for (row = 0; row<m_rows; row++)
			{
				for (col = 0; col<m_cols; col++)
				{
					(*this)(row, col) -= other(row, col);
				}
			}

			return *this;
		}

		const Matrix Matrix::operator*(const Matrix& other) const
		{
			if (m_cols != other.m_rows)
			{
				throw MismatchedDimensions();
			}

			Matrix product(m_rows, other.m_cols);

			for (unsigned int row = 0; row<m_rows; row++)
			{
				for (unsigned int col = 0; col<other.m_cols; col++)
				{
					double total = 0.0;

					unsigned int i;
					for (i = 0; i<m_cols; i++)
					{
						total += (*this)(row, i) * other(i, col);
					}

					product(row, col) = total;
				}
			}

			return product;
		}

		Matrix& Matrix::operator*=(const Matrix& other)
		{
			*this = (*this)*other;

			return *this;
		}

		const Matrix Matrix::operator*(const double d) const
		{
			Matrix m(*this);

			unsigned int row, col;
			for (row = 0; row<m_rows; row++)
			{
				for (col = 0; col<m_cols; col++)
				{
					m(row, col) *= d;
				}
			}

			return m;
		}

		Matrix& Matrix::operator*=(const double d)
		{
			*this = (*this)*d;

			return *this;
		}

		const Matrix Matrix::operator/(const Matrix& other) const
		{
			return (*this)*(other.transpose());
		}

		Matrix& Matrix::operator/=(const Matrix& other)
		{
			*this = (*this) / other;

			return *this;
		}

		const Matrix Matrix::operator/(const double d) const
		{
			Matrix m(*this);

			unsigned int row, col;
			for (row = 0; row<m_rows; row++)
			{
				for (col = 0; col<m_cols; col++)
				{
					m(row, col) /= d;
				}
			}

			return m;
		}

		Matrix& Matrix::operator/=(double d)
		{
			*this = (*this) / d;

			return *this;
		}

		Matrix& Matrix::operator=(const Matrix& other)
		{
			resize(other.m_rows, other.m_cols);
			m_dir = other.m_dir;
			m_data = other.m_data;

			return *this;
		}

		void Matrix::populate(std::vector<double> values)
		{
			assert(values.size() == m_rows * m_cols);

			for (unsigned int row = 0; row<m_rows; row++)
			{
				for (unsigned int col = 0; col<m_cols; col++)
				{
					(*this)(row, col) = values[row*m_cols + col];
				}
			}
		}

		Matrix Matrix::transpose() const
		{
			Matrix t(*this);

			switch (m_dir)
			{
			case ROW_COL: t.m_dir = COL_ROW; break;
			case COL_ROW: t.m_dir = ROW_COL; break;
			}

			t.m_rows = m_cols;
			t.m_cols = m_rows;

			return t;
		}

		Matrix Matrix::inverse() const
		{
			if (m_rows != m_cols)
			{
				throw NonInvertibleMatrix();
			}

			double det = determinant();

			if (det == 0.0)
			{
				throw SingularMatrix();
			}

			return cofactorMatrix().transpose() / det;
		}

		double Matrix::determinant() const
		{
			if (m_rows != m_cols)
			{
				throw NonSquareMatrix();
			}

			double a, b, c, d, e, f, g, h, i;
			switch (m_rows)
			{
				//Special cases (for speed)
			case 1:
				a = (*this)(0, 0);
				return a;
				break;
			case 2:
				a = (*this)(0, 0);
				b = (*this)(0, 1);
				c = (*this)(1, 0);
				d = (*this)(1, 1);
				return a * d - b * c;
				break;
			case 3:
				a = (*this)(0, 0);
				b = (*this)(0, 1);
				c = (*this)(0, 2);
				d = (*this)(1, 0);
				e = (*this)(1, 1);
				f = (*this)(1, 2);
				g = (*this)(2, 0);
				h = (*this)(2, 1);
				i = (*this)(2, 2);
				return a * e*i + b * f*g + c * d*h - g * e*c - h * f*a - i * d*b;
				break;

				//General case
			default:
				double total = 0.0;
				for (unsigned int row = 0; row<m_rows; row++)
				{
					total += (*this)(row, 0) * cofactor(row, 0);
				}

				return total;
				break;
			}
		}

		double Matrix::cofactor(unsigned int row, unsigned int col) const
		{
			if (m_rows != m_cols)
			{
				throw NonSquareMatrix();
			}

			if (row >= m_rows || col >= m_cols)
			{
				throw SubscriptOutOfRange();
			}

			Matrix minor(m_rows - 1, m_cols - 1);

			unsigned int minorRow, minorCol;
			for (minorRow = 0; minorRow<minor.rows(); minorRow++)
			{
				for (minorCol = 0; minorCol<minor.cols(); minorCol++)
				{
					unsigned int thisRow = (minorRow<row) ? minorRow : minorRow + 1;
					unsigned int thisCol = (minorCol<col) ? minorCol : minorCol + 1;
					minor(minorRow, minorCol) = (*this)(thisRow, thisCol);
				}
			}
			double det = minor.determinant();

			bool evenSum = (row + col) % 2 == 0;
			return evenSum ? det : -det;
		}

		Matrix Matrix::cofactorMatrix() const
		{
			Matrix cf(m_rows, m_cols);

			unsigned int row, col;
			for (row = 0; row<m_rows; row++)
			{
				for (col = 0; col<m_cols; col++)
				{
					cf(row, col) = cofactor(row, col);
				}
			}

			return cf;
		}

		//TODO:  May not be efficient for large matrices
		Matrix Matrix::solve(Matrix values)
		{
			if (m_rows != values.rows() || values.cols() != 1)
			{
				throw MismatchedDimensions();
			}

			if (m_rows != m_cols)
			{
				throw NonSquareMatrix();
			}

			double det;
			if ((det = determinant()) == 0.0f)
			{
				throw NonInvertibleMatrix();
			}

			Matrix results(values.rows(), 1);

			for (unsigned int row = 0; row<values.rows(); row++)
			{
				Matrix curValueMatrix(*this);

				for (unsigned int i = 0; i<values.rows(); i++)
				{
					curValueMatrix(i, 1) = values(i, 1);
				}
				results(row, 1) = curValueMatrix.determinant() / det;
			}

			return results;
		}
	}
}