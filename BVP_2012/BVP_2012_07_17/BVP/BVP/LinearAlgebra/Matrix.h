#pragma once
#include <array>
#include <stdlib.h>
#include <time.h>
#include "../Utils/AuxUtils.h"

namespace LinAlg
{
	/// <summary>
	/// Simple implementation of a matrix class
	/// </summary>
	template <class R, int RowDim, int ColDim>
	class Matrix
	{
		static_assert(RowDim > 0 && ColDim > 0, "Invalid matrix dimensions");
	private:
		std::array<std::array<R, ColDim>, RowDim> data{};
	public:

		/// <summary>
		/// Constructor from 2d array
		/// </summary>
		/// <param name="d"></param>
		Matrix(const std::array<std::array<R, ColDim>, RowDim>& d) : data{d}
		{}

		/// <summary>
		/// Default constructor
		/// </summary>
		Matrix() = default;

		/// <summary>
		/// Sub-script operator
		/// </summary>
		std::array<R, ColDim>& operator [](const int i)
		{
			return data[i];
		}

		/// <summary>
		/// Sub-script operator (const)
		/// </summary>
		const std::array<R, ColDim>& operator [](const int i) const
		{
			return data[i];
		}

		/// <summary>
		/// += operator
		/// </summary>
		Matrix<R, RowDim, ColDim>& operator +=(const Matrix<R, RowDim, ColDim>& rhs)
		{
			for (int i = 0; i < RowDim; i++)
				for (int j = 0; j < ColDim; j++)
					data[i][j] += rhs.data[i][j];

			return *this;
		}

		/// <summary>
		/// -= operator
		/// </summary>
		Matrix<R, RowDim, ColDim>& operator -=(const Matrix<R, RowDim, ColDim>& rhs)
		{
			for (int i = 0; i < RowDim; i++)
				for (int j = 0; j < ColDim; j++)
					data[i][j] -= rhs.data[i][j];

			return *this;
		}

		/// <summary>
		/// *= operator (matrix-scalar in place multiplication)
		/// </summary>
		Matrix<R, RowDim, ColDim>& operator *=(const R& rhs)
		{
			for (int i = 0; i < RowDim; i++)
				for (int j = 0; j < ColDim; j++)
					data[i][j] *= rhs;

			return *this;
		}

		/// <summary>
		/// /= operator (matrix-scalar in place division)
		/// </summary>
		Matrix<R, RowDim, ColDim>& operator /=(const R& rhs)
		{
			return *this *= R(1)/rhs;
		}

		/// <summary>
		/// Returns determinant of the matrix
		/// </summary>
		R Determinant() const
		{
			static_assert(RowDim == ColDim && RowDim == 2, "Only 2x2 matrices are supported");
			return data[0][0] * data[1][1] - data[1][0] * data[0][1];
		}

		/// <summary>
		/// Returns inverse matrix
		/// </summary>
		Matrix<R, RowDim, ColDim> Inverse() const
		{
			static_assert(RowDim == ColDim && RowDim == 2, "Only 2x2 matrices are supported");

			const auto one_over_determinant = R(1) / Determinant();

			Matrix<R, RowDim, ColDim> result;

			result[0][0] =  one_over_determinant * data[1][1];
			result[0][1] = -one_over_determinant * data[0][1];
			result[1][0] = -one_over_determinant * data[1][0];
			result[1][1] =  one_over_determinant * data[0][0];

			return result;
		}

		/// <summary>
		/// Returns transposed matrix
		/// </summary>
		Matrix<R, ColDim, RowDim> Transpose() const
		{
			Matrix<R, ColDim, RowDim> result;

			for (int i = 0; i < RowDim; i++)
				for (int j = 0; j < ColDim; j++)
					result[j][i] = data[i][j];

			return result;
		}

		/// <summary>
		/// Randomize;
		/// </summary>
		static void Randomize()
		{
			srand(time(NULL));
		}

		/// <summary>
		/// Returns randomly generated matrix of the current dimensions
		/// </summary>
		static Matrix<R, RowDim, ColDim> Random()
		{
			Matrix<R, RowDim, ColDim> result{};

			for (int i = 0; i < RowDim; i++)
				for (int j = 0; j < ColDim; j++)
					result[i][j] = ((R)rand()) / (R)RAND_MAX;

			return result;
		}

		/// <summary>
		/// Euclidean notm of the matrix
		/// </summary>
		R Norm() const
		{
			R norm{};

			for (int i = 0; i < RowDim; i++)
				for (int j = 0; j < ColDim; j++)
					norm += data[i][j] * data[i][j];

			return auxutils::Sqrt(norm);

		}

		/// <summary>
		/// Returns identity matrix
		/// </summary>
		static Matrix<R, RowDim, ColDim> Identity()
		{
			static_assert(RowDim == ColDim, "Identity matrix is square!");

			Matrix<R, RowDim, ColDim> identity{};

			for (int i = 0; i < RowDim; i++)
				identity[i][i] = R(1);

			return identity;
		}
	};

	/// <summary>
	/// + operation
	/// </summary>
	template <class R, int RowDim, int ColDim>
	Matrix<R, RowDim, ColDim> operator +(Matrix<R, RowDim, ColDim> lhs, const Matrix<R, RowDim, ColDim>& rhs)
	{
		return lhs += rhs;
	}

	/// <summary>
	/// - operation
	/// </summary>
	template <class R, int RowDim, int ColDim>
	Matrix<R, RowDim, ColDim> operator -(Matrix<R, RowDim, ColDim> lhs, const Matrix<R, RowDim, ColDim>& rhs)
	{
		return lhs -= rhs;
	}

	/// <summary>
	/// matrix-scalar multiplication
	/// </summary>
	template <class R, int RowDim, int ColDim>
	Matrix<R, RowDim, ColDim> operator *(Matrix<R, RowDim, ColDim> lhs, const R& rhs)
	{
		return lhs *= rhs;
	}

	/// <summary>
	/// scalar-matrix multiplication
	/// </summary>
	template <class R, int RowDim, int ColDim>
	Matrix<R, RowDim, ColDim> operator *(const R& lhs, Matrix<R, RowDim, ColDim> rhs)
	{
		return rhs *= lhs;
	}

	/// <summary>
	/// matrix-scalar division
	/// </summary>
	template <class R, int RowDim, int ColDim>
	Matrix<R, RowDim, ColDim> operator /(Matrix<R, RowDim, ColDim> lhs, const R& rhs)
	{
		return lhs /= rhs;
	}

	/// <summary>
	/// scalar-matrix division
	/// </summary>
	template <class R, int RowDim, int ColDim>
	Matrix<R, RowDim, ColDim> operator /(const R& lhs, Matrix<R, RowDim, ColDim> rhs)
	{
		return rhs /= lhs;
	}

	/// <summary>
	/// matrix-matrix multiplication
	/// </summary>
	template <class R, int RowDim, int ColDim, int ColDim2>
	Matrix<R, RowDim, ColDim2> operator *(const Matrix<R, RowDim, ColDim>& lhs, const Matrix<R, ColDim, ColDim2>& rhs)
	{
		Matrix<R, RowDim, ColDim2> result{};
		for (int i = 0; i < RowDim; i++)
			for (int j = 0; j < ColDim2; j++)
				for (int k = 0; k < ColDim; k++)
					result[i][j] += lhs[i][k] * rhs[k][j];

		return result;
	}

	/// <summary>
	/// unary "-" operator
	/// </summary>
	template <class R, int RowDim, int ColDim>
	Matrix<R, RowDim, ColDim> operator -(Matrix<R, RowDim, ColDim> m)
	{
		for (int i = 0; i < RowDim; i++)
			for (int j = 0; j < ColDim; j++)
				m[i][j] = -m[i][j];

		return m;
	}
}
