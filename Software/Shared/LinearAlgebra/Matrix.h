#pragma once
#include <array>
#include <stdlib.h>
#include <time.h>
#include <limits>
#include "../Utils/Utils.h"

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
		std::array<int, RowDim> row_map{};

		/// <summary>
		/// Initializes row map
		/// </summary>
		void initialize_row_map()
		{
			for (int row_id = 0; row_id < RowDim; row_id++)
				row_map[row_id] = row_id;
		}

	public:

		/// <summary>
		/// Constructor from 2d array
		/// </summary>
		/// <param name="d"></param>
		Matrix(const std::array<std::array<R, ColDim>, RowDim>& d) : data{d}
		{
			initialize_row_map();
		}

		/// <summary>
		/// Default constructor
		/// </summary>
		Matrix()
		{
			initialize_row_map();
		}

		/// <summary>
		/// Sub-script operator
		/// </summary>
		std::array<R, ColDim>& operator [](const int i)
		{
			return data[row_map[i]];
		}

		/// <summary>
		/// Sub-script operator (const)
		/// </summary>
		const std::array<R, ColDim>& operator [](const int i) const
		{
			return data[row_map[i]];
		}

		/// <summary>
		/// Swap rows with the given indices
		/// It is responsibility of the caller to ensure that the indices are within the acceptable range
		/// </summary>
		void swap_rows(const int row_id_0, const int row_id_1)
		{
			const auto temp = row_map[row_id_0];
			row_map[row_id_0] = row_map[row_id_1];
			row_map[row_id_1] = temp;
		}

		/// <summary>
		/// Multiplies the row with index `row_to_add_id` by the given factor and adds it to the row with index `row_to_add_to_id`
		/// It is responsibility of the caller to ensure that the indices are within the acceptable range
		/// </summary>
		void add_rows(const int row_to_add_id, const int row_to_add_to_id, const R& factor)
		{
			const auto& row_to_add = (*this)[row_to_add_id];
			auto& row_to_add_to = (*this)[row_to_add_to_id];
			for (int col_id = 0; col_id < ColDim; col_id++)
				row_to_add_to[col_id] += row_to_add[col_id] * factor;
		}

		/// <summary>
		/// Multiplies given row by the given factor
		/// </summary>
		void multiply_row(const int row_id, const R& factor)
		{
			auto& row = (*this)[row_id];

			for (int col_id = 0; col_id < ColDim; col_id++)
				row[col_id] *= factor;
		}

		/// <summary>
		/// Tries to calculate inverse matrix for the current one
		/// throws an exception if the current matrix is not invertable
		/// Can be calle on square matrices only
		/// </summary>
		Matrix<R, RowDim, ColDim> Inverse() const;

		/// <summary>
		/// += operator
		/// </summary>
		Matrix<R, RowDim, ColDim>& operator +=(const Matrix<R, RowDim, ColDim>& rhs)
		{
			for (int i = 0; i < RowDim; i++)
				for (int j = 0; j < ColDim; j++)
					(*this)[i][j] += rhs[i][j];

			return *this;
		}

		/// <summary>
		/// -= operator
		/// </summary>
		Matrix<R, RowDim, ColDim>& operator -=(const Matrix<R, RowDim, ColDim>& rhs)
		{
			for (int i = 0; i < RowDim; i++)
				for (int j = 0; j < ColDim; j++)
					(*this)[i][j] -= rhs[i][j];

			return *this;
		}

		/// <summary>
		/// *= operator (matrix-scalar in place multiplication)
		/// </summary>
		Matrix<R, RowDim, ColDim>& operator *=(const R& rhs)
		{
			for (int i = 0; i < RowDim; i++)
				for (int j = 0; j < ColDim; j++)
					(*this)[i][j] *= rhs;

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
		/// Returns transposed matrix
		/// </summary>
		Matrix<R, ColDim, RowDim> Transpose() const
		{
			Matrix<R, ColDim, RowDim> result;

			for (int i = 0; i < RowDim; i++)
				for (int j = 0; j < ColDim; j++)
					result[j][i] = (*this)[i][j];

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
					norm += (*this)[i][j] * (*this)[i][j];

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

		/// <summary>
		/// "Equal" operator
		/// </summary>
		bool operator == (const Matrix<R, RowDim, ColDim>& matr) const
		{
			for (int i = 0; i < RowDim; i++)
				for (int j = 0; j < ColDim; j++)
					if ((*this)[i][j] != matr[i][j])
						return false;

			return true;
		}

		/// <summary>
		/// "Not equal" opertor
		/// </summary>
		bool operator != (const Matrix<R, RowDim, ColDim>& matr) const
		{
			return !((*this) == matr);
		}

		/// <summary>
		/// Returns `true` if at least one element of the matrix is "not a number"
		/// </summary>
		/// <returns></returns>
		bool is_nan() const
		{
			return (*this) != (*this);
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

	/// <summary>
	/// Returns matrix inverse for the given one
	/// Or throws exception if the argyument matix is not invertible 
	/// </summary>
	template <class R, int Dim>
	Matrix<R, Dim, Dim> inverse(const Matrix<R, Dim, Dim>& m)
	{
		static_assert(Dim != 1 && Dim != 2, "More efficient implementation must be considered");

		static const R EPSILON_WEAK = 100 * std::numeric_limits<R>::epsilon();

		Matrix<R, Dim, Dim> result{};

		auto temp = m;
		result = Matrix<R, Dim, Dim>::Identity();

		for (int col_id = 0; col_id < Dim - 1; col_id++)
		{
			int best_row_id = col_id;
			auto best_element_abs = auxutils::Abs(temp[col_id][col_id]);
			//in the given column we search for the biggest element below the diagonal
			for (int row_id = col_id + 1; row_id < Dim; row_id++)
			{
				const auto candidate_abs = auxutils::Abs(temp[row_id][col_id]);
				if (candidate_abs > best_element_abs)
				{
					best_row_id = row_id;
					best_element_abs = candidate_abs;
				}
			}

			if (best_row_id != col_id)
			{
				temp.swap_rows(best_row_id, col_id);
				result.swap_rows(best_row_id, col_id);
			}

			if (best_element_abs < EPSILON_WEAK)
				continue;

			const auto diagonal_element_inverted = R(1) / temp[col_id][col_id];
			for (int row_id = col_id + 1; row_id < Dim; row_id++)
			{
				const auto factor = -temp[row_id][col_id] * diagonal_element_inverted;
				temp.add_rows(col_id, row_id, factor);
				result.add_rows(col_id, row_id, factor);
			}
		}

		for (int i = 0; i < Dim; i++)
			if (auxutils::Abs(temp[i][i]) < EPSILON_WEAK)
				throw std::exception("Singular matrix");

		for (int col_id = 0; col_id < Dim; col_id++)
		{
			const auto factor = R(1) / temp[col_id][col_id];
			temp.multiply_row(col_id, factor);
			result.multiply_row(col_id, factor);

			for (int row_id = col_id - 1; row_id >= 0; row_id--)
			{
				const auto factor = -temp[row_id][col_id];
				temp.add_rows(col_id, row_id, factor);
				result.add_rows(col_id, row_id, factor);
			}
		}

		return result;
	}

	/// <summary>
	/// Returns matrix inverse for the given one
	/// Or throws exception if the argyument matix is not invertible 
	/// Partial specialization for the case of matrices 2x2
	/// </summary>
	template <class R>
	Matrix<R, 2, 2> inverse(const Matrix<R, 2, 2>& m)
	{
		static const R EPSILON_WEAK = 100 * std::numeric_limits<R>::epsilon();

		Matrix<R, 2, 2> result{};

		const auto det = m[0][0] * m[1][1] - m[1][0] * m[0][1];

		if (auxutils::Abs(det) < EPSILON_WEAK)
			throw std::exception("Singular matrix");

		const auto one_over_determinant = R(1) / det;

		result[0][0] = one_over_determinant * m[1][1];
		result[0][1] = -one_over_determinant * m[0][1];
		result[1][0] = -one_over_determinant * m[1][0];
		result[1][1] = one_over_determinant * m[0][0];

		return result;
	}

	/// <summary>
	/// Returns matrix inverse for the given one
	/// Or throws exception if the argyument matix is not invertible 
	/// Partial specialization for the case of matrices 1x1
	/// </summary>
	template <class R>
	Matrix<R, 1, 1> inverse(const Matrix<R, 1, 1>& m)
	{
		static const R EPSILON_WEAK = 100 * std::numeric_limits<R>::epsilon();

		Matrix<R, 1, 1> result{};

		if (auxutils::Abs(m[0][0]) < EPSILON_WEAK)
			throw std::exception("Singular matrix");

		result[0][0] = R(1)/ m[0][0];

		return result;
	}

	/// <summary>
	/// Tries to calculate inverse matrix for the current one
	/// throws an exception if the current matrix is not invertable
	/// Can be calle on square matrices only
	/// </summary>
	template <class R, int RowDim, int ColDim>
	Matrix<R, RowDim, ColDim> Matrix<R, RowDim, ColDim>::Inverse() const
	{
		return inverse(*this);
	}

}
