#pragma once
#include <vector>
#include <array>

template <class T, int Size>
class JordanBlock 
{
private:
	std::array<T, Size> _data;

	static_assert(Size >= 1, "Invalid size");

	//REturns inverse block (if defined)
	JordanBlock<T, Size> GetInverse() const
	{
		JordanBlock<T, Size> inverse;

		inverse._data[0] = T(1)/_data[0];

		for (int i = 1; i< Size; i++)
		{
			T b = T(0);
			for (int j = 1; j <= i; j++)
				b += _data[j]*inverse._data[i - j];

			inverse._data[i] = -b*inverse._data[0];
		}

		return inverse;
	}

public:

	//Default constructor
	JordanBlock()
	{
		_data.fill(T(0));
	}

	//Constructor
	JordanBlock(const JordanBlock<T, Size>& source)
	{
		_data = source._data;
	}

	//Constructor
	JordanBlock(const T& value)
	{
		_data.fill(T(0));
		_data[0] = value;
	}

	JordanBlock(const std::array<T, Size>& source)
	{
		_data = source;
	}

	//Subscript operator
	inline T operator[](const int index) const
	{
		return _data[index];
	}

	JordanBlock<T, Size>& operator =(const JordanBlock<T, Size>& anotherBlock)
	{
		_data = anotherBlock._data;
		return this;
	}

	void operator +=(const JordanBlock<T, Size>& anotherBlock)
	{
		for (int i = 0; i< Size; i++)
			_data[i] += anotherBlock._data[i];
	}

	JordanBlock<T, Size> operator +(const JordanBlock<T, Size>& anotherBlock) const
	{
		JordanBlock<T,Size> result(*this);
		result+=anotherBlock;
		return result;
	}

	void operator +=(const T& scalar)
	{
		_data[0] += scalar;
	}

	JordanBlock<T, Size> operator +(const T& scalar) const
	{
		JordanBlock<T,Size> result(*this);
		result += scalar;
		return result;
	}

	void operator -=(const JordanBlock<T, Size>& anotherBlock)
	{
		for (int i = 0; i< Size; i++)
			_data[i] -= anotherBlock._data[i];
	}

	JordanBlock<T, Size> operator -(const JordanBlock<T, Size>& anotherBlock) const
	{
		JordanBlock<T, Size> result(*this);
		result -= anotherBlock;
		return result;
	}

	void operator -=(const T& scalar)
	{
		*this += -scalar;
	}

	JordanBlock<T, Size> operator -(const T& scalar) const
	{
		JordanBlock<T, Size> result(*this);
		result -= scalar;
		return result;
	}

	JordanBlock<T, Size> operator -() const
	{
		JordanBlock<T, Size> result;
		result -= *this;
		return result;
	}

	JordanBlock<T, Size> operator *(const JordanBlock<T, Size>& anotherBlock) const
	{
		JordanBlock<T, Size> result;

		for (int i = 0; i< Size; i++)
			for (int j = 0; j <= i; j++)
				result._data[i] += _data[j]*anotherBlock._data[i - j];

		return result;
	}

	void operator *=(const JordanBlock<T, Size>& anotherBlock)
	{
		*this = *this*anotherBlock;
	}

	void operator *=(const T& factor)
	{
		for (int i = 0; i<Size; i++)
			_data[i] *= factor;
	}

	JordanBlock<T, Size> operator *(const T& factor)
	{
		JordanBlock<T, Size> result(*this);
		result *= factor;
		return result;
	}

	JordanBlock<T, Size> operator /(const JordanBlock<T, Size>& anotherBlock) const
	{
		return *this*anotherBlock.GetInverse();
	}

	void operator /=(const JordanBlock<T, Size>& anotherBlock)
	{
		*this = *this*anotherBlock.GetInverse();
	}

	void operator /=(const T& denominator)
	{
		T factor = T(1)/denominator;
		*this *= factor;
	}

	JordanBlock<T, Size> operator /(const T& denominator) const
	{
		JordanBlock<T, Size> result(*this);
		result /= denominator;
		return result;
	}

    template <class T1, class T2, int Size1>
	friend JordanBlock<T1, Size1> operator /(const T2& factor, const JordanBlock<T1, Size1>& anotherBlock);
};

    template <class T, int Size>
	JordanBlock<T, Size> operator +(const T& scalar, const JordanBlock<T, Size>& anotherBlock)
	{
		return anotherBlock + scalar;
	}

    template <class T, int Size>
	JordanBlock<T, Size> operator +(const int& scalar, const JordanBlock<T, Size>& anotherBlock)
	{
		return anotherBlock + scalar;
	}

    template <class T, int Size>
	JordanBlock<T, Size> operator -(const T& scalar, const JordanBlock<T, Size>& anotherBlock)
	{
		return - anotherBlock + scalar;
	}

    template <class T, int Size>
	JordanBlock<T, Size> operator -(const int& scalar, const JordanBlock<T, Size>& anotherBlock)
	{
		return -anotherBlock + scalar;
	}

    template <class T, int Size>
	JordanBlock<T, Size> operator *(const T& factor, const JordanBlock<T, Size>& anotherBlock)
	{
		return anotherBlock*factor;
	}

    template <class T, int Size>
	JordanBlock<T, Size> operator *(const int& factor, const JordanBlock<T, Size>& anotherBlock)
	{
		return anotherBlock*factor;
	}

    template <class T, class T1, int Size>
	JordanBlock<T, Size> operator /(const T1& factor, const JordanBlock<T, Size>& anotherBlock)
	{
		return anotherBlock.GetInverse()*factor;
	}
