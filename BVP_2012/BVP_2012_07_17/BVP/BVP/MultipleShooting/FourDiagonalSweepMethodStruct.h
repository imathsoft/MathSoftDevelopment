#ifndef GUARD_FOUR_DIAGONAL_SWEEP_METHOD_STRUCT
#define GUARD_FOUR_DIAGONAL_SWEEP_METHOD_STRUCTd


template <class T>
class FourDiagonalSweepMethodStruct
{
private:
	static const int _colCount = 5;
	std::vector<std::array<T, _colCount>> _body;
	int _rowCount;
public:
	///Constructor
	FourDiagonalSweepMethodStruct(int rowCount)
	{
		_rowCount = max(rowCount, 0);
		_body.resize(_rowCount);
	}

	///REturns count of columns
	int RowCount()
	{
		return _rowCount;
	}

	///Destructor
	~FourDiagonalSweepMethodStruct()
	{
	}

	///Subscript operatorl
	std::array<T, _colCount>& operator[](int rowIndex)
	{			
		return _body.at(rowIndex);
	}

	///Method to save the struct to file
	void SaveToFile(char* filename, std::streamsize pecision = 15)
	{
		 ofstream file;
		 file.precision(pecision);
		 file.open (filename);
		 for (int i = 0; i < _rowCount; i++)
		 {
			 for (int j = 0; j < _colCount; j++)
			     file << _body[i][j] << " "; 
			 file << endl;
		 }
         file.close();
	}
};

#endif