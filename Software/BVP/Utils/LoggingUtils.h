#pragma once
#include <fstream>
#include "AuxUtils.h"

/// <summary>
/// Fucntionality to log data in a text file as a table (space, comma, or whatever separated values, that can be loaded in Excel)
/// </summary>
template <int ColNum>
class table_logger
{
private:
	std::ofstream _file;
	std::string _separator = " ";
public:

	/// <summary>
	/// Constructor
	/// </summary>
	/// <param name="file_name">Name of the file to save the log</param>
	/// <param name="head_line">Head line of the table to be logged</param>
	table_logger(const std::string file_name, const std::array<std::string, ColNum>& head_line, std::string separator = " ", const bool append = false)
	{
		_file.open(file_name, append ? std::ios::app : std::ios::out);
		_separator = separator;

		for (int caption_id = 0; caption_id < ColNum; caption_id++)
			if (caption_id != ColNum - 1) 
				_file << head_line[caption_id] << separator;
			else
				_file << head_line[caption_id] << std::endl;
	}

	/// <summary>
	/// Writes a line with the given data to the underlying file
	/// </summary>
	template<class T>
	void write_line(const std::array<T, ColNum>& data_line)
	{
		if (!_file.is_open())
			throw std::exception("Attempt to write to a closed file");

		for (int item_id = 0; item_id < ColNum; item_id++)
			if (item_id != ColNum - 1)
				_file << auxutils::ToString(data_line[item_id]) << _separator;
			else
				_file << auxutils::ToString(data_line[item_id]) << std::endl;
	}

	/// <summary>
	/// Closes the underlying file
	/// </summary>
	void close()
	{
		_file.close();
	}

	/// <summary>
	/// Destructor
	/// </summary>
	~table_logger()
	{
		_file.close();
	}
};
