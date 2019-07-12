#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstddef>
#include <exception>

class LUMatrixException : public std::exception
	{
	public:
		LUMatrixException(const int& e, const int i = 0, const int j = 0) {_errno = e; _i = i; _j = j;}
		const char* what()
		{
			switch (_errno)
			{
				case 1: return "bad matrix dimensions";
				case 2: return "index out of bounds";
				case 3: return "row index out of bounds";
				case 4: return "linearly dependent row detected";
				case 5: return "bad input";
				default: return "unknown error";
			}
		}
		int i() const {return _i;}
		int j() const {return _j;}
	private:
		int _errno, _i, _j;
	};

class LUMatrix
{
	friend std::ostream& operator<<(std::ostream& os, const LUMatrix& obj);
	friend std::istream& operator>>(std::istream& is, LUMatrix& obj) throw(LUMatrixException);
	
public:
	
	LUMatrix(const int& n)
	{
		try
		{
			if (n < 1)
			{
				throw LUMatrixException(1, n);
			}
			_dim = n;
			_swaps = 0;
			_matrix = new double*[_dim];
			for (int i = 0; i < _dim; ++i)
			{
				_matrix[i] = new double[_dim];
			}
			_order = new int[_dim];
		}
		catch (LUMatrixException& e)
		{
			std::cerr << "LUMatrix::LUMatrix(const int&) : " << e.what() << ' ' << e.i() << '\n';
			_dim = 0;
			_swaps = 0;
			_matrix = nullptr;
			_order = nullptr;
		}
	}
	
	~LUMatrix()
	{
		if (_matrix)
		{
			for (int i = 0; i < _dim; ++i)
			{
				delete[] _matrix[i];
			}
			delete[] _matrix;
		}
		if (_order)
		{
			delete[] _order;
		}
	}
	
	void decompose(const double& tolerance)
	{
	try
	{
		_tolerance = (tolerance < 0.0 ? 1.0/1024 : tolerance);
		for (int i = 0; i < _dim; ++i)
		{
			_order[i] = i;
		}
		for (int pivot = 0, swap, temp; pivot < _dim; ++pivot)
		{
			// swap
			swap = index_best(pivot);
			if (swap != pivot)
			{
				rswap(pivot, swap);
				temp = _order[pivot];
				_order[pivot] = _order[swap];
				_order[swap] = temp;
				++_swaps;
			}
			// elimination
			for (int row_elim = pivot + 1; row_elim < _dim; ++row_elim)
			{
				double c = _matrix[row_elim][pivot] / _matrix[pivot][pivot];
				for (int col = pivot + 1; col < _dim; ++col)
				{
					_matrix[row_elim][col] -= (c * _matrix[pivot][col]);
				}
				_matrix[row_elim][pivot] = c;
			}
		}
	}
	catch (LUMatrixException& e)
	{
		std::cerr << "void LUMatrix::decompose(const double&) : " << e.what() << '\n';
	}
	}
	
	void display(std::ostream& os)
	{
		const int w = 14;
		os << "Matrix L\n";
		for (int i = 0, j; i < _dim; ++i)
		{
			for (j = 0; j < i; ++j)
			{
				os << std::setw(w) << _matrix[i][j];
			}
			os << std::setw(w) << '1';
			for (++j; j < _dim; ++j)
			{
				os << std::setw(w) << '0';
			}
			os << '\n';
		}
		os << "Matrix U\n";
		for (int i = 0, j; i < _dim; ++i)
		{
			for (j = 0; j < i; ++j)
			{
				os << std::setw(w) << '0';
			}
			for (; j < _dim; ++j)
			{
				os << std::setw(w) << _matrix[i][j];
			}
			os << '\n';
		}
		os << "Swap vector ";
		for (int i = 0; i < _dim; ++i)
		{
			os << (_order[i] + 1) << ' ';
		}
		os << '\n';
		os << "swaps: " << _swaps << '\n';
		return;
	}
	
	double& at(const int& i, const int& j)
	{
		try
		{
			if ((i < 1) || (j < 1) || (i >= _dim) || (j >= _dim))
			{
				throw LUMatrixException(2, i, j);
			}
		}
		catch (LUMatrixException& e)
		{
			std::cerr << "double& LUMatrix::at(const int&, const int&) : " << e.what() << ' ' << e.i() << ',' << e.j() << '\n';
		}
	}
	
	double* operator[](const int& row)
	{
		try
		{
			if ((row < 0) || (row >= _dim)) throw LUMatrixException(3, row);
			return _matrix[row];
		}
		catch (LUMatrixException& e)
		{
			std::cerr << "double* LUMatrix::operator[](const int&) : " << e.what() << ' ' << e.i() << '\n';
			return nullptr;
		}
	}
	
private:
	
	void rswap(const int& i, const int& j)
	{
		double* temp = _matrix[i];
		_matrix[i] = _matrix[j];
		_matrix[j] = temp;
		return;
	}
	
	int index_best(const int& index_pivot) throw(LUMatrixException)
	{
		int index = index_pivot;
		double best = 0.0;
		double row_max;
		for (int i = index_pivot; i < _dim; ++i)
		{
			row_max = std::abs(_matrix[i][index_pivot]);
			for (int j = index_pivot; j < _dim; ++j)
			{
				if (std::abs(_matrix[i][j]) > row_max)
				{
					row_max = std::abs(_matrix[i][j]);
				}
			}
			if (row_max < _tolerance)
			{
				throw LUMatrixException(4);
			}
			if (std::abs(_matrix[i][index_pivot]) / row_max > best)
			{
				best = std::abs(_matrix[i][index_pivot] / row_max);
				index = i;
			}
		}
		return index;
	}
	
	int _dim;
	int _swaps;
	int* _order;
	double _tolerance;
	double** _matrix;
};

std::ostream& operator<<(std::ostream& os, const LUMatrix& obj)
{
	const int w = 14;
	for (int i = 0; i < obj._dim; ++i)
	{
		for (int j = 0; j < obj._dim; ++j)
		{
			os << std::setw(w) << obj._matrix[i][j] << ' ';
		}
		os << '\n';
	}
	os << "swaps: " << obj._swaps << '\n';
	return os;
}

std::istream& operator>>(std::istream& is, LUMatrix& obj) throw(LUMatrixException)
{
	for (int i = 0; i < obj._dim; ++i)
	{
		for (int j = 0; j < obj._dim; ++j)
		{
			std::cout << "(" << i+1 << ',' << j+1 << ") = ";
			is >> obj._matrix[i][j];
			if (std::cin.fail())
			{
				throw LUMatrixException(5);
			}
			std::cout << '\n';
		}
	}
	return is;
}

int main(int argc, char** argv)
{
	try
	{
		int n = 0;
		{
			std::cout << "n = ";
			std::cin >> n;
			std::cin.clear();
			std::cin.ignore();
		} while ((n < 3) || (n > 1000000));
		LUMatrix matrix = LUMatrix(n);
		std::cin >> matrix;
		matrix.decompose(1e-6);
		matrix.display(std::cout);
	}
	catch (LUMatrixException& e)
	{
		std::cerr << e.what() << '\n';
	}
	return 0;;
}
