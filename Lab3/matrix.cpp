#include "matrix.h"
#include <iomanip>

Matrix::Matrix(uint dim): n(dim) 
{
	// орти в порядку X0, X1, X2, X3...
	for (uint i = 0; i < n; i++) {
		matrix.push_back(vector<double>(n));
		ortArr.push_back(i);
		barr.push_back(0);
	}
}

void Matrix::inpMatrix()
{
	clearMatrix();
	double val;
	vector<double> temp;
	for (uint i = 0; i < n; i++) {
		for (uint j = 0; j < n; j++) {
			if (!(cin >> val)) exit(200);
			temp.push_back(val);
		}
		matrix.push_back(temp);
		temp.clear();
	}
}

void Matrix::outMatrix()const
{
	for (uint i = 0; i < n; i++) 
		for (uint j = 0; j < n; j++) {
			cout << setw(7) << matrix[i].at(j);
			if (j == n - 1) cout << " | " << setw(7) << barr[i] << "\n";
			else cout << " ";
		}
}

void Matrix::inpBarr()
{
	double val;
	barr.clear();
	for (uint i = 0; i < n; i++) {
		if (!(cin >> val)) exit(123);
		barr.push_back(val);
	}
}

double Matrix::maxElemRow(uint i)
{
	double max = -DBL_MAX + 1;
	for (uint j = i; j < n; j++) {
		if (matrix[i][j] >= max)
			max = matrix[i][j];
	}
	return max;
}

uint Matrix::columnIndByElem(double elem, uint row, uint startColumn)
{
	for (uint j = startColumn; j < n; j++) {
		if (matrix[row][j] == elem)
			return j;
	}
	exit(600);
}

void Matrix::swapColumns(uint j1, uint j2)
{
	vector<double> temp;
	for (uint i = 0; i < n; i++) 
		temp.push_back(matrix[i][j1]);
	for (uint i = 0; i < n; i++)
		matrix[i][j1] = matrix[i][j2];
	for (uint i = 0; i < n; i++)
		matrix[i][j2] = temp[i];
	temp.clear();
}

void Matrix::swapOrts(uint j1, uint j2)
{
	double temp = ortArr[j1];
	ortArr[j1] = ortArr[j2];
	ortArr[j2] = temp;
}

double Matrix::diagComp()
{
	double res = 1;
	for (uint i = 0; i < n; i++)
		res *= matrix[i][i];
	return res;
}

void Matrix::setDiadOne()
{
	for (uint i = 0; i < n; i++)
		for (uint j = 0; j < n; j++)
			(i == j) ? matrix[i][j] = 1 : matrix[i][j] = 0;
}

void Matrix::setMatVal(uint i, uint j, double val)
{
	if(i >= n || j >= n) 
		exit(800);
	matrix[i][j] = val;
}

double Matrix::getNorm()
{
	double max = -DBL_MAX;
	double sum = 0;
	for(uint j = 0; j < n; j++){
		for(uint i = 0; i < n; i++) {
			sum += abs(matrix[i][j]);
		}
		if (sum >= max)
			max = sum;
		sum = 0;
	}
	return max;
}

void Matrix::clearMatrix()
{
	matrix.clear();
}

istream& operator>>(istream& in, Matrix& g)
{
	uint n = g.getN();
	double val;
	cout << "Enter " << n * n << " doubles, " << n << " per line, divided by spaces\n";
	g.inpMatrix();
	cout << "Enter " << n << " doubles inline(b column)\n";
	g.inpBarr();
	return in;
}

ostream& operator<<(ostream& out, const Matrix& g)
{
	cout << "RESULT:\n";
	g.outMatrix();
	return out;
}
