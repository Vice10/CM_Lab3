#pragma once
#include <istream>
#include <ostream>
#include <vector>
#include <iostream>
#include <algorithm>
using namespace std;
typedef unsigned int uint;

class Matrix
{
public:
	Matrix(uint dim);
	Matrix(const Matrix& g) = delete;
	Matrix(Matrix&& g) = delete;
	Matrix& operator=(const Matrix& g) = delete;
	uint getN() { return n; };
	void inpMatrix();
	void outMatrix()const;
	void inpBarr();
	double maxElemRow(uint i);
	uint columnIndByElem(double elem, uint row, uint startColumn);
	void swapColumns(uint i, uint j);
	void swapOrts(uint i, uint j);
	double diagComp();
	void setDiadOne();
	void setMatVal(uint i, uint j, double val);
	double getNorm();
	void setXelem(int i, int j, double val) {
		if (i < n && j < n)  matrix[i][j] = val; else exit(300);
	}
	void setMatElems(vector<vector<double>> mat) { matrix = mat; }
	void setBvec(vector<double> bvec) { barr = bvec; }

protected:
	uint n;
	vector<vector<double>> matrix;
	vector<double> barr;
	void clearMatrix();
	vector<uint> ortArr;
};

istream& operator >> (istream& in, Matrix& g);
ostream& operator << (ostream& out, const Matrix& g);