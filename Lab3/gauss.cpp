#include "gauss.h"
using namespace std;

Gauss::Gauss(uint dim) : Matrix(dim) { xsol.resize(dim); }

Gauss::Gauss(const Gauss& c) : Matrix(c.n)
{
	for (uint i = 0; i < n; i++) {
		matrix.push_back(vector<double>(n));
	}

	for (uint i = 0; i < n; i++) {
		for (uint j = 0; j < n; j++) {
			matrix[i][j] = c.matrix[i][j];
		}
		barr.push_back(c.barr[i]);
		ortArr[i] = c.ortArr[i];
	}
}

void Gauss::direct()
{
	Matrix diagMat{ n };
	diagMat.setDiadOne();

	vector<double> lmj(n);
	for (uint i = 0; i < n - 1; i++) {
		double max = maxElemRow(i); // get max elem in row
		uint maxCol = columnIndByElem(max, i, i); // get column index of max elem
		swapColumns(i, maxCol);
		diagMat.swapColumns(i, maxCol);
		swapOrts(i, maxCol);
		diagMat.swapOrts(i, maxCol);
		for (uint m = i + 1; m < n; m++) {
			lmj[m] = matrix[m][i] / matrix[i][i];
		}
		for (uint j = i; j < n; j++) {
			for (uint m2 = i + 1; m2 < n; m2++) {
				matrix[m2][j] = matrix[m2][j] - matrix[i][j] * lmj[m2];
				corVal(matrix[m2][j]);
			}
		}
		for (uint m2 = i + 1; m2 < n; m2++) {
			barr[m2] = barr[m2] - barr[i] * lmj[m2];
			corVal(barr[m2]);
		}
	}
	lmj.clear();
}

void Gauss::reverse()
{
	vector<double> xval(n);
	double sum = 0;
	xval[n - 1] = barr[n - 1] / matrix[n - 1][n - 1];
	for (int i = n - 2; i >= 0; i--) {
		for(int j = n - 1; j > i; j--){
			sum += matrix[i][j] * xval[j];
		}
		xval[i] = (barr[i] - sum) / matrix[i][i];
		sum = 0;
	}

	for (uint i = 0; i < n; i++)
		xsol[ortArr[i]] = xval[i];
	xval.clear();
}

void Gauss::dispSol()
{
	for (uint i = 0; i < n; i++)
		cout << "y" << i + 1 << " = " << setw(7) << xsol[i] << "\n";
	cout << "y" << 10 << " = " << 1 << "\n";
}

double Gauss::getDet()const
{
	Gauss g{*(this)};
	g.direct();
	return g.diagComp();
}

Gauss Gauss::getInv(double det)
{
	Gauss mat{ n };
	for (uint i = 0; i < 3; i++)
	{
		for (uint j = 0; j < 3; j++)
		{
			if (i == j) {
				mat.setMatVal(i, j, ((matrix[(i + 1) % 3][(j + 1) % 3] *
					matrix[(i + 2) % 3][(j + 2) % 3]) - (matrix[(i + 1) % 3][(j + 2) % 3] *
						matrix[(i + 2) % 3][(j + 1) % 3])) / det);
			}
			else
			{
				mat.setMatVal(j, i, ((matrix[(i + 1) % 3][(j + 1) % 3] *
					matrix[(i + 2) % 3][(j + 2) % 3]) - (matrix[(i + 1) % 3][(j + 2) % 3] *
						matrix[(i + 2) % 3][(j + 1) % 3])) / det);
			}
			/*cout << ((matrix[(i + 1) % 3][(j + 1) % 3] *
				matrix[(i + 2) % 3][(j + 2) % 3]) - (matrix[(i + 1) % 3][(j + 2) % 3] *
					matrix[(i + 2) % 3][(j + 1) % 3])) / det << "\t";*/
		}
		cout << "\n";
	}
	return mat;
}

void Gauss::corVal(double& val)
{
	if (abs(val) < 1e-10)
		val = 0;
}
