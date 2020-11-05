#pragma once
#include <iomanip>
#include "matrix.h"


class Gauss :
	public Matrix
{
public:
	Gauss(uint dim);
	Gauss(const Gauss&);
	void direct();
	void reverse();
	void dispSol();
	double getDet()const;
	Gauss getInv(double d);
private:
	vector<double> xsol;
	void corVal(double& val);
};

