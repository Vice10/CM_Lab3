#include <iostream>
#include <vector>
#include "matrix.h"
#include "gauss.h"
using namespace std;

const double h = 0.1;

double p(double xi) { return 1.0; }
double q(double xi) { return -1 / xi; }
double a(double xi) { return 1 / (h * h) - 1 / (2 * h); }
double b(double xi) { return -2 / (h * h) - 1 / xi; }
double c(double xi) {return 1 / (h * h) + 1 / (2 * h); }
double f(double xi) { return 8 * xi * xi - 8 * xi + 1.5; }

int main()
{
    vector<vector<double>> mat;
    vector<double> bvec;
    vector<double> tmp;
    cout << b(h) << "y1 + " << c(h) << "y2 = d1* = " << f(h) - a(h) * 0 << "\n";
    tmp.push_back(b(h)); tmp.push_back(c(h)); tmp.resize(9);
    mat.push_back(tmp); tmp.clear();
    bvec.push_back(f(h) - a(h) * 0);

    for (int i = 2; i <= 8; i++) {
        cout << a(i * h) << "y" << i - 1 << " + " << b(i * h) << "y" << i << " + " << c(i * h) << "y" << i + 1
            << " = d" << i << " = " << f(i*h) <<"\n";
        tmp.resize(9);
        tmp[(i - 1) - 1] = a(i * h); tmp[(i) - 1] = b(i * h); tmp[(i + 1) - 1] = c(i * h);
        mat.push_back(tmp);
        tmp.clear();
        bvec.push_back(f(i * h));
    }
    cout << a(h*9) << "y8 + " << b(h*9) <<"y9" << " = d9* = " << f(h*9) - c(h*9)<<"\n"; // y10 = 1
    tmp.resize(9);
    tmp[8 - 1] = a(h * 9); tmp[9 - 1] = b(h * 9);
    mat.push_back(tmp);
    bvec.push_back(f(h * 9) - c(h * 9));


    int n = 9;
    Gauss matr{ 9 };
    matr.setMatElems(mat);
    matr.setBvec(bvec);
    cout << "\n\n" << matr;
    matr.direct();
    cout << "\n\n" << matr << "\n\n";
    matr.reverse();
    matr.dispSol();
}