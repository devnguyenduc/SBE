#ifndef BLOCH_H
#define BLOCH_H

#include <complex>
using namespace std;

//parameters
#define PI 3.14159265358979
#define N 100
#define var_epsilon 2*300/N
#define hbar 658.5
#define Er 4.2
#define dephasing_time 100.0
#define Chi 0.1
#define phase sqrt(Er)*var_epsilon/PI

//declaration of g(n,l) function
double g(int n, int l);
//declaration of energy function E
complex<double> E(int n, double t, complex<double> nsolutions[]);
//declaration of Rabi frequency function
complex<double> Omega(int n, double t, complex<double> nsolutions[]);
//declaration of pn
complex<double> pn(int n, double t, complex<double> nsolutions[]);
//declaration of fn
double fn(int n, double t, complex<double> nsolutions[]);
//declaration of monitoring
complex<double> solve(int index, double t, complex<double> nsolutions[]);

#endif
