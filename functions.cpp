#include "bloch.h"
#include <iostream>
#include <math.h>
#include <complex>
using namespace std;
const complex<double> i_(0.0,1.0);
//define function g
double g(int n, int l)
{
	if (n==l)
		return 0;
	else
		return (1/(sqrt(n*var_epsilon))*(log(sqrt(n)+sqrt(l))-log(abs(sqrt(n)-sqrt(l)))));
}
//define function E
complex<double> E(int n, double t, complex<double> nsolutions[])
{
	complex<double> sum{0};
	cout << N / 2 << endl;
	for (int l=0; l<(N/2); l++)
	{
		sum += g(n/2,l+1)*2*nsolutions[2*l];
	}
	return phase*sum;

}
//define function Omega
complex<double> Omega(int n, double t, complex<double> nsolutions[])
{
	complex<double> sum{0};
	for(int l=0;l<(N/2);l++)
	{
		sum+= g(n/2,l+1)*nsolutions[2*l+1];
	}
	cout << "Sum" << sum << endl;
	double time=sqrt(PI)*Chi*exp(-t*t/(25*25))/(2*25);
	return time+phase*sum/hbar;
}
//define function pn
complex<double> pn(int n, double t, complex<double> nsolutions[])
{
	return (-i_/hbar*((n+1)/2.0*var_epsilon-30.0-E(n,t,nsolutions))*nsolutions[n]+i_*(1.0-2.0*nsolutions[n-1])*Omega(n,t,nsolutions)-nsolutions[n]/dephasing_time);
}
//define function fn
double fn(int n, double t, complex<double> nsolutions[])
{
	return (-2.0*imag(Omega(n,t,nsolutions)*conj(nsolutions[n+1])));
}
//define function solve
complex<double> solve(int index, double t, complex<double> nsolutions[])
{
	if(index % 2 !=0)
		return pn(index,t,nsolutions);
	else
		return fn(index,t,nsolutions);
}
