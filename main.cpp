//#include "bloch.h"
#include <iostream>
#include "class.h"
#include <complex>
complex<double> *zeros(int num){
	complex<double> *x = new complex<double>[num];
	for(int i = 0;i < num;i++){
		x[i] = 0;
	}
	return x;
}
int main()
{
	//number of interation 
	int loop = 288;
	//create an object subDev of SBE class
	SBE *subDev=new SBE();
	//generate object rk4 of RungKutta class
	RungeKutta4 *rk4=new RungeKutta4(-75,501,loop,subDev);
	//create an array of 100 zeros
	complex<double> *test = zeros(100);
	//std::cout << test[5] << '\n';
	//call method add_solution
	rk4->add_solution(test);
	//call method run of class RungeKutta
	rk4->run();
	delete rk4;
	return 0;
}
