#include <iostream>
#include <complex>
#include "bloch.h"

using namespace std;
complex<double> *addVector(int dim, complex<double> *vec1, complex<double> *vec2)
{
	complex<double> *vec=new complex<double>;
	for(int i=0;i<dim;i++)
	{
		vec[i]=vec1[i]+vec2[i];	
	}
	return vec;
}
complex<double> *matMul(int dim, double factor, complex<double> *vec) 
{
	complex<double> *temp=new complex<double>;
	for(int i=0;i<dim;i++)
	{
		temp[i]=factor*vec[i];	
	}
	return temp;
}
class SBE
{
public:
	SBE(){}
	complex<double> solve(int j, double t, complex<double> *nsolutions){
		return solve(j,t,nsolutions);
	}
	~SBE(){}
};

class RungeKutta4
{
public:
	SBE *dev = new SBE(); 
	double start,end;
	int loop;
	double h;
	complex<double> *solution;
	complex<double> **result;
	RungeKutta4(double start, double end, int loop, SBE *dev)
	{
		this->dev = dev;
		this->start=start;
		this->end=end;
		this->loop=loop;
		this->h=(end-start)/loop;
		this->result = new complex<double>*[this->loop];
		for(int i = 0; i < this->loop; i++){
			this->result[i] = new complex<double>[N];
		}
	}
	~RungeKutta4(){ 
		delete this->dev; 
		delete this->solution;
		delete this->result;       
	}
	void add_solution(complex<double> *solution)
	{
		this->solution = solution;
	}
	void run()
	{
		complex<double> *k1 = new complex<double>[N];
		complex<double> *k2 = new complex<double>[N];
		complex<double> *k3 = new complex<double>[N];
		complex<double> *k4 = new complex<double>[N];
		
		for(int j=0;j<N;j++)
		{
			this->result[0][j]= this->solution[j];
		}
		for(int i=1;i<this->loop;i++)
		{
			for(int j=0;j<N;j++)
			{
				k1[j]=this->h*this->dev->solve(j,this->start+this->h*i,this->result[i-1]);
			}
			for(int j=0;j<N;j++)
			{
				k2[j]=this->h*this->dev->solve(j,this->start+this->h*(i+1/2),addVector(N,this->result[i-1],matMul(N,0.5,k1)));
			};
			
			for(int j=0;j<N;j++)
			{
				k3[j]=this->h*this->dev->solve(j,this->start+this->h*(i+1/2),addVector(N,this->result[i-1],matMul(N,0.5,k2)));
			};
			
			for(int j=0;j<N;j++)
			{
				k4[j]=this->h*this->dev->solve(j,this->start+this->h*(i+1),addVector(N,this->result[i-1],k3));
			};
			for(int j=0;j<N;j++)
			{
				this->result[i][j]=result[i-1][j]+(k1[j]+2.0*k2[j]+2.0*k3[j]+k4[j])/6.0;
			}
		}
	delete[] k1;
	delete[] k2;
	delete[] k3;
	delete[] k4;

	}
};
