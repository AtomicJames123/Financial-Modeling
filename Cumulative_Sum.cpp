#include<iostream>
#include<tgmath.h> 

using namespace std;

double cum_norm(double x)
{
const double root = sqrt(0.5);
return 0.5*(1.0 + erf(x*root));
}

int main() {
	double input = -.224;
	double output = 0;
	output = cum_norm(input);
	
	cout << output << endl;
	
	return 0;
}
