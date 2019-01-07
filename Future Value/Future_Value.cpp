#include<iostream>
#include<tgmath.h> 

using namespace std;


double future_value(double F0, double t0, double t1, double r)	{
	double r_decimal = 0.01*r;
	double F1 = F0*exp(r_decimal*(t1-t0));
	return F1;
}

int main() {
	int CashFlow;
	int CurrentTime;
	int FutureTime;
	int Rate;
	
	CashFlow = 1000;
	CurrentTime = 10;
	FutureTime = 15;
	Rate = 5;
	
	int value;
	
	while (0 <= CashFlow) {
		value = future_value(CashFlow,CurrentTime,FutureTime,Rate);
		cout << value << endl;
		CashFlow = CashFlow - 100;
	}
	
	return 0;
}
