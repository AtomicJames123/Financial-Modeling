#include<iostream>
#include<tgmath.h> 

using namespace std;

double e = 2.71828182845904523536;

int df_and_r(double F0, double F1, double t0, double t1, double & df, double & r)	{
	
	if (t1-t0 == 0.0) {
		df = 0;
		r = 0;
		return -1;
	}
	
	if ((F0 <= 0.0) || (F1 <= 0.0)) {
		df = 0;
		r = 0;
		return -2;
	}
	
	df = F0 / F1; // DiscountFactor equation
	
	double temp1 = log(df);
	double temp2 = (t1 - t0);
	double temp3 = -100.0;
	double temp4 = (temp1 / temp2);
	
	r = temp4 * temp3;
	
	return 0;
} // Function for Question 1.2

int main() {
	double CashFlow;
	double FutureCashFlow;
	double CurrentTime;
	double FutureTime;
	double DiscountFactor;
	double Rate;
	
	CashFlow = 100;
	FutureCashFlow = 325;
	CurrentTime = 10;
	FutureTime = 20;
	
	int value; // For testing purposes
	
	value = df_and_r(CashFlow, FutureCashFlow, CurrentTime, FutureTime, DiscountFactor, Rate);
	
	cout << "The discount factor is: " << DiscountFactor << endl; // DiscountFactor value
	cout << "The rate is: " << Rate << endl;
	
	return 0;
}
