#include<iostream>
#include<tgmath.h> 
#include <cmath> 

using namespace std;

double F = 100.0;
double c = 8.0;
int n = 3;

double B_market = 100.0;
double tol = 0.0001;
int max_iter = 100;

void price_from_yield(double F, double c, double y, int n, double & B) {
	double y_decimal = 0.01*y; // correct y_decimal
	B = 0;
	
	double temp = .5 * c;
	double temp2 = 1 + (.5 * y_decimal);
	
	for (int i = 1; i <= n - 1; i++) {
		B = B + temp / pow(temp2,i);
	} // Calculates sums from 1 until n - 1
	
	double temp3 = ((F + temp) / (pow(temp2,n))); // Calculates the final sum for the power n
	
	B = B + temp3;
		
}

int yield_from_price_Bisection(double F, double c, int n, double B_market, double tol, int max_iter, double & y, int & num_iter)	{
	
	num_iter = 0;
	double By_low = 0.0;
	double By_high = 0.0;
	
	double y_low = 0.0;
	double y_high = 100.0; // End of step 2
	
	price_from_yield(F,c,y_low,n,By_low); // Step 3
	double LowTemp = By_low - B_market;
	
	if ((abs(LowTemp)) <= tol) {
		return 0; // success
	} // Step 4
	
	if (By_low < B_market) {
		y = 0;
		return 1; // fail
	} // Step 5
	
	price_from_yield(F,c,y_high,n,By_high); // Step 6
	double HighTemp = By_high - B_market;
	
	if ((abs(HighTemp)) <= tol) {
		y = y_high; // y = yhigh
		return 0; // success
	} // Step 7
	
	if (By_high > B_market) {
		y = 0;
		return 1; // fail
	} // Step 8
	
	cout << "Testing for Step 9: " << By_low << " > " << B_market << " > " << By_high << endl; // Step 9
	
	double B;
	double temp = y_low + y_high;
	
	for (int i = 0; i < max_iter; i++) {
		
		double temp = y_low + y_high;
		y = temp / 2.0; 
		price_from_yield(F,c,y,n,B); //  End of Step 12
		
		if (abs(B - B_market) <= tol) {
			return 0;
		} // Step 13
		
		else if (B > B_market) {
			y_low = y;
		} // Step 14
		
		else if (B < B_market) {
			y_high = y;
		} // Step 15
		
		if (y_high - y_low <= tol) {
			return 0;
		} // Step 16
		
	} // Step 10 -> 17 
	
	y = 0; // Step 18
	
	return 1; // Fail is iteration does not have a good enough answer for y, End of Step 18
} 

int main() {
	
	double YieldOutput = 0.0; // Store the yield output from the function
	int numberIteration;
	int value;
	
	value = yield_from_price_Bisection(F, c, n, B_market, tol, max_iter, YieldOutput, numberIteration);
	
	cout << "This is the yield: " << YieldOutput << endl;
	cout << "NumberIteration: " << numberIteration << endl;
	cout << "The value of the function is: " << value << endl;
	
	return 0;
}
