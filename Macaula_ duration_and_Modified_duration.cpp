#include<iostream>
#include<tgmath.h> 

using namespace std;

void price_from_yield1(double F, double c, double y, int n, double & B) {
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

void price_from_yield2(double F, double c, double y, int n, double & B, double & D_mac, double & D_mod) {
	double y_decimal = 0.01*y; // correct y_decimal
	B = 0.0;
	
	double temp1 = .5 * c;
	double temp2 = 1 + (.5 * y_decimal);
	
	for (int i = 1; i <= n - 1; i++) {
		B = B + temp1 / pow(temp2,i);
	} // Calculates sums from 1 until n - 1
	
	double temp3 = ((F + temp1) / (pow(temp2,n))); // Calculates the final sum for the power n
	
	B = B + temp3; // End of old Function
	
	
	D_mac = 0.0; // New Calculation
	
	double D_macIteration = 0.0; // Stores D_mac sum from 1 to n - 1
	
	for (int j = 1; j <= n - 1; ++j) {
		double time_j = 0.5*j;
		D_macIteration = D_macIteration + (j / 2) * (temp1 / pow(temp2,j));
	} // Calculate sums from 1 until n - 1
	
	double Final_D_mac_Sum = (n / 2.0) * ((F + temp1) / pow(temp2,n)); // Calculates final n sum
	
	double D_mac_temp = Final_D_mac_Sum + D_macIteration; // Adds the n - 1 sums to the final sum thus having in total n sums
	
	D_mac = (1 / B) * D_mac_temp;
	
	D_mod = D_mac / temp2;
	
}

int main() {
	double B;
	double B2;
	double D_mac;
	double D_mod;
	
	price_from_yield1(100,0,5,3,B); // Old Function
	
	cout << B << endl;
	
	price_from_yield2(100,0,5,15,B2,D_mac,D_mod); // New Function
	
	cout << B2 << endl;
	
	cout << D_mac << endl;
	cout << D_mod << endl;
	
	return 0;
}
