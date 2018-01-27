#include<iostream>
#include <vector>
#include <cmath>
#include <stdlib.h>
#include <iomanip> 
#include <fstream> 

using namespace std;

//James Goodman
// CS 365 FINAL EXAM/HOMEWORK 8
// 12/20/2017

class Derivative
{
public:
	virtual ~Derivative() {}
	
	virtual double TerminalPayoff(double S) const { return 0; }
	virtual int ValuationTests(double S, double & V) const { return 0; }
	
	// data
	double r;
	double q;
	double sigma;
	double T;
	
protected:
	Derivative() { r = 0; q = 0; sigma = 0; T = 0; }
};

class Option : public Derivative
{
public:
	Option() { K = 0; isCall = false; isAmerican = false; }
	virtual ~Option() {}
	virtual double TerminalPayoff(double S) const;
	virtual int ValuationTests(double S, double &V) const;
	
	// data
	double K;
	bool isCall;
	bool isAmerican;
};

double Option::TerminalPayoff(double S) const
{
	double intrinsic = 0.0;
	
	if (isCall) {
		if (S > K) {
			intrinsic = S - K;
		}
	}
		
	if (!isCall) {
		if (S < K) {
			intrinsic = K - S;
		}
	}
		
	return intrinsic;
// *** RETURN TERMINAL PAYOFF FOR PUT OR CALL OPTION ***
}

int Option::ValuationTests(double S, double &V) const
{
	double payoff = 0;
	
	if (isAmerican) {
		if (isCall) {
			payoff = (S - K);
			if (payoff > V) {
				V = payoff;
			} 
		}
				
		if (!isCall) {
			payoff = (K - S);
			if (payoff > V) {
				V = payoff;
			} 
		}
	}
	
	return 0;
// *** TEST IF THE VALUE OF V SHOULD BE UPDATED TO THE INTRINSIC VALUE ***
}

class BinomialModel
{
public:
	BinomialModel(int n);
	~BinomialModel();
	int FairValue(int n, Derivative * p_derivative, double S, double t0, double & V);
	int ImpliedVolatility(int n, Derivative * p_derivative, double S, double t0,
						  double target, double & implied_vol, int & num_iter);

private:
	// methods
	void Clear();
	int Allocate(int n);
	int ImpliedVolatilityPrivate(int n, Derivative * p_derivative, double S, double t0,
						  double target, double & implied_vol, int & num_iter);

	// data
	int n_tree;
	double **stock_nodes;
	double **value_nodes;
	const double tol = 0.0001;
	const int max_iter = 100;
};

BinomialModel::BinomialModel(int n)
{
	n_tree = 0;
	stock_nodes = NULL;
	value_nodes = NULL;
	Allocate(n);
}

BinomialModel::~BinomialModel()
{
	Clear();
}

void BinomialModel::Clear()
{	
	delete [] stock_nodes;
	delete [] value_nodes;
	
// *** WRITE THE FUNCTION TO RELEASE ALLOCATED MEMORY ***
}

int BinomialModel::Allocate(int n)
{
	if (n <= n_tree) return 0;
	
	// deallocate old tree
	Clear();
	
	// allocate memory
	n_tree = n;
	//cout << n_tree << endl;
	
	// allocate memory
	int i = 0;
	int j = 0;
	stock_nodes = new double*[n_tree+1];
	value_nodes = new double*[n_tree+1];
	
	double *S_tmp;
	double *V_tmp;
	
	for (i = 0; i <= n_tree; ++i) {
		stock_nodes[i] = new double[n_tree+1];
		value_nodes[i] = new double[n_tree+1];
		S_tmp = stock_nodes[i];
		V_tmp = value_nodes[i];
		
		for (j = 0; j <= n_tree; ++j) {
			S_tmp[j] = 0;
			V_tmp[j] = 0;
		}
	}

// *** WRITE THE FUNCTION TO ALLOCATE NEW MEMORY ***
}

int BinomialModel::FairValue(int n, Derivative * p_derivative, double S, double t0, double & V) {
	
	int rc = 0;
	V = 0.0;
		
	// validation checks
	if ((n < 1) || (S <= 0) || (p_derivative == NULL)) {
		return 1;
	}
	
	// declaration of local variables (I use S_tmp and V_tmp)
	double *S_tmp;
	double *V_tmp;
	int i;
	int j;
	
	// calculate parameters
	double dt = (p_derivative->T-t0)/double(n);
	double df = exp(-p_derivative->r*dt);
	double growth = exp((p_derivative->r-p_derivative->q)*dt);
	double u = exp(p_derivative->sigma*sqrt(dt));
	double d = 1.0/u;
	double p_prob = (growth - d)/(u-d);
	double q_prob = 1.0 - p_prob;
	
	// more validation checks
	if (p_prob < 0.0) {
		return 1;
	}
	
	if (p_prob > 1.0) {
		return 1;
	}
	
	// allocate memory if required (call Allocate(n))
	Allocate(n);
	
	// set up stock prices in tree
	S_tmp = stock_nodes[0];
	S_tmp[0] = S;
	
	for (i = 1; i <= n; ++i) {
		double * prev = stock_nodes[i-1];
		S_tmp = stock_nodes[i];
		S_tmp[0] = prev[0] * d;
		for (j = 1; j <= n; ++j) {
			S_tmp[j] = S_tmp[j-1]*u*u;
		}
	}
	
	// set terminal payoff (call virtual function in derivative class to calculate payoff)
	i = n;
	S_tmp = stock_nodes[i];
	V_tmp = value_nodes[i];
	
	for (j = 0; j <= n; ++j) {
		V_tmp[j] = p_derivative->TerminalPayoff(S_tmp[j]);
	}

	// valuation loop (call virtual function in derivative class for valuation tests)
	for (i = n-1; i >= 0; --i) {
		S_tmp = stock_nodes[i];
		V_tmp = value_nodes[i];
		double * V_next = value_nodes[i+1];
		for (j = 0; j <= i; ++j) {
			V_tmp[j] = df*(p_prob*V_next[j+1] + q_prob*V_next[j]);
			p_derivative->ValuationTests(S_tmp[j], V_tmp[j]); // VALUATION TESTS
		}
	}

	// derivative fair value
	V_tmp = value_nodes[0];
	V = V_tmp[0];
	
	// deallocate memory (if necessary)
	//Clear();
	
	return 0;

}

int BinomialModel::ImpliedVolatilityPrivate(int n, Derivative * p_derivative,
											double S, double t0, double target, double & implied_vol, int & num_iter) {
													
	implied_vol = 0;
	num_iter = 0;
	
	double sigma_low = 0.01; // 1%
	double sigma_high = 3.0; // 300%
	double FV_low = 0;
	double FV_high = 0;
	double FV = 0;
										
	p_derivative->sigma = sigma_low;
	FairValue(n, p_derivative, S, t0, FV_low);
	double diff_FV_low = FV_low - target;
	
	if (fabs(diff_FV_low) <= tol) {
		implied_vol = p_derivative->sigma;
		return 0;
	}
	
	p_derivative->sigma = sigma_high;
	FairValue(n,p_derivative,S,t0,FV_high);
	double diff_FV_high = FV_high - target;
	
	if (fabs(diff_FV_high) <= tol) {
		implied_vol = p_derivative->sigma;
		return 0;
	}
	
	if (diff_FV_low * diff_FV_high > 0) {
		implied_vol = 0;
		return 1; // fail
	}
	
	int i;
	double diff_FV;
	
	for (i = 0; i < max_iter; i++) {
		p_derivative->sigma = 0.5*(sigma_low + sigma_high);
		FairValue(n, p_derivative, S, t0, FV);
		diff_FV = FV - target;
		
		if (fabs(diff_FV) <= tol) {
			implied_vol = p_derivative->sigma;
			num_iter = i;
			return 0;
		}
		
		if ((diff_FV_low * diff_FV) > 0) {
			sigma_low = p_derivative->sigma;
		}
		
		else {
			sigma_high = p_derivative->sigma;
		}
		
		if (fabs(sigma_low - sigma_high) <= tol) {
			implied_vol = p_derivative->sigma;
			num_iter = i;
			return 0;
		}
		
	}
	
	num_iter = max_iter;
	implied_vol = 0;
	return 1;

}

int BinomialModel::ImpliedVolatility(int n, Derivative * p_derivative, double S, double t0,
double target, double & implied_vol, int & num_iter)
{
	int rc = 0;
	const double saved_vol = p_derivative->sigma;
	rc = ImpliedVolatilityPrivate(n, p_derivative, S, t0, target, implied_vol, num_iter);
	p_derivative->sigma = saved_vol;
	return rc;
}

int binomial_test()
{
	int rc = 0;
	// output file
	//std::ofstream ofs("output.txt");
	double S = 100;
	double K = 100;
	double r = 0.1;
	double q = 0.1;
	double sigma = 0.5;
	double T = 0.4;
	double t0 = 0;
	Option Eur_put;
	Eur_put.r = r;
	Eur_put.q = q;
	Eur_put.sigma = sigma;
	Eur_put.T = T;
	Eur_put.K = K;
	Eur_put.isCall = false;
	Eur_put.isAmerican = false;
	Option Am_put;
	Am_put.r = r;
	Am_put.q = q;
	Am_put.sigma = sigma;
	Am_put.T = T;
	Am_put.K = K;
	Am_put.isCall = false;
	Am_put.isAmerican = true;
	Option Eur_call;
	Eur_call.r = r;
	Eur_call.q = q;
	Eur_call.sigma = sigma;
	Eur_call.T = T;
	Eur_call.K = K;
	Eur_call.isCall = true;
	Eur_call.isAmerican = false;
	Option Am_call;
	Am_call.r = r;
	Am_call.q = q;
	Am_call.sigma = sigma;
	Am_call.T = T;
	Am_call.K = K;
	Am_call.isCall = true;
	Am_call.isAmerican = true;
	double FV_Am_put = 0;
	double FV_Eur_put = 0;
	double FV_Am_call = 0;
	double FV_Eur_call = 0;
	int n = 4;
	BinomialModel binom(n);
	
	/*)double dS = 0.1;
	int imax = 2000;
	int i;
	for (i = 1; i <= imax; ++i) {
		S = i*dS;
		rc = binom.FairValue(n, &Am_put, S, t0, FV_Am_put);
		rc = binom.FairValue(n, &Eur_put, S, t0, FV_Eur_put);
		rc = binom.FairValue(n, &Am_call, S, t0, FV_Am_call);
		rc = binom.FairValue(n, &Eur_call, S, t0, FV_Eur_call);
		ofs << std::setw(16) << S << " ";
		ofs << std::setw(16) << FV_Am_put << " ";
		ofs << std::setw(16) << FV_Eur_put << " ";
		ofs << std::setw(16) << FV_Am_call << " ";
		ofs << std::setw(16) << FV_Eur_call << " ";
		ofs << std::endl;
	}**/
	
	rc = binom.FairValue(n, &Am_put, S, t0, FV_Am_put);
	rc = binom.FairValue(n, &Eur_put, S, t0, FV_Eur_put);
	rc = binom.FairValue(n, &Am_call, S, t0, FV_Am_call);
	rc = binom.FairValue(n, &Eur_call, S, t0, FV_Eur_call);

	//cout << S << endl;
	//cout << FV_Am_put << endl;
	//cout << FV_Eur_put << endl;
	//cout << FV_Am_call << endl;
	//cout << FV_Eur_call << endl;
	//cout << endl;
	
	return 0;
}

int binomial_testhw9ImpliedVolatility() {
	
	double S = 105;
	double K = 100;
	double r = 0.05;
	double q = 0.02;
	double T = 0.5;
	double t0 = 0;
	double outputvol = 0;
	double outputvol2 = 0;
	int numIter = 0;
	
	int rc = 0;
	
	Option Eur_put;
	Eur_put.r = r;
	Eur_put.q = q;
	Eur_put.T = T;
	Eur_put.K = K;
	Eur_put.isCall = false;
	Eur_put.isAmerican = false;
	
	Option Am_put;
	Am_put.r = r;
	Am_put.q = q;
	Am_put.T = T;
	Am_put.K = K;
	Am_put.isCall = false;
	Am_put.isAmerican = true;
	
	Option Eur_call;
	Eur_call.r = r;
	Eur_call.q = q;
	Eur_call.T = T;
	Eur_call.K = K;
	Eur_call.isCall = true;
	Eur_call.isAmerican = false;
	
	Option Am_call;
	Am_call.r = r;
	Am_call.q = q;
	Am_call.T = T;
	Am_call.K = K;
	Am_call.isCall = true;
	Am_call.isAmerican = true;
	
	double FV_Am_put = 9;
	double FV_Eur_put = 9;
	double FV_Am_call = 9;
	double FV_Eur_call = 9;
	int n = 100;
	
	BinomialModel binom(n);

	rc = binom.ImpliedVolatility(n,&Eur_call,S,t0,FV_Eur_call,outputvol,numIter);
	
	cout << "The output vol of a European Call with Market price 9 dollars " << outputvol << endl;
	
	cout << "NumerIter: " <<  numIter << endl;
	
	rc = binom.ImpliedVolatility(n,&Am_call,S,t0,FV_Am_call,outputvol2,numIter);
	
	cout << "The output vol of a American put with Market price 9 dollars " << outputvol2 << endl;
	
	cout << "NumerIter: " <<  numIter << endl;
	
	return 0;
}



int main() {
	int value;
	int value2;
	//value = binomial_test();
	value2 = binomial_testhw9ImpliedVolatility();
	return 0;
}
