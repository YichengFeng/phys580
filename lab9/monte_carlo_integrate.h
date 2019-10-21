#ifndef MONTECARLOINTEGRATE_H
#define MONTECARLOINTEGRATE_H

#include <iostream>
#include <vector>
#include <cmath>
#include <time.h>

using namespace std;


class MonteCarloIntegrate
{
private:
	// input
	int _dim;
	vector<double> _x1;
	vector<double> _x2;
	double (*_f)(const vector<double> &x);
	int _trial;
	int _alg;
	int _n;

	// middle
	double _vol;

	// output
	double _output;
	double _error;

public:
	MonteCarloIntegrate();
	MonteCarloIntegrate(int dim, vector<double> x1, vector<double> x2, double (*f)(const vector<double> &x));

	void set_x1(vector<double> x1);
	void set_x2(vector<double> x2);
	void set_f(double (*f)(const vector<double> &x));
	void set_alg(int alg);
	void set_n(int n);
	void set_dim(int dim);
	void set_trial(int trial);

	vector<double> get_x1(void) const;
	vector<double> get_x2(void) const;
	int get_alg(void) const;
	int get_n(void) const;
	int get_dim(void) const;
	int get_trial(void) const;
	double get_output(void) const;
	double get_error(void) const;

	bool check() const;
	double cal_once(); 
	double cal(); 
};

inline void MonteCarloIntegrate::set_x1(vector<double> x1) { _x1 = x1; }
inline void MonteCarloIntegrate::set_x2(vector<double> x2) { _x2 = x2; }
inline void MonteCarloIntegrate::set_f(double (*f)(const vector<double> &x)) { _f = f; }
inline void MonteCarloIntegrate::set_alg(int alg) { _alg = alg; }
inline void MonteCarloIntegrate::set_dim(int dim) { _dim = dim; }
inline void MonteCarloIntegrate::set_trial(int trial) { _trial = trial; }
inline void MonteCarloIntegrate::set_n(int n) { _n = n; }

inline vector<double> MonteCarloIntegrate::get_x1(void) const { return _x1; }
inline vector<double> MonteCarloIntegrate::get_x2(void) const { return _x2; }
inline int MonteCarloIntegrate::get_alg(void) const { return _alg; }
inline int MonteCarloIntegrate::get_dim(void) const { return _dim; }
inline int MonteCarloIntegrate::get_trial(void) const { return _trial; }
inline int MonteCarloIntegrate::get_n(void) const { return _n; }
inline double MonteCarloIntegrate::get_output(void) const { return _output; }
inline double MonteCarloIntegrate::get_error(void) const { return _error; }


#endif
