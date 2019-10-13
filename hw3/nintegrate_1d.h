#ifndef NINTEGRATE1D_H
#define NINTEGRATE1D_H

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;


class NIntegrate1D
{
private:
	// input
	double _x1;
	double _x2;
	double (*_f)(double x);
	int _alg;
	int _n;

	// middle
	double _sgn;

	// output
	double _output;

public:
	NIntegrate1D();
	NIntegrate1D(double x1, double x2, double (*f)(double x));

	void set_x1(double x1);
	void set_x2(double x2);
	void set_f(double (*f)(double x));
	void set_alg(int alg);
	void set_n(int n);

	double get_x1(void) const;
	double get_x2(void) const;
	int get_alg(void) const;
	int get_n(void) const;
	double get_output(void) const;

	void cal_trapzoidal();
	void cal_Simpson();
	void cal_Romberg();
	double cal(); 
};

inline void NIntegrate1D::set_x1(double x1) { _x1 = x1; }
inline void NIntegrate1D::set_x2(double x2) { _x2 = x2; }
inline void NIntegrate1D::set_f(double (*f)(double x)) { _f = f; }
inline void NIntegrate1D::set_alg(int alg) { _alg = alg; }
inline void NIntegrate1D::set_n(int n) { _n = n; }

inline double NIntegrate1D::get_x1(void) const { return _x1; }
inline double NIntegrate1D::get_x2(void) const { return _x2; }
inline int NIntegrate1D::get_alg(void) const { return _alg; }
inline int NIntegrate1D::get_n(void) const { return _n; }
inline double NIntegrate1D::get_output(void) const { return _output; }


#endif
