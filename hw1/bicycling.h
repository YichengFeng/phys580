#ifndef BICYCLING_H
#define BICYCLING_H

#include <iostream>
#include <vector>
#include "runge_kutta.h"

using namespace std;


class Bicycling
{
private:
	const int _n_eqns = 1;
	int _n_stps;

	double _dt;
	double _t_start;
	double _v_start;
	vector<double> _x_start;

	vector<double> _t;
	vector< vector<double> > _x;

	static double _C;
	static double _P;
	static double _A;
	static double _m;
	static double _rho;

	static double f_v(double t, const vector<double> &x);
	static bool stop(double t, const vector<double> &x);

public:
	Bicycling();
	~Bicycling();

	void set_dt(double dt);
	void set_t_start(double t_start);
	void set_v_start(double v_start);

	double get_dt() const;
	double get_t_start() const;
	double get_v_start() const;

	static void set_C(double C);
	static void set_P(double P);
	static void set_A(double A);
	static void set_m(double m);
	static void set_rho(double rho);

	static double get_C();
	static double get_P();
	static double get_A();
	static double get_m();
	static double get_rho();

	bool check() const;
	void cal();

	vector<double> get_t() const;
	vector< vector<double> > get_x() const;
	int get_n_stps() const;
};

inline void Bicycling::set_dt(double dt) { _dt = dt; }
inline void Bicycling::set_t_start(double t_start) { _t_start = t_start; }
inline void Bicycling::set_v_start(double v_start) { _v_start = v_start; _x_start[0] = _v_start; }

inline double Bicycling::get_dt() const { return _dt; }
inline double Bicycling::get_t_start() const { return _t_start; }
inline double Bicycling::get_v_start() const { return _v_start; }

inline void Bicycling::set_C(double C) { _C = C; };
inline void Bicycling::set_P(double P) { _P = P; };
inline void Bicycling::set_A(double A) { _A = A; };
inline void Bicycling::set_m(double m) { _m = m; };
inline void Bicycling::set_rho(double rho) { _rho = rho; };

inline double Bicycling::get_C() { return _C; }
inline double Bicycling::get_P() { return _P; }
inline double Bicycling::get_A() { return _A; }
inline double Bicycling::get_m() { return _m; }
inline double Bicycling::get_rho() { return _rho; }

inline vector<double> Bicycling::get_t() const { return _t; }
inline vector< vector<double> > Bicycling::get_x() const { return _x; }
inline int Bicycling::get_n_stps() const { return _n_stps; }

#endif
