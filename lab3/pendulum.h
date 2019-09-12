#ifndef PENDULUM_H
#define PENDULUM_H

#include <iostream>
#include <vector>
#include <cmath>
#include "runge_kutta.h"

using namespace std;


class Pendulum
{
private:
	const int _n_eqns = 2;

	// inputs
	// initial condition
	int _alg; // 0: Euler Cromer; 1: Euler; 2: RK2; 4: RK4;
	double _t_start;
	double _theta_start;
	double _omega_start;

	vector<double> _x_start;

	double _mode;

	// physics condition
	static double _g;
	static double _l;
	static double _q;
	static double _F; // driving force: magnitude
	static double _O; // driving force: frequency
	static double _P; // driving force: phase
	static double _n_periods;

	// calcution setting
	double _dt;
	double _n_stps;

	// middle
	double _energy_start;

	// output
	vector<double> _t;
	vector<double> _theta;
	vector<double> _omega;
	vector<double> _energy;
	vector< vector<double> > _x;

	vector< vector<double> > _x_simplified_exact;
	vector< vector<double> > _x_simplified_error;
	vector<double> _energy_exact;
	vector<double> _energy_error;

	static double f_theta(double t, const vector<double> &x);
	static double f_omega_simplified(double t, const vector<double> &x);
	static double f_omega_physics(double t, const vector<double> &x);
	static double f_energy_simplified(double t, const vector<double> &x);
	static double f_energy_physics(double t, const vector<double> &x);

	// the exact solution for the simplified ODE
	static double f_omega_simplified_exact(double t, const vector<double> &x);
	static double f_theta_simplified_exact(double t, const vector<double> &x);
	static double f_energy_simplified_exact(double t, const vector<double> &x);

	static bool stop(double t, const vector<double> &x);

public:
	Pendulum();
	~Pendulum();

	void set_mode(int mode);
	void set_alg(int alg);
	void set_dt(double dt);
	void set_t_start(double t_start);
	void set_theta_start(double theta_start);
	void set_omega_start(double omega_start);

	static void set_g(double g);
	static void set_l(double l);
	static void set_q(double q);
	static void set_F(double F);
	static void set_O(double O);
	static void set_P(double P);
	static void set_n_periods(double n_periods);

	static double get_g();
	static double get_l();
	static double get_q();
	static double get_F();
	static double get_O();
	static double get_P();
	static double get_n_periods();

	bool check() const;
	void cal();
	double cal_period() const;
	double cal_analytical_period() const;

	vector<double> get_t() const;
	vector<double> get_energy() const;
	vector<double> get_energy_exact() const;
	vector<double> get_energy_error() const;
	vector< vector<double> > get_x() const;
	vector< vector<double> > get_x_simplified_exact() const;
	vector< vector<double> > get_x_simplified_error() const;
	int get_n_eqns() const;
	int get_n_stps() const;
	double get_dt() const;
};


inline void Pendulum::set_mode(int mode) { _mode = mode; }
inline void Pendulum::set_alg(int alg) { _alg = alg; }
inline void Pendulum::set_dt(double dt) { _dt = dt; }
inline void Pendulum::set_t_start(double t_start) { _t_start = t_start; }
inline void Pendulum::set_theta_start(double theta_start) { _theta_start = theta_start; _x_start[0] = theta_start; }
inline void Pendulum::set_omega_start(double omega_start) { _omega_start = omega_start; _x_start[1] = omega_start; }

inline void Pendulum::set_g(double g) { _g = g; }
inline void Pendulum::set_l(double l) { _l = l; }
inline void Pendulum::set_q(double q) { _q = q; }
inline void Pendulum::set_F(double F) { _F = F; }
inline void Pendulum::set_O(double O) { _O = O; }
inline void Pendulum::set_P(double P) { _P = P; }
inline void Pendulum::set_n_periods(double n_periods) { _n_periods = n_periods; }

inline double Pendulum::get_g() { return _g; }
inline double Pendulum::get_l() { return _l; }
inline double Pendulum::get_q() { return _q; }
inline double Pendulum::get_F() { return _F; }
inline double Pendulum::get_O() { return _O; }
inline double Pendulum::get_P() { return _P; }
inline double Pendulum::get_n_periods() { return _n_periods; }

inline vector<double> Pendulum::get_t() const { return _t; }
inline vector<double> Pendulum::get_energy() const { return _energy; }
inline vector<double> Pendulum::get_energy_exact() const { return _energy_exact; }
inline vector<double> Pendulum::get_energy_error() const { return _energy_error; }
inline vector< vector<double> > Pendulum::get_x() const { return _x; }
inline vector< vector<double> > Pendulum::get_x_simplified_exact() const { return _x_simplified_exact; }
inline vector< vector<double> > Pendulum::get_x_simplified_error() const { return _x_simplified_error; }
inline int Pendulum::get_n_eqns() const { return _n_eqns; }
inline int Pendulum::get_n_stps() const { return _n_stps; }
inline double Pendulum::get_dt() const { return _dt; }

#endif
