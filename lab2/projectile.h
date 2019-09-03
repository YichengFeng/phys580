#ifndef PROJECTILE_H
#define PROJECTILE_H

#include <iostream>
#include <vector>
#include <cmath>
#include "runge_kutta.h"

using namespace std;


class Projectile
{
private:
	const int _n_eqns = 4;

	int _mode;

	// inputs
	// initial condition
	double _t_start;
	double _x_start;
	double _y_start;
	double _vx_start;
	double _vy_start;
	int _rk_order;

	vector<double> _c_start;

	// physics condition
	static double _g;
	static double _B2_m;
	static double _Y0;
	static double _T0;
	static double _a;
	static double _gamma;

	// calculation setting
	double _dt;
	double _n_stps;

	// output
	vector<double> _t;
	vector< vector<double> > _x;

	double _range;

	static double f_x(double t, const vector<double> &x);
	static double f_y(double t, const vector<double> &x);
	static double f_vx_no_drag(double t, const vector<double> &x);
	static double f_vy_no_drag(double t, const vector<double> &x);
	static double f_vx_constant_drag(double t, const vector<double> &x);
	static double f_vy_constant_drag(double t, const vector<double> &x);
	static double f_vx_isothermal_drag(double t, const vector<double> &x);
	static double f_vy_isothermal_drag(double t, const vector<double> &x);
	static double f_vx_adiabatic_drag(double t, const vector<double> &x);
	static double f_vy_adiabatic_drag(double t, const vector<double> &x);

	static bool stop(double t, const vector<double> &x);

public:
	Projectile();
	~Projectile();

	void set_mode(int mode);
	void set_rk_order(int rk_order);
	void set_dt(double dt);
	void set_t_start(double t_start);
	void set_x_start(double x_start);
	void set_y_start(double y_start);
	void set_vx_start(double vx_start);
	void set_vy_start(double vy_start);
	void set_v_theta_start(double v_start, double theta_start);

	static void set_g(double g);
	static void set_B2_m(double B2_m);
	static void set_Y0(double Y0);
	static void set_a(double a);
	static void set_T0(double T0);
	static void set_gamma(double gamma);

	static double get_g();
	static double get_B2_m();
	static double get_Y0();
	static double get_a();
	static double get_T0();
	static double get_gamma();

	bool check() const;
	void cal();
	double cal_range();
	double search_theta_for_max_range(double v);

	vector<double> get_t() const;
	vector< vector<double> > get_x() const;
	int get_n_eqns() const;
	int get_n_stps() const;
	double get_dt() const;
	double get_range() const;
};


inline void Projectile::set_mode(int mode) { _mode = mode; }
inline void Projectile::set_rk_order(int rk_order) { _rk_order = rk_order; }
inline void Projectile::set_dt(double dt) { _dt = dt; }
inline void Projectile::set_t_start(double t_start) { _t_start = t_start; }
inline void Projectile::set_x_start(double x_start) { _x_start = x_start; _c_start[0] = x_start; }
inline void Projectile::set_y_start(double y_start) { _y_start = y_start; _c_start[1] = y_start; }
inline void Projectile::set_vx_start(double vx_start) { _vx_start = vx_start; _c_start[2] = vx_start; }
inline void Projectile::set_vy_start(double vy_start) { _vy_start = vy_start; _c_start[3] = vy_start; }
inline void Projectile::set_v_theta_start(double v_start, double theta_start) {
	_vx_start = v_start*cos(theta_start*M_PI/180.0);
	_vy_start = v_start*sin(theta_start*M_PI/180.0);
	_c_start[2] = _vx_start;
	_c_start[3] = _vy_start;
}

inline void Projectile::set_g(double g) { _g = g; }
inline void Projectile::set_B2_m(double B2_m) { _B2_m = B2_m; }
inline void Projectile::set_Y0(double Y0) { _Y0 = Y0; }
inline void Projectile::set_a(double a) { _a = a; }
inline void Projectile::set_T0(double T0) { _T0 = T0; }
inline void Projectile::set_gamma(double gamma) { _gamma = gamma; }

inline double Projectile::get_g() { return _g; }
inline double Projectile::get_B2_m() { return _B2_m; }
inline double Projectile::get_Y0() { return _Y0; }
inline double Projectile::get_a() { return _a; }
inline double Projectile::get_T0() { return _T0; }
inline double Projectile::get_gamma() { return _gamma; }

inline vector<double> Projectile::get_t() const { return _t; }
inline vector< vector<double> > Projectile::get_x() const { return _x; }
inline int Projectile::get_n_stps() const { return _n_stps; }
inline double Projectile::get_dt() const { return _dt; }
inline double Projectile::get_range() const { return _range; }


#endif
