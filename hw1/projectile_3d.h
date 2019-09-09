#ifndef PROJECTILE_3D_H
#define PROJECTILE_3D_H

#include <iostream>
#include <vector>
#include <cmath>
#include "runge_kutta.h"

using namespace std;


class Projectile3D
{
private:
	const int _n_eqns = 6;

	int _mode;

	// inputs
	// initial condition
	double _t_start;
	double _x_start;
	double _y_start;
	double _z_start;
	double _vx_start;
	double _vy_start;
	double _vz_start;
	int _rk_order;

	vector<double> _c_start;

	// physics condition
	static double _g;
	static double _B2_m;
	static double _Y0;
	static double _T0;
	static double _T_ref;
	static double _a;
	static double _gamma;
	static double _wx;
	static double _wy;
	static double _wz;

	// calculation setting
	double _dt;
	double _n_stps;

	// output
	vector<double> _t;
	vector< vector<double> > _x;

	double _range;

	static double f_x(double t, const vector<double> &x);
	static double f_y(double t, const vector<double> &x);
	static double f_z(double t, const vector<double> &x);
	static double f_vx_no_drag(double t, const vector<double> &x);
	static double f_vy_no_drag(double t, const vector<double> &x);
	static double f_vz_no_drag(double t, const vector<double> &x);
	static double f_vx_constant_drag(double t, const vector<double> &x);
	static double f_vy_constant_drag(double t, const vector<double> &x);
	static double f_vz_constant_drag(double t, const vector<double> &x);
	static double f_vx_isothermal_drag(double t, const vector<double> &x);
	static double f_vy_isothermal_drag(double t, const vector<double> &x);
	static double f_vz_isothermal_drag(double t, const vector<double> &x);
	static double f_vx_adiabatic_drag(double t, const vector<double> &x);
	static double f_vy_adiabatic_drag(double t, const vector<double> &x);
	static double f_vz_adiabatic_drag(double t, const vector<double> &x);

	static bool stop(double t, const vector<double> &x);

public:
	Projectile3D();
	~Projectile3D();

	void set_mode(int mode);
	void set_rk_order(int rk_order);
	void set_dt(double dt);
	void set_t_start(double t_start);
	void set_x_start(double x_start);
	void set_y_start(double y_start);
	void set_z_start(double z_start);
	void set_vx_start(double vx_start);
	void set_vy_start(double vy_start);
	void set_vz_start(double vz_start);
	void set_v_theta_start(double v_start, double theta_start);

	static void set_g(double g);
	static void set_B2_m(double B2_m);
	static void set_Y0(double Y0);
	static void set_a(double a);
	static void set_T0(double T0);
	static void set_T_ref(double T_ref);
	static void set_gamma(double gamma);
	static void set_omega_cartesian(double wx, double wy, double wz);
	static void set_omega_polar(double w, double la, double phi);

	static double get_g();
	static double get_B2_m();
	static double get_Y0();
	static double get_a();
	static double get_T0();
	static double get_T_ref();
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


inline void Projectile3D::set_mode(int mode) { _mode = mode; }
inline void Projectile3D::set_rk_order(int rk_order) { _rk_order = rk_order; }
inline void Projectile3D::set_dt(double dt) { _dt = dt; }
inline void Projectile3D::set_t_start(double t_start) { _t_start = t_start; }
inline void Projectile3D::set_x_start(double x_start) { _x_start = x_start; _c_start[0] = x_start; }
inline void Projectile3D::set_y_start(double y_start) { _y_start = y_start; _c_start[1] = y_start; }
inline void Projectile3D::set_z_start(double z_start) { _z_start = z_start; _c_start[2] = z_start; }
inline void Projectile3D::set_vx_start(double vx_start) { _vx_start = vx_start; _c_start[3] = vx_start; }
inline void Projectile3D::set_vy_start(double vy_start) { _vy_start = vy_start; _c_start[4] = vy_start; }
inline void Projectile3D::set_vz_start(double vz_start) { _vz_start = vz_start; _c_start[5] = vz_start; }
inline void Projectile3D::set_v_theta_start(double v_start, double theta_start) {
	_vx_start = v_start*cos(theta_start*M_PI/180.0);
	_vy_start = v_start*sin(theta_start*M_PI/180.0);
	_c_start[3] = _vx_start;
	_c_start[4] = _vy_start;
}

inline void Projectile3D::set_g(double g) { _g = g; }
inline void Projectile3D::set_B2_m(double B2_m) { _B2_m = B2_m; }
inline void Projectile3D::set_Y0(double Y0) { _Y0 = Y0; }
inline void Projectile3D::set_a(double a) { _a = a; }
inline void Projectile3D::set_T0(double T0) { _T0 = T0; }
inline void Projectile3D::set_T_ref(double T_ref) { _T_ref = T_ref; }
inline void Projectile3D::set_gamma(double gamma) { _gamma = gamma; }
inline void Projectile3D::set_omega_cartesian(double wx, double wy, double wz) { _wx = wx; _wy = wy; _wz = wz; }
inline void Projectile3D::set_omega_polar(double w, double la, double phi) {
	_wx = +w*cos(la*M_PI/180)*sin(phi*M_PI/180);
	_wy = +w*sin(la*M_PI/180);
	_wz = -w*cos(la*M_PI/180)*cos(phi*M_PI/180);
}

inline double Projectile3D::get_g() { return _g; }
inline double Projectile3D::get_B2_m() { return _B2_m; }
inline double Projectile3D::get_Y0() { return _Y0; }
inline double Projectile3D::get_a() { return _a; }
inline double Projectile3D::get_T0() { return _T0; }
inline double Projectile3D::get_T_ref() { return _T_ref; }
inline double Projectile3D::get_gamma() { return _gamma; }

inline vector<double> Projectile3D::get_t() const { return _t; }
inline vector< vector<double> > Projectile3D::get_x() const { return _x; }
inline int Projectile3D::get_n_stps() const { return _n_stps; }
inline double Projectile3D::get_dt() const { return _dt; }
inline double Projectile3D::get_range() const { return _range; }


#endif
