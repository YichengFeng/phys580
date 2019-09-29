#ifndef PLANETSPINORBIT_H
#define PLANETSPINORBIT_H

#include <iostream>
#include <vector>
#include "runge_kutta.h"

using namespace std;


class PlanetSpinOrbit
{
private:
	const int _n_eqns = 6;

	// inputs
	// initial condition
	double _t_start;
	double _x_start;
	double _y_start;
	double _theta_start;
	double _vx_start;
	double _vy_start;
	double _omega_start;

	vector<double> _c_start;

	// calculation setting
	int _alg;
	double _dt;

	// physics condition
	static double _GM;
	static double _t_end;
	double _a;
	double _e;
	double _mp;
	double _energy;

	// output
	vector<double> _t;
	vector< vector<double> > _x; // x[0]: x; x[1]: y; x[2]: theta; x[3]: vx; x[4]: vy; x[5]: omega
	int _n_stps;
	double _period;

	static double f_x(double t, const vector<double> &x);
	static double f_y(double t, const vector<double> &x);
	static double f_theta(double t, const vector<double> &x);
	static double f_vx(double t, const vector<double> &x);
	static double f_vy(double t, const vector<double> &x);
	static double f_omega(double t, const vector<double> &x);

	static bool stop(double t, const vector<double> &x);

public:
	PlanetSpinOrbit();
	PlanetSpinOrbit(double a, double e, double mp);
	~PlanetSpinOrbit();

	void increase_v_start(double r);

	void set_t_start(double t_start);
	void set_x_start(double x_start);
	void set_y_start(double y_start);
	void set_theta_start(double theta_start);
	void set_vx_start(double vx_start);
	void set_vy_start(double vy_start);
	void set_omega_start(double omega_start);
	void set_alg(int alg);
	void set_dt(double dt);

	static void set_GM(double GM);
	static void set_t_end(double t_end);

	double get_t_start() const;
	double get_x_start() const;
	double get_y_start() const;
	double get_theta_start() const;
	double get_vx_start() const;
	double get_vy_start() const;
	double get_omega_start() const;
	int get_alg() const;
	double get_dt() const;
	double get_energy() const;

	static double get_GM();
	static double get_t_end();

	bool check() const;
	void cal();
	double cal_period();

	vector<double> get_t() const;
	vector< vector<double> > get_x() const;
	int get_n_stps() const;
	double get_period() const;

};


inline void PlanetSpinOrbit::set_t_start(double t_start) { _t_start = t_start; }
inline void PlanetSpinOrbit::set_x_start(double x_start) { _x_start = x_start; _c_start[0] = x_start; }
inline void PlanetSpinOrbit::set_y_start(double y_start) { _y_start = y_start; _c_start[1] = y_start; }
inline void PlanetSpinOrbit::set_theta_start(double theta_start) { _theta_start = theta_start; _c_start[2] = theta_start; }
inline void PlanetSpinOrbit::set_vx_start(double vx_start) { _vx_start = vx_start; _c_start[3] = vx_start; }
inline void PlanetSpinOrbit::set_vy_start(double vy_start) { _vy_start = vy_start; _c_start[4] = vy_start; }
inline void PlanetSpinOrbit::set_omega_start(double omega_start) { _omega_start = omega_start; _c_start[5] = omega_start; }
inline void PlanetSpinOrbit::set_alg(int alg) { _alg = alg; }
inline void PlanetSpinOrbit::set_dt(double dt) { _dt = dt; }

inline void PlanetSpinOrbit::set_GM(double GM) { _GM = GM; }
inline void PlanetSpinOrbit::set_t_end(double t_end) { _t_end = t_end; }

inline double PlanetSpinOrbit::get_t_start() const { return _t_start; }
inline double PlanetSpinOrbit::get_x_start() const { return _x_start; }
inline double PlanetSpinOrbit::get_y_start() const { return _y_start; }
inline double PlanetSpinOrbit::get_theta_start() const { return _theta_start; }
inline double PlanetSpinOrbit::get_vx_start() const { return _vx_start; }
inline double PlanetSpinOrbit::get_vy_start() const { return _vy_start; }
inline double PlanetSpinOrbit::get_omega_start() const { return _omega_start; }
inline int PlanetSpinOrbit::get_alg() const { return _alg; }
inline double PlanetSpinOrbit::get_dt() const { return _dt; }
inline double PlanetSpinOrbit::get_energy() const { return _energy; }

inline double PlanetSpinOrbit::get_GM() { return _GM; }
inline double PlanetSpinOrbit::get_t_end() { return _t_end; }

inline vector<double> PlanetSpinOrbit::get_t() const { return _t; }
inline vector< vector<double> > PlanetSpinOrbit:: get_x() const { return _x; }
inline int PlanetSpinOrbit::get_n_stps() const { return _n_stps; }
inline double PlanetSpinOrbit::get_period() const { return _period; }


#endif
