#ifndef PLANETORBIT_H
#define PLANETORBIT_H

#include <iostream>
#include <vector>
#include "runge_kutta.h"

using namespace std;


class PlanetOrbit
{
private:
	const int _n_eqns = 4;

	// inputs
	// initial condition
	double _t_start;
	double _x_start;
	double _y_start;
	double _vx_start;
	double _vy_start;

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
	vector< vector<double> > _x; // x[0]: x; x[1]: y; x[2]: vx; x[3]: vy
	int _n_stps;
	double _period;

	static double f_x(double t, const vector<double> &x);
	static double f_y(double t, const vector<double> &x);
	static double f_vx(double t, const vector<double> &x);
	static double f_vy(double t, const vector<double> &x);

	static bool stop(double t, const vector<double> &x);

public:
	PlanetOrbit();
	PlanetOrbit(double a, double e, double mp);
	~PlanetOrbit();

	void increase_v_start(double r);

	void set_t_start(double t_start);
	void set_x_start(double x_start);
	void set_y_start(double y_start);
	void set_vx_start(double vx_start);
	void set_vy_start(double vy_start);
	void set_alg(int alg);
	void set_dt(double dt);

	static void set_GM(double GM);
	static void set_t_end(double t_end);

	double get_t_start() const;
	double get_x_start() const;
	double get_y_start() const;
	double get_vx_start() const;
	double get_vy_start() const;
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


inline void PlanetOrbit::set_t_start(double t_start) { _t_start = t_start; }
inline void PlanetOrbit::set_x_start(double x_start) { _x_start = x_start; _c_start[0] = x_start; }
inline void PlanetOrbit::set_y_start(double y_start) { _y_start = y_start; _c_start[1] = y_start; }
inline void PlanetOrbit::set_vx_start(double vx_start) { _vx_start = vx_start; _c_start[2] = vx_start; }
inline void PlanetOrbit::set_vy_start(double vy_start) { _vy_start = vy_start; _c_start[3] = vy_start; }
inline void PlanetOrbit::set_alg(int alg) { _alg = alg; }
inline void PlanetOrbit::set_dt(double dt) { _dt = dt; }

inline void PlanetOrbit::set_GM(double GM) { _GM = GM; }
inline void PlanetOrbit::set_t_end(double t_end) { _t_end = t_end; }

inline double PlanetOrbit::get_t_start() const { return _t_start; }
inline double PlanetOrbit::get_x_start() const { return _x_start; }
inline double PlanetOrbit::get_y_start() const { return _y_start; }
inline double PlanetOrbit::get_vx_start() const { return _vx_start; }
inline double PlanetOrbit::get_vy_start() const { return _vy_start; }
inline int PlanetOrbit::get_alg() const { return _alg; }
inline double PlanetOrbit::get_dt() const { return _dt; }
inline double PlanetOrbit::get_energy() const { return _energy; }

inline double PlanetOrbit::get_GM() { return _GM; }
inline double PlanetOrbit::get_t_end() { return _t_end; }

inline vector<double> PlanetOrbit::get_t() const { return _t; }
inline vector< vector<double> > PlanetOrbit:: get_x() const { return _x; }
inline int PlanetOrbit::get_n_stps() const { return _n_stps; }
inline double PlanetOrbit::get_period() const { return _period; }


#endif
