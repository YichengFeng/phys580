#ifndef RELATIVITYPRECESSION_H
#define RELATIVITYPRECESSION_H

#include <iostream>
#include <vector>
#include "runge_kutta.h"

using namespace std;


class RelativityPrecession
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
	static double _alpha;
	static double _t_end;
	double _a;
	double _e;
	double _mp;
	double _energy;

	// middle
	vector<double> _rr;

	// output
	vector<double> _t;
	vector< vector<double> > _x; // x[0]: x; x[1]: y; x[2]: vx; x[3]: vy
	int _n_stps;
	double _period;
	vector<double> _perihelion_t;
	vector<double> _perihelion_x;
	vector<double> _perihelion_y;
	vector<double> _perihelion_theta;

	static double f_x(double t, const vector<double> &x);
	static double f_y(double t, const vector<double> &x);
	static double f_vx(double t, const vector<double> &x);
	static double f_vy(double t, const vector<double> &x);

	static bool stop(double t, const vector<double> &x);

public:
	RelativityPrecession();
	RelativityPrecession(double a, double e, double mp);
	~RelativityPrecession();

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
	static void set_alpha(double alpha);

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
	static double get_alpha();

	bool check() const;
	void cal();
	double cal_period();
	void cal_perihelion();

	vector<double> get_t() const;
	vector< vector<double> > get_x() const;
	int get_n_stps() const;
	double get_period() const;

	vector<double> get_perihelion_t() const;
	vector<double> get_perihelion_x() const;
	vector<double> get_perihelion_y() const;
	vector<double> get_perihelion_theta() const;
};


inline void RelativityPrecession::set_t_start(double t_start) { _t_start = t_start; }
inline void RelativityPrecession::set_x_start(double x_start) { _x_start = x_start; _c_start[0] = x_start; }
inline void RelativityPrecession::set_y_start(double y_start) { _y_start = y_start; _c_start[1] = y_start; }
inline void RelativityPrecession::set_vx_start(double vx_start) { _vx_start = vx_start; _c_start[2] = vx_start; }
inline void RelativityPrecession::set_vy_start(double vy_start) { _vy_start = vy_start; _c_start[3] = vy_start; }
inline void RelativityPrecession::set_alg(int alg) { _alg = alg; }
inline void RelativityPrecession::set_dt(double dt) { _dt = dt; }

inline void RelativityPrecession::set_GM(double GM) { _GM = GM; }
inline void RelativityPrecession::set_t_end(double t_end) { _t_end = t_end; }
inline void RelativityPrecession::set_alpha(double alpha) { _alpha = alpha; }

inline double RelativityPrecession::get_t_start() const { return _t_start; }
inline double RelativityPrecession::get_x_start() const { return _x_start; }
inline double RelativityPrecession::get_y_start() const { return _y_start; }
inline double RelativityPrecession::get_vx_start() const { return _vx_start; }
inline double RelativityPrecession::get_vy_start() const { return _vy_start; }
inline int RelativityPrecession::get_alg() const { return _alg; }
inline double RelativityPrecession::get_dt() const { return _dt; }
inline double RelativityPrecession::get_energy() const { return _energy; }

inline double RelativityPrecession::get_GM() { return _GM; }
inline double RelativityPrecession::get_t_end() { return _t_end; }
inline double RelativityPrecession::get_alpha() { return _alpha; }

inline vector<double> RelativityPrecession::get_t() const { return _t; }
inline vector< vector<double> > RelativityPrecession::get_x() const { return _x; }
inline int RelativityPrecession::get_n_stps() const { return _n_stps; }
inline double RelativityPrecession::get_period() const { return _period; }

inline vector<double> RelativityPrecession::get_perihelion_t() const { return _perihelion_t; }
inline vector<double> RelativityPrecession::get_perihelion_x() const { return _perihelion_x; }
inline vector<double> RelativityPrecession::get_perihelion_y() const { return _perihelion_y; }
inline vector<double> RelativityPrecession::get_perihelion_theta() const { return _perihelion_theta; }

#endif
