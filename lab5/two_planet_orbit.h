#ifndef TWOPLANETORBIT_H
#define TWOPLANETORBIT_H

#include <iostream>
#include <vector>
#include "runge_kutta.h"

using namespace std;


class TwoPlanetOrbit
{
private:
	const int _n_eqns = 8;

	// inputs
	// initial condition
	double _t_start;
	double _x1_start;
	double _y1_start;
	double _x2_start;
	double _y2_start;
	double _vx1_start;
	double _vy1_start;
	double _vx2_start;
	double _vy2_start;

	vector<double> _c_start;

	// calculation setting
	int _alg;
	double _dt;

	// physics condition
	static double _GM;
	static double _t_end;
	static double _m1;
	static double _m2;
	double _a1;
	double _a2;
	double _e1;
	double _e2;

	// output
	vector<double> _t;
	vector< vector<double> > _x;
	int _n_stps;
	double _period;

	static double f_x1(double t, const vector<double> &x);
	static double f_y1(double t, const vector<double> &x);
	static double f_x2(double t, const vector<double> &x);
	static double f_y2(double t, const vector<double> &x);
	static double f_vx1(double t, const vector<double> &x);
	static double f_vy1(double t, const vector<double> &x);
	static double f_vx2(double t, const vector<double> &x);
	static double f_vy2(double t, const vector<double> &x);

	static bool stop(double t, const vector<double> &x);

public:
	TwoPlanetOrbit();
	TwoPlanetOrbit(double a1, double e1, double m1, double a2, double e2, double m2);
	~TwoPlanetOrbit();

	void set_t_start(double t_start);
	void set_x1_start(double x1_start);
	void set_y1_start(double y1_start);
	void set_x2_start(double x2_start);
	void set_y2_start(double y2_start);
	void set_vx1_start(double vx1_start);
	void set_vy1_start(double vy1_start);
	void set_vx2_start(double vx2_start);
	void set_vy2_start(double vy2_start);
	void set_alg(int alg);
	void set_dt(double dt);

	static void set_GM(double GM);
	static void set_t_end(double t_end);
	static void set_m1(double m1);
	static void set_m2(double m2);

	double get_t_start() const;
	double get_x1_start() const;
	double get_y1_start() const;
	double get_x2_start() const;
	double get_y2_start() const;
	double get_vx1_start() const;
	double get_vy1_start() const;
	double get_vx2_start() const;
	double get_vy2_start() const;
	int get_alg() const;
	double get_dt() const;

	static double get_GM();
	static double get_t_end();
	static double get_m1();
	static double get_m2();

	bool check() const;
	void cal();
	double cal_period();

	vector<double> get_t() const;
	vector< vector<double> > get_x() const;
	int get_n_stps() const;
	double get_period() const;

};


inline void TwoPlanetOrbit::set_t_start(double t_start) { _t_start = t_start; }
inline void TwoPlanetOrbit::set_x1_start(double x1_start) { _x1_start = x1_start; _c_start[0] = x1_start; }
inline void TwoPlanetOrbit::set_y1_start(double y1_start) { _y1_start = y1_start; _c_start[1] = y1_start; }
inline void TwoPlanetOrbit::set_x2_start(double x2_start) { _x2_start = x2_start; _c_start[2] = x2_start; }
inline void TwoPlanetOrbit::set_y2_start(double y2_start) { _y2_start = y2_start; _c_start[3] = y2_start; }
inline void TwoPlanetOrbit::set_vx1_start(double vx1_start) { _vx1_start = vx1_start; _c_start[4] = vx1_start; }
inline void TwoPlanetOrbit::set_vy1_start(double vy1_start) { _vy1_start = vy1_start; _c_start[5] = vy1_start; }
inline void TwoPlanetOrbit::set_vx2_start(double vx2_start) { _vx2_start = vx2_start; _c_start[6] = vx2_start; }
inline void TwoPlanetOrbit::set_vy2_start(double vy2_start) { _vy2_start = vy2_start; _c_start[7] = vy2_start; }
inline void TwoPlanetOrbit::set_alg(int alg) { _alg = alg; }
inline void TwoPlanetOrbit::set_dt(double dt) { _dt = dt; }

inline void TwoPlanetOrbit::set_GM(double GM) { _GM = GM; }
inline void TwoPlanetOrbit::set_t_end(double t_end) { _t_end = t_end; }

inline double TwoPlanetOrbit::get_t_start() const { return _t_start; }
inline double TwoPlanetOrbit::get_x1_start() const { return _x1_start; }
inline double TwoPlanetOrbit::get_y1_start() const { return _y1_start; }
inline double TwoPlanetOrbit::get_x2_start() const { return _x2_start; }
inline double TwoPlanetOrbit::get_y2_start() const { return _y2_start; }
inline double TwoPlanetOrbit::get_vx1_start() const { return _vx1_start; }
inline double TwoPlanetOrbit::get_vy1_start() const { return _vy1_start; }
inline double TwoPlanetOrbit::get_vx2_start() const { return _vx2_start; }
inline double TwoPlanetOrbit::get_vy2_start() const { return _vy2_start; }
inline int TwoPlanetOrbit::get_alg() const { return _alg; }
inline double TwoPlanetOrbit::get_dt() const { return _dt; }

inline double TwoPlanetOrbit::get_GM() { return _GM; }
inline double TwoPlanetOrbit::get_t_end() { return _t_end; }

inline vector<double> TwoPlanetOrbit::get_t() const { return _t; }
inline vector< vector<double> > TwoPlanetOrbit:: get_x() const { return _x; }
inline int TwoPlanetOrbit::get_n_stps() const { return _n_stps; }
inline double TwoPlanetOrbit::get_period() const { return _period; }


#endif
