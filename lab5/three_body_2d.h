#ifndef THREEBODY2D_H
#define THREEBODY2D_H

#include <iostream>
#include <vector>
#include "runge_kutta.h"

using namespace std;


class ThreeBody2D
{
private:
	const int _n_eqns = 12;

	// inputs
	// initial condition
	double _t_start;
	double _x1_start;
	double _y1_start;
	double _x2_start;
	double _y2_start;
	double _x3_start;
	double _y3_start;
	double _vx1_start;
	double _vy1_start;
	double _vx2_start;
	double _vy2_start;
	double _vx3_start;
	double _vy3_start;

	vector<double> _c_start;

	// calculation setting
	int _alg;
	double _dt;

	// physics condition
	static double _GM;
	static double _t_end;
	static double _m1;
	static double _m2;
	static double _m3;
	double _a1;
	double _a2;
	double _a3;
	double _e1;
	double _e2;
	double _e3;

	// output
	vector<double> _t;
	vector< vector<double> > _x; 
	int _n_stps;
	double _period;

	static double f_x1(double t, const vector<double> &x);
	static double f_y1(double t, const vector<double> &x);
	static double f_x2(double t, const vector<double> &x);
	static double f_y2(double t, const vector<double> &x);
	static double f_x3(double t, const vector<double> &x);
	static double f_y3(double t, const vector<double> &x);
	static double f_vx1(double t, const vector<double> &x);
	static double f_vy1(double t, const vector<double> &x);
	static double f_vx2(double t, const vector<double> &x);
	static double f_vy2(double t, const vector<double> &x);
	static double f_vx3(double t, const vector<double> &x);
	static double f_vy3(double t, const vector<double> &x);

	static bool stop(double t, const vector<double> &x);

public:
	ThreeBody2D();
	ThreeBody2D(double a1, double e1, double m1, double a2, double e2, double m2, double a3, double e3, double m3);
	~ThreeBody2D();

	void set_t_start(double t_start);
	void set_x1_start(double x1_start);
	void set_y1_start(double y1_start);
	void set_x2_start(double x2_start);
	void set_y2_start(double y2_start);
	void set_x3_start(double x3_start);
	void set_y3_start(double y3_start);
	void set_vx1_start(double vx1_start);
	void set_vy1_start(double vy1_start);
	void set_vx2_start(double vx2_start);
	void set_vy2_start(double vy2_start);
	void set_vx3_start(double vx3_start);
	void set_vy3_start(double vy3_start);
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
	double get_x3_start() const;
	double get_y3_start() const;
	double get_vx1_start() const;
	double get_vy1_start() const;
	double get_vx2_start() const;
	double get_vy2_start() const;
	double get_vx3_start() const;
	double get_vy3_start() const;
	int get_alg() const;
	double get_dt() const;

	static double get_GM();
	static double get_t_end();
	static double get_m1();
	static double get_m2();
	static double get_m3();

	bool check() const;
	void cal();
	double cal_period();

	vector<double> get_t() const;
	vector< vector<double> > get_x() const;
	int get_n_stps() const;
	double get_period() const;

};


inline void ThreeBody2D::set_t_start(double t_start) { _t_start = t_start; }
inline void ThreeBody2D::set_x1_start(double x1_start) { _x1_start = x1_start; _c_start[0] = x1_start; }
inline void ThreeBody2D::set_y1_start(double y1_start) { _y1_start = y1_start; _c_start[1] = y1_start; }
inline void ThreeBody2D::set_x2_start(double x2_start) { _x2_start = x2_start; _c_start[2] = x2_start; }
inline void ThreeBody2D::set_y2_start(double y2_start) { _y2_start = y2_start; _c_start[3] = y2_start; }
inline void ThreeBody2D::set_x3_start(double x3_start) { _x3_start = x3_start; _c_start[4] = x3_start; }
inline void ThreeBody2D::set_y3_start(double y3_start) { _y3_start = y3_start; _c_start[5] = y3_start; }
inline void ThreeBody2D::set_vx1_start(double vx1_start) { _vx1_start = vx1_start; _c_start[6] = vx1_start; }
inline void ThreeBody2D::set_vy1_start(double vy1_start) { _vy1_start = vy1_start; _c_start[7] = vy1_start; }
inline void ThreeBody2D::set_vx2_start(double vx2_start) { _vx2_start = vx2_start; _c_start[8] = vx2_start; }
inline void ThreeBody2D::set_vy2_start(double vy2_start) { _vy2_start = vy2_start; _c_start[9] = vy2_start; }
inline void ThreeBody2D::set_vx3_start(double vx3_start) { _vx3_start = vx3_start; _c_start[10] = vx3_start; }
inline void ThreeBody2D::set_vy3_start(double vy3_start) { _vy3_start = vy3_start; _c_start[11] = vy3_start; }
inline void ThreeBody2D::set_alg(int alg) { _alg = alg; }
inline void ThreeBody2D::set_dt(double dt) { _dt = dt; }

inline void ThreeBody2D::set_GM(double GM) { _GM = GM; }
inline void ThreeBody2D::set_t_end(double t_end) { _t_end = t_end; }

inline double ThreeBody2D::get_t_start() const { return _t_start; }
inline double ThreeBody2D::get_x1_start() const { return _x1_start; }
inline double ThreeBody2D::get_y1_start() const { return _y1_start; }
inline double ThreeBody2D::get_x2_start() const { return _x2_start; }
inline double ThreeBody2D::get_y2_start() const { return _y2_start; }
inline double ThreeBody2D::get_x3_start() const { return _x3_start; }
inline double ThreeBody2D::get_y3_start() const { return _y3_start; }
inline double ThreeBody2D::get_vx1_start() const { return _vx1_start; }
inline double ThreeBody2D::get_vy1_start() const { return _vy1_start; }
inline double ThreeBody2D::get_vx2_start() const { return _vx2_start; }
inline double ThreeBody2D::get_vy2_start() const { return _vy2_start; }
inline double ThreeBody2D::get_vx3_start() const { return _vx3_start; }
inline double ThreeBody2D::get_vy3_start() const { return _vy3_start; }
inline int ThreeBody2D::get_alg() const { return _alg; }
inline double ThreeBody2D::get_dt() const { return _dt; }

inline double ThreeBody2D::get_GM() { return _GM; }
inline double ThreeBody2D::get_t_end() { return _t_end; }

inline vector<double> ThreeBody2D::get_t() const { return _t; }
inline vector< vector<double> > ThreeBody2D:: get_x() const { return _x; }
inline int ThreeBody2D::get_n_stps() const { return _n_stps; }
inline double ThreeBody2D::get_period() const { return _period; }


#endif
