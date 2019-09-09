#ifndef RUNGEKUTTA_H
#define RUNGEKUTTA_H

#include <iostream>
#include <vector>

using namespace std;


class RungeKutta
{
private:
	// inputs
	int _n_eqns;
	double _dt;
	double _t_start;
	vector<double> _x_start;
	vector<double (*)(double t, const vector<double> &x)> _fv;
	bool (*_stop)(double t, const vector<double> &x);

	// middle
	int _n_stps;

	// output
	vector<double> _t;
	vector< vector<double> > _x;

public:
	RungeKutta();
	RungeKutta(int n_eqns, double dt, double t_start, const vector<double> &x_start);
	RungeKutta(int n_eqns, double dt, double t_start, const vector<double> &x_start, vector<double (*)(double t, const vector<double> &x)> fv, bool (*stop)(double t, const vector<double> &x));

	static bool def_stop(double t, const vector<double> &x);

	void set_n_eqns(int n_eqns);
	void set_dt(double dt);
	void set_t_start(double t_start);
	void set_x_start(const vector<double> &x_start);
	void set_fv(vector<double (*)(double t, const vector<double> &x)> fv);
	void load_f(double (*f)(double t, const vector<double> &x));
	void set_stop(bool (*stop)(double t, const vector<double> &x));

	bool check() const;

	void cal_rk1();
	void cal_rk2();
	void cal_rk4();

	int get_n_stps() const;
	vector<double> get_t() const;
	vector< vector<double> > get_x() const;
};


inline void RungeKutta::set_n_eqns(int n_eqns) { _n_eqns = n_eqns; }
inline void RungeKutta::set_dt(double dt) { _dt = dt; }
inline void RungeKutta::set_t_start(double t_start) { _t_start = t_start; }
inline void RungeKutta::set_x_start(const vector<double> &x_start) { for(int i=0; i<x_start.size(); i++) { _x_start.push_back(x_start[i]); } }
inline void RungeKutta::set_fv(vector<double (*)(double t, const vector<double> &x)> fv) { for(int i=0; i<fv.size(); i++) { _fv.push_back(fv[i]); } }
inline void RungeKutta::load_f(double (*f)(double t, const vector<double> &x)) { _fv.push_back(f); }
inline void RungeKutta::set_stop(bool (*stop)(double t, const vector<double> &x)) { _stop = stop; }

inline int RungeKutta::get_n_stps() const { return _n_stps; }
inline vector<double> RungeKutta::get_t() const { return _t; }
inline vector< vector<double> > RungeKutta::get_x() const { return _x; }


#endif
