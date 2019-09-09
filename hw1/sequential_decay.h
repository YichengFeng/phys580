#ifndef SEQUENTIAL_DECAY_H
#define SEQUENTIAL_DECAY_H

#include <iostream>
#include <vector>

#include "runge_kutta.h"

using namespace std;


class SequentialDecay
{
private:
	// inputs
	const int _n_eqns = 2;
	double _dt;
	double _t_start;
	double _N_A_start;
	double _N_B_start;
	vector<double> _x_start;

	// physics parameter
	static double _tau_A_B;

	// middle
	int _n_stps;
	int _n_stps_exact;

	// output
	vector<double> _t;
	vector< vector<double> > _x;
	vector<double> _N_A;
	vector<double> _N_B;

	vector<double> _t_exact;
	vector< vector<double> > _x_exact;
	vector<double> _N_A_exact;
	vector<double> _N_B_exact;

	vector< vector<double> > _x_error;
	vector<double> _N_A_error;
	vector<double> _N_B_error;

	static bool stop(double t, const vector<double> &x);
	static double f_N_A(double t, const vector<double> &x);
	static double f_N_B(double t, const vector<double> &x);

	double f_N_A_exact(double t);
	double f_N_B_exact(double t);

public:
	SequentialDecay();
	~SequentialDecay();

	void set_t_start(double t_start);
	void set_N_A_start(double N_A_start);
	void set_N_B_start(double N_B_start);
	void set_dt(double dt);

	double get_t_start() const;
	double get_N_A_start() const;
	double get_N_B_start() const;
	double get_dt() const;

	static void set_tau_A_B(double tau_A_B);
	static double get_tau_A_B();

	bool check() const;
	void cal();
	void cal_exact();
	void cal_error();

	vector<double> get_t() const;
	vector< vector<double> > get_x() const;
	vector<double> get_N_A() const;
	vector<double> get_N_B() const;
	int get_n_stps() const;

	vector<double> get_t_exact() const;
	vector< vector<double> > get_x_exact() const;
	vector<double> get_N_A_exact() const;
	vector<double> get_N_B_exact() const;
	int get_n_stps_exact() const;

	vector< vector<double> > get_x_error() const;
	vector<double> get_N_A_error() const;
	vector<double> get_N_B_error() const;
};


inline void SequentialDecay::set_t_start(double t_start) { _t_start = t_start; }
inline void SequentialDecay::set_N_A_start(double N_A_start) { _N_A_start = N_A_start; _x_start[0] = _N_A_start; }
inline void SequentialDecay::set_N_B_start(double N_B_start) { _N_B_start = N_B_start; _x_start[1] = _N_B_start; }
inline void SequentialDecay::set_dt(double dt) { _dt = dt; }

inline double SequentialDecay::get_t_start() const { return _t_start; }
inline double SequentialDecay::get_N_A_start() const { return _N_A_start; }
inline double SequentialDecay::get_N_B_start() const { return _N_B_start; }
inline double SequentialDecay::get_dt() const { return _dt; }

inline void SequentialDecay::set_tau_A_B(double tau_A_B) { _tau_A_B = tau_A_B; }
inline double SequentialDecay::get_tau_A_B() { return _tau_A_B; }

inline vector<double> SequentialDecay::get_t() const { return _t; }
inline vector< vector<double> > SequentialDecay::get_x() const { return _x; }
inline vector<double> SequentialDecay::get_N_A() const { return _N_A; }
inline vector<double> SequentialDecay::get_N_B() const { return _N_B; }
inline int SequentialDecay::get_n_stps() const { return _n_stps; }

inline vector<double> SequentialDecay::get_t_exact() const { return _t_exact; }
inline vector< vector<double> > SequentialDecay::get_x_exact() const { return _x_exact; }
inline vector<double> SequentialDecay::get_N_A_exact() const { return _N_A_exact; }
inline vector<double> SequentialDecay::get_N_B_exact() const { return _N_B_exact; }
inline int SequentialDecay::get_n_stps_exact() const { return _n_stps_exact; }

inline vector< vector<double> > SequentialDecay::get_x_error() const { return _x_error; }
inline vector<double> SequentialDecay::get_N_A_error() const { return _N_A_error; }
inline vector<double> SequentialDecay::get_N_B_error() const { return _N_B_error; }


#endif	
