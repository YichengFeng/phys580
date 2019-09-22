#ifndef ANHARMONICS_H
#define ANHARMONICS_H

#include <iostream>
#include <vector>
#include <cmath>
#include "runge_kutta.h"

using namespace std;


class Anharmonics
{
private:
	const int _n_eqns = 2;

	// inputs
	// initial condition
	double _t_start;
	double _theta_start;
	double _omega_start;

	vector<double> _x_start;

	// parameters for the ODE
	static double _alpha;
	static double _k;	
	static double _t_end;

	// calculation setting
	int _alg; // 0: Euler Cromer; 1: Euler; 2: RK2; 4: RK4;
	double _dt;
	double _n_stps;

	// output
	vector<double> _t;
	vector<double> _theta;
	vector<double> _omega;
	vector< vector<double> > _x;

	// ODEs
	static double f_theta(double t, const vector<double> &x);
	static double f_omega(double t, const vector<double> &x);

	static bool stop(double t, const vector<double> &x);

public:
	Anharmonics();
	~Anharmonics();

	void set_t_start(double t_start);
	void set_theta_start(double theta_start);
	void set_omega_start(double omega_start);
	void set_alg(int alg);
	void set_dt(double dt);

	static void set_alpha(double alpha);
	static void set_k(double k);
	static void set_t_end(double t_end);

	double get_t_start() const;
	double get_theta_start() const;
	double get_omega_start() const;
	double get_dt() const;
	int get_n_stps() const;

	static double get_alpha();
	static double get_k();
	static double get_t_end();

	bool check() const;
	void cal();
	double cal_period() const;

	vector<double> get_t() const;
	vector< vector<double> > get_x() const;
};


inline void Anharmonics::set_t_start(double t_start) { _t_start = t_start; }
inline void Anharmonics::set_theta_start(double theta_start) { _theta_start = theta_start; _x_start[0] = theta_start; }
inline void Anharmonics::set_omega_start(double omega_start) { _omega_start = omega_start; _x_start[1] = omega_start; }
inline void Anharmonics::set_alg(int alg) { _alg = alg; }
inline void Anharmonics::set_dt(double dt) { _dt = dt; }

inline void Anharmonics::set_alpha(double alpha) { _alpha = alpha; }
inline void Anharmonics::set_k(double k) { _k = k; }
inline void Anharmonics::set_t_end(double t_end) { _t_end = t_end; }

inline double Anharmonics::get_t_start() const { return _t_start; }
inline double Anharmonics::get_theta_start() const { return _theta_start; }
inline double Anharmonics::get_omega_start() const { return _omega_start; }
inline double Anharmonics::get_dt() const { return _dt; }
inline int Anharmonics::get_n_stps() const { return _n_stps; }

inline double Anharmonics::get_alpha() { return _alpha; }
inline double Anharmonics::get_k() { return _k; }
inline double Anharmonics::get_t_end() { return _t_end; }

inline vector<double> Anharmonics::get_t() const { return _t; }
inline vector< vector<double> > Anharmonics::get_x() const { return _x; }


#endif
