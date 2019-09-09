#include <iostream>
#include <cmath>
#include "sequential_decay.h"

using namespace std;


//-------------------------------------------------------------------------//

double SequentialDecay::_tau_A_B = 1;

//-------------------------------------------------------------------------//

bool SequentialDecay::stop(double t, const vector<double> &x) {
	return t>7;
}

//-------------------------------------------------------------------------//

double SequentialDecay::f_N_A(double t, const vector<double> &x) {
	return -x[0]/_tau_A_B;
}

//-------------------------------------------------------------------------//

double SequentialDecay::f_N_B(double t, const vector<double> &x) {
	return x[0]/_tau_A_B-x[1];
}

//-------------------------------------------------------------------------//

double SequentialDecay::f_N_A_exact(double t) {
	return _N_A_start*exp(-t/_tau_A_B);
}

//-------------------------------------------------------------------------//

double SequentialDecay::f_N_B_exact(double t) {
	if(_tau_A_B == 1) {
		return (_N_B_start + _N_A_start*t/_tau_A_B)*exp(-t);
	} else {
		return _N_A_start/(1-_tau_A_B)*(exp(-t)-exp(-t/_tau_A_B))+_N_B_start*exp(-t);
	}
}

//-------------------------------------------------------------------------//

SequentialDecay::SequentialDecay() {
	_dt = 0.02;
	_t_start = 0;
	_N_A_start = 1000;
	_N_B_start = 0;
	_x_start.push_back(_N_A_start);
	_x_start.push_back(_N_B_start);
}

//-------------------------------------------------------------------------//

SequentialDecay::~SequentialDecay() {
	
}

//-------------------------------------------------------------------------//

bool SequentialDecay::check() const {
	if(_x_start.size() != _n_eqns) {
		cout << "ERROR: equation number not match!" << endl;
		return false;
	}

	return true;
}

//-------------------------------------------------------------------------//

void SequentialDecay::cal() {
	if(!check()) {
		cout << "ERROR: check() not passed! no calculation!" << endl;
	}

	RungeKutta rk(_n_eqns, _dt, _t_start, _x_start);
	rk.set_stop(stop);
	rk.load_f(f_N_A);
	rk.load_f(f_N_B);

	rk.cal_rk1();

	_t = rk.get_t();
	_x = rk.get_x();

	_N_A = _x[0];
	_N_B = _x[1];

	_n_stps = rk.get_n_stps();
}

//-------------------------------------------------------------------------//

void SequentialDecay::cal_exact() {
	double t = _t_start;
	vector<double> x;
	x.push_back(_N_A_start);
	x.push_back(_N_B_start);

	_n_stps_exact = 0;

	_t_exact.clear();
	_N_A_exact.clear();
	_N_B_exact.clear();
	_x_exact.clear();

	while(!stop(t, x)) {
		_t_exact.push_back(t);
		_N_A_exact.push_back(x[0]);
		_N_B_exact.push_back(x[1]);
		_n_stps_exact ++;

		t += _dt;
		x[0] = f_N_A_exact(t);
		x[1] = f_N_B_exact(t);
	}

	_t_exact.push_back(t);
	_N_A_exact.push_back(x[0]);
	_N_B_exact.push_back(x[1]);
	_n_stps_exact ++;

	_x_exact.push_back(_N_A_exact);
	_x_exact.push_back(_N_B_exact);
}

//-------------------------------------------------------------------------//

void SequentialDecay::cal_error() {
	if(_N_A.size() != _N_A_exact.size() || _N_B.size() != _N_B_exact.size()) {
		cout << "ERROR: cal_error(): vector size not matched!" << endl;
		cout << _N_A.size() << " != " << _N_A_exact.size() << endl;
		return;
	}

	for(int i=0; i<_N_A.size(); i++) {
		if( _t[i] != _t_exact[i]) {
			cout << "ERROR: time not match!" << endl;
		}
		//double error = fabs(_N_A[i] - _N_A_exact[i])/_N_A_exact[i];
		double error = (_N_A[i] - _N_A_exact[i])/_N_A_exact[i];
		_N_A_error.push_back(error);
	}

	for(int i=0; i<_N_B.size(); i++) {
		//double error = fabs(_N_B[i] - _N_B_exact[i])/_N_B_exact[i];
		double error = (_N_B[i] - _N_B_exact[i])/_N_B_exact[i];
		_N_B_error.push_back(error);
	}

	_x_error.push_back(_N_A_error);
	_x_error.push_back(_N_B_error);
}
