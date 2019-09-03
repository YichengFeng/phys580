#include <iostream>
#include <vector>

#include "runge_kutta.h"

using namespace std;


//-------------------------------------------------------------------------//

bool RungeKutta::def_stop(double t, const vector<double> &x) {
	return true;
}

//-------------------------------------------------------------------------//

RungeKutta::RungeKutta() {
	_n_eqns = 0;
	_dt = 1;
	_t_start = 0;
	_stop = def_stop;

	_n_stps = 0;
}

//-------------------------------------------------------------------------//

RungeKutta::RungeKutta(int n_eqns, double dt, double t_start, const vector<double> &x_start) {
	_n_eqns = n_eqns;
	_dt = dt;
	_t_start = t_start;
	for(int i=0; i<x_start.size(); i++) {
		_x_start.push_back(x_start[i]);
		vector<double> xtmp;
		_x.push_back(xtmp);
	}
	_stop = def_stop;

	_n_stps = 0;
}

//-------------------------------------------------------------------------//

RungeKutta::RungeKutta(int n_eqns, double dt, double t_start, const vector<double> &x_start, vector<double (*)(double t, const vector<double> &x)> fv, bool (*stop)(double t, const vector<double> &x)) {
	_n_eqns = n_eqns;
	_dt = dt;
	_t_start = t_start;
	for(int i=0; i<x_start.size(); i++) {
		_x_start.push_back(x_start[i]);
		vector<double> xtmp;
		_x.push_back(xtmp);
	}
	for(int i=0; i<fv.size(); i++) {
		_fv.push_back(fv[i]);
	}
	_stop = stop;

	if(n_eqns != fv.size()) {
		cout << "ERROR: equation number not match!" << endl;
	}
	if(n_eqns != x_start.size()) {
		cout << "ERROR: initial condition number not match!" << endl;
	}

	_n_stps = 0;
}

//-------------------------------------------------------------------------//

bool RungeKutta::check() const {
	if(_n_eqns == 0) {
		cout << "WARNING: zero input equation!" << endl;
		return false;
	}
	if(_n_eqns != _fv.size()) {
		cout << "ERROR: equation number not match!" << endl;
		return false;
	}
	if(_n_eqns != _x_start.size()) {
		cout << "ERROR: initial condition number not match!" << endl;
		return false;
	}
	for(int i=0; i<_fv.size(); i++) {
		if(_fv[i] == nullptr) {
			cout << "ERROR: invalid input function!" << endl;
			return false;
		}
	}

	return true;
}

//-------------------------------------------------------------------------//

void RungeKutta::cal_rk1() {
	if(!check()) {
		cout << "WARNING: check() not passed, no calculation!" << endl;
		return;
	}
	//cout << "check() passed" << endl;

	_n_stps = 0;

	double t = _t_start;
	vector<double> x = _x_start;

	vector<double> F1(_n_eqns);

	vector<double> x1(_n_eqns);

	double t1;

	int nfv = _fv.size(); // nfv = _n_eqns
	while(!_stop(t, x)) {
		// Fill into output
		_t.push_back(t);
		for(int ifv=0; ifv<nfv; ifv++) {
			_x[ifv].push_back(x[ifv]);
		}
		_n_stps ++;

		// F1
		t1 = t;
		for(int ifv=0; ifv<nfv; ifv++) {
			x1[ifv] = x[ifv];
		}
		for(int ifv=0; ifv<nfv; ifv++) {
			F1[ifv] = _fv[ifv](t1, x1);
		}

		t += _dt;
		for(int ifv=0; ifv<nfv; ifv++) {
			x[ifv] += F1[ifv]*_dt;
		}
	}

	_t.push_back(t);
	for(int ifv=0; ifv<nfv; ifv++) {
		_x[ifv].push_back(x[ifv]);
	}
	_n_stps ++;
}

//-------------------------------------------------------------------------//

void RungeKutta::cal_rk2() {
	if(!check()) {
		cout << "WARNING: check() not passed, no calculation!" << endl;
		return;
	}
	//cout << "check() passed" << endl;

	_n_stps = 0;

	double t = _t_start;
	vector<double> x = _x_start;

	vector<double> F1(_n_eqns);
	vector<double> F2(_n_eqns);

	vector<double> x1(_n_eqns);
	vector<double> x2(_n_eqns);

	double t1, t2;

	int nfv = _fv.size(); // nfv = _n_eqns
	while(!_stop(t, x)) {
		// Fill into output
		_t.push_back(t);
		for(int ifv=0; ifv<nfv; ifv++) {
			_x[ifv].push_back(x[ifv]);
		}
		_n_stps ++;

		// F1
		t1 = t;
		for(int ifv=0; ifv<nfv; ifv++) {
			x1[ifv] = x[ifv];
		}
		for(int ifv=0; ifv<nfv; ifv++) {
			F1[ifv] = _fv[ifv](t1, x1);
		}
		// F2
		t2 = t + 0.5*_dt;
		for(int ifv=0; ifv<nfv; ifv++) {
			x2[ifv] = x[ifv] + F1[ifv]*0.5*_dt;
		}
		for(int ifv=0; ifv<nfv; ifv++) {
			F2[ifv] = _fv[ifv](t2, x2);
		}

		t += _dt;
		for(int ifv=0; ifv<nfv; ifv++) {
			x[ifv] += F2[ifv]*_dt;
		}
	}

	_t.push_back(t);
	for(int ifv=0; ifv<nfv; ifv++) {
		_x[ifv].push_back(x[ifv]);
	}
	_n_stps ++;
}

//-------------------------------------------------------------------------//

void RungeKutta::cal_rk4() {
	if(!check()) {
		cout << "WARNING: check() not passed, no calculation!" << endl;
		return;
	}
	//cout << "check() passed" << endl;

	_n_stps = 0;

	double t = _t_start;
	vector<double> x = _x_start;

	vector<double> F1(_n_eqns);
	vector<double> F2(_n_eqns);
	vector<double> F3(_n_eqns);
	vector<double> F4(_n_eqns);

	vector<double> x1(_n_eqns);
	vector<double> x2(_n_eqns);
	vector<double> x3(_n_eqns);
	vector<double> x4(_n_eqns);

	double t1, t2, t3, t4;

	int nfv = _fv.size(); // nfv = _n_eqns
	while(!_stop(t, x)) {
		// Fill into output
		_t.push_back(t);
		for(int ifv=0; ifv<nfv; ifv++) {
			_x[ifv].push_back(x[ifv]);
		}
		_n_stps ++;

		// F1
		t1 = t;
		for(int ifv=0; ifv<nfv; ifv++) {
			x1[ifv] = x[ifv];
		}
		for(int ifv=0; ifv<nfv; ifv++) {
			F1[ifv] = _fv[ifv](t1, x1);
		}
		// F2
		t2 = t + 0.5*_dt;
		for(int ifv=0; ifv<nfv; ifv++) {
			x2[ifv] = x[ifv] + F1[ifv]*0.5*_dt;
		}
		for(int ifv=0; ifv<nfv; ifv++) {
			F2[ifv] = _fv[ifv](t2, x2);
		}
		// F3
		t3 = t + 0.5*_dt;
		for(int ifv=0; ifv<nfv; ifv++) {
			x3[ifv] = x[ifv] + F2[ifv]*0.5*_dt;
		}
		for(int ifv=0; ifv<nfv; ifv++) {
			F3[ifv] = _fv[ifv](t3, x3);
		}
		// F4
		t4 = t + _dt;
		for(int ifv=0; ifv<nfv; ifv++) {
			x4[ifv] = x[ifv] + F3[ifv]*_dt;
		}
		for(int ifv=0; ifv<nfv; ifv++) {
			F4[ifv] = _fv[ifv](t4, x4);
		}

		t += _dt;
		for(int ifv=0; ifv<nfv; ifv++) {
			x[ifv] += (F1[ifv]/6 + F2[ifv]/3 + F3[ifv]/3 + F4[ifv]/6)*_dt;
		}
	}

	_t.push_back(t);
	for(int ifv=0; ifv<nfv; ifv++) {
		_x[ifv].push_back(x[ifv]);
	}
	_n_stps ++;
}

