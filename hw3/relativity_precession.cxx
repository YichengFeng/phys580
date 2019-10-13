#include <iostream>
#include <vector>
#include <cmath>

#include "relativity_precession.h"

using namespace std;


//-------------------------------------------------------------------------//

double RelativityPrecession::_GM = 4*M_PI*M_PI;
double RelativityPrecession::_t_end = 5;
double RelativityPrecession::_alpha = 0;

//-------------------------------------------------------------------------//

RelativityPrecession::RelativityPrecession() {
	_a = 0.38665;
	_e = 0.206;
	_mp = 3.30e23/1.99e30;
	_GM = 4*M_PI*M_PI;
	_t_end = 5;
	_alpha = 0;

	_t_start = 0;
	_x_start = _a*(1-_e);
	_y_start = 0;
	_vx_start = 0;
	_vy_start = sqrt(_GM*(1+_e)/(1-_e)/_a);
	_c_start.push_back(_x_start);
	_c_start.push_back(_y_start);
	_c_start.push_back(_vx_start);
	_c_start.push_back(_vy_start);

	_energy = 0.5*_mp*(_vy_start*_vy_start + _vx_start*_vx_start);
	_energy += -_GM*_mp/sqrt(_x_start*_x_start + _y_start*_y_start);

	_alg = 0;
	_dt = 0.001;
}

//-------------------------------------------------------------------------//

RelativityPrecession::RelativityPrecession(double a, double e, double mp) {
	_a = a;
	_e = e;
	_mp = mp;
	_GM = 4*M_PI*M_PI;
	_t_end = 5;
	_alpha = 0;

	_t_start = 0;
	_x_start = _a*(1-_e);
	_y_start = 0;
	_vx_start = 0;
	_vy_start = sqrt(_GM*(1+_e)/(1-_e)/_a);
	_c_start.push_back(_x_start);
	_c_start.push_back(_y_start);
	_c_start.push_back(_vx_start);
	_c_start.push_back(_vy_start);

	_energy = 0.5*_mp*(_vy_start*_vy_start + _vx_start*_vx_start);
	_energy += -_GM*_mp/sqrt(_x_start*_x_start + _y_start*_y_start);

	_alg = 0;
	_dt = 0.001;
}

//-------------------------------------------------------------------------//

RelativityPrecession::~RelativityPrecession() {

}

//-------------------------------------------------------------------------//

void RelativityPrecession::increase_v_start(double r) {
	_t_end = 1;
	_vy_start = sqrt(2*_GM/_x_start)*r;
	_c_start[3] = _vy_start*r;
}

//-------------------------------------------------------------------------//

double RelativityPrecession::f_x(double t, const vector<double> &x) {
	return x[2];
}

//-------------------------------------------------------------------------//

double RelativityPrecession::f_y(double t, const vector<double> &x) {
	return x[3];
}

//-------------------------------------------------------------------------//

double RelativityPrecession::f_vx(double t, const vector<double> &x) {
	double r = sqrt(x[0]*x[0]+x[1]*x[1]);
	return -_GM*x[0]/(r*r*r)*(1+_alpha/(r*r));
}

//-------------------------------------------------------------------------//

double RelativityPrecession::f_vy(double t, const vector<double> &x) {
	double r = sqrt(x[0]*x[0]+x[1]*x[1]);
	return -_GM*x[1]/(r*r*r)*(1+_alpha/(r*r));
}

//-------------------------------------------------------------------------//

bool RelativityPrecession::stop(double t, const vector<double> &x) {
	return t>_t_end;
}

//-------------------------------------------------------------------------//

bool RelativityPrecession::check() const {
	// size of ODE set
	if(_c_start.size() != _n_eqns) {
		cout << "ERROR: equation number not match!" << endl;
		return false;
	}
	// _alg: 0: Euler-Cromer; 1: Euler; 2: RK2; 4: RK4
	if(_alg!=0 && _alg!=1 && _alg!=2 && _alg!=4) {
		cout << "ERROR: numerical method selection fails!" << endl;
		return false;
	}

	return true;
}

//-------------------------------------------------------------------------//

void RelativityPrecession::cal() {
	if(!check()) return;

	RungeKutta rk(_n_eqns, _dt, _t_start, _c_start);
	rk.load_f(f_x);
	rk.load_f(f_y);
	rk.load_f(f_vx);
	rk.load_f(f_vy);
	rk.set_stop(stop);

	// _alg: 0: Euler-Cromer; 1: Euler; 2: RK2; 4: RK4
	if(_alg == 0) {
		rk.cal_Euler_Cromer();
	} else if(_alg == 1) {
		rk.cal_rk1();
	} else if(_alg == 2) {
		rk.cal_rk2();
	} else if(_alg == 4) {
		rk.cal_rk4();
	} else {
		cout << "ERROR: numerical method selection fails!" << endl;
	}

	_t = rk.get_t();
	_x = rk.get_x();
	_n_stps = rk.get_n_stps();
}

//-------------------------------------------------------------------------//

double RelativityPrecession::cal_period() {
	double period = -1;
	vector<double> t0s;
	double t0;
	for(int i=0; i<_n_stps-1; i++) {
		if(_x[0][i]*_x[0][i+1] > 0) continue;
		t0 = (_x[0][i+1]*_t[i]-_x[0][i]*_t[i+1])/(_x[0][i+1]-_x[0][i]);
		t0s.push_back(t0);
	}

	if(t0s.size()>=3) {
		period = t0s[2]-t0s[0];
	}

	_period = period;

	return period;
}

//-------------------------------------------------------------------------//

void RelativityPrecession::cal_perihelion() {
	for(int i=0; i<_n_stps; i++) {
		_rr.push_back(_x[0][i]*_x[0][i]+_x[1][i]*_x[1][i]);
	}
	_perihelion_t.push_back(_t[0]);
	_perihelion_x.push_back(_x[0][0]);
	_perihelion_y.push_back(_x[1][0]);
	_perihelion_theta.push_back(0);

	for(int i=1; i<_n_stps-1; i++) {
		if(_rr[i]<=_rr[i-1] && _rr[i]<_rr[i+1]) {
			_perihelion_t.push_back(_t[i]);
			_perihelion_x.push_back(_x[0][i]);
			_perihelion_y.push_back(_x[1][i]);
			double tmp_theta = atan2(_x[1][i], _x[0][i]);
			_perihelion_theta.push_back(tmp_theta);
		}
	}
}

//-------------------------------------------------------------------------//

