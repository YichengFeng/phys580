#include <iostream>
#include <vector>
#include <cmath>
#include "pendulum.h"

using namespace std;


//-------------------------------------------------------------------------//

double Pendulum::_g = 9.8;
double Pendulum::_l = 9.8;
double Pendulum::_q = 0;
double Pendulum::_F = 0;
double Pendulum::_O = 0;
double Pendulum::_P = 0;
double Pendulum::_n_periods = 2;

//-------------------------------------------------------------------------//

Pendulum::Pendulum() {
	_dt = 0.02;
	_mode = 0;
	_alg = 1;
	_t_start = 0;
	_theta_start = 0;
	_omega_start = 0;
	_x_start.push_back(_theta_start);
	_x_start.push_back(_omega_start);

	_g = 9.8;
	_l = 9.8;
	_q = 0;
	_F = 0;
	_O = 0;
	_P = 0;
	_n_periods = 2;
}

//-------------------------------------------------------------------------//

Pendulum::~Pendulum() {

}

//-------------------------------------------------------------------------//

// x[0]: theta; x[1]: omega
double Pendulum::f_theta(double t, const vector<double> &x) {
	// dtheta/dt = omega
	return x[1];
}

//-------------------------------------------------------------------------//

double Pendulum::f_omega_simplified(double t, const vector<double> &x) {
	return -_g/_l*x[0] - _q*x[1] + _F*sin(_O*t+_P);
}

//-------------------------------------------------------------------------//

double Pendulum::f_energy_simplified(double t, const vector<double> &x) {
	// 2E/m/l^2
	return _g/_l*x[0]*x[0]+x[1]*x[1];
}

//-------------------------------------------------------------------------//

double Pendulum::f_omega_physics(double t, const vector<double> &x) {
	return -_g/_l*sin(x[0]) - _q*x[1] + _F*sin(_O*t+_P);
}

//-------------------------------------------------------------------------//

double Pendulum::f_energy_physics(double t, const vector<double> &x) {
	// 2E/m/l^2
	return _g/_l*(1-cos(x[0]))+x[1]*x[1];
}

//-------------------------------------------------------------------------//

// the exact solution of the simplified ODE
double Pendulum::f_theta_simplified_exact(double t, const vector<double> &x_start) {
	double Omg = sqrt(_g/_l);
	double Amp = sqrt(x_start[0]*x_start[0]+x_start[1]*x_start[1]/Omg/Omg);
	double Phi = atan2(x_start[0], x_start[1]/Omg);
	return Amp*sin(Omg*t+Phi);
}

//-------------------------------------------------------------------------//

// the exact solution of the simplified ODE
double Pendulum::f_omega_simplified_exact(double t, const vector<double> &x_start) {
	double Omg = sqrt(_g/_l);
	double Amp = sqrt(x_start[0]*x_start[0]+x_start[1]*x_start[1]/Omg/Omg);
	double Phi = atan2(x_start[0], x_start[1]/Omg);
	return Amp*Omg*cos(Omg*t+Phi);
}

//-------------------------------------------------------------------------//

bool Pendulum::stop(double t, const vector<double> &x) {
	return t > _n_periods*2*M_PI;
}

//-------------------------------------------------------------------------//

bool Pendulum::check() const {
	if(_x_start.size() != _n_eqns) {
		cout << "ERROR: equation number not match!" << endl;
		return false;
	}
	// _mode: 0: simplified; 1: physics
	if(_mode!=0 && _mode!=1) {
		cout << "ERROR: ODE set not specified!" << endl;
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

void Pendulum::cal() {
	if(!check()) return;

	RungeKutta rk(_n_eqns, _dt, _t_start, _x_start);
	rk.load_f(f_theta);
	// _mode: 0: simplified; 1: physics
	double (*f_omega)(double t, const vector<double> &x);
	double (*f_energy)(double t, const vector<double> &x);
	if(_mode == 0) {
		f_omega = f_omega_simplified;
		f_energy = f_energy_simplified;
	} else if(_mode == 1) {
		f_omega = f_omega_physics;
		f_energy = f_energy_physics;
	} else {
		cout << "ERROR: ODE set not specified!" << endl;
	}
	rk.load_f(f_omega);
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

	_energy.clear();
	_energy_exact.clear();
	_energy_error.clear();
	for(int i=0; i<_n_stps; i++) {
		double t = _t[i];
		vector<double> x;
		for(int j=0; j<_n_eqns; j++) {
			x.push_back(_x[j][i]);
		}
		double energy = f_energy(t, x);
		_energy.push_back(energy);
		_energy_exact.push_back(_energy[0]);
		_energy_error.push_back(energy-_energy[0]);
	}

	_x_simplified_exact.clear();
	_x_simplified_error.clear();
	vector<double> theta_simplified_exact;
	vector<double> theta_simplified_error;
	vector<double> omega_simplified_exact;
	vector<double> omega_simplified_error;
	if(_mode == 0) {
		for(int i=0; i<_n_stps; i++) {
			double tmp_theta = f_theta_simplified_exact(_t[i], _x_start);
			theta_simplified_exact.push_back(tmp_theta);
			theta_simplified_error.push_back(_x[0][i]-tmp_theta);
			double tmp_omega = f_omega_simplified_exact(_t[i], _x_start);
			omega_simplified_exact.push_back(tmp_omega);
			omega_simplified_error.push_back(_x[1][i]-tmp_omega);
		}
		_x_simplified_exact.push_back(theta_simplified_exact);
		_x_simplified_exact.push_back(omega_simplified_exact);
		_x_simplified_error.push_back(theta_simplified_error);
		_x_simplified_error.push_back(omega_simplified_error);
	}
}

//-------------------------------------------------------------------------//

double Pendulum::cal_period() const {
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

	return period;
}

//-------------------------------------------------------------------------//

double Pendulum::cal_analytical_period() const {
	double period = -1;
	double energy = _energy[0];
	double theta_max = acos(1-energy*_l/_g);
	if(_mode == 0) {
		period = 2*M_PI*sqrt(_l/_g);
	} else if(_mode == 1) {
		period = 2*M_PI*sqrt(_l/_g);
		period *= 1 + pow(theta_max,2)/16 + pow(theta_max,4)*11/3072;
		//period *= 1 + pow(theta_max,2)/16 + pow(theta_max,4)*11/3072 + pow(theta_max,6)*173/737280;
	}

	return period;
}

