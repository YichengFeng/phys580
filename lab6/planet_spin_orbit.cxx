#include <iostream>
#include <vector>
#include <cmath>

#include "planet_spin_orbit.h"

using namespace std;


//-------------------------------------------------------------------------//

double PlanetSpinOrbit::_GM = 4*M_PI*M_PI;
double PlanetSpinOrbit::_t_end = 5;

//-------------------------------------------------------------------------//

PlanetSpinOrbit::PlanetSpinOrbit() {
	_a = 1/(1-0.123);
	_e = 0;
	_mp = 5.97e24/1.99e30;
	_GM = 4*M_PI*M_PI;
	_t_end = 5;

	_t_start = 0;
	_x_start = _a*(1-_e);
	_y_start = 0;
	_theta_start = 0;
	_vx_start = 0;
	//_vy_start = sqrt(_GM*(1+_e)/(1-_e)/_a*(1+_mp));
	_vy_start = sqrt(_GM*(1+_e)/(1-_e)/_a);
	_omega_start = 0;
	_c_start.push_back(_x_start);
	_c_start.push_back(_y_start);
	_c_start.push_back(_theta_start);
	_c_start.push_back(_vx_start);
	_c_start.push_back(_vy_start);
	_c_start.push_back(_omega_start);

	_energy = 0.5*_mp*(_vy_start*_vy_start + _vx_start*_vx_start);
	_energy += -_GM*_mp/sqrt(_x_start*_x_start + _y_start*_y_start);

	_alg = 0;
	_dt = 0.001;
}

//-------------------------------------------------------------------------//

PlanetSpinOrbit::PlanetSpinOrbit(double a, double e, double mp) {
	_a = a;
	_e = e;
	_mp = mp;
	_GM = 4*M_PI*M_PI;
	_t_end = 5;

	_t_start = 0;
	_x_start = _a*(1-_e);
	_y_start = 0;
	_theta_start = 0;
	_vx_start = 0;
	//_vy_start = sqrt(_GM*(1+_e)/(1-_e)/_a*(1+_mp));
	_vy_start = sqrt(_GM*(1+_e)/(1-_e)/_a);
	_omega_start = 0;
	_c_start.push_back(_x_start);
	_c_start.push_back(_y_start);
	_c_start.push_back(_theta_start);
	_c_start.push_back(_vx_start);
	_c_start.push_back(_vy_start);
	_c_start.push_back(_omega_start);

	_energy = 0.5*_mp*(_vy_start*_vy_start + _vx_start*_vx_start);
	_energy += -_GM*_mp/sqrt(_x_start*_x_start + _y_start*_y_start);

	_alg = 0;
	_dt = 0.001;
}

//-------------------------------------------------------------------------//

PlanetSpinOrbit::~PlanetSpinOrbit() {

}

//-------------------------------------------------------------------------//

void PlanetSpinOrbit::increase_v_start(double r) {
	_t_end = 1;
	_vy_start = sqrt(2*_GM/_x_start)*r;
	_c_start[3] = _vy_start*r;
}

//-------------------------------------------------------------------------//

double PlanetSpinOrbit::f_x(double t, const vector<double> &x) {
	return x[3];
}

//-------------------------------------------------------------------------//

double PlanetSpinOrbit::f_y(double t, const vector<double> &x) {
	return x[4];
}

//-------------------------------------------------------------------------//

double PlanetSpinOrbit::f_theta(double t, const vector<double> &x) {
	return x[5];
}

//-------------------------------------------------------------------------//

double PlanetSpinOrbit::f_vx(double t, const vector<double> &x) {
	return -_GM*x[0]/pow(x[0]*x[0]+x[1]*x[1], 1.5);
}

//-------------------------------------------------------------------------//

double PlanetSpinOrbit::f_vy(double t, const vector<double> &x) {
	return -_GM*x[1]/pow(x[0]*x[0]+x[1]*x[1], 1.5);
}

//-------------------------------------------------------------------------//

double PlanetSpinOrbit::f_omega(double t, const vector<double> &x) {
	return -3*_GM/pow(x[0]*x[0]+x[1]*x[1], 2.5)*(x[0]*cos(x[2])+x[1]*sin(x[2]))*(x[0]*sin(x[2])-x[1]*cos(x[2]));
}

//-------------------------------------------------------------------------//

bool PlanetSpinOrbit::stop(double t, const vector<double> &x) {
	return t>_t_end;
}

//-------------------------------------------------------------------------//

bool PlanetSpinOrbit::check() const {
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

void PlanetSpinOrbit::cal() {
	if(!check()) return;

	RungeKutta rk(_n_eqns, _dt, _t_start, _c_start);
	rk.load_f(f_x);
	rk.load_f(f_y);
	rk.load_f(f_theta);
	rk.load_f(f_vx);
	rk.load_f(f_vy);
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
}

//-------------------------------------------------------------------------//

double PlanetSpinOrbit::cal_period() {
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

