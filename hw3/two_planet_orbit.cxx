#include <iostream>
#include <vector>
#include <cmath>

#include "two_planet_orbit.h"

using namespace std;


//-------------------------------------------------------------------------//

double TwoPlanetOrbit::_GM = 4*M_PI*M_PI;
double TwoPlanetOrbit::_t_end = 5;
double TwoPlanetOrbit::_m1 = 3.003e-6;
double TwoPlanetOrbit::_m2 = 9.546e-4;

//-------------------------------------------------------------------------//

TwoPlanetOrbit::TwoPlanetOrbit() {
	// planet 1: Earth
	// planet 2: Moon

	_a1 = 1;
	_e1 = 0;//0.017;
	_m1 = 5.97e24/1.99e30;
	_a2 = 1+2.570e-3;
	_e2 = 0.;
	_m2 = 3.69397e-8;;

	_GM = 4*M_PI*M_PI;
	_t_end = 2;

	_t_start = 0;
	_x1_start = _a1*(1-_e1);
	_y1_start = 0;
	_x2_start = _a2*(1-_e2);
	_y2_start = 0;
	_vx1_start = 0;
	_vy1_start = sqrt(_GM*(1+_e1)/(1-_e1)/_a1);
	_vx2_start = 0;
	_vy2_start = _vy1_start*_a2/_a1+2*M_PI/29.53*365.2422*2.570e-3;
	_c_start.push_back(_x1_start);
	_c_start.push_back(_y1_start);
	_c_start.push_back(_x2_start);
	_c_start.push_back(_y2_start);
	_c_start.push_back(_vx1_start);
	_c_start.push_back(_vy1_start);
	_c_start.push_back(_vx2_start);
	_c_start.push_back(_vy2_start);

	_alg = 4;
	_dt = 0.0001;
}

//-------------------------------------------------------------------------//

TwoPlanetOrbit::TwoPlanetOrbit(double a1, double e1, double m1, double a2, double e2, double m2) {
	// planet 1: Earth
	// planet 2: Jupiter

	_a1 = a1;
	_e1 = e1;
	_m1 = m1;
	_a2 = a2;
	_e2 = e2;
	_m2 = m2;

	_GM = 4*M_PI*M_PI;
	_t_end = 5;

	_t_start = 0;
	_x1_start = _a1*(1-_e1);
	_y1_start = 0;
	_x2_start = _a2*(1-_e2);
	_y2_start = 0;
	_vx1_start = 0;
	_vy1_start = sqrt(_GM*(1+_e1)/(1-_e1)/_a1*(1+_m1));
	_vx2_start = 0;
	_vy2_start = sqrt(_GM*(1+_e2)/(1-_e2)/_a2*(1+_m2));
	_c_start.push_back(_x1_start);
	_c_start.push_back(_y1_start);
	_c_start.push_back(_x2_start);
	_c_start.push_back(_y2_start);
	_c_start.push_back(_vx1_start);
	_c_start.push_back(_vy1_start);
	_c_start.push_back(_vx2_start);
	_c_start.push_back(_vy2_start);

	_alg = 0;
	_dt = 0.001;
}

//-------------------------------------------------------------------------//

TwoPlanetOrbit::~TwoPlanetOrbit() {

}

//-------------------------------------------------------------------------//

double TwoPlanetOrbit::f_x1(double t, const vector<double> &x) {
	return x[4];
}

//-------------------------------------------------------------------------//

double TwoPlanetOrbit::f_y1(double t, const vector<double> &x) {
	return x[5];
}

//-------------------------------------------------------------------------//

double TwoPlanetOrbit::f_x2(double t, const vector<double> &x) {
	return x[6];
}

//-------------------------------------------------------------------------//

double TwoPlanetOrbit::f_y2(double t, const vector<double> &x) {
	return x[7];
}

//-------------------------------------------------------------------------//

double TwoPlanetOrbit::f_vx1(double t, const vector<double> &x) {
	double rEJ = sqrt((x[0]-x[2])*(x[0]-x[2])+(x[1]-x[3])*(x[1]-x[3]));
	return -_GM*x[0]/pow(x[0]*x[0]+x[1]*x[1], 1.5) - _GM*_m2*(x[0]-x[2])/rEJ/rEJ/rEJ;
}

//-------------------------------------------------------------------------//

double TwoPlanetOrbit::f_vy1(double t, const vector<double> &x) {
	double rEJ = sqrt((x[0]-x[2])*(x[0]-x[2])+(x[1]-x[3])*(x[1]-x[3]));
	return -_GM*x[1]/pow(x[0]*x[0]+x[1]*x[1], 1.5) - _GM*_m2*(x[1]-x[3])/rEJ/rEJ/rEJ;
}

//-------------------------------------------------------------------------//

double TwoPlanetOrbit::f_vx2(double t, const vector<double> &x) {
	double rEJ = sqrt((x[0]-x[2])*(x[0]-x[2])+(x[1]-x[3])*(x[1]-x[3]));
	return -_GM*x[2]/pow(x[2]*x[2]+x[3]*x[3], 1.5) - _GM*_m1*(x[2]-x[0])/rEJ/rEJ/rEJ;
}

//-------------------------------------------------------------------------//

double TwoPlanetOrbit::f_vy2(double t, const vector<double> &x) {
	double rEJ = sqrt((x[0]-x[2])*(x[0]-x[2])+(x[1]-x[3])*(x[1]-x[3]));
	return -_GM*x[3]/pow(x[2]*x[2]+x[3]*x[3], 1.5) - _GM*_m1*(x[3]-x[1])/rEJ/rEJ/rEJ;
}

//-------------------------------------------------------------------------//

bool TwoPlanetOrbit::stop(double t, const vector<double> &x) {
	return t>_t_end;
}

//-------------------------------------------------------------------------//

bool TwoPlanetOrbit::check() const {
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

void TwoPlanetOrbit::cal() {
	if(!check()) return;

	RungeKutta rk(_n_eqns, _dt, _t_start, _c_start);
	rk.load_f(f_x1);
	rk.load_f(f_y1);
	rk.load_f(f_x2);
	rk.load_f(f_y2);
	rk.load_f(f_vx1);
	rk.load_f(f_vy1);
	rk.load_f(f_vx2);
	rk.load_f(f_vy2);
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

double TwoPlanetOrbit::cal_period() {
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

