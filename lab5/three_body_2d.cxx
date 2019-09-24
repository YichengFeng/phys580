#include <iostream>
#include <vector>
#include <cmath>

#include "three_body_2d.h"

using namespace std;


//-------------------------------------------------------------------------//

double ThreeBody2D::_GM = 4*M_PI*M_PI;
double ThreeBody2D::_t_end = 5;
double ThreeBody2D::_m1 = 3.003e-6;
double ThreeBody2D::_m2 = 9.546e-4;
double ThreeBody2D::_m3 = 1;

//-------------------------------------------------------------------------//

ThreeBody2D::ThreeBody2D() {
	// body 1: Earth
	// body 2: Jupiter
	// body 3: Sun

	_a1 = 1;
	_e1 = 0.017;
	_m1 = 5.97e24/1.99e30;
	_a2 = 5.203;
	_e2 = 0.048;
	_m2 = 1.90e27/1.99e30;
	_a3 = 0;
	_e3 = 0;
	_m3 = 1;

	_GM = 4*M_PI*M_PI;
	_t_end = 5;

	_t_start = 0;
	_x1_start = _a1*(1-_e1);
	_y1_start = 0;
	_x2_start = _a2*(1-_e2);
	_y2_start = 0;
	_x3_start = 0;
	_y3_start = 0;
	_vx1_start = 0;
	_vy1_start = sqrt(_GM*(1+_e1)/(1-_e1)/_a1);
	_vx2_start = 0;
	_vy2_start = sqrt(_GM*(1+_e2)/(1-_e2)/_a2);
	_vx3_start = 0;
	_vy3_start = 0;

	// automatically go to the frame of mass center
	double cx = (_m1*_x1_start + _m2*_x2_start + _m3*_x3_start)/(_m1 + _m2 + _m3);
	double cy = (_m1*_y1_start + _m2*_y2_start + _m3*_y3_start)/(_m1 + _m2 + _m3);
	double cvx = (_m1*_vx1_start + _m2*_vx2_start + _m3*_vx3_start)/(_m1 + _m2 + _m3);
	double cvy = (_m1*_vy1_start + _m2*_vy2_start + _m3*_vy3_start)/(_m1 + _m2 + _m3);

	_x1_start -= cx;
	_y1_start -= cy;
	_x2_start -= cx;
	_y2_start -= cy;
	_x3_start -= cx;
	_y3_start -= cy;
	_vx1_start -= cvx; 
	_vy1_start -= cvy;
	_vx2_start -= cvx;
	_vy2_start -= cvy;
	_vx3_start -= cvx;
	_vy3_start -= cvy;

	_c_start.push_back(_x1_start);
	_c_start.push_back(_y1_start);
	_c_start.push_back(_x2_start);
	_c_start.push_back(_y2_start);
	_c_start.push_back(_x3_start);
	_c_start.push_back(_y3_start);
	_c_start.push_back(_vx1_start);
	_c_start.push_back(_vy1_start);
	_c_start.push_back(_vx2_start);
	_c_start.push_back(_vy2_start);
	_c_start.push_back(_vx3_start);
	_c_start.push_back(_vy3_start);

	_alg = 0;
	_dt = 0.001;
}

//-------------------------------------------------------------------------//

ThreeBody2D::ThreeBody2D(double a1, double e1, double m1, double a2, double e2, double m2, double a3, double e3, double m3) {
	// body 1: Earth
	// body 2: Jupiter
	// body 3: Sun

	_a1 = a1;
	_e1 = e1;
	_m1 = m1;
	_a2 = a2;
	_e2 = e2;
	_m2 = m2;
	_a3 = a3;
	_e3 = e3;
	_m3 = m3;

	_GM = 4*M_PI*M_PI;
	_t_end = 5;

	_t_start = 0;
	_x1_start = _a1*(1-_e1);
	_y1_start = 0;
	_x2_start = _a2*(1-_e2);
	_y2_start = 0;
	_x3_start = 0;
	_y3_start = 0;
	_vx1_start = 0;
	_vy1_start = sqrt(_GM*(1+_e1)/(1-_e1)/_a1);
	_vx2_start = 0;
	_vy2_start = sqrt(_GM*(1+_e2)/(1-_e2)/_a2);
	_vx3_start = 0;
	_vy3_start = 0;

	// automatically go to the frame of mass center
	double cx = (_m1*_x1_start + _m2*_x2_start + _m3*_x3_start)/(_m1 + _m2 + _m3);
	double cy = (_m1*_y1_start + _m2*_y2_start + _m3*_y3_start)/(_m1 + _m2 + _m3);
	double cvx = (_m1*_vx1_start + _m2*_vx2_start + _m3*_vx3_start)/(_m1 + _m2 + _m3);
	double cvy = (_m1*_vy1_start + _m2*_vy2_start + _m3*_vy3_start)/(_m1 + _m2 + _m3);

	_x1_start -= cx;
	_y1_start -= cy;
	_x2_start -= cx;
	_y2_start -= cy;
	_x3_start -= cx;
	_y3_start -= cy;
	_vx1_start -= cvx; 
	_vy1_start -= cvy;
	_vx2_start -= cvx;
	_vy2_start -= cvy;
	_vx3_start -= cvx;
	_vy3_start -= cvy;

	_c_start.push_back(_x1_start);
	_c_start.push_back(_y1_start);
	_c_start.push_back(_x2_start);
	_c_start.push_back(_y2_start);
	_c_start.push_back(_x3_start);
	_c_start.push_back(_y3_start);
	_c_start.push_back(_vx1_start);
	_c_start.push_back(_vy1_start);
	_c_start.push_back(_vx2_start);
	_c_start.push_back(_vy2_start);
	_c_start.push_back(_vx3_start);
	_c_start.push_back(_vy3_start);

	_alg = 0;
	_dt = 0.001;
}

//-------------------------------------------------------------------------//

ThreeBody2D::~ThreeBody2D() {

}

//-------------------------------------------------------------------------//

double ThreeBody2D::f_x1(double t, const vector<double> &x) {
	return x[6];
}

//-------------------------------------------------------------------------//

double ThreeBody2D::f_y1(double t, const vector<double> &x) {
	return x[7];
}

//-------------------------------------------------------------------------//

double ThreeBody2D::f_x2(double t, const vector<double> &x) {
	return x[8];
}

//-------------------------------------------------------------------------//

double ThreeBody2D::f_y2(double t, const vector<double> &x) {
	return x[9];
}

//-------------------------------------------------------------------------//

double ThreeBody2D::f_x3(double t, const vector<double> &x) {
	return x[10];
}

//-------------------------------------------------------------------------//

double ThreeBody2D::f_y3(double t, const vector<double> &x) {
	return x[11];
}

//-------------------------------------------------------------------------//

double ThreeBody2D::f_vx1(double t, const vector<double> &x) {
	double rEJ = sqrt((x[0]-x[2])*(x[0]-x[2])+(x[1]-x[3])*(x[1]-x[3]));
	double rES = sqrt((x[0]-x[4])*(x[0]-x[4])+(x[1]-x[5])*(x[1]-x[5]));
	double rJS = sqrt((x[2]-x[4])*(x[2]-x[4])+(x[3]-x[5])*(x[3]-x[5]));
	return -_GM*_m2*(x[0]-x[2])/rEJ/rEJ/rEJ -_GM*_m3*(x[0]-x[4])/rES/rES/rES;
}

//-------------------------------------------------------------------------//

double ThreeBody2D::f_vy1(double t, const vector<double> &x) {
	double rEJ = sqrt((x[0]-x[2])*(x[0]-x[2])+(x[1]-x[3])*(x[1]-x[3]));
	double rES = sqrt((x[0]-x[4])*(x[0]-x[4])+(x[1]-x[5])*(x[1]-x[5]));
	double rJS = sqrt((x[2]-x[4])*(x[2]-x[4])+(x[3]-x[5])*(x[3]-x[5]));
	return -_GM*_m2*(x[1]-x[3])/rEJ/rEJ/rEJ -_GM*_m3*(x[1]-x[5])/rES/rES/rES;
}

//-------------------------------------------------------------------------//

double ThreeBody2D::f_vx2(double t, const vector<double> &x) {
	double rEJ = sqrt((x[0]-x[2])*(x[0]-x[2])+(x[1]-x[3])*(x[1]-x[3]));
	double rES = sqrt((x[0]-x[4])*(x[0]-x[4])+(x[1]-x[5])*(x[1]-x[5]));
	double rJS = sqrt((x[2]-x[4])*(x[2]-x[4])+(x[3]-x[5])*(x[3]-x[5]));
	return -_GM*_m3*(x[2]-x[4])/rJS/rJS/rJS -_GM*_m1*(x[2]-x[0])/rEJ/rEJ/rEJ;
}

//-------------------------------------------------------------------------//

double ThreeBody2D::f_vy2(double t, const vector<double> &x) {
	double rEJ = sqrt((x[0]-x[2])*(x[0]-x[2])+(x[1]-x[3])*(x[1]-x[3]));
	double rES = sqrt((x[0]-x[4])*(x[0]-x[4])+(x[1]-x[5])*(x[1]-x[5]));
	double rJS = sqrt((x[2]-x[4])*(x[2]-x[4])+(x[3]-x[5])*(x[3]-x[5]));
	return -_GM*_m3*(x[3]-x[5])/rJS/rJS/rJS -_GM*_m1*(x[3]-x[1])/rEJ/rEJ/rEJ;
}

//-------------------------------------------------------------------------//

double ThreeBody2D::f_vx3(double t, const vector<double> &x) {
	double rEJ = sqrt((x[0]-x[2])*(x[0]-x[2])+(x[1]-x[3])*(x[1]-x[3]));
	double rES = sqrt((x[0]-x[4])*(x[0]-x[4])+(x[1]-x[5])*(x[1]-x[5]));
	double rJS = sqrt((x[2]-x[4])*(x[2]-x[4])+(x[3]-x[5])*(x[3]-x[5]));
	return -_GM*_m1*(x[4]-x[0])/rES/rES/rES -_GM*_m2*(x[4]-x[2])/rJS/rJS/rJS;
}

//-------------------------------------------------------------------------//

double ThreeBody2D::f_vy3(double t, const vector<double> &x) {
	double rEJ = sqrt((x[0]-x[2])*(x[0]-x[2])+(x[1]-x[3])*(x[1]-x[3]));
	double rES = sqrt((x[0]-x[4])*(x[0]-x[4])+(x[1]-x[5])*(x[1]-x[5]));
	double rJS = sqrt((x[2]-x[4])*(x[2]-x[4])+(x[3]-x[5])*(x[3]-x[5]));
	return -_GM*_m1*(x[5]-x[1])/rES/rES/rES -_GM*_m2*(x[5]-x[3])/rJS/rJS/rJS;
}

//-------------------------------------------------------------------------//

bool ThreeBody2D::stop(double t, const vector<double> &x) {
	return t>_t_end;
}

//-------------------------------------------------------------------------//

bool ThreeBody2D::check() const {
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

void ThreeBody2D::cal() {
	if(!check()) return;

	RungeKutta rk(_n_eqns, _dt, _t_start, _c_start);
	rk.load_f(f_x1);
	rk.load_f(f_y1);
	rk.load_f(f_x2);
	rk.load_f(f_y2);
	rk.load_f(f_x3);
	rk.load_f(f_y3);
	rk.load_f(f_vx1);
	rk.load_f(f_vy1);
	rk.load_f(f_vx2);
	rk.load_f(f_vy2);
	rk.load_f(f_vx3);
	rk.load_f(f_vy3);
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

double ThreeBody2D::cal_period() {
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

