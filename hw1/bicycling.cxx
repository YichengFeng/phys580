#include <iostream>
#include <vector>
#include "bicycling.h"

using namespace std;


//-------------------------------------------------------------------------//

double Bicycling::_C = 0.5;
double Bicycling::_P = 400;
double Bicycling::_A = 0.33;
double Bicycling::_m = 70;
double Bicycling::_rho = 1.292;

//-------------------------------------------------------------------------//

Bicycling::Bicycling() {
	_C = 0.5;
	_P = 400;
	_A = 0.33;
	_m = 70;
	_rho = 1.292;

	_dt = 0.1;
	_t_start = 0;
	_v_start = 4;
	_x_start.push_back(_v_start);
}

//-------------------------------------------------------------------------//

Bicycling::~Bicycling() {

}

//-------------------------------------------------------------------------//

double Bicycling::f_v(double t, const vector<double> &x) {
	return _P/_m/x[0] - 0.5/_m*_C*_rho*_A*x[0]*x[0];
}

//-------------------------------------------------------------------------//

bool Bicycling::stop(double t, const vector<double> &x) {
	return t>100;
}

//-------------------------------------------------------------------------//

bool Bicycling::check() const {
	if(_n_eqns != _x_start.size()) {
		cout << "ERROR: equation number not match!" << endl;
		return false;
	}

	return true;
}

//-------------------------------------------------------------------------//

void Bicycling::cal() {
	if(!check()) {
		cout << "ERROR: check() not passed, no calculation!" << endl;
		return;
	}

	RungeKutta rk(_n_eqns, _dt, _t_start, _x_start);
	rk.set_stop(stop);
	rk.load_f(f_v);

	rk.cal_rk1();

	_t = rk.get_t();
	_x = rk.get_x();

	_n_stps = rk.get_n_stps();
}

