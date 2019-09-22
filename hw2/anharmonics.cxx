#include <iostream>
#include <vector>
#include "anharmonics.h"

using namespace std;


//-------------------------------------------------------------------------//

double Anharmonics::_alpha = 1;
double Anharmonics::_k = 1;	
double Anharmonics::_t_end = 10;	

//-------------------------------------------------------------------------//

Anharmonics::Anharmonics() {
	_t_start = 0;
	_theta_start = 0.2;
	_omega_start = 0;
	_alg = 1; // 0: Euler Cromer; 1: Euler; 2: RK2; 4: RK4;
	_dt = 0.02;

	_x_start.push_back(_theta_start);
	_x_start.push_back(_omega_start);

	_alpha = 1;
	_k = 1;
	_t_end = 10;
}

//-------------------------------------------------------------------------//

Anharmonics::~Anharmonics() {

}

//-------------------------------------------------------------------------//

double Anharmonics::f_theta(double t, const vector<double> &x) {
	return x[1];
}

//-------------------------------------------------------------------------//

double Anharmonics::f_omega(double t, const vector<double> &x) {
	int sign = 0;
	if(x[0]>0) sign = 1;
	if(x[0]<0) sign = -1;
	return -_k*sign*pow(fabs(x[0]), _alpha);
}

//-------------------------------------------------------------------------//

bool Anharmonics::stop(double t, const vector<double> &x) {
	return t>_t_end;
}

//-------------------------------------------------------------------------//

bool Anharmonics::check() const {
	
	if(_x_start.size() != _n_eqns) {
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

void Anharmonics::cal() {

	if(!check()) return;

	RungeKutta rk(_n_eqns, _dt, _t_start, _x_start);
	rk.load_f(f_theta);
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

double Anharmonics::cal_period() const {
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

