#include <iostream>
#include <vector>
#include <cmath>
#include "qm_pow4.h"

using namespace std;


//-------------------------------------------------------------------------//

QMPow4::QMPow4() {

	QMPow4(0.01, 0, 0.1);
}

//-------------------------------------------------------------------------//

QMPow4::QMPow4(double dx, double E, double dE) {

	_dx = dx;
	_E = E;
	_dE = dE;

	_x_match = 0.7;
	_x_left = -3;
	_x_right = 3;
	_y_min = -2;
	_y_max = 2;
	_r = 1.0;//1.0/(3.0*3.0*3.0*3.0);
	_devdiff = 999;

	_matchL = round((_x_match-_x_left)/_dx);
	_matchR = round((_x_right-_x_match)/_dx);

	_fr = 0.01;
	_T = 0.1;
	_n_total = 0;

	_iend = round((_x_right-_x_left)/_dx) + 3;
	for(int i=0; i<_iend; i++) {
		_x.push_back((i-1)*_dx + _x_left);
		_psi.push_back(0);
		_psi_old.push_back(0);
	}
	double delta = 2.0;
	int i_delta = round((+1 - _x_left)/_dx) + 1;
	int i_start = round((-1 - _x_left)/_dx) + 1;
	for(int i=i_start; i<=i_delta; i++) {
		_psi[i] = 1.0/sqrt(delta);
	}

	_dpsi = _fr/sqrt(delta);
	var_psi_normalize();
	_varE = var_energy();
}

//-------------------------------------------------------------------------//

void QMPow4::cal_potential_table() {

	_potential_table.clear();
	for(int i=0; i<_x.size(); i++) {
		double V = potential(_x[i]);
		if(isnan(V) || isinf(V)) V = 0;
		_potential_table.push_back(V*_r);
	}
}

//-------------------------------------------------------------------------//

bool QMPow4::check() const {

	if(_x_left>_x_right) {
		cout << "ERROR: range invalid!" << endl;
		return false;
	}

	return true;
}

//-------------------------------------------------------------------------//

double QMPow4::potential(double x) {

	return x*x*x*x;
}

//-------------------------------------------------------------------------//

void QMPow4::cal_match_once() {

	if(!check()) return;

	_xL.clear();
	_xR.clear();
	_psiL.clear();
	_psiR.clear();
	_psi.clear();

	_psiL.push_back(0);
	_psiL.push_back(1e-2*_dx);
	_psiR.push_back(0);
	_psiR.push_back(1e-2*_dx);

	_xL.push_back(_x_left-_dx);
	_xL.push_back(_x_left);
	_xR.push_back(_x_right+_dx);
	_xR.push_back(_x_right);

	double psimatch;

	for(int i=1; i<_matchL+16; i++) {
		double x = _x_left + (i-1)*_dx;
		double tmpL = 2*_psiL[i] - _psiL[i-1] - 2*(_E - potential(x))*_dx*_dx*_psiL[i];
		_psiL.push_back(tmpL);
		_xL.push_back(x+_dx);
	}
	psimatch = _psiL[_matchL+1];
	for(int i=0; i<_psiL.size(); i++) {
		_psiL[i] /= psimatch;
	}

	for(int i=1; i<_matchR+16; i++) {
		double x = _x_right - (i-1)*_dx;
		double tmpR = 2*_psiR[i] - _psiR[i-1] - 2*(_E - potential(x))*_dx*_dx*_psiR[i];
		_psiR.push_back(tmpR);
		_xR.push_back(x-_dx);
	}
	psimatch = _psiR[_matchR+1];
	for(int i=0; i<_psiR.size(); i++) {
		_psiR[i] /= psimatch;
	}
}

//-------------------------------------------------------------------------//

void QMPow4::adjust_match_once() {

	cal_match_once();

	double devL = (_psiL[_matchL+2] - _psiL[_matchL])/(2*_dx);
	double devR = (_psiR[_matchR] - _psiR[_matchR+2])/(2*_dx);
	double tmpdiff = devL - devR;

	if(tmpdiff==0) {
		_dE = 0;
	} else if(tmpdiff*_devdiff<0) {
		_dE = -0.5*_dE;
	} else if(tmpdiff>0 && tmpdiff>_devdiff) {
		_dE = -_dE;
	} else if(tmpdiff<0 && tmpdiff<_devdiff) {
		_dE = -_dE;
	} else {
		_dE = _dE;
	}

	_E += _dE;
	_devdiff = tmpdiff;
}

//-------------------------------------------------------------------------//

void QMPow4::adjust_match_until(double ee) {

	while(fabs(_dE)>ee) {
		adjust_match_once();
	}
}

//-------------------------------------------------------------------------//

bool QMPow4::cal_var_once() {

	int n = 1 + floor(1.0*rand()/RAND_MAX*(_iend-2));
	int m = n + floor(1.0*rand()/RAND_MAX*(_iend-1-n));
	double p = 1.0*rand()/RAND_MAX;
	for(int i=n; i<=m; i++) {
		_psi_old[i] = _psi[i];
		_psi[i] += 2*(p-0.5)*_dpsi;
	}
	_newE = var_energy();
	double deltaE = _newE - _varE;
	bool keep;
	if(_newE >= _varE) {
		keep = false;
		if(_T > 0) {
			p = 1.0*rand()/RAND_MAX;
			if(p <= exp(-deltaE/_T)) keep = true;
		}
	} else {
		keep = true;
	}

	if(keep) {
		_varE = _newE;
		var_psi_normalize();
	} else {
		for(int i=n; i<=m; i++) {
			_psi[i] = _psi_old[i];
		}
	}

	return keep;
}

//-------------------------------------------------------------------------//

int QMPow4::cal_var_attempts(int n_attempts) {

	int n_move = 0;
	for(int i_attempts=0; i_attempts<n_attempts; i_attempts++) {
		if(cal_var_once()) n_move ++;
	}
	_n_total += n_attempts;
	_n_moves += n_move;

	return n_move;
}

//-------------------------------------------------------------------------//

void QMPow4::var_psi_normalize() {

	double sum = 0;
	for(int i=0; i<_iend; i++) {
		sum += _dx*_psi[i]*_psi[i];
	}
	for(int i=0; i<_iend; i++) {
		_psi[i] /= sum;
	}
}

//-------------------------------------------------------------------------//

double QMPow4::var_energy() {

	double energy = 0;
	double sum = 0;
	for(int i=1; i<_iend-1; i++) {
		double x = _x_left + (i-1)*_dx;
		double V = potential(x);
		energy += _dx*V*_psi[i]*_psi[i] - 0.5*_psi[i]*(_psi[i+1]+_psi[i-1]-2*_psi[i])/_dx;
		sum += _psi[i]*_psi[i]*_dx;
	}
	energy /= sum;
	return energy;
}

//-------------------------------------------------------------------------//
