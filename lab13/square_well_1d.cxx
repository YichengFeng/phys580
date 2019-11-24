#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include "square_well_1d.h"

using namespace std;


//-------------------------------------------------------------------------//

SquareWell1D::SquareWell1D() {

	SquareWell1D(1, 0.01, 1e3, 0, 0, 2, 1.1, 0.02);
}

//-------------------------------------------------------------------------//

SquareWell1D::SquareWell1D(int parity, double dx, double Vmax, double a, double V2, double b, double E, double dE) {

	_parity = parity;
	_dx = dx;
	_Vmax = Vmax;
	_E = E;
	_dE = dE;
	_a = a;
	_V2 = V2;
	_b = b;

	_psiAt1 = 999;
	_iter = -1;

	_xmin = -_L*1.3;
	_xmax = +_L*1.3;
	_ymin = -b;
	_ymax = +b;

	_r = _ymax/max(fabs(_Vmax), fabs(_V2));
	_nhalf = _xmax/_dx;
}

//-------------------------------------------------------------------------//

bool SquareWell1D::check() const {

	if(_parity!=1 && _parity!=-1) {
		cout << "ERROR: parity invalid!" << endl;
		return false;
	}
	if(_a<0) {
		cout << "ERROR: barrier edge invalid!" << endl;
		return false;
	}
	if(_b<0) {
		cout << "ERROR: vertical scale invalid!" << endl;
		return false;
	}
	if(_L<0) {
		cout << "ERROR: well edge invalid!" << endl;
		return false;
	}

	return true;
}

//-------------------------------------------------------------------------//

double SquareWell1D::potential(double x) {

	if(fabs(x)<=_a) {
		return _V2;
	} else if(fabs(x)<_L) {
		return 0;
	} else {
		return _Vmax;
	}
}

//-------------------------------------------------------------------------//

void SquareWell1D::cal_once() {

	if(!check()) return;

	_xL.clear();
	_xR.clear();
	_psiL.clear();
	_psiR.clear();
	_x.clear();
	_psi.clear();

	if(_parity==1) {
		if(_a!=0) {
			_psiR.push_back(0.1);
			_psiR.push_back(0.1);
			_psiL.push_back(0.1);
			_psiL.push_back(0.1);
		} else {
			_psiR.push_back(1);
			_psiR.push_back(1);
			_psiL.push_back(1);
			_psiL.push_back(1);
		}
	} else {
		_psiR.push_back(-_dx);
		_psiR.push_back(0);
		_psiL.push_back(-_dx);
		_psiL.push_back(0);
	}
	_xR.push_back(-_dx);
	_xR.push_back(0);
	_xL.push_back(_dx);
	_xL.push_back(0);

	int i = 0;
	while(true) {
		i++;
		double x = (i-1)*_dx;
		double tmpR = 2*_psiR[i] - _psiR[i-1] - 2*(_E - potential(x))*_dx*_dx*_psiR[i];
		_xR.push_back(i*_dx);
		_xL.push_back(-i*_dx);
		_psiR.push_back(tmpR);
		_psiL.push_back(_parity*tmpR);

		if(fabs(_psiR[i+1])>_b) break;
	}

	int n = _psiL.size();
	for(int i=0; i<n; i++) {
		_x.push_back(_xL[n-1-i]);
		_psi.push_back(_psiL[n-1-i]);
	}
	for(int i=2; i<n; i++) {
		_x.push_back(_xR[i]);
		_psi.push_back(_psiR[i]);
	}
}

//-------------------------------------------------------------------------//

void SquareWell1D::adjust_once() {

	cal_once();
	double psiAt1Now = _psiR[floor(_L/_dx+1)];
	if(psiAt1Now==0) {
		_dE = 0;
	} else if(psiAt1Now*_psiAt1<0) {
		_dE = -0.5*_dE;
	} else if(psiAt1Now>0 && psiAt1Now>_psiAt1) {
		_dE = -_dE;
	} else if(psiAt1Now<0 && psiAt1Now<_psiAt1) {
		_dE = -_dE;
	} else {
		_dE = _dE;
	}

	_E += _dE;
	_psiAt1 = psiAt1Now;
}

//-------------------------------------------------------------------------//

void SquareWell1D::adjust_until(double ee) {

	int n = 0;
	while(fabs(_dE)>ee) {
		adjust_once();
		if(n == _iter) break;
		n ++;
	}
}

//-------------------------------------------------------------------------//

void SquareWell1D::cal_potential_table() {

	_potential_table.clear();
	for(int i=0; i<_x.size(); i++) {
		_potential_table.push_back(potential(_x[i])*_r);
	}
}

//-------------------------------------------------------------------------//

