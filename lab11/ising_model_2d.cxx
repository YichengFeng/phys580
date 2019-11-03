#include <iostream>
#include <vector>
#include <cmath>
#include "ising_model_2d.h"

using namespace std;


//-------------------------------------------------------------------------//

IsingModel2D::IsingModel2D() {

	_N = 10;
	_T = 2;
	_H = 0;
	_MCS = 20;
	_c = 2; // (1: all up; -1: all down; 2: random)

	_now = 0;

	_spin.clear();
	for(int i=0; i<_N; i++) {
		vector<int> tmp;
		for(int j=0; j<_N; j++) {
			int one_spin;
			if(_c == 1) {
				one_spin = 1;
			} else if(_c == -1) {
				one_spin = -1;
			} else {
				one_spin = 1.0*rand()/RAND_MAX<0.5?-1:1;
			}
			tmp.push_back(one_spin);
		}
		_spin.push_back(tmp);
	}

	_m = 0;
	_E = 0;
	for(int i=0; i<_N; i++) {
		int il = i-1<  0? i-1+_N:i-1;
		int ir = i+1>=_N? i+1-_N:i+1;
		for(int j=0; j<_N; j++) {
			int jl = j-1<  0? j-1+_N:j-1;
			int jr = j+1>=_N? j+1-_N:j+1;

			_m += _spin[i][j];
			double sum = _spin[il][j] + _spin[ir][j] + _spin[i][jl] + _spin[i][jr];
			_E -= _spin[i][j]*(_H + 0.5*sum);
		}
	}
}

//-------------------------------------------------------------------------//

IsingModel2D::IsingModel2D(int N, double T, double H, int MCS, int c) {

	_N = N;
	_T = T;
	_H = H;
	_MCS = MCS;
	_c = c; // (1: all up; -1: all down; 2: random)

	_now = 0;

	_spin.clear();
	for(int i=0; i<_N; i++) {
		vector<int> tmp;
		for(int j=0; j<_N; j++) {
			int one_spin;
			if(_c == 1) {
				one_spin = 1;
			} else if(_c == -1) {
				one_spin = -1;
			} else {
				one_spin = 1.0*rand()/RAND_MAX<0.5?-1:1;
			}
			tmp.push_back(one_spin);
		}
		_spin.push_back(tmp);
	}

	_m = 0;
	_E = 0;
	for(int i=0; i<_N; i++) {
		int il = i-1<  0? i-1+_N:i-1;
		int ir = i+1>=_N? i+1-_N:i+1;
		for(int j=0; j<_N; j++) {
			int jl = j-1<  0? j-1+_N:j-1;
			int jr = j+1>=_N? j+1-_N:j+1;

			_m += _spin[i][j];
			double sum = _spin[il][j] + _spin[ir][j] + _spin[i][jl] + _spin[i][jr];
			_E -= _spin[i][j]*(_H + 0.5*sum);
		}
	}
}

//-------------------------------------------------------------------------//

IsingModel2D::~IsingModel2D() {

}

//-------------------------------------------------------------------------//

bool IsingModel2D::check() const {

	if(_c!=-1 && _c!=1 && _c!=2) {
		cout << "ERROR: initial spin configuration invalid!" << endl;
		return false;
	}

	return true;
}

//-------------------------------------------------------------------------//

void IsingModel2D::cal_once() {

	for(int i=0; i<_N; i++) {
		int il = i-1<  0? i-1+_N:i-1;
		int ir = i+1>=_N? i+1-_N:i+1;
		for(int j=0; j<_N; j++) {
			int jl = j-1<  0? j-1+_N:j-1;
			int jr = j+1>=_N? j+1-_N:j+1;

			double sum = _spin[il][j] + _spin[ir][j] + _spin[i][jl] + _spin[i][jr];
			double DeltaE = 2*_spin[i][j]*(sum + _H);
			if(DeltaE<0 || exp(-DeltaE/_T)>1.0*rand()/RAND_MAX) {
				_spin[i][j] = -_spin[i][j];
				_m += 2*_spin[i][j];
				_E += DeltaE;
			}
		}
	}
}

//-------------------------------------------------------------------------//

void IsingModel2D::cal_until(int t) {

	if(!check()) return;

	for(int i=_now; i<t; i++) {
		cal_once();
	}

	_now = t;
}

//-------------------------------------------------------------------------//

