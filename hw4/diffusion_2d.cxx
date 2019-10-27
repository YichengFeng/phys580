#include "diffusion_2d.h"
#include <iostream>
#include <vector>
#include <string>
#include <cmath>

using namespace std;


//-------------------------------------------------------------------------//

Diffusion2D::Diffusion2D() {

	_name = "default";

	_M = 100;
	_N = 2*_M + 1;
	_P = 100;
	_P0 = _P;
	_D = 3;
	_S = 0;
	_L = 5;

	_mod = 0;
	_mod_sample = 0;
	_n_t = 6;
	_t_max = 1000;
	_t = 0;

	for(int i=0; i<_P; i++) {
		_x.push_back(0);
		_y.push_back(0);
	}
}

//-------------------------------------------------------------------------//

Diffusion2D::Diffusion2D(string name, int M, int P, int D, int n_t, double t_max) {

	_name = name;

	_M = M;
	_N = 2*_M + 1;
	_P = P;
	_P0 = _P;
	_D = D;
	_S = 0;
	_L = 5;

	_mod = 0;
	_mod_sample = 0;
	_n_t = n_t;
	_t_max = t_max;
	_t = 0;

	for(int i=0; i<_P; i++) {
		_x.push_back(0);
		_y.push_back(0);
	}
}

//-------------------------------------------------------------------------//

Diffusion2D::~Diffusion2D() {

}

//-------------------------------------------------------------------------//

void Diffusion2D::cal_until(double t_stop) {

	double delta[4][2] = {{1,0}, {-1,0}, {0,1}, {0,-1}};

	for(int i=int(_t); i<t_stop; i++) {
		int pid = int(1.0*_P*rand()/RAND_MAX);
		int idx = int(4.0*rand()/RAND_MAX);
		_x[pid] += delta[idx][0];
		_y[pid] += delta[idx][1];
		// periodic condition
		if(_x[pid]<-_M) _x[pid] += _N;
		if(_x[pid]> _M) _x[pid] -= _N;
		if(_y[pid]<-_M) _y[pid] += _N;
		if(_y[pid]> _M) _y[pid] -= _N;
	}

	_t = t_stop;

	cal_entropy();
}

//-------------------------------------------------------------------------//

void Diffusion2D::cal_leak_until(double t_stop) {

	double delta[4][2] = {{1,0}, {-1,0}, {0,1}, {0,-1}};

	for(int i=int(_t); i<t_stop; i++) {
		if(_P == 0) break;
		int pid = int(1.0*_P0*rand()/RAND_MAX);
		if(pid>=_P) continue;
		int idx = int(4.0*rand()/RAND_MAX);
		_x[pid] += delta[idx][0];
		_y[pid] += delta[idx][1];
		if(_x[pid]>=-_L && _x[pid]<=_L && _y[pid]>=_M) {
			_x.erase(_x.begin()+pid);
			_y.erase(_y.begin()+pid);
			_P -= 1;
			continue;
		}
		// reflection condition
		if(_x[pid]<-_M) _x[pid] = -2*_M - _x[pid];
		if(_x[pid]> _M) _x[pid] =  2*_M - _x[pid];
		if(_y[pid]<-_M) _y[pid] = -2*_M - _y[pid];
		if(_y[pid]> _M) _y[pid] =  2*_M - _y[pid];
	}

	_t = t_stop;

	cal_entropy();
}

//-------------------------------------------------------------------------//

bool Diffusion2D::check() const {

	if(_mod!=0 && _mod!=1) {
		cout << "ERROR: invalid mode!" << endl;
		return false;
	}
	if(_mod_sample!=0 && _mod_sample!=1) {
		cout << "ERROR: invalid sample mode!" << endl;
		return false;
	}

	return true;
}

//-------------------------------------------------------------------------//

void Diffusion2D::cal_all() {

	if(!check()) return;

	double pt;
	double t_stop;
	if(_mod_sample == 0) {
		pt = pow(_t_max, 1.0/_n_t);
		t_stop = 1.0;
	}
	if(_mod_sample == 1) {
		pt = _t_max/_n_t;
		t_stop = 0;
	}
	for(int i=0; i<_n_t; i++) {
		if(_mod_sample == 0) t_stop *= pt;
		if(_mod_sample == 1) t_stop += pt;
		if(_mod == 0) cal_until(t_stop);
		if(_mod == 1) cal_leak_until(t_stop);
		_t_sample.push_back(floor(t_stop));
		_x_sample.push_back(_x);
		_y_sample.push_back(_y);
		_S_sample.push_back(_S);
	}
}

//-------------------------------------------------------------------------//

double Diffusion2D::cal_entropy() {

	int dx = _N/_D;
	int dy = _N/_D;

	vector< vector<double> > ndist;
	for(int i=0; i<_D; i++) {
		vector<double> tmp;
		for(int j=0; j<_D; j++) {
			tmp.push_back(0);
		}
		ndist.push_back(tmp);
	}

	for(int i=0; i<_P; i++) {
		int ix = int((_x[i]+_M)/dx);
		if(ix<0) ix = 0;
		if(ix>=_D) ix = _D-1;
		int iy = int((_y[i]+_M)/dy);
		if(iy<0) iy = 0;
		if(iy>=_D) iy = _D-1;

		ndist[ix][iy] += 1.0;
	}

	double tmpS = 0;
	for(int ix=0; ix<_D; ix++) {
		for(int iy=0; iy<_D; iy++) {
			double ntmp = ndist[ix][iy];
			if(ntmp == 0) {
				tmpS -= 0;
			} else {
				tmpS -= 1.0*ntmp/_P*log(1.0*ntmp/_P);
			}
		}
	}

	_S = tmpS;
	return _S;
}

//-------------------------------------------------------------------------//

