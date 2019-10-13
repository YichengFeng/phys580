#include <iostream>
#include <vector>
#include <cmath>
#include "capacitor_2d.h"

using namespace std;


//-------------------------------------------------------------------------//

Capacitor2D::Capacitor2D() {
	_a = 0.5; // plate width
	_b = 0.5; // gap
	_d = 0.1;
	_acc = 1e-4;
	_tmp_acc = 1;
	_alg = 0;
	_periodic = 0;
	_n_iter = 0;

	_M = int(0.5/_d+0.5);
	_d = 1.0/(2*_M);
	_N = 2*_M+1; // linear grid size
	_aM = int(_a*_M + 0.5);
	_bM = int(_b*_M + 0.5);

	_alpha = 2.0/(1.0+M_PI/_N);

	for(int i=0; i<_N; i++) {
		vector<double> tmp;
		for(int j=0; j<_N; j++) {
			tmp.push_back(0);
		}
		_V.push_back(tmp);
		_rho.push_back(tmp);
	}

	for(int i=_M+1-_aM; i<=_M+1+_aM; i++) {
		_V[i][_M+1+_bM] = -1;
		_V[i][_M+1-_bM] = 1;
	}
}

//-------------------------------------------------------------------------//

Capacitor2D::Capacitor2D(double a, double b, double d, double acc) {
	_a = a;
	_b = b;
	_d = d;
	_acc = acc;
	_tmp_acc = 1;
	_alg = 0;
	_periodic = 0;
	_n_iter = 0;

	_M = int(0.5/_d+0.5);
	_d = 1.0/(2*_M);
	_N = 2*_M+1; // linear grid size
	_aM = int(_a*_M + 0.5);
	_bM = int(_b*_M + 0.5);

	_alpha = 2.0/(1.0+M_PI/_N);

	for(int i=0; i<_N; i++) {
		vector<double> tmp;
		for(int j=0; j<_N; j++) {
			tmp.push_back(0);
		}
		_V.push_back(tmp);
		_rho.push_back(tmp);
	}

	for(int i=_M+1-_aM; i<=_M+1+_aM; i++) {
		_V[i][_M+1+_bM] = -1;
		_V[i][_M+1-_bM] = 1;
	}
}

//-------------------------------------------------------------------------//

Capacitor2D::~Capacitor2D() {

}

//-------------------------------------------------------------------------//

void Capacitor2D::set_rho(double x1, double x2, double y1, double y2, double r) {
	int i1 = int((x1+0.5)/_d+0.5);
	int i2 = int((x2+0.5)/_d+0.5);
	int ni = fabs(i2-i1) + 1;
	int j1 = int((y1+0.5)/_d+0.5);
	int j2 = int((y2+0.5)/_d+0.5);
	int nj = fabs(j2-j1) + 1;
	for(int j=j1; j<=j2; j++) {
		for(int i=i1; i<=i2; i++) {
			_rho[j][i] = r/ni/nj;
		}
	}
}

//-------------------------------------------------------------------------//

bool Capacitor2D::check() const {

	bool isPass = true;

	int n1 = _V.size();
	int n2 = _V[0].size();
	if(n1 != n2) isPass = false;
	for(int i=0; i<n1; i++) {
		if(_V[i].size() !=  n2) isPass = false;
	}
	if(_alg!=0 && _alg!=1 && _alg!=2) {
		cout << "ERROR: algorithm not specified!" << endl;
		isPass = false;
	}
	if(_periodic!=0 && _periodic!=1) {
		cout << "ERROR: boundary condition not specified!" << endl;
		isPass = false;
	}

	return isPass;
}

//-------------------------------------------------------------------------//

void Capacitor2D::cal_once_Jacobi() {

	if(!check()) return;

	vector< vector<double> > Vnew = _V;
	_tmp_acc = 0;

	for(int j=1; j<_N-1; j++) {
		for(int i=1; i<_N-1; i++) {
			if(_V[i][j]==1 || _V[i][j]==-1) continue;
			Vnew[i][j] = (_V[i-1][j] + _V[i+1][j] + _V[i][j-1] + _V[i][j+1] + _rho[i][j])*0.25;
			_tmp_acc += fabs(Vnew[i][j] - _V[i][j]);
		}
	}

	_tmp_acc = _tmp_acc/_N/_N;

	_V = Vnew;

	_n_iter ++;
}

//-------------------------------------------------------------------------//

void Capacitor2D::cal_once_Jacobi_periodic() {

	if(!check()) return;

	vector< vector<double> > Vnew = _V;
	_tmp_acc = 0;

	for(int j=0; j<_N; j++) {
		for(int i=0; i<_N; i++) {
			if(_V[i][j]==1 || _V[i][j]==-1) continue;
			int jl = j-1; if(jl < 0) jl += _N;
			int jr = j+1; if(jl>=_N) jr -= _N;
			int il = i-1; if(il < 0) il += _N;
			int ir = i+1; if(ir>=_N) ir -= _N;
			Vnew[i][j] = (_V[il][j] + _V[ir][j] + _V[i][jl] + _V[i][jr] + _rho[i][j])*0.25;
			_tmp_acc += fabs(Vnew[i][j] - _V[i][j]);
		}
	}

	_tmp_acc = _tmp_acc/_N/_N;

	_V = Vnew;

	_n_iter ++;
}

//-------------------------------------------------------------------------//

void Capacitor2D::cal_once_GaussSeidel() {

	if(!check()) return;

	vector< vector<double> > Vold = _V;
	_tmp_acc = 0;

	for(int j=1; j<_N-1; j++) {
		for(int i=1; i<_N-1; i++) {
			if(_V[i][j]==1 || _V[i][j]==-1) continue;
			_V[i][j] = (_V[i-1][j] + _V[i+1][j] + _V[i][j-1] + _V[i][j+1] + _rho[i][j])*0.25;
			_tmp_acc += fabs(Vold[i][j] - _V[i][j]);
		}
	}

	_tmp_acc = _tmp_acc/_N/_N;

	_n_iter ++;
}

//-------------------------------------------------------------------------//

void Capacitor2D::cal_once_GaussSeidel_periodic() {

	if(!check()) return;

	vector< vector<double> > Vold = _V;
	_tmp_acc = 0;

	for(int j=0; j<_N; j++) {
		for(int i=0; i<_N; i++) {
			if(_V[i][j]==1 || _V[i][j]==-1) continue;
			int jl = j-1; if(jl < 0) jl += _N;
			int jr = j+1; if(jl>=_N) jr -= _N;
			int il = i-1; if(il < 0) il += _N;
			int ir = i+1; if(ir>=_N) ir -= _N;
			_V[i][j] = (_V[il][j] + _V[ir][j] + _V[i][jl] + _V[i][jr] + _rho[i][j])*0.25;
			_tmp_acc += fabs(Vold[i][j] - _V[i][j]);
		}
	}

	_tmp_acc = _tmp_acc/_N/_N;

	_n_iter ++;
}

//-------------------------------------------------------------------------//

void Capacitor2D::cal_once_SOR() {

	if(!check()) return;

	vector< vector<double> > Vold = _V;
	_tmp_acc = 0;

	for(int j=1; j<_N-1; j++) {
		for(int i=1; i<_N-1; i++) {
			if(_V[i][j]==1 || _V[i][j]==-1) continue;
			_V[i][j] = (_V[i-1][j] + _V[i+1][j] + _V[i][j-1] + _V[i][j+1] + _rho[i][j])*0.25;
		}
	}

	for(int j=1; j<_N-1; j++) {
		for(int i=1; i<_N-1; i++) {
			if(_V[i][j]==1 || _V[i][j]==-1) continue;
			_V[i][j] = _alpha*_V[i][j] + (1-_alpha)*Vold[i][j];
			_tmp_acc += fabs(Vold[i][j] - _V[i][j]);
		}
	}

	_tmp_acc = _tmp_acc/_N/_N;

	_n_iter ++;
}

//-------------------------------------------------------------------------//

void Capacitor2D::cal_once_SOR_periodic() {

	if(!check()) return;

	vector< vector<double> > Vold = _V;
	_tmp_acc = 0;

	for(int j=0; j<_N; j++) {
		for(int i=0; i<_N; i++) {
			if(_V[i][j]==1 || _V[i][j]==-1) continue;
			int jl = j-1; if(jl < 0) jl += _N;
			int jr = j+1; if(jl>=_N) jr -= _N;
			int il = i-1; if(il < 0) il += _N;
			int ir = i+1; if(ir>=_N) ir -= _N;
			_V[i][j] = (_V[il][j] + _V[ir][j] + _V[i][jl] + _V[i][jr] + _rho[i][j])*0.25;
		}
	}

	for(int j=1; j<_N-1; j++) {
		for(int i=1; i<_N-1; i++) {
			if(_V[i][j]==1 || _V[i][j]==-1) continue;
			_V[i][j] = _alpha*_V[i][j] + (1-_alpha)*Vold[i][j];
			_tmp_acc += fabs(Vold[i][j] - _V[i][j]);
		}
	}

	_tmp_acc = _tmp_acc/_N/_N;

	_n_iter ++;
}

//-------------------------------------------------------------------------//

void Capacitor2D::cal() {

	if(!check()) return;

	if(_alg == 0) {
		if(_periodic == 0){
			while(_tmp_acc > _acc) {
				cal_once_Jacobi();
			}
		}
		if(_periodic == 1){
			while(_tmp_acc > _acc) {
				cal_once_Jacobi_periodic();
			}
		}
	}
	if(_alg == 1) {
		if(_periodic == 0){
			while(_tmp_acc > _acc) {
				cal_once_GaussSeidel();
			}
		}
		if(_periodic == 1){
			while(_tmp_acc > _acc) {
				cal_once_GaussSeidel_periodic();
			}
		}
	}
	if(_alg == 2) {
		if(_periodic == 0){
			//while(_tmp_acc > _acc && _tmp_acc<1000) {
			while(_tmp_acc > _acc) {
				cal_once_SOR();
			}
		}
		if(_periodic == 1){
			//while(_tmp_acc > _acc && _tmp_acc<1000) {
			while(_tmp_acc > _acc) {
				cal_once_SOR_periodic();
			}
		}
	}
}
