#include <iostream>
#include <vector>
#include <cmath>
#include "charge_dist_3d.h"

using namespace std;


//-------------------------------------------------------------------------//

ChargeDist3D::ChargeDist3D() {
	_d = 0.1;
	_acc = 1e-4;
	_tmp_acc = 1;
	_alg = 0;
	_periodic = 0;
	_n_iter = 0;

	_M = int(0.5/_d+0.5);
	_d = 1.0/(2*_M);
	_N = 2*_M+1; // linear grid size

	_alpha = 2.0/(1.0+M_PI/_N);

	for(int k=0; k<_N; k++) {
		vector< vector<double> > tmp2d;
		for(int j=0; j<_N; j++) {
			vector<double> tmp;
			for(int i=0; i<_N; i++) {
				tmp.push_back(0);
			}
			tmp2d.push_back(tmp);
		}
		_V.push_back(tmp2d);
		_rho.push_back(tmp2d);
	}
}

//-------------------------------------------------------------------------//

ChargeDist3D::ChargeDist3D(double d, double acc) {
	_d = d;
	_acc = acc;
	_tmp_acc = 1;
	_alg = 0;
	_periodic = 0;
	_n_iter = 0;

	_M = int(0.5/_d+0.5);
	_d = 1.0/(2*_M);
	_N = 2*_M+1; // linear grid size

	_alpha = 2.0/(1.0+M_PI/_N);

	for(int k=0; k<_N; k++) {
		vector< vector<double> > tmp2d;
		for(int j=0; j<_N; j++) {
			vector<double> tmp;
			for(int i=0; i<_N; i++) {
				tmp.push_back(0);
			}
			tmp2d.push_back(tmp);
		}
		_V.push_back(tmp2d);
		_rho.push_back(tmp2d);
	}
}

//-------------------------------------------------------------------------//

ChargeDist3D::~ChargeDist3D() {

}

//-------------------------------------------------------------------------//

void ChargeDist3D::set_rho(double x1, double x2, double y1, double y2, double z1, double z2, double r) {
	int i1 = int((x1+0.5)/_d+0.5);
	int i2 = int((x2+0.5)/_d+0.5);
	int ni = fabs(i2-i1) + 1;
	int j1 = int((y1+0.5)/_d+0.5);
	int j2 = int((y2+0.5)/_d+0.5);
	int nj = fabs(j2-j1) + 1;
	int k1 = int((z1+0.5)/_d+0.5);
	int k2 = int((z2+0.5)/_d+0.5);
	int nk = fabs(k2-k1) + 1;
	for(int k=k1; k<=k2; k++) {
		for(int j=j1; j<=j2; j++) {
			for(int i=i1; i<=i2; i++) {
				_rho[k][j][i] = r/ni/nj/nk;
			}
		}
	}
}

//-------------------------------------------------------------------------//

bool ChargeDist3D::check() const {

	bool isPass = true;

	int n1 = _V.size();
	int n2 = _V[0].size();
	int n3 = _V[0][0].size();
	if(n1 != n2) isPass = false;
	if(n1 != n3) isPass = false;
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

void ChargeDist3D::cal_once_Jacobi() {

	if(!check()) return;

	vector< vector< vector<double> > > Vnew = _V;
	_tmp_acc = 0;

	for(int i=1; i<_N-1; i++) {
		for(int j=1; j<_N-1; j++) {
			for(int k=1; k<_N-1; k++) {
				Vnew[i][j][k] = (_V[i-1][j][k] + _V[i+1][j][k] + _V[i][j-1][k] + _V[i][j+1][k] + _V[i][j][k-1] + _V[i][j][k+1] + _rho[i][j][k])/6;
				_tmp_acc += fabs(Vnew[i][j][k] - _V[i][j][k]);
			}
		}
	}

	_tmp_acc = _tmp_acc/_N/_N/_N;

	_V = Vnew;

	_n_iter ++;
}

//-------------------------------------------------------------------------//

void ChargeDist3D::cal_once_Jacobi_periodic() {

	if(!check()) return;

	vector< vector< vector<double> > > Vnew = _V;
	_tmp_acc = 0;

	for(int i=0; i<_N; i++) {
		for(int j=0; j<_N; j++) {
			for(int k=0; k<_N; k++) {
				int il = i-1; if(il < 0) il += _N;
				int ir = i+1; if(ir>=_N) ir -= _N;
				int jl = j-1; if(jl < 0) jl += _N;
				int jr = j+1; if(jl>=_N) jr -= _N;
				int kl = k-1; if(kl < 0) kl += _N;
				int kr = k+1; if(kr>=_N) kr -= _N;
				Vnew[i][j][k] = (_V[il][j][k] + _V[ir][j][k] + _V[i][jl][k] + _V[i][jr][k] + _V[i][j][kl] + _V[i][j][kr] + _rho[i][j][k])/6;
				_tmp_acc += fabs(Vnew[i][j][k] - _V[i][j][k]);
			}
		}
	}

	_tmp_acc = _tmp_acc/_N/_N/_N;

	_V = Vnew;

	_n_iter ++;
}

//-------------------------------------------------------------------------//

void ChargeDist3D::cal_once_GaussSeidel() {

	if(!check()) return;

	vector< vector< vector<double> > > Vold = _V;
	_tmp_acc = 0;

	for(int i=1; i<_N-1; i++) {
		for(int j=1; j<_N-1; j++) {
			for(int k=1; k<_N-1; k++) {
				_V[i][j][k] = (_V[i-1][j][k] + _V[i+1][j][k] + _V[i][j-1][k] + _V[i][j+1][k] + _V[i][j][k-1] + _V[i][j][k+1] + _rho[i][j][k])/6;
				_tmp_acc += fabs(Vold[i][j][k] - _V[i][j][k]);
			}
		}
	}

	_tmp_acc = _tmp_acc/_N/_N/_N;

	_n_iter ++;
}

//-------------------------------------------------------------------------//

void ChargeDist3D::cal_once_GaussSeidel_periodic() {

	if(!check()) return;

	vector< vector< vector<double> > > Vold = _V;
	_tmp_acc = 0;

	for(int i=0; i<_N; i++) {
		for(int j=0; j<_N; j++) {
			for(int k=0; k<_N; k++) {
				int il = i-1; if(il < 0) il += _N;
				int ir = i+1; if(ir>=_N) ir -= _N;
				int jl = j-1; if(jl < 0) jl += _N;
				int jr = j+1; if(jl>=_N) jr -= _N;
				int kl = k-1; if(kl < 0) kl += _N;
				int kr = k+1; if(kr>=_N) kr -= _N;
				_V[i][j][k] = (_V[il][j][k] + _V[ir][j][k] + _V[i][jl][k] + _V[i][jr][k] + _V[i][j][kl] + _V[i][j][kr] + _rho[i][j][k])/6;
				_tmp_acc += fabs(Vold[i][j][k] - _V[i][j][k]);
			}
		}
	}

	_tmp_acc = _tmp_acc/_N/_N/_N;

	_n_iter ++;
}

//-------------------------------------------------------------------------//

void ChargeDist3D::cal_once_SOR() {

	if(!check()) return;

	vector< vector< vector<double> > > Vold = _V;
	_tmp_acc = 0;

	for(int i=1; i<_N-1; i++) {
		for(int j=1; j<_N-1; j++) {
			for(int k=1; k<_N-1; k++) {
				_V[i][j][k] = (_V[i-1][j][k] + _V[i+1][j][k] + _V[i][j-1][k] + _V[i][j+1][k] + _V[i][j][k-1] + _V[i][j][k+1] + _rho[i][j][k])/6;
			}
		}
	}

	for(int i=1; i<_N-1; i++) {
		for(int j=1; j<_N-1; j++) {
			for(int k=1; k<_N-1; k++) {
				_V[i][j][k] = _alpha*_V[i][j][k] + (1-_alpha)*Vold[i][j][k];
				_tmp_acc += fabs(Vold[i][j][k] - _V[i][j][k]);
			}
		}
	}

	_tmp_acc = _tmp_acc/_N/_N/_N;

	_n_iter ++;
}

//-------------------------------------------------------------------------//

void ChargeDist3D::cal_once_SOR_periodic() {

	if(!check()) return;

	vector< vector< vector<double> > > Vold = _V;
	_tmp_acc = 0;

	for(int i=0; i<_N; i++) {
		for(int j=0; j<_N; j++) {
			for(int k=0; k<_N; k++) {
				int il = i-1; if(il < 0) il += _N;
				int ir = i+1; if(ir>=_N) ir -= _N;
				int jl = j-1; if(jl < 0) jl += _N;
				int jr = j+1; if(jl>=_N) jr -= _N;
				int kl = k-1; if(kl < 0) kl += _N;
				int kr = k+1; if(kr>=_N) kr -= _N;
				_V[i][j][k] = (_V[il][j][k] + _V[ir][j][k] + _V[i][jl][k] + _V[i][jr][k] + _V[i][j][kl] + _V[i][j][kr] + _rho[i][j][k])/6;
			}
		}
	}

	for(int i=1; i<_N-1; i++) {
		for(int j=1; j<_N-1; j++) {
			for(int k=1; k<_N-1; k++) {
				_V[i][j][k] = _alpha*_V[i][j][k] + (1-_alpha)*Vold[i][j][k];
				_tmp_acc += fabs(Vold[i][j][k] - _V[i][j][k]);
			}
		}
	}

	_tmp_acc = _tmp_acc/_N/_N/_N;

	_n_iter ++;
}

//-------------------------------------------------------------------------//

void ChargeDist3D::cal() {

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
