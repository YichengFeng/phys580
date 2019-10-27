#include "percolation_2d.h"
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;


//-------------------------------------------------------------------------//

Percolation2D::Percolation2D() {

	_p = 0.593;
	_L = 100;
	_trial = 50;

	for(int i=0; i<_L; i++) {
		vector<bool> tmp_occ;
		vector<int>  tmp_cid;
		for(int j=0; j<_L; j++) {
			tmp_occ.push_back(false);
			tmp_cid.push_back(0);
		}
		_occupied.push_back(tmp_occ);
		_cluster_id.push_back(tmp_cid);
	}
}

//-------------------------------------------------------------------------//

Percolation2D::Percolation2D(double p, int L, int trial) {

	_p = p;
	_L = L;
	_trial = trial;

	for(int i=0; i<_L; i++) {
		vector<bool> tmp_occ;
		vector<int>  tmp_cid;
		for(int j=0; j<_L; j++) {
			tmp_occ.push_back(false);
			tmp_cid.push_back(0);
		}
		_occupied.push_back(tmp_occ);
		_cluster_id.push_back(tmp_cid);
	}
}

//-------------------------------------------------------------------------//

void Percolation2D::cal_perculation() {

	for(int i=0; i<_L; i++) {
		for(int j=0; j<_L; j++) {
			if(1.0*rand()/RAND_MAX<_p) {
				_occupied[i][j] = true;
			} else {
				_occupied[i][j] = false;
			}
		}
	}
}

//-------------------------------------------------------------------------//

void Percolation2D::cal_DFS(int ix, int iy, int cid) {

	if(ix<0 || ix>=_L) return;
	if(iy<0 || iy>=_L) return;

	if(_cluster_id[iy][ix] == 0 && _occupied[iy][ix]) {
		_cluster_id[iy][ix] = cid;
	} else {
		return;
	}

	cal_DFS(ix+1, iy, cid);
	cal_DFS(ix-1, iy, cid);
	cal_DFS(ix, iy+1, cid);
	cal_DFS(ix, iy-1, cid);
}

//-------------------------------------------------------------------------//

void Percolation2D::cal_cluster() {

	for(int i=0; i<_L; i++) {
		for(int j=0; j<_L; j++) {
			_cluster_id[i][j] = 0;
		}
	}

	int cid = 0;

	for(int iy=0; iy<_L; iy++) {
		for(int ix=0; ix<_L; ix++) {
			if(_cluster_id[iy][ix]==0 && _occupied[iy][ix]) {
				cid ++;
				cal_DFS(ix, iy, cid);
			}
		}
	}

	_cluster_size.clear();
	_cluster_x.clear();
	_cluster_y.clear();
	_spanning_cluster_id.clear();

	for(int i=0; i<=cid; i++) {
		_cluster_size.push_back(0);
		_cluster_x.push_back(vector<int>());
		_cluster_y.push_back(vector<int>());
	}

	for(int iy=0; iy<_L; iy++) {
		for(int ix=0; ix<_L; ix++) {
			int idx = _cluster_id[iy][ix];
			_cluster_size[idx] ++;
			_cluster_x[idx].push_back(ix);
			_cluster_y[idx].push_back(iy);
		}
	}

	int max = 0;
	int max_id = 1;
	for(int i=1; i<=cid; i++) {
		if(_cluster_size[i]>max) {
			max = _cluster_size[i];
			max_id = i;
		}
	}
	for(int i=1; i<=cid; i++) {
		if(_cluster_size[i]==max) {
			_spanning_cluster_id.push_back(i);
		}
	}

	_n_cluster = cid;
}

//-------------------------------------------------------------------------//

bool Percolation2D::check() {

	int n1 = 0;
	for(int iy=0; iy<_L; iy++) {
		for(int ix=0; ix<_L; ix++) {
			if(_occupied[iy][ix]) n1++;
		}
	}

	int n2 = 0;
	for(int i=1; i<_cluster_size.size(); i++) {
		n2 += _cluster_size[i];
	}

	_n_occupied = n1;

	if(n1 != n2) {
		cout << "ERROR: occupied number not match!" << endl;
		return false;
	}

	return true;
}

//-------------------------------------------------------------------------//

void Percolation2D::cal_PS() {

	if(!check()) return;
	
	//_P = 1.0*_spanning_cluster_id.size()*_cluster_size[_spanning_cluster_id[0]]/_n_occupied;
	_P = 1.0*_cluster_size[_spanning_cluster_id[0]]/_n_occupied;

	double sum2 = 0;
	double tmpi = 0;
	for(int i=1; i<_cluster_size.size(); i++) {
		if(tmpi < _spanning_cluster_id.size() && i == _spanning_cluster_id[tmpi]){
			tmpi ++;
			continue;
		}
		sum2 += 1.0*_cluster_size[i]*_cluster_size[i];
	}
	_S = sum2/_L/_L;
}

//-------------------------------------------------------------------------//

void Percolation2D::cal_one_trial() {

	cal_perculation();
	cal_cluster();
	cal_PS();
}

//-------------------------------------------------------------------------//

void Percolation2D::cal_trials() {

	double P_sum = 0;
	double S_sum = 0;

	for(int i=0; i<_trial; i++) {
		cal_one_trial();
		P_sum += _P;
		S_sum += _S;
		_vP.push_back(_P);
		_vS.push_back(_S);
	}

	_P_average = P_sum/_trial;
	_S_average = S_sum/_trial;
}

//-------------------------------------------------------------------------//


