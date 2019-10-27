#include <time.h>
#include "random_walk_3d.h"

using namespace std;


//-------------------------------------------------------------------------//

RandomWalk3D::RandomWalk3D() {

	//srand(time(nullptr));

	_n_step = 1000;
	_n_walk = 100;
	_mod = 0;
	_d = 1;
	_min_step = 0.1;

	for(int i=0; i<_n_step; i++) {
		_r2.push_back(0);
		_r4.push_back(0);
		_r2_std.push_back(0);
		_t.push_back(i);
		_good_walk.push_back(0);
	}
}

//-------------------------------------------------------------------------//

RandomWalk3D::RandomWalk3D(int n_step, int n_walk, int mod) {

	//srand(time(nullptr));

	_n_step = n_step;
	_n_walk = n_walk;
	_mod = mod;
	_d = 1;
	_min_step = 0.1;

	for(int i=0; i<_n_step; i++) {
		_r2.push_back(0);
		_r4.push_back(0);
		_r2_std.push_back(0);
		_t.push_back(i);
		_good_walk.push_back(0);
	}
}

//-------------------------------------------------------------------------//

RandomWalk3D::~RandomWalk3D() {

}

//-------------------------------------------------------------------------//

void RandomWalk3D::cal_rw_1() {

	_x.clear();
	_y.clear();
	_z.clear();

	_x.push_back(0);
	_y.push_back(0);
	_z.push_back(0);

	double delta[6][3] = {{1,0,0}, {-1,0,0}, {0,1,0}, {0,-1,0}, {0,0,1}, {0,0,-1}};

	for(int i=1; i<_n_step; i++) {

		double tmprand = 6.0*rand()/RAND_MAX;

		int idx = int(tmprand);

		_x.push_back(_x[i-1]+delta[idx][0]);
		_y.push_back(_y[i-1]+delta[idx][1]);
		_z.push_back(_z[i-1]+delta[idx][2]);
	}
}

//-------------------------------------------------------------------------//

void RandomWalk3D::cal_rw_d() {

	_x.clear();
	_y.clear();
	_z.clear();

	_x.push_back(0);
	_y.push_back(0);
	_z.push_back(0);

	double delta[6][3] = {{1,0,0}, {-1,0,0}, {0,1,0}, {0,-1,0}, {0,0,1}, {0,0,-1}};

	for(int i=1; i<_n_step; i++) {

		double d = _d*1.0*rand()/RAND_MAX;

		double tmprand = 6.0*rand()/RAND_MAX;

		int idx = int(tmprand);

		_x.push_back(_x[i-1]+delta[idx][0]*d);
		_y.push_back(_y[i-1]+delta[idx][1]*d);
		_z.push_back(_z[i-1]+delta[idx][2]*d);
	}
}

//-------------------------------------------------------------------------//

bool RandomWalk3D::cal_saw() {

	vector< vector< vector<bool> > > occupied;
	for(int i=0; i<_n_step*2+1; i++) {
		vector< vector<bool> > tmp2;
		for(int j=0; j<_n_step*2+1; j++) {
			vector<bool> tmp1;
			for(int k=0; k<_n_step*2+1; k++) {
				tmp1.push_back(false);
			}
			tmp2.push_back(tmp1);
		}
		occupied.push_back(tmp2);
	}

	_x.clear();
	_y.clear();
	_z.clear();

	_x.push_back(0);
	_y.push_back(0);
	_z.push_back(0);
	occupied[_n_step][_n_step][_n_step] = true;

	_x.push_back(1);
	_y.push_back(0);
	_z.push_back(0);
	occupied[_n_step+1][_n_step][_n_step] = true;

	int dir = 2;
	int nodir = 1;
	double delta[6][3] = {{1,0,0}, {-1,0,0}, {0,1,0}, {0,-1,0}, {0,0,1}, {0,0,-1}};
	int nx;
	int ny;
	int nz;
	int idx;

	for(int i=2; i<_n_step; i++) {

		idx = int(5.0*rand()/RAND_MAX);
		if(idx>=nodir) idx++;

		nx = int(_x[i-1]+delta[idx][0]);
		ny = int(_y[i-1]+delta[idx][1]);
		nz = int(_z[i-1]+delta[idx][2]);

		dir = idx;
		if(idx==0) nodir = 1;
		if(idx==1) nodir = 0;
		if(idx==2) nodir = 3;
		if(idx==3) nodir = 2;
		if(idx==4) nodir = 5;
		if(idx==5) nodir = 4;

		if(occupied[_n_step+nx][_n_step+ny][_n_step+nz]) {
			break;
		} else {
			occupied[_n_step+nx][_n_step+ny][_n_step+nz] = true;
			_x.push_back(1.0*nx);
			_y.push_back(1.0*ny);
			_z.push_back(1.0*nz);
		}
	}

	if(_x.size()<_min_step*_n_step) return false;

	return true;
}

//-------------------------------------------------------------------------//

void RandomWalk3D::cal() {

	for(int i_walk=0; i_walk<_n_walk; i_walk++) {

		if(_mod == 0) cal_rw_1();
		if(_mod == 1) cal_rw_d();
		if(_mod == 2) {
			if(!cal_saw()) continue;
		}

		for(int i_step=0; i_step<_x.size(); i_step++) {
			double r2tmp = _x[i_step]*_x[i_step] + _y[i_step]*_y[i_step] + _z[i_step]*_z[i_step];
			_r2[i_step] += r2tmp;
			_r4[i_step] += r2tmp*r2tmp;
			_good_walk[i_step] += 1.0;
		}
	}
	for(int i_step=0; i_step<_n_step; i_step++) {
		if(_good_walk[i_step] == 0) continue;
		_r2[i_step] /= _good_walk[i_step];
		_r4[i_step] /= _good_walk[i_step];
		_r2_std[i_step] = sqrt(_r4[i_step]-_r2[i_step]*_r2[i_step]);
	}
}
