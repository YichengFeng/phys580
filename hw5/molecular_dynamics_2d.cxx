#include <iostream>
#include <vector>
#include <cmath>
#include "molecular_dynamics_2d.h"

using namespace std;


//-------------------------------------------------------------------------//

MolecularDynamics2D::MolecularDynamics2D() {

	MolecularDynamics2D(25, 5, 0.01, 0.05, 0);
}

//-------------------------------------------------------------------------//

MolecularDynamics2D::MolecularDynamics2D(int N, double L, double dt, double dmax, double vmax) {

	_N = N;
	_rec = 10;
	_L = L;
	_dt = dt;
	_dmax = dmax;
	_vmax = vmax;
	_t_now = 0;

	_tag_displace_now = 0;
	_tag_distance_now = 0;

	_x.clear();
	_y.clear();
	_vx.clear();
	_vy.clear();

	for(int i=0; i<_N; i++) {
		_x_old.push_back(0);
		_x_now.push_back(0);
		_x_new.push_back(0);
		_y_old.push_back(0);
		_y_now.push_back(0);
		_y_new.push_back(0);

		_vx_now.push_back(0);
		_vy_now.push_back(0);
	}

	double sqN = sqrt(_N);
	double grid = sqN==floor(sqN)?_L/sqN:_L/(sqN+1);

	int n = 0;
	double i = 0;
	while(i < _L) {
		double j = 0;
		while(j < _L) {
			if(n >= _N) break;
			_x_now[n] = i + 0.5*grid + _dmax*(1.0*rand()/RAND_MAX-0.5)*grid*sqrt(2);
			_y_now[n] = j + 0.5*grid + _dmax*(1.0*rand()/RAND_MAX-0.5)*grid*sqrt(2);
			_vx_now[n] = _vmax*(1.0*rand()/RAND_MAX-0.5)*sqrt(2);
			_vy_now[n] = _vmax*(1.0*rand()/RAND_MAX-0.5)*sqrt(2);
			_x_old[n] = _x_now[n] - _vx_now[n]*_dt;
			_y_old[n] = _y_now[n] - _vy_now[n]*_dt;

			n ++;
			j += grid;
		}
		i += grid;
	}

	_x_start = _x_old;
	_y_start = _y_old;

	cout << "initialization completed" << endl;
}

//-------------------------------------------------------------------------//

MolecularDynamics2D::~MolecularDynamics2D() {

}

//-------------------------------------------------------------------------//

bool MolecularDynamics2D::check() const {

	if(_x_old.size() != _N) return false;
	if(_x_now.size() != _N) return false;
	if(_x_new.size() != _N) return false;
	if(_y_old.size() != _N) return false;
	if(_y_now.size() != _N) return false;
	if(_y_new.size() != _N) return false;

	if(_vx_now.size() != _N) return false;
	if(_vy_now.size() != _N) return false;

	return true;
}

//-------------------------------------------------------------------------//

void MolecularDynamics2D::cal_once() {

	for(int i=0; i<_N; i++) {
		double fx = 0;
		double fy = 0;

		for(int j=0; j<_N; j++) {
			if(j == i) continue;

			double dx = _x_now[i] - _x_now[j];
			double dy = _y_now[i] - _y_now[j];

			if(fabs(dx) > 0.5*_L) dx -= dx>0?_L:-_L;
			if(fabs(dy) > 0.5*_L) dy -= dy>0?_L:-_L;

			double r = sqrt(dx*dx + dy*dy);
			if(r < 3) {
				double fij = 24.0*(2.0/pow(r,13) - 1.0/pow(r,7));
				fx += fij*dx/r;
				fy += fij*dy/r;
			}
		}

		_x_new[i] = 2*_x_now[i] - _x_old[i] + fx*_dt*_dt;
		_y_new[i] = 2*_y_now[i] - _y_old[i] + fy*_dt*_dt;
		_vx_now[i] = (_x_new[i] - _x_old[i])/(2.0*_dt);
		_vy_now[i] = (_y_new[i] - _y_old[i])/(2.0*_dt);

		if(i == 0) {
			double dx = _x_new[i] - _x_start[i];
			double dy = _y_new[i] - _y_start[i];
			if(fabs(dx) > 0.5*_L) dx -= dx>0?_L:-_L;
			if(fabs(dy) > 0.5*_L) dy -= dy>0?_L:-_L;
			_tag_displace_now = sqrt(dx*dx + dy*dy);
		} else if(i == 1) {
			double dx = _x_new[i] - _x_new[0];
			double dy = _y_new[i] - _y_new[0];
			if(fabs(dx) > 0.5*_L) dx -= dx>0?_L:-_L;
			if(fabs(dy) > 0.5*_L) dy -= dy>0?_L:-_L;
			_tag_distance_now = sqrt(dx*dx + dy*dy);
		}

		if(_x_new[i]<0) {
			_x_new[i] += _L;
			_x_now[i] += _L;
		} else if(_x_new[i]>_L) {
			_x_new[i] -= _L;
			_x_now[i] -= _L;
		}
		if(_y_new[i]<0) {
			_y_new[i] += _L;
			_y_now[i] += _L;
		} else if(_y_new[i]>_L) {
			_y_new[i] -= _L;
			_y_now[i] -= _L;
		}
	}

	_x_old = _x_now;
	_x_now = _x_new;
	_y_old = _y_now;
	_y_now = _y_new;
}

//-------------------------------------------------------------------------//

void MolecularDynamics2D::cal_ET() {

	double Ekin = 0;
	double Epot = 0;

	for(int i=0; i<_N; i++) {
		Ekin += 0.5*(_vx_now[i]*_vx_now[i] + _vy_now[i]*_vy_now[i]);

		for(int j=i+1; j<_N; j++) {

			double dx = _x_old[j] - _x_old[i];
			double dy = _y_old[j] - _y_old[i];

			if(fabs(dx) > 0.5*_L) dx -= dx>0?_L:-_L;
			if(fabs(dy) > 0.5*_L) dy -= dy>0?_L:-_L;

			double r = sqrt(dx*dx + dy*dy);
			if(r < 3) {
				double invr6 = 1.0/pow(r,6);
				double invr12 = invr6*invr6;
				Epot += 4*(invr12 - invr6);
			}
		}
	}

	_E_now = Ekin + Epot;
	_T_now = Ekin/_N;
}

//-------------------------------------------------------------------------//

void MolecularDynamics2D::cal_until(double t_end) {

	if(!check()) {
		cout << "ERROR: check() not pass!" << endl;
		return;
	}
	cout << "check passed" << endl;

	_t.clear();
	_x.clear();
	_y.clear();
	_vx.clear();
	_vy.clear();
	_tag_displace.clear();
	_tag_distance.clear();

	_E.clear();
	_T.clear();

	int n = 0;

	while(_t_now < t_end) {

		if(n%_rec == 0) {
			_t.push_back(_t_now);
			_x.push_back(_x_now);
			_y.push_back(_y_now);
			_vx.push_back(_vx_now);
			_vy.push_back(_vy_now);
			_tag_displace.push_back(_tag_displace_now);
			_tag_distance.push_back(_tag_distance_now);

			cal_ET();
			_E.push_back(_E_now);
			_T.push_back(_T_now);
		}

		cal_once();

		n ++;
		_t_now += _dt;
	}
}

//-------------------------------------------------------------------------//

void MolecularDynamics2D::change_velocity(double factor) {

	for(int i=0; i<_N; i++) {
		_x_old[i] = _x_now[i] - factor*(_x_now[i] - _x_old[i]);
		_y_old[i] = _y_now[i] - factor*(_y_now[i] - _y_old[i]);
	}
}
