#include "loop_magnetic_field.h"

using namespace std;


//-------------------------------------------------------------------------//

LoopMagneticField::LoopMagneticField() {

	_n_loop = 1;
	_n_phi = 20;
	_R = 0.5;
	_pitch = 0;
	_d = 0.05;
	_mod = 0;

	_M = int(1.0/_d+0.5);
	_d = 1.0/_M;
	_N = 2*_M+1;

	for(int i=0; i<_n_loop; i++) {
		_z_loop.push_back(-0.5*_pitch*(_n_loop-1) + i*_pitch);
		_I_loop.push_back(1.0);
	}
}

//-------------------------------------------------------------------------//

LoopMagneticField::LoopMagneticField(int n_phi, double R, int n_loop, double pitch, double d) {

	_n_loop = n_loop;
	_n_phi = n_phi;
	_R = R;
	_pitch = pitch;
	_d = d;
	_mod = 0;

	_M = int(1.0/_d+0.5);
	_d = 1.0/_M;
	_N = 2*_M+1;

	for(int i=0; i<_n_loop; i++) {
		_z_loop.push_back(-0.5*_pitch*(_n_loop-1) + i*_pitch);
		_I_loop.push_back(1.0);
	}
}

//-------------------------------------------------------------------------//

LoopMagneticField::~LoopMagneticField() {

}

//-------------------------------------------------------------------------//

vector<double> LoopMagneticField::cal_one_grid(double x, double y, double z) {

	double bx = 0;
	double by = 0;
	double bz = 0;

	double dphi = 2*M_PI/_n_phi;

	for(int i_loop=0; i_loop<_n_loop; i_loop++) {
		for(int i=0; i<_n_phi; i++) {
			double phi = i*dphi;
			double dlx = -_R*dphi*sin(phi);
			double dly =  _R*dphi*cos(phi);
			double dlz = 0;
			double rx = x - _R*cos(phi);
			double ry = y - _R*sin(phi);
			double rz = z - _z_loop[i_loop];
			double r = sqrt(rx*rx + ry*ry + rz*rz);
			if(r>_R*1e-4) {
				bx += _I_loop[i_loop]*(dly*rz - dlz*ry)/(r*r*r);
				by += _I_loop[i_loop]*(dlz*rx - dlx*rz)/(r*r*r);
				bz += _I_loop[i_loop]*(dlx*ry - dly*rx)/(r*r*r);
			}
		}
	}

	vector<double> b(3);
	b[0] = bx;
	b[1] = by;
	b[2] = bz;

	return b;
}

//-------------------------------------------------------------------------//

vector<double> LoopMagneticField::cal_one_grid_helical(double x, double y, double z) {

	double bx = 0;
	double by = 0;
	double bz = 0;

	double dphi = 2*M_PI/_n_phi*_n_loop;

	double z_low = -0.5*_pitch*_n_loop;

	for(int i=0; i<_n_phi; i++) {
		double phi = i*dphi;
		double dlx = -_R*dphi*sin(phi);
		double dly =  _R*dphi*cos(phi);
		double dlz =  _pitch/2.0/M_PI;
		double rx = x - _R*cos(phi);
		double ry = y - _R*sin(phi);
		double rz = z - (z_low+_pitch/2.0/M_PI*phi);
		double r = sqrt(rx*rx + ry*ry + rz*rz);
		if(r>_R*1e-4) {
			bx += _I_loop[0]*(dly*rz - dlz*ry)/(r*r*r);
			by += _I_loop[0]*(dlz*rx - dlx*rz)/(r*r*r);
			bz += _I_loop[0]*(dlx*ry - dly*rx)/(r*r*r);
		}
	}

	vector<double> b(3);
	b[0] = bx;
	b[1] = by;
	b[2] = bz;

	return b;
}

//-------------------------------------------------------------------------//

void LoopMagneticField::cal() {

	for(int k=0; k<_N; k++) {
		double z = _d*(k-_M);
		vector<double> tmpbx;
		vector<double> tmpby;
		vector<double> tmpbz;
		for(int i=0; i<_N; i++) {
			double x = _d*(i-_M);
			vector<double> b;
			if(_mod==0) b = cal_one_grid(x,0,z);
			if(_mod==1) b = cal_one_grid_helical(x,0,z);
			tmpbx.push_back(b[0]);
			tmpby.push_back(b[1]);
			tmpbz.push_back(b[2]);
		}
		_bx.push_back(tmpbx);
		_by.push_back(tmpby);
		_bz.push_back(tmpbz);
	}
}
