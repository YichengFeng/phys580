#include <iostream>
#include "projectile_3d.h"

using namespace std;


//-------------------------------------------------------------------------//

// static variables
double Projectile3D::_g    = 9.8;
double Projectile3D::_B2_m = 4e-5;
double Projectile3D::_Y0   = 1e4;
double Projectile3D::_a    = 6.5e-3;
double Projectile3D::_T0   = 300;
double Projectile3D::_T_ref= 300;
double Projectile3D::_gamma= 1.4;
double Projectile3D::_wx   = 0;
double Projectile3D::_wy   = 0;
double Projectile3D::_wz   = 0;

//-------------------------------------------------------------------------//

Projectile3D::Projectile3D() {
	_mode = 0;
	_dt = 0.02;
	_x_start = 0;
	_y_start = 0;
	_z_start = 0;
	_rk_order = 1;

	for(int i=0; i<_n_eqns; i++) {
		_c_start.push_back(0);
	}
}

//-------------------------------------------------------------------------//

Projectile3D::~Projectile3D() {

}

//-------------------------------------------------------------------------//

// x[0]: x; x[1]: y; x[2]: z; x[3]: vx; x[4]: vy; x[5]: vz
double Projectile3D::f_x(double t, const vector<double> &x) {
	// dx/dt = vx
	return x[3];
}

//-------------------------------------------------------------------------//

double Projectile3D::f_y(double t, const vector<double> &x) {
	// dy/dt = vy
	return x[4];
}

//-------------------------------------------------------------------------//

double Projectile3D::f_z(double t, const vector<double> &x) {
	// dz/dt = vz
	return x[5];
}

//-------------------------------------------------------------------------//

double Projectile3D::f_vx_no_drag(double t, const vector<double> &x) {
	double drag = 0;
	return -drag*x[3]-2*(_wy*x[5]-_wz*x[4]);
}

//-------------------------------------------------------------------------//

double Projectile3D::f_vy_no_drag(double t, const vector<double> &x) {
	double drag = 0;
	return -drag*x[4]-_g-2*(_wz*x[3]-_wx*x[5]);
}

//-------------------------------------------------------------------------//

double Projectile3D::f_vz_no_drag(double t, const vector<double> &x) {
	double drag = 0;
	return -drag*x[5]-2*(_wx*x[4]-_wy*x[3]);
}

//-------------------------------------------------------------------------//

double Projectile3D::f_vx_constant_drag(double t, const vector<double> &x) {
	double drag = _B2_m*sqrt(x[3]*x[3]+x[4]*x[4]+x[5]*x[5]);
	return -drag*x[3]-2*(_wy*x[5]-_wz*x[4]);
}

//-------------------------------------------------------------------------//

double Projectile3D::f_vy_constant_drag(double t, const vector<double> &x) {
	double drag = _B2_m*sqrt(x[3]*x[3]+x[4]*x[4]+x[5]*x[5]);
	return -drag*x[4]-_g-2*(_wz*x[3]-_wx*x[5]);
}

//-------------------------------------------------------------------------//

double Projectile3D::f_vz_constant_drag(double t, const vector<double> &x) {
	double drag = _B2_m*sqrt(x[3]*x[3]+x[4]*x[4]+x[5]*x[5]);
	return -drag*x[5]-2*(_wx*x[4]-_wy*x[3]);
}

//-------------------------------------------------------------------------//

double Projectile3D::f_vx_isothermal_drag(double t, const vector<double> &x) {
	double drag = _B2_m*sqrt(x[3]*x[3]+x[4]*x[4]+x[5]*x[5])*exp(-x[1]/_Y0);
	return -drag*x[3]-2*(_wy*x[5]-_wz*x[4]);
}

//-------------------------------------------------------------------------//

double Projectile3D::f_vy_isothermal_drag(double t, const vector<double> &x) {
	double drag = _B2_m*sqrt(x[3]*x[3]+x[4]*x[4]+x[5]*x[5])*exp(-x[1]/_Y0);
	return -drag*x[4]-_g-2*(_wz*x[3]-_wx*x[5]);
}

//-------------------------------------------------------------------------//

double Projectile3D::f_vz_isothermal_drag(double t, const vector<double> &x) {
	double drag = _B2_m*sqrt(x[3]*x[3]+x[4]*x[4]+x[5]*x[5])*exp(-x[1]/_Y0);
	return -drag*x[5]-2*(_wx*x[4]-_wy*x[3]);
}

//-------------------------------------------------------------------------//

double Projectile3D::f_vx_adiabatic_drag(double t, const vector<double> &x) {
	double alpha = 1/(_gamma-1);
	double drag = _B2_m*pow(_T0/_T_ref, alpha)*sqrt(x[3]*x[3]+x[4]*x[4]+x[5]*x[5])*pow(1-_a*x[1]/_T0, alpha);
	return -drag*x[3]-2*(_wy*x[5]-_wz*x[4]);
}

//-------------------------------------------------------------------------//

double Projectile3D::f_vy_adiabatic_drag(double t, const vector<double> &x) {
	double alpha = 1/(_gamma-1);
	double drag = _B2_m*pow(_T0/_T_ref, alpha)*sqrt(x[3]*x[3]+x[4]*x[4]+x[5]*x[5])*pow(1-_a*x[1]/_T0, alpha);
	return -drag*x[4]-_g-2*(_wz*x[3]-_wx*x[5]);
}

//-------------------------------------------------------------------------//

double Projectile3D::f_vz_adiabatic_drag(double t, const vector<double> &x) {
	double alpha = 1/(_gamma-1);
	double drag = _B2_m*pow(_T0/_T_ref, alpha)*sqrt(x[3]*x[3]+x[4]*x[4]+x[5]*x[5])*pow(1-_a*x[1]/_T0, alpha);
	return -drag*x[5]-2*(_wx*x[4]-_wy*x[3]);
}

//-------------------------------------------------------------------------//

bool Projectile3D::stop(double t, const vector<double> &x) {
	return x[1]<0;
}

//-------------------------------------------------------------------------//

bool Projectile3D::check() const {
	if(_mode!=0 && _mode!=1 && _mode!=2 && _mode!=3) {
		cout << "ERROR: invalid mode!" << endl;
		cout << "0: no drag; 1: constant drag; 2: isothermal drag; 3: adiabatic drag." << endl;
		return false;
	}

	if(_rk_order!=1 && _rk_order!=2 && _rk_order!=4) {
		cout << "ERROR: invalid order of RK method!" << endl;
		cout << "1: RK1 (Euler); 2: RK2; 4: RK4." << endl;
		return false;
	}

	return true;
}

//-------------------------------------------------------------------------//

void Projectile3D::cal() {
	if(!check()) return;

	RungeKutta rk(_n_eqns, _dt, _t_start, _c_start);
	rk.load_f(f_x);
	rk.load_f(f_y);
	rk.load_f(f_z);
	if(_mode == 0) {
		rk.load_f(f_vx_no_drag);
		rk.load_f(f_vy_no_drag);
		rk.load_f(f_vz_no_drag);
	} else if(_mode == 1) {
		rk.load_f(f_vx_constant_drag);
		rk.load_f(f_vy_constant_drag);
		rk.load_f(f_vz_constant_drag);
	} else if(_mode == 2) {
		rk.load_f(f_vx_isothermal_drag);
		rk.load_f(f_vy_isothermal_drag);
		rk.load_f(f_vz_isothermal_drag);
	} else if(_mode == 3) {
		rk.load_f(f_vx_adiabatic_drag);
		rk.load_f(f_vy_adiabatic_drag);
		rk.load_f(f_vz_adiabatic_drag);
	} else {
		cout << "ERROR: invalid mode!" << endl;
		cout << "0: no drag; 1: constant drag; 2: isothermal drag; 3: adiabatic drag." << endl;
	}

	rk.set_stop(stop);

	if(_rk_order == 1) {
		rk.cal_rk1();
	} else if(_rk_order == 2) {
		rk.cal_rk2();
	} else if(_rk_order == 4) {
		rk.cal_rk4();
	} else {
		cout << "ERROR: invalid order of RK method!" << endl;
		cout << "1: RK1 (Euler); 2: RK2; 4: RK4." << endl;
	}

	_t = rk.get_t();
	_x = rk.get_x();
	_n_stps = rk.get_n_stps();
}

//-------------------------------------------------------------------------//

double Projectile3D::cal_range() {
	int n = _n_stps - 1;
	//cout << "y[n-1] = " << _x[1][n-1] << "; y[n] = " << _x[1][n] << endl;
	_range = ( _x[1][n]*_x[0][n-1] - _x[1][n-1]*_x[0][n] )/( _x[1][n] - _x[1][n-1] );
	return _range;
}

//-------------------------------------------------------------------------//

double Projectile3D::search_theta_for_max_range(double v = 700) {
	double dtheta = 32;
	double theta = 45;
	double theta_tmp_left;
	double theta_tmp_right;
	//double v = 700;
	double range = -1;
	double range_tmp_left;
	double range_tmp_right;

	do {
		theta_tmp_left = theta - dtheta;
		_c_start[3] = v*cos(theta_tmp_left*M_PI/180);
		_c_start[4] = v*sin(theta_tmp_left*M_PI/180);
		cal();
		range_tmp_left = cal_range();

		theta_tmp_right = theta + dtheta;
		_c_start[3] = v*cos(theta_tmp_right*M_PI/180);
		_c_start[4] = v*sin(theta_tmp_right*M_PI/180);
		cal();
		range_tmp_right = cal_range();

		if(range_tmp_left <= range && range_tmp_right <= range) {
			dtheta /= 2.0;
			continue;
		}
		if(range_tmp_left > range) {
			range = range_tmp_left;
			theta = theta_tmp_left;
		}
		if(range_tmp_right > range) {
			range = range_tmp_right;
			theta = theta_tmp_right;
		}
	} while(dtheta>0.05);

	cout << "left uncertainty: " << range-range_tmp_left << endl;
	cout << "right uncertainty: " << range-range_tmp_right << endl;

	return theta;
}

//-------------------------------------------------------------------------//

