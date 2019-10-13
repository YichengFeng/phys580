#include "nintegrate_1d.h"

using namespace std;


//-------------------------------------------------------------------------//

NIntegrate1D::NIntegrate1D() {

}

//-------------------------------------------------------------------------//

NIntegrate1D::NIntegrate1D(double x1, double x2, double (*f)(double x)) {

	_x1 = x1;
	_x2 = x2;
	_f  = f;
	_n  = 100;
	_alg = 0;
}

//-------------------------------------------------------------------------//

void NIntegrate1D::cal_trapzoidal() {

	double dx = (_x2 - _x1)/_n;
	double wt;

	_output = 0;

	for(int i=0; i<=_n; i++) {
		if(i==0 || i==_n) {
			wt = 1;
		} else {
			wt = 2;
		}

		double x = _x1 + i*dx;
		_output += wt*_f(x);
	}

	_output *= 0.5*dx;
}

//-------------------------------------------------------------------------//

void NIntegrate1D::cal_Simpson() {

	if(_n%2 == 1) _n ++;

	double dx = (_x2 - _x1)/_n;
	double wt;

	_output = 0;

	for(int i=0; i<=_n; i++) {
		if(i==0 || i==_n) {
			wt = 1;
		} else {
			wt = i%2==0?2:4;
		}

		double x = _x1 + i*dx;
		_output += wt*_f(x);
	}

	_output *= dx/3.0;
}

//-------------------------------------------------------------------------//

void NIntegrate1D::cal_Romberg() {

	_n = int(pow(2, log2(1.0*_n)) + 1);

}

//-------------------------------------------------------------------------//

double NIntegrate1D::cal() {

	if(_alg == 0) {
		cal_trapzoidal();
	}
	if(_alg == 1) {
		cal_Simpson();
	}
	if(_alg == 2) {
		cal_Romberg();
	}

	return _output;
}
