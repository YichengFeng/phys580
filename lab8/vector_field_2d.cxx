#include <iostream>
#include "vector_field_2d.h"

using namespace std;


//-------------------------------------------------------------------------//

TVectorField2D::TVectorField2D(TString name, int nx, double xl, double xr, int ny, double yl, double yr, const vector< vector<double> > &vx, const vector< vector<double> > &vy) {

	_name = name;

	_nx = nx;
	_xl = xl;
	_xr = xr;
	_ny = ny;
	_yl = yl;
	_yr = yr;
	_vx = vx;
	_vy = vy;

	_isCal = false;

	_frame = new TH2D("frame_"+_name, "frame_"+_name, 1,_xl,_xr, 1,_yl,_yr);

}

//-------------------------------------------------------------------------//

TVectorField2D::TVectorField2D(TString name, double xl, double xr, double yl, double yr, const vector< vector<double> > &vx, const vector< vector<double> > &vy) {

	_name = name;

	_nx = vx[0].size();
	_xl = xl;
	_xr = xr;
	_ny = vx.size();
	_yl = yl;
	_yr = yr;
	_vx = vx;
	_vy = vy;

	_isCal = false;

	_frame = new TH2D("frame_"+_name, "frame_"+_name, 1,_xl,_xr, 1,_yl,_yr);

}

//-------------------------------------------------------------------------//

TVectorField2D::~TVectorField2D() {

}

//-------------------------------------------------------------------------//

void TVectorField2D::FindMax() {

	_xmax = 0;
	_ymax = 0;

	for(int j=0; j<_ny; j++) {
		for(int i=0; i<_nx; i++) {
			if(_xmax < fabs(_vx[j][i])) _xmax = fabs(_vx[j][i]);
			if(_ymax < fabs(_vy[j][i])) _ymax = fabs(_vy[j][i]);
		}
	}

	_max = _xmax>_ymax?_xmax:_ymax;
}

//-------------------------------------------------------------------------//

bool TVectorField2D::Check() {

	bool isPass = true;

	if( _max == 0 ) {
		isPass = false;
		cout << "Zero vector field" << endl;
		return false;
	}
	if( _xmax == 0 && _ymax != 0) {
		_xmax = 1;
	}
	if( _xmax != 0 && _ymax == 0) {
		_ymax = 1;
	}

	return isPass;
}

//-------------------------------------------------------------------------//

void TVectorField2D::Cal() {

	FindMax();

	if(!Check()) {
		return;
	}

	double dx = (_xr - _xl)/_nx;
	double dy = (_yr - _yl)/_ny;
	double dr = sqrt(dx*dx + dy*dy);

	for(int j=0; j<_ny; j++) {
		vector<TArrow*> vf_tmp;
		for(int i=0; i<_nx; i++) {
			double axc = _xl + (0.5+i)*dx;
			double ayc = _yl + (0.5+j)*dy;
			double arx = _vx[j][i]/_xmax;
			double ary = _vy[j][i]/_ymax;
			double ar  = sqrt(arx*arx*dx*dx + ary*ary*dy*dy);

			double axl = axc - 0.45*arx*dx;
			double axr = axc + 0.45*arx*dx;
			double ayl = ayc - 0.45*ary*dy;
			double ayr = ayc + 0.45*ary*dy;

			TArrow *a = new TArrow(axl,ayl,axr,ayr, 0.1*ar);
			vf_tmp.push_back(a);
		}
		_vf.push_back(vf_tmp);
	}

	_isCal = true;

}

//-------------------------------------------------------------------------//

void TVectorField2D::Draw() {

	if(!_isCal) {
		Cal();
	}

	_frame->Draw();
	for(int j=0; j<_ny; j++) {
		for(int i=0; i<_nx; i++) {
			_vf[j][i]->Draw();
		}
	}

}
