#ifndef VECTORFIELD2D_H
#define VECTORFIELD2D_H

#include <TCanvas.h>
#include <TH2.h>
#include <TArrow.h>
#include <TAxis.h>
#include <TString.h>

using namespace std;


class TVectorField2D
{
private:
	// input
	TString _name;

	int _nx;
	double _xl;
	double _xr;
	int _ny;
	double _yl;
	double _yr;

	vector< vector<double> > _vx;
	vector< vector<double> > _vy;

	// middle
	double _xmax;
	double _ymax;
	double _max;
	bool _isCal;

	// output
	TH2D *_frame;
	vector< vector<TArrow*> > _vf;

	bool Check();
	void FindMax();
	void Cal();

public:
	TVectorField2D(TString name, int nx, double xl, double xr, int ny, double yl, double yr, const vector< vector<double> > &vx, const vector< vector<double> > &vy);
	TVectorField2D(TString name, double xl, double xr, double yl, double yr, const vector< vector<double> > &vx, const vector< vector<double> > &vy);
	~TVectorField2D();

	TAxis *GetXaxis();
	TAxis *GetYaxis();
	void Draw();

};


inline TAxis *TVectorField2D::GetXaxis() { return _frame->GetXaxis(); }
inline TAxis *TVectorField2D::GetYaxis() { return _frame->GetYaxis(); }


#endif
