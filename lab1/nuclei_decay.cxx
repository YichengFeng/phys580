#include <iostream>
#include <TROOT.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMath.h>
#include <TAxis.h>

#include "nuclei_decay.h"

using namespace std;


//-------------------------------------------------------------------------//
NucleiDecay::NucleiDecay() {
	_dt = 0.05;
	_t_start = 0;
	_t_end = 5;
	_N0 = 1000;
}


//-------------------------------------------------------------------------//
NucleiDecay::NucleiDecay(double dt) {
	_dt = dt;
	_t_start = 0;
	_t_end = 5;
	_N0 = 1000;
}


//-------------------------------------------------------------------------//
NucleiDecay::NucleiDecay(double dt, double t_start, double t_end, double N0) {
	_dt = dt;
	_t_start = t_start;
	_t_end = t_end;
	_N0 = N0;
}


//-------------------------------------------------------------------------//
NucleiDecay::~NucleiDecay() {

}


//-------------------------------------------------------------------------//
double NucleiDecay::fxt(double x, double t) {
	return -x;
}


//-------------------------------------------------------------------------//
void NucleiDecay::cal_g_exact() {

	int nt = (int)((_t_end - _t_start)/_dt)+1;
	double 	*N_exact = new double[nt];
	double 	*t = new double[nt];

	N_exact[0] = _N0;
	t[0] = _t_start;

	for(int i=1; i<nt; i++) {
		t[i] = t[i-1] + _dt;
		N_exact[i] = _N0*TMath::Exp(-t[i]);
	}

	_g_exact = new TGraph(nt, t, N_exact);

	return;
}


//-------------------------------------------------------------------------//
void NucleiDecay::cal_g_euler() {

	int nt = (int)((_t_end - _t_start)/_dt)+1;
	double 	*N_euler = new double[nt];
	double 	*t = new double[nt];

	N_euler[0] = _N0;
	t[0] = _t_start;

	for(int i=1; i<nt; i++) {
		t[i] = t[i-1] + _dt;
		N_euler[i] = N_euler[i-1] + fxt(N_euler[i-1], t[i-1])*_dt;
	}

	_g_euler = new TGraph(nt, t, N_euler);

	return;
}


//-------------------------------------------------------------------------//
TGraph* NucleiDecay::cal_g_error(TGraph *g) {

	if(!g || !_g_exact) {
		cout << "g or _g_exact empty!" << endl;
		return nullptr;
	}

	int nt = g->GetN();
	double *t = g->GetX();
	double *N = g->GetY();

	double *t_exact = _g_exact->GetX();
	double *N_exact = _g_exact->GetY();

	double *N_error = new double[nt];

	for(int i=0; i<nt; i++) {
		if(t[i] != t_exact[i]) {
			cout << "time not match!" << endl;
			return nullptr;
		}
	}

	for(int i=0; i<nt; i++) {
		N_error[i] = fabs( N[i] - N_exact[i] );
	}

	TGraph *g_error = new TGraph(nt, t, N_error);

	return g_error;
}


//-------------------------------------------------------------------------//
TGraph* NucleiDecay::cal_g_relative_error(TGraph *g) {

	if(!g || !_g_exact) {
		cout << "g or _g_exact empty!" << endl;
		return nullptr;
	}

	int nt = g->GetN();
	double *t = g->GetX();
	double *N = g->GetY();

	double *t_exact = _g_exact->GetX();
	double *N_exact = _g_exact->GetY();

	double *N_relative_error = new double[nt];

	for(int i=0; i<nt; i++) {
		if(t[i] != t_exact[i]) {
			cout << "time not match!" << endl;
			return nullptr;
		}
	}

	for(int i=0; i<nt; i++) {
		N_relative_error[i] = fabs( N[i] - N_exact[i] )/N_exact[i];
	}

	TGraph *g_relative_error = new TGraph(nt, t, N_relative_error);

	return g_relative_error;
}


//-------------------------------------------------------------------------//
void NucleiDecay::cal_g_rk2() {

	int nt = (int)((_t_end - _t_start)/_dt)+1;
	double 	*N_rk2 = new double[nt];
	double 	*t = new double[nt];

	N_rk2[0] = _N0;
	t[0] = _t_start;

	for(int i=1; i<nt; i++) {
		t[i] = t[i-1] + _dt;

		double F1 = fxt(N_rk2[i-1], t[i-1]);
		double F2 = fxt(N_rk2[i-1]+F1*0.5*_dt, t[i-1]+0.5*_dt);
		N_rk2[i] = N_rk2[i-1] + F2*_dt;
	}

	_g_rk2 = new TGraph(nt, t, N_rk2);

	return;
}


//-------------------------------------------------------------------------//
void NucleiDecay::cal_g_rk4() {

	int nt = (int)((_t_end - _t_start)/_dt)+1;
	double 	*N_rk4 = new double[nt];
	double 	*t = new double[nt];

	N_rk4[0] = _N0;
	t[0] = _t_start;

	for(int i=1; i<nt; i++) {
		t[i] = t[i-1] + _dt;

		double F1 = fxt(N_rk4[i-1], t[i-1]);
		double F2 = fxt(N_rk4[i-1]+F1*0.5*_dt, t[i-1]+0.5*_dt);
		double F3 = fxt(N_rk4[i-1]+F2*0.5*_dt, t[i-1]+0.5*_dt);
		double F4 = fxt(N_rk4[i-1]+F3*_dt, t[i-1]+_dt);
		N_rk4[i] = N_rk4[i-1] + ( F1/6 + F2/3 + F3/3 + F4/6 )*_dt;
	}

	_g_rk4 = new TGraph(nt, t, N_rk4);

	return;
}


//-------------------------------------------------------------------------//
void NucleiDecay::cal_g() {

	cal_g_exact();
	cal_g_euler();
	cal_g_rk2();
	cal_g_rk4();
	_g_error = cal_g_error(_g_euler);
	_g_relative_error = cal_g_relative_error(_g_euler);
	_g_rk2_error = cal_g_error(_g_rk2);
	_g_rk2_relative_error = cal_g_relative_error(_g_rk2);
	_g_rk4_error = cal_g_error(_g_rk4);
	_g_rk4_relative_error = cal_g_relative_error(_g_rk4);

	return;
}
