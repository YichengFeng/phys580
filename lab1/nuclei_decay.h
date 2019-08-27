#ifndef NUCLEIDECAY_H
#define NUCLEIDECAY_H

#include <iostream>
#include <TROOT.h>
#include <TGraph.h>

using namespace std;


class NucleiDecay
{
private:
	// i/o
	double _dt;
	double _t_start;
	double _t_end;
	double _N0;

	TGraph *_g_euler;
	TGraph *_g_exact;
	TGraph *_g_error;
	TGraph *_g_relative_error;
	TGraph *_g_rk2;
	TGraph *_g_rk2_error;
	TGraph *_g_rk2_relative_error;
	TGraph *_g_rk4;
	TGraph *_g_rk4_error;
	TGraph *_g_rk4_relative_error;

public:
	NucleiDecay();
	NucleiDecay(double dt);
	NucleiDecay(double dt, double t_start, double t_end, double N0);

	~NucleiDecay();

	double fxt(double x, double t);

	void set_dt(double dt);
	void set_t_start(double t_start);
	void set_t_end(double t_end);
	void set_N0(double N0);

	double get_dt() const;
	double get_t_start() const;
	double get_t_end() const;
	double get_N0() const;

	void cal_g_exact();
	void cal_g_euler();
	void cal_g_rk2();
	void cal_g_rk4();
	TGraph* cal_g_error(TGraph *g);
	TGraph* cal_g_relative_error(TGraph *g);
	void cal_g();

	TGraph *get_g_euler() const;
	TGraph *get_g_exact() const;
	TGraph *get_g_error() const;
	TGraph *get_g_relative_error() const;
	TGraph *get_g_rk2() const;
	TGraph *get_g_rk2_error() const;
	TGraph *get_g_rk2_relative_error() const;
	TGraph *get_g_rk4() const;
	TGraph *get_g_rk4_error() const;
	TGraph *get_g_rk4_relative_error() const;
};


inline void NucleiDecay::set_dt(double dt) { _dt = dt; }
inline void NucleiDecay::set_t_start(double t_start) { _t_start = t_start; }
inline void NucleiDecay::set_t_end(double t_end) { _t_end = t_end; }
inline void NucleiDecay::set_N0(double N0) { _N0 = N0; }

inline double NucleiDecay::get_dt() const { return _dt; }
inline double NucleiDecay::get_t_start() const { return _t_start; }
inline double NucleiDecay::get_t_end() const { return _t_end; }
inline double NucleiDecay::get_N0() const { return _N0; }

inline TGraph* NucleiDecay::get_g_euler() const { return _g_euler; }
inline TGraph* NucleiDecay::get_g_exact() const { return _g_exact; }
inline TGraph* NucleiDecay::get_g_error() const { return _g_error; }
inline TGraph* NucleiDecay::get_g_relative_error() const { return _g_relative_error; }
inline TGraph* NucleiDecay::get_g_rk2() const { return _g_rk2; }
inline TGraph* NucleiDecay::get_g_rk2_error() const { return _g_rk2_error; }
inline TGraph* NucleiDecay::get_g_rk2_relative_error() const { return _g_rk2_relative_error; }
inline TGraph* NucleiDecay::get_g_rk4() const { return _g_rk4; }
inline TGraph* NucleiDecay::get_g_rk4_error() const { return _g_rk4_error; }
inline TGraph* NucleiDecay::get_g_rk4_relative_error() const { return _g_rk4_relative_error; }


#endif
