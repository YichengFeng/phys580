#ifndef SQUAREWELL1D_H
#define SQUAREWELL1D_H

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;


class SquareWell1D
{
private:
	const double _L = 1;

	// input
	int _parity; // even: +1; odd: -1;
	double _dx;
	double _Vmax; // the height of the well, should be large
	double _E;
	double _dE;
	int _iter;
	double _boundary; // boundary where psi = 0
	// barrier
	double _a; // barrier edge, 0 means no barrier
	double _V2;
	double _b; // y range scaling [-b,b]

	// middle
	double _r;
	double _xmin;
	double _xmax;
	double _ymin;
	double _ymax;
	int _nhalf;

	double _psiAt1;

	vector<double> _xL;
	vector<double> _xR;
	vector<double> _psiL;
	vector<double> _psiR;

	// output
	vector<double> _x;
	vector<double> _psi;
	vector<double> _potential_table;

public:
	SquareWell1D();
	SquareWell1D(int parity, double dx, double Vmax, double a, double V2, double b, double E, double dE);

	void set_parity(int parity);
	void set_iter(int iter);
	void set_boundary(double boundary);
	void set_dx(double dx);
	void set_Vmax(double Vmax);
	void set_E(double E);
	void set_dE(double dE);
	void set_a(double a);
	void set_V2(double V2);
	void set_b(double b);

	int get_parity() const;
	double get_boundary() const;
	double get_dx() const;
	double get_Vmax() const;
	double get_E() const;
	double get_dE() const;
	double get_a() const;
	double get_V2() const;
	double get_b() const;

	vector<double> get_x() const;
	vector<double> get_psi() const;
	vector<double> get_potential_table() const;

	bool check() const;
	double potential(double x);
	void cal_once();
	void adjust_once();
	void adjust_until(double ee);
	void cal_potential_table();
};


inline void SquareWell1D::set_parity(int parity) { _parity = parity; }
inline void SquareWell1D::set_iter(int iter) { _iter = iter; }
inline void SquareWell1D::set_boundary(double boundary) { _boundary = boundary; }
inline void SquareWell1D::set_dx(double dx) { _dx = dx; }
inline void SquareWell1D::set_Vmax(double Vmax) { _Vmax = Vmax; }
inline void SquareWell1D::set_E(double E) { _E = E; }
inline void SquareWell1D::set_dE(double dE) { _dE = dE; }
inline void SquareWell1D::set_a(double a) { _a = a; }
inline void SquareWell1D::set_V2(double V2) { _V2 = V2; }
inline void SquareWell1D::set_b(double b) { _b = b; }

inline int SquareWell1D::get_parity() const { return _parity; }
inline double SquareWell1D::get_boundary() const { return _boundary; }
inline double SquareWell1D::get_dx() const { return _dx; }
inline double SquareWell1D::get_Vmax() const { return _Vmax; }
inline double SquareWell1D::get_E() const { return _E; }
inline double SquareWell1D::get_dE() const { return _dE; }
inline double SquareWell1D::get_a() const { return _a; }
inline double SquareWell1D::get_V2() const { return _V2; }
inline double SquareWell1D::get_b() const { return _b; }

inline vector<double> SquareWell1D::get_x() const { return _x; }
inline vector<double> SquareWell1D::get_psi() const { return _psi; }
inline vector<double> SquareWell1D::get_potential_table() const { return _potential_table; }


#endif
