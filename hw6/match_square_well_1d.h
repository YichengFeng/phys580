#ifndef MATCHSQUAREWELL1D_H
#define MATCHSQUAREWELL1D_H

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;


class MatchSquareWell1D
{
private:
	// input
	double _dx;
	double _E;
	double _dE;
	double _fr;
	double _T;
	double _Vmax;
	double _a;
	double _V2;
	int _Nmax;

	// middle
	double _L;
	double _x_match;
	double _x_left;
	double _x_right;
	double _y_min;
	double _y_max;
	double _r;
	double _devdiff;

	// output
	vector<double> _psiL;
	vector<double> _psiR;
	int _matchL;
	int _matchR;
	vector<double> _xL;
	vector<double> _xR;

	// variation
	double _newE;
	double _varE;
	double _dpsi;
	int _n_total;
	int _n_moves;
	int _iend;
	vector<double> _x;
	vector<double> _potential_table;
	vector<double> _psi;
	vector<double> _psi_old;

public:
	MatchSquareWell1D();
	MatchSquareWell1D(double Vmax, double a, double V2, double dx, double E, double dE);

	void set_dx(double dx);
	void set_E(double E);
	void set_dE(double dE);
	void set_fr(double fr);
	void set_T(double T);
	void set_Nmax(int Nmax);

	double get_dx() const;
	double get_E() const;
	double get_dE() const;
	double get_varE() const;

	vector<double> get_psiL() const;
	vector<double> get_psiR() const;
	vector<double> get_xL() const;
	vector<double> get_xR() const;

	vector<double> get_x() const;
	vector<double> get_potential_table() const;
	vector<double> get_psi() const;

	double potential(double x);
	void cal_potential_table();

	bool check() const;
	void cal_match_once();
	void adjust_match_once();
	void adjust_match_until(double ee);

	void var_psi_normalize();
	double var_energy();
	bool cal_var_once();
	int cal_var_attempts(int n_attempts);
};


inline void MatchSquareWell1D::set_dx(double dx) { _dx = dx; }
inline void MatchSquareWell1D::set_E(double E) { _E = E; }
inline void MatchSquareWell1D::set_dE(double dE) { _dE = dE; }
inline void MatchSquareWell1D::set_fr(double fr) { _fr = fr; }
inline void MatchSquareWell1D::set_T(double T) { _T = T; }
inline void MatchSquareWell1D::set_Nmax(int Nmax) { _Nmax = Nmax; }

inline double MatchSquareWell1D::get_dx() const { return _dx; }
inline double MatchSquareWell1D::get_E() const { return _E; }
inline double MatchSquareWell1D::get_dE() const { return _dE; }
inline double MatchSquareWell1D::get_varE() const { return _varE; }

inline vector<double> MatchSquareWell1D::get_psiL() const { return _psiL; } 
inline vector<double> MatchSquareWell1D::get_psiR() const { return _psiR; } 
inline vector<double> MatchSquareWell1D::get_xL() const { return _xL; } 
inline vector<double> MatchSquareWell1D::get_xR() const { return _xR; } 

inline vector<double> MatchSquareWell1D::get_x() const { return _x; } 
inline vector<double> MatchSquareWell1D::get_potential_table() const { return _potential_table; } 
inline vector<double> MatchSquareWell1D::get_psi() const { return _psi; } 


#endif
