#ifndef ISINGMODEL2D_H
#define ISINGMODEL2D_H

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;


class IsingModel2D
{
private:
	// input
	int _N; // lattice size
	double _T; // temperature (J/k_B)
	double _H; // magnetic field (J/mu)
	int _MCS; // Monte Carlo steps per spin
	int _c; // initial spin configuration (1: all up; -1: all down; 2: random)

	// middle
	int _now;

	// output
	double _m; // magnetization per spin <s>
	double _E; // energy per spin
	vector< vector<int> > _spin; // spin 2d distribution

public:
	IsingModel2D();
	IsingModel2D(int N, double T, double H, int MCS, int c);
	~IsingModel2D();

	void set_N(int N);
	void set_T(double T);
	void set_H(double H);
	void set_MCS(int MCS);
	void set_c(int c);

	int get_N() const;
	double get_T() const;
	double get_H() const;
	int get_MCS() const;
	int get_c() const;

	double get_m() const;
	double get_E() const;
	vector< vector<int> > get_spin() const;

	bool check() const;
	void cal_once();
	void cal_until(int t);
	void cal_all();

};


inline void IsingModel2D::set_N(int N) { _N = N; }
inline void IsingModel2D::set_T(double T) { _T = T; }
inline void IsingModel2D::set_H(double H) { _H = H; }
inline void IsingModel2D::set_MCS(int MCS) { _MCS = MCS; }
inline void IsingModel2D::set_c(int c) { _c = c; }
 
inline int IsingModel2D::get_N() const { return _N; }
inline double IsingModel2D::get_T() const { return _T; }
inline double IsingModel2D::get_H() const { return _H; }
inline int IsingModel2D::get_MCS() const { return _MCS; }
inline int IsingModel2D::get_c() const { return _c; }

inline double IsingModel2D::get_m() const { return _m; }
inline double IsingModel2D::get_E() const { return _E; }
inline vector< vector<int> > IsingModel2D::get_spin() const { return _spin; }

#endif
