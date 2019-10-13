#ifndef CAPACITOR2D_H
#define CAPACITOR2D_H

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;


class Capacitor2D
{
private:
	// input
	double _a; // plates width
	double _b; // gap between the two plates
	double _d; // grid size
	double _acc; // average error
	double _alpha; // SOR
	int _alg; // 0: Jacobi; 1: Gauss-Seidel; 2: SOR
	int _periodic;

	// middle
	int _M;
	int _N;
	int _aM; // plate width in grids
	int _bM; // plate gap in grids
	vector< vector<double> > _rho; // charge density

	// output
	vector< vector<double> > _V;
	double _tmp_acc; // actual average
	int _n_iter;

public:
	Capacitor2D();
	Capacitor2D(double a, double b, double d, double acc);
	~Capacitor2D();

	void set_a(double a);
	void set_b(double b);
	void set_d(double d);
	void set_acc(double acc);
	void set_alpha(double alpha);
	void set_alg(int alg);
	void set_periodic(int periodic);
	void set_rho(double x1, double x2, double y1, double y2, double r);

	double get_a() const;
	double get_b() const;
	double get_d() const;
	double get_acc() const;
	double get_alpha() const;
	double get_tmp_acc() const;
	int get_alg() const;
	int get_periodic() const;
	int get_n_iter() const;
	vector< vector<double> > get_V() const;

	bool check() const;
	// calculations: Jacobi, Gauss-Seidel, Succussive Over R. (SOR)
	void cal_once_Jacobi();
	void cal_once_Jacobi_periodic();
	void cal_once_GaussSeidel();
	void cal_once_GaussSeidel_periodic();
	void cal_once_SOR();
	void cal_once_SOR_periodic();
	void cal();

};


inline void Capacitor2D::set_a(double a) { _a = a; }
inline void Capacitor2D::set_b(double b) { _b = b; }
inline void Capacitor2D::set_d(double d) { _d = d; }
inline void Capacitor2D::set_acc(double acc) { _acc = acc; }
inline void Capacitor2D::set_alpha(double alpha) { _alpha = alpha; }
inline void Capacitor2D::set_alg(int alg) { _alg = alg; }
inline void Capacitor2D::set_periodic(int periodic) { _periodic = periodic; }

inline double Capacitor2D::get_a() const { return _a; }
inline double Capacitor2D::get_b() const { return _b; }
inline double Capacitor2D::get_d() const { return _d; }
inline double Capacitor2D::get_acc() const { return _acc; }
inline double Capacitor2D::get_alpha() const { return _alpha; }
inline double Capacitor2D::get_tmp_acc() const { return _tmp_acc; }
inline int Capacitor2D::get_alg() const { return _alg; }
inline int Capacitor2D::get_periodic() const { return _periodic; }
inline int Capacitor2D::get_n_iter() const { return _n_iter; }

inline vector< vector<double> > Capacitor2D::get_V() const { return _V; }


#endif
