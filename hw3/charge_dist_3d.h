#ifndef CHARGEDIST3D_H
#define CHARGEDIST3D_H

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;


class ChargeDist3D
{
private:
	// input
	double _d; // grid size
	double _acc; // average error
	double _alpha; // SOR
	int _alg; // 0: Jacobi; 1: Gauss-Seidel; 2: SOR
	int _periodic;

	// middle
	int _M;
	int _N;
	vector< vector< vector<double> > > _rho; // charge density

	// output
	vector< vector< vector<double> > > _V;
	double _tmp_acc; // actual average
	int _n_iter;

public:
	ChargeDist3D();
	ChargeDist3D(double d, double acc);
	~ChargeDist3D();

	void set_d(double d);
	void set_acc(double acc);
	void set_alpha(double alpha);
	void set_alg(int alg);
	void set_periodic(int periodic);
	void set_rho(double x1, double x2, double y1, double y2, double z1, double z2, double r);

	double get_d() const;
	double get_acc() const;
	double get_alpha() const;
	double get_tmp_acc() const;
	int get_alg() const;
	int get_periodic() const;
	int get_n_iter() const;
	vector< vector< vector<double> > > get_V() const;

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


inline void ChargeDist3D::set_d(double d) { _d = d; }
inline void ChargeDist3D::set_acc(double acc) { _acc = acc; }
inline void ChargeDist3D::set_alpha(double alpha) { _alpha = alpha; }
inline void ChargeDist3D::set_alg(int alg) { _alg = alg; }
inline void ChargeDist3D::set_periodic(int periodic) { _periodic = periodic; }

inline double ChargeDist3D::get_d() const { return _d; }
inline double ChargeDist3D::get_acc() const { return _acc; }
inline double ChargeDist3D::get_alpha() const { return _alpha; }
inline double ChargeDist3D::get_tmp_acc() const { return _tmp_acc; }
inline int ChargeDist3D::get_alg() const { return _alg; }
inline int ChargeDist3D::get_periodic() const { return _periodic; }
inline int ChargeDist3D::get_n_iter() const { return _n_iter; }

inline vector< vector< vector<double> > > ChargeDist3D::get_V() const { return _V; }


#endif
