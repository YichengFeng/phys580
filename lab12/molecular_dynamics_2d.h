#ifndef MOLECULARDYNAMICS2D_H
#define MOLECULARDYNAMICS2D_H

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;


class MolecularDynamics2D
{
private:
	int _N;
	int _rec;
	double _L;
	double _dt;
	double _dmax;
	double _vmax;

	double _t_now;
	vector<double> _x_start;
	vector<double> _x_old;
	vector<double> _x_now;
	vector<double> _x_new;
	vector<double> _y_start;
	vector<double> _y_old;
	vector<double> _y_now;
	vector<double> _y_new;
	vector<double> _vx_now;
	vector<double> _vy_now;

	double _tag_displace_now;
	double _tag_distance_now;
	vector<double> _tag_displace;
	vector<double> _tag_distance;

	double _T_now;
	double _E_now;
	vector<double> _T;
	vector<double> _E;

	vector<double> _t;
	vector< vector<double> > _x;
	vector< vector<double> > _y;
	vector< vector<double> > _vx;
	vector< vector<double> > _vy;

public:
	MolecularDynamics2D();
	MolecularDynamics2D(int N, double L, double dt, double dmax, double vmax);
	~MolecularDynamics2D();

	void set_rec(int rec);

	int get_N() const;
	int get_rec() const;
	double get_L() const;
	double get_dmax() const;
	double get_vmax() const;

	vector<double> get_T() const;
	vector<double> get_E() const;
	vector<double> get_tag_displace() const;
	vector<double> get_tag_distance() const;
	vector<double> get_t() const;
	vector< vector<double> > get_x() const;
	vector< vector<double> > get_y() const;
	vector< vector<double> > get_vx() const;
	vector< vector<double> > get_vy() const;

	bool check() const;

	void cal_once();
	void cal_until(double t_end);
	void cal_ET();

	void change_velocity(double factor);
};


inline void MolecularDynamics2D::set_rec(int rec) { _rec = rec; }

inline int MolecularDynamics2D::get_N() const { return _N; }
inline int MolecularDynamics2D::get_rec() const { return _rec; }
inline double MolecularDynamics2D::get_L() const { return _L; }
inline double MolecularDynamics2D::get_dmax() const { return _dmax; }
inline double MolecularDynamics2D::get_vmax() const { return _vmax; }

inline vector<double> MolecularDynamics2D::get_T() const { return _T; }
inline vector<double> MolecularDynamics2D::get_E() const { return _E; }
inline vector<double> MolecularDynamics2D::get_tag_displace() const { return _tag_displace; }
inline vector<double> MolecularDynamics2D::get_tag_distance() const { return _tag_distance; }
inline vector<double> MolecularDynamics2D::get_t() const { return _t; }

inline vector< vector<double> > MolecularDynamics2D::get_x() const { return _x; }
inline vector< vector<double> > MolecularDynamics2D::get_y() const { return _y; }
inline vector< vector<double> > MolecularDynamics2D::get_vx() const { return _vx; }
inline vector< vector<double> > MolecularDynamics2D::get_vy() const { return _vy; }


#endif
