#ifndef LOOPMAGNETICFIELD_H
#define LOOPMAGNETICFIELD_H

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;


class LoopMagneticField
{
private:
	// input
	int _n_loop;
	int _n_phi;
	double _R;
	double _pitch;
	double _d;
	int _mod; // 0: discrete loops; 1: helical coil

	// middle
	int _M; // half of the grid number in each cordinate.
	int _N;
	vector<double> _z_loop;
	vector<double> _I_loop;

	// output
	vector< vector<double> > _bx;
	vector< vector<double> > _by;
	vector< vector<double> > _bz;

public:
	LoopMagneticField();
	LoopMagneticField(int n_phi, double R, int n_loop, double pitch, double d);
	~LoopMagneticField();

	void set_n_loop(int n_loop);
	void set_n_phi(int n_phi);
	void set_R(double r);
	void set_pitch(double pitch);
	void set_d(double d);
	void set_M(int M);
	void set_mod(int mod);
	void set_I_loop(int i, double I);

	int get_n_loop() const;
	int get_n_phi() const;
	double get_R() const;
	double get_pitch() const;
	double get_d() const;
	int get_M() const;
	int get_mod() const;
	vector< vector<double> > get_bx() const;
	vector< vector<double> > get_by() const;
	vector< vector<double> > get_bz() const;

	bool check() const;
	vector<double> cal_one_grid(double x, double y, double z);
	vector<double> cal_one_grid_helical(double x, double y, double z);
	void cal();
};


inline void LoopMagneticField::set_n_loop(int n_loop) { _n_loop = n_loop; }
inline void LoopMagneticField::set_n_phi(int n_phi) { _n_phi = n_phi; }
inline void LoopMagneticField::set_R(double R) { _R = R; }
inline void LoopMagneticField::set_pitch(double pitch) { _pitch = pitch; }
inline void LoopMagneticField::set_d(double d) { _d = d; }
inline void LoopMagneticField::set_M(int M) { _M = M; }
inline void LoopMagneticField::set_mod(int mod) { _mod = mod; }
inline void LoopMagneticField::set_I_loop(int i, double I) { _I_loop[i] = I; }

inline int LoopMagneticField::get_n_loop() const { return _n_loop; }
inline int LoopMagneticField::get_n_phi() const { return _n_phi; }
inline double LoopMagneticField::get_R() const { return _R; }
inline double LoopMagneticField::get_pitch() const { return _pitch; }
inline double LoopMagneticField::get_d() const { return _d; }
inline int LoopMagneticField::get_M() const { return _M; }
inline int LoopMagneticField::get_mod() const { return _mod; }
inline vector< vector<double> > LoopMagneticField::get_bx() const { return _bx; }
inline vector< vector<double> > LoopMagneticField::get_by() const { return _by; }
inline vector< vector<double> > LoopMagneticField::get_bz() const { return _bz; }

#endif
