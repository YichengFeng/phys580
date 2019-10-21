#ifndef RANDOMWALK2D_H
#define RANDOMWALK2D_H

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;


class RandomWalk2D
{
private:
	int _n_step;
	int _n_walk;
	int _mod;
	double _d;
	double _min_step;

	vector<double> _t;
	vector<double> _x;
	vector<double> _y;
	vector<double> _r2;
	vector<double> _r4;
	vector<double> _r2_std;
	vector<double> _good_walk;

public:
	RandomWalk2D();
	RandomWalk2D(int n_step, int n_walk, int mod);
	~RandomWalk2D();

	void set_n_step(int n_step);
	void set_n_walk(int n_walk);
	void set_mod(int mod);
	void set_d(double d);
	void set_min_step(double min_step);

	int get_n_step() const;
	int get_n_walk() const;
	int get_mod() const;
	double get_d() const;
	double get_min_step() const;
	vector<double> get_t() const;
	vector<double> get_x() const;
	vector<double> get_y() const;
	vector<double> get_r2() const;
	vector<double> get_r4() const;
	vector<double> get_r2_std() const;
	vector<double> get_good_walk() const;

	void cal_rw_1();
	void cal_rw_d();
	bool cal_saw();

	void cal();
};

inline void RandomWalk2D::set_n_step(int n_step) { _n_step = n_step; }
inline void RandomWalk2D::set_n_walk(int n_walk) { _n_walk = n_walk; }
inline void RandomWalk2D::set_mod(int mod) { _mod = mod; }
inline void RandomWalk2D::set_d(double d) { _d = d; }
inline void RandomWalk2D::set_min_step(double min_step) { _min_step = min_step; }

inline int RandomWalk2D::get_n_step() const { return _n_step; }
inline int RandomWalk2D::get_n_walk() const { return _n_walk; }
inline int RandomWalk2D::get_mod() const { return _mod; }
inline double RandomWalk2D::get_d() const { return _d; }
inline double RandomWalk2D::get_min_step() const { return _min_step; }

inline vector<double> RandomWalk2D::get_t() const { return _t; }
inline vector<double> RandomWalk2D::get_x() const { return _x; }
inline vector<double> RandomWalk2D::get_y() const { return _y; }
inline vector<double> RandomWalk2D::get_r2() const { return _r2; }
inline vector<double> RandomWalk2D::get_r4() const { return _r4; }
inline vector<double> RandomWalk2D::get_r2_std() const { return _r2_std; }
inline vector<double> RandomWalk2D::get_good_walk() const { return _good_walk; }


#endif
