#ifndef DIFFUSION2D_H
#define DIFFUSION2D_H

#include <iostream>
#include <vector>
#include <cmath>
#include <string>

using namespace std;


class Diffusion2D
{
private:
	string _name;

	int _M; // half edge
	int _N; // full edge
	int _P; // number of particles
	int _P0; // initial number of particles
	int _D; // edge division for entropy
	int _L; // half leak length
	double _S; // entropy

	int _mod; // 0:normal; 1:leak
	int _mod_sample; // 0:equal times; 1:equal diff;

	int _n_t; // number of sample time
	double _t_max; // total evolution time
	double _t;
	vector<double> _x;
	vector<double> _y;

	vector<double> _t_sample;
	vector< vector<double> > _x_sample;
	vector< vector<double> > _y_sample;
	vector<double> _S_sample;

public:
	Diffusion2D();
	Diffusion2D(string name, int M, int P, int D, int n_t, double t_max);
	~Diffusion2D();

	void set_name(string name);
	void set_M(int M);
	void set_P(int P);
	void set_D(int D);
	void set_L(int L);
	void set_mod(int mod);
	void set_mod_sample(int mod_sample);
	void set_n_t(int n_t);
	void set_t_max(double t_max);

	string get_name() const;
	int get_M() const;
	int get_N() const;
	int get_P() const;
	int get_D() const;
	int get_L() const;
	int get_mod() const;
	int get_mod_sample() const;
	int get_n_t() const;
	double get_t_max() const;
	double get_S() const;

	vector<double> get_t_sample() const;
	vector< vector<double> > get_x_sample() const;
	vector< vector<double> > get_y_sample() const;
	vector<double> get_S_sample() const;

	double get_t() const;
	vector<double> get_x() const;
	vector<double> get_y() const;

	bool check() const;
	void cal_all();
	void cal_until(double t_stop);
	void cal_leak_until(double t_stop);
	double cal_entropy();
};


inline void Diffusion2D::set_name(string name) { _name = name; }
inline void Diffusion2D::set_M(int M) { _M = M; }
inline void Diffusion2D::set_P(int P) { _P = P; }
inline void Diffusion2D::set_D(int D) { _D = D; }
inline void Diffusion2D::set_L(int L) { _L = L; }
inline void Diffusion2D::set_mod(int mod) { _mod = mod; }
inline void Diffusion2D::set_mod_sample(int mod_sample) { _mod_sample = mod_sample; }
inline void Diffusion2D::set_n_t(int n_t) { _n_t = n_t; }
inline void Diffusion2D::set_t_max(double t_max) { _t_max = t_max; }

inline string Diffusion2D::get_name() const { return _name; }
inline int Diffusion2D::get_M() const { return _M; }
inline int Diffusion2D::get_N() const { return _N; }
inline int Diffusion2D::get_P() const { return _P; }
inline int Diffusion2D::get_D() const { return _D; }
inline int Diffusion2D::get_L() const { return _L; }
inline int Diffusion2D::get_mod() const { return _mod; }
inline int Diffusion2D::get_mod_sample() const { return _mod_sample; }
inline int Diffusion2D::get_n_t() const { return _n_t; }
inline double Diffusion2D::get_t_max() const { return _t_max; }
inline double Diffusion2D::get_S() const { return _S; }

inline vector<double> Diffusion2D::get_t_sample() const { return _t_sample; }
inline vector< vector<double> > Diffusion2D::get_x_sample() const { return _x_sample; }
inline vector< vector<double> > Diffusion2D::get_y_sample() const { return _y_sample; }
inline vector<double> Diffusion2D::get_S_sample() const { return _S_sample; }

inline double Diffusion2D::get_t() const { return _t; }
inline vector<double> Diffusion2D::get_x() const { return _x; }
inline vector<double> Diffusion2D::get_y() const { return _y; }


#endif
