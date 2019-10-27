#ifndef PERCOLATION2D_H
#define PERCOLATION2D_H

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;


class Percolation2D
{
private:
	double _p;
	int _L; // lattice edge
	int _trial; // number of trials

	// for each site
	int _n_occupied;
	vector< vector<bool> > _occupied;
	vector< vector<int> > _cluster_id;
	// for each cluster
	int _n_cluster;
	vector<int> _cluster_size;
	vector< vector<int> > _cluster_x;
	vector< vector<int> > _cluster_y;
	vector<int> _spanning_cluster_id;
	// for each trial
	double _P;
	double _S;
	vector<double> _vP;
	vector<double> _vS;

	// for average
	double _P_average;
	double _S_average;

public:
	Percolation2D();
	Percolation2D(double p, int L, int trial);

	void set_p(double p);
	void set_L(int L);
	void set_trial(int trial);

	double get_p() const;
	int get_L() const;
	int get_trial() const;

	// for each site
	int get_n_occupied() const;
	vector< vector<bool> > get_occupied() const;
	vector< vector<int> > get_cluster_id() const;
	// for each cluster
	int get_n_cluster() const;
	vector<int> get_cluster_size() const;
	vector< vector<int> > get_cluster_x() const;
	vector< vector<int> > get_cluster_y() const;
	vector<int> get_spanning_cluster_id() const;
	// for average
	double get_P_average() const;
	double get_S_average() const;

	void cal_perculation();
	void cal_DFS(int ix, int iy, int cid);
	void cal_cluster();
	bool check();
	void cal_PS();
	void cal_one_trial();
	void cal_trials();
};

inline void Percolation2D::set_p(double p) { _p = p; }
inline void Percolation2D::set_L(int L) { _L = L; }
inline void Percolation2D::set_trial(int trial) { _trial = trial; }

inline double Percolation2D::get_p() const { return _p; }
inline int Percolation2D::get_L() const { return _L; }
inline int Percolation2D::get_trial() const { return _trial; }

// for each site
inline int Percolation2D::get_n_occupied() const { return _n_occupied; }
inline vector< vector<bool> > Percolation2D::get_occupied() const { return _occupied; }
inline vector< vector<int> > Percolation2D::get_cluster_id() const { return _cluster_id; }
// for each cluster
inline int Percolation2D::get_n_cluster() const { return _n_cluster; }
inline vector<int> Percolation2D::get_cluster_size() const { return _cluster_size; }
inline vector< vector<int> > Percolation2D::get_cluster_x() const { return _cluster_x; }
inline vector< vector<int> > Percolation2D::get_cluster_y() const { return _cluster_y; }
inline vector<int> Percolation2D::get_spanning_cluster_id() const { return _spanning_cluster_id; }
// for average
inline double Percolation2D::get_P_average() const { return _P_average; }
inline double Percolation2D::get_S_average() const { return _S_average; }

#endif
