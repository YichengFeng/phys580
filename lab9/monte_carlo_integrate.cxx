#include "monte_carlo_integrate.h"

using namespace std;


//-------------------------------------------------------------------------//

MonteCarloIntegrate::MonteCarloIntegrate() {

}

//-------------------------------------------------------------------------//

MonteCarloIntegrate::MonteCarloIntegrate(int dim, vector<double> x1, vector<double> x2, double (*f)(const vector<double> &x)) {

	_x1 = x1;
	_x2 = x2;
	_f  = f;
	_n  = 1000;
	_alg = 0;
	_dim = dim;
	_trial = 1000;

	_vol = 1;
	for(int i=0; i<_dim; i++) {
		_vol *= _x2[i]-_x1[i];
	}

	srand(time(nullptr));
}

//-------------------------------------------------------------------------//

bool MonteCarloIntegrate::check() const {

	if(_x1.size() != _dim || _x2.size() != _dim) {
		cout << "ERROR: dimension not match!" << endl;
		return false;
	}

	return true;
}

//-------------------------------------------------------------------------//

double MonteCarloIntegrate::cal_once() {

	if(!check()) return 0;

	double tmpy = 0;
	double sumy = 0;
	double sumyy = 0;

	for(int i=0; i<_n; i++) {
		vector<double> x_tmp;
		for(int j=0; j<_dim; j++) {
			//srand(time(nullptr));
			double tmp = _x1[j] + (_x2[j]-_x1[j])*rand()/RAND_MAX;
			x_tmp.push_back(tmp);
		}

		tmpy = _f(x_tmp);

		sumy += tmpy;
	}

	return sumy*_vol/_n;
}

//-------------------------------------------------------------------------//

double MonteCarloIntegrate::cal() {

	_output = 0;
	_error = 0;

	double sum = 0;
	double sum2 = 0;
	double tmpy = 0;

	for(int i=0; i<_trial; i++) {
		tmpy = cal_once();
		sum += tmpy;
		sum2 += tmpy*tmpy;
	}

	_output = sum/_trial;
	double sd = sqrt(sum2/_trial - sum*sum/_trial/_trial);
	_error = sd/sqrt(1.0*_trial);

	return _output;
}
