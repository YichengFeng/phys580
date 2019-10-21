#include <iostream>
#include <vector>
#include <cmath>
#include "monte_carlo_integrate.h"
#include "random_walk_2d.h"

#include <TGraph.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TF1.h>
#include <TLegend.h>

using namespace std;


//-------------------------------------------------------------------------//

double f1(const vector<double> &x) {

	return sqrt(4-x[0]*x[0]);
}

//-------------------------------------------------------------------------//

double f2(const vector<double> &x) {

	return 1.0/sqrt(4-x[0]*x[0]);
}

//-------------------------------------------------------------------------//

double f3(const vector<double> &x) {

	return exp(-x[0])*log(x[0]);
}

//-------------------------------------------------------------------------//

vector<double> monte_carlo_integrate_print(double (*f)(const vector<double> &x), double xi, double xf, int trial) {

	vector<double> x1;
	vector<double> x2;
	x1.push_back(xi);
	x2.push_back(xf);

	MonteCarloIntegrate mci(1, x1, x2, f);
	mci.set_n(1000); // per trial
	mci.set_trial(trial); // number of trials
	mci.cal();
	double value = mci.get_output();
	double error = mci.get_error();

	cout << trial << "   " << value << " +- " << error << endl;

	vector<double> output(2);
	output[0] = value;
	output[1] = error;

	return output;
}

//-------------------------------------------------------------------------//

void monte_carlo_pi(void) {

	const int nn = 6;
	double trial[nn] = {5, 10, 100, 500, 1000, 20000};
	double value[nn];
	double error[nn];
	double logN[nn];
	double logerror[nn];
	for(int i=0; i<nn; i++) {
		vector<double> vtmp = monte_carlo_integrate_print(f1, 0, 2, int(trial[i]));
		value[i] = vtmp[0];
		error[i] = vtmp[1];
		logN[i] = log10(1000.0*trial[i]);
		logerror[i] = log10(error[i]);
	}
	TGraph *g_log_N_error = new TGraph(nn, logN, logerror);
	g_log_N_error->GetXaxis()->SetTitle("log10(N)");
	g_log_N_error->GetYaxis()->SetTitle("log10(error)");
	g_log_N_error->SetMarkerStyle(20);
	TCanvas *c_log_N_error = new TCanvas("log_N_error", "log_N_error");
	g_log_N_error->Draw("PLA");
	TF1 *func = new TF1("func", "[0]+[1]*x", 0,8);
	func->SetLineStyle(7);
	g_log_N_error->Fit(func);
	c_log_N_error->SaveAs("./plot/log_N_error.pdf");
}

//-------------------------------------------------------------------------//

void random_walk_2d_plot(int mod) {

	gStyle->SetOptFit(0);

	int n_step = 100;
	int n_walk = 1000;

	TString str_mod[3] = {"fixed", "continuous", "saw"};
	TString str_adj;
	TString str_tmp;

	str_adj = str_mod[mod];

	RandomWalk2D rw2d(n_step, n_walk, mod);
	rw2d.cal();
	vector<double> t = rw2d.get_t();
	vector<double> r2 = rw2d.get_r2();
	vector<double> r2_std = rw2d.get_r2_std();

	vector<double> t_log;
	vector<double> r2_log;
	for(int i=0; i<n_step; i++) {
		if(t[i]>0 && r2[i]>0) {
			t_log.push_back(log10(t[i]));
			r2_log.push_back(log10(r2[i]));
		}
	}

	double p0;
	double p0e;
	double p1;
	double p1e;

	str_tmp = "r2_"+str_adj;
	TGraph *g_r2 = new TGraph(n_step, &t[0], &r2[0]);
	TCanvas *c_r2 = new TCanvas("c_"+str_tmp, "c_"+str_tmp);
	g_r2->GetXaxis()->SetTitle("t (step)");
	g_r2->GetYaxis()->SetTitle("<r^{2}>");
	g_r2->SetLineWidth(2);
	g_r2->SetMarkerStyle(20);
	g_r2->Draw("PLA");
	TF1 *func_r2 = new TF1("func_r2", "[0]+[1]*x", 0,1.0*n_step);
	func_r2->SetLineStyle(7);
	g_r2->Fit(func_r2);
	p0 = func_r2->GetParameter(0);
	p0e= func_r2->GetParError(0);
	p1 = func_r2->GetParameter(1);
	p1e= func_r2->GetParError(1);
	TLegend *l_r2 = new TLegend(0.1,0.7,0.6,0.9);
	l_r2->SetFillStyle(0);
	l_r2->SetBorderSize(0);
	l_r2->AddEntry((TObject*)0, Form("slope = %.4f #pm %.4f",p1,p1e), "");
	l_r2->AddEntry((TObject*)0, Form("intercept = %.4f #pm %.4f",p0,p0e), "");
	l_r2->Draw("same");
	c_r2->SaveAs("./plot/"+str_tmp+".pdf");

	str_tmp = "r2_log_"+str_adj;
	TGraph *g_r2_log = new TGraph(t_log.size(), &t_log[0], &r2_log[0]);
	TCanvas *c_r2_log = new TCanvas("c_"+str_tmp, "c_"+str_tmp);
	g_r2_log->GetXaxis()->SetTitle("log10(t) (step)");
	g_r2_log->GetYaxis()->SetTitle("log10(<r^{2}>)");
	g_r2_log->SetLineWidth(2);
	g_r2_log->SetMarkerStyle(20);
	g_r2_log->Draw("PLA");
	TF1 *func_r2_log = new TF1("func_r2_log", "[0]+[1]*x", 0,1.4);
	func_r2_log->SetLineStyle(7);
	g_r2_log->Fit(func_r2_log, "R");
	p0 = func_r2_log->GetParameter(0);
	p0e= func_r2_log->GetParError(0);
	p1 = func_r2_log->GetParameter(1); 
	p1e= func_r2_log->GetParError(1);
	TLegend *l_r2_log = new TLegend(0.1,0.7,0.6,0.9);
	l_r2_log->SetFillStyle(0);
	l_r2_log->SetBorderSize(0);
	l_r2_log->AddEntry((TObject*)0, Form("slope = %.4f #pm %.4f",p1,p1e), "");
	l_r2_log->AddEntry((TObject*)0, Form("intercept = %.4f #pm %.4f",p0,p0e), "");
	l_r2_log->Draw("same");
	c_r2_log->SaveAs("./plot/"+str_tmp+".pdf");

	str_tmp = "r2_std_"+str_adj;
	TGraph *g_r2_std = new TGraph(n_step, &t[0], &r2_std[0]);
	TCanvas *c_r2_std = new TCanvas("c_"+str_tmp, "c_"+str_tmp);
	g_r2_std->GetXaxis()->SetTitle("t (step)");
	g_r2_std->GetYaxis()->SetTitle("#sigma(r^{2})");
	g_r2_std->SetLineWidth(2);
	g_r2_std->SetMarkerStyle(20);
	g_r2_std->Draw("PLA");
	TF1 *func_r2_std = new TF1("func_r2_std", "[0]*x", 0,1.0*n_step);
	func_r2_std->SetLineStyle(7);
	c_r2_std->SaveAs("./plot/"+str_tmp+".pdf");
}

//-------------------------------------------------------------------------//

int main() {

	gStyle->SetTitleOffset(1.4,"y");
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);
	gStyle->SetOptFit();

	// (1)
	//monte_carlo_pi();
	//double ee = 1e-8;
	//monte_carlo_integrate_print(f2, -2+ee, 2-ee, 10000);
	//monte_carlo_integrate_print(f3, 0+ee, 30, 10000);

	// (2) (3)
	random_walk_2d_plot(0);
	random_walk_2d_plot(1);
	random_walk_2d_plot(2);

}
