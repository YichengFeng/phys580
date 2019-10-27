#include <iostream>
#include <vector>
#include <cmath>
#include "random_walk_2d.h"
#include "random_walk_3d.h"
#include "diffusion_2d.h"
#include "percolation_2d.h"

#include <TGraph.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TF1.h>
#include <TLegend.h>
#include <TH2.h>

using namespace std;


//-------------------------------------------------------------------------//

void random_walk_2d_plot(int mod) {

	int n_step = 100;
	int n_walk = 1000;

	TString str_mod[3] = {"fixed", "continuous", "saw"};
	TString str_adj;
	TString str_tmp;

	str_adj = str_mod[mod] + "_2d";

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

void random_walk_3d_plot(int mod) {

	int n_step = 50;
	int n_walk = 200;

	TString str_mod[3] = {"fixed", "continuous", "saw"};
	TString str_adj;
	TString str_tmp;

	str_adj = str_mod[mod] + "_3d";

	RandomWalk3D rw3d(n_step, n_walk, mod);
	rw3d.cal();
	vector<double> t = rw3d.get_t();
	vector<double> r2 = rw3d.get_r2();
	vector<double> r2_std = rw3d.get_r2_std();

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

void diffusion_2d_plot(int mod, int M, int P, double t_max) {

	TString str_mod[2] = {"close", "leak"};
	TString str_adj = str_mod[mod] + "_";
	TString str_tmp;

	//int M = 100;
	//int P = 300;
	//double t_max = 1e8;
	int D = 5;
	int n_t = 128;

	Diffusion2D df2d((string)str_adj, M, P, D, n_t, t_max);
	df2d.set_mod(mod);
	df2d.set_mod_sample(mod);
	df2d.cal_all();
	vector<double> t_sample = df2d.get_t_sample();
	vector<double> S_sample = df2d.get_S_sample();
	vector<double> x_end = df2d.get_x();
	vector<double> y_end = df2d.get_y();
	vector< vector<double> > x_sample = df2d.get_x_sample();
	vector< vector<double> > y_sample = df2d.get_y_sample();
	vector<double> n_sample;
	vector<double> log_n_sample;
	for(int i=0; i<x_sample.size(); i++) {
		double n_tmp = 1.0*x_sample[i].size();
		n_sample.push_back(n_tmp);
		if(n_tmp>0) log_n_sample.push_back(log(n_tmp));
	}

	str_tmp = str_adj + "entropy";
	TCanvas *c_entropy = new TCanvas("c_"+str_tmp, "c_"+str_tmp);
	TGraph *g_entropy = new TGraph(t_sample.size(), &t_sample[0], &S_sample[0]);
	g_entropy->GetXaxis()->SetTitle("time");
	g_entropy->GetYaxis()->SetTitle("entropy");
	g_entropy->Draw("PLA");
	c_entropy->SetLogx();
	c_entropy->SaveAs("./plot/"+str_tmp+".pdf");

	str_tmp = str_adj + "number";
	TCanvas *c_number = new TCanvas("c_"+str_tmp, "c_"+str_tmp);
	//TGraph *g_number = new TGraph(t_sample.size(), &t_sample[0], &n_sample[0]);
	TGraph *g_number = new TGraph(log_n_sample.size(), &t_sample[0], &log_n_sample[0]);
	g_number->GetXaxis()->SetTitle("time");
	g_number->GetYaxis()->SetTitle("ln(n)");
	g_number->Draw("PLA");
	TF1 *func_number = new TF1("func_number", "[0]+[1]*x", 0, t_max);
	func_number->SetLineStyle(7);
	func_number->SetLineColor(kRed);
	g_number->Fit(func_number);
	TLegend *l_number = new TLegend(0.7,0.6,0.9,0.9);
	l_number->AddEntry(g_number, "simulation", "l");
	l_number->AddEntry(func_number, "fit", "l");
	l_number->AddEntry((TObject*)0, Form("#tau = %.2e", -1/func_number->GetParameter(1)), "");
	l_number->Draw("same");
	c_number->SaveAs("./plot/"+str_tmp+".pdf");

	TH2D *hframe = new TH2D("frame", "frame", 1,-1.0*M,1.0*M, 1,-1.0*M,1.0*M);

	str_tmp = str_adj + "xy_end";
	TCanvas *cxy_end = new TCanvas("c_"+str_tmp, "c_"+str_tmp, 600,600);
	TGraph *gxy_end = new TGraph(x_end.size(), &x_end[0], &y_end[0]);
	gxy_end->GetXaxis()->SetTitle("x");
	gxy_end->GetYaxis()->SetTitle("y");
	gxy_end->SetMarkerStyle(20);
	hframe->Draw();
	gxy_end->Draw("P same");
	cxy_end->SaveAs("./plot/"+str_tmp+".pdf");

	TGraph *g_xy[8];
	TLegend *l_xy[8];
	TCanvas *c_xy[8];
	for(int i=0; i<8; i++) {
		int j = i*16;
		str_tmp = str_adj + Form("xy_t%d",i+1);
		g_xy[i] = new TGraph(x_sample[j].size(), &x_sample[j][0], &y_sample[j][0]);
		g_xy[i]->SetMarkerStyle(20);
		g_xy[i]->GetXaxis()->SetTitle("x");
		g_xy[i]->GetYaxis()->SetTitle("y");
		l_xy[i] = new TLegend(0.1,0.9,0.9,0.96);
		l_xy[i]->SetFillStyle(0);
		l_xy[i]->SetBorderSize(0);
		l_xy[i]->SetNColumns(2);
		l_xy[i]->AddEntry((TObject*)0, Form("t = %.2e",t_sample[j]), "");
		l_xy[i]->AddEntry((TObject*)0, Form("n = %d", (int)x_sample[j].size()), "");
		c_xy[i] = new TCanvas("c_"+str_tmp, "c_"+str_tmp, 600,600);
		hframe->Draw();
		g_xy[i]->Draw("P");
		l_xy[i]->Draw("same");
		c_xy[i]->SaveAs("./plot/"+str_tmp+".pdf");
	}
}

//-------------------------------------------------------------------------//

void percolation_2d_dimension() {

	TString str_adj = "Percolation";
	TString str_tmp;

	const int nL = 7;
	double Ls[nL] = {50, 100, 200, 400, 600, 800, 1000};
	double L2P[nL];
	double logLs[nL];
	double logL2P[nL];
	for(int i=0; i<nL; i++) {
		Percolation2D pc2d(0.593, int(Ls[i]), 50);
		pc2d.cal_trials();
		L2P[i] = pc2d.get_P_average()*Ls[i]*Ls[i];
		cout << L2P[i] << endl;
		logLs[i] = log(Ls[i]);
		logL2P[i] = log(L2P[i]);
	}

	str_tmp = str_adj + "_L_L2P";
	TCanvas *c_L_L2P = new TCanvas(str_tmp, str_tmp);
	TGraph *g_L_L2P = new TGraph(nL, Ls, L2P);
	g_L_L2P->GetXaxis()->SetTitle("L");
	g_L_L2P->GetYaxis()->SetTitle("L^{2}P(p_{c})");
	g_L_L2P->SetMarkerStyle(20);
	g_L_L2P->Draw("PLA");
	c_L_L2P->SaveAs("./plot/"+str_tmp+".pdf");

	str_tmp = str_adj + "_logL_logL2P";
	TCanvas *c_logL_logL2P = new TCanvas(str_tmp, str_tmp);
	TGraph *g_logL_logL2P = new TGraph(nL, logLs, logL2P);
	g_logL_logL2P->GetXaxis()->SetTitle("ln(L)");
	g_logL_logL2P->GetYaxis()->SetTitle("ln(L^{2}P(p_{c}))");
	g_logL_logL2P->SetMarkerStyle(20);
	g_logL_logL2P->Draw("PLA");
	TF1 *func_logL_logL2P = new TF1("func_logL_logL2P", "[0]+[1]*x");
	func_logL_logL2P->SetLineColor(kRed);
	func_logL_logL2P->SetLineStyle(7);
	g_logL_logL2P->Fit(func_logL_logL2P);
	TLegend *l_logL_logL2P = new TLegend(0.1,0.9,0.9,0.96);
	l_logL_logL2P->SetFillStyle(0);
	l_logL_logL2P->SetBorderSize(0);
	l_logL_logL2P->SetNColumns(2);
	l_logL_logL2P->AddEntry((TObject*)0, Form("slope = %.5f",func_logL_logL2P->GetParameter(1)), "");
	l_logL_logL2P->AddEntry((TObject*)0, Form("intercept = %.5f",func_logL_logL2P->GetParameter(0)), "");
	l_logL_logL2P->Draw("same");
	c_logL_logL2P->SaveAs("./plot/"+str_tmp+".pdf");
}

//-------------------------------------------------------------------------//

int main() {

	gStyle->SetTitleOffset(1.4,"y");
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);

	random_walk_2d_plot(2);
	random_walk_3d_plot(2);
	diffusion_2d_plot(0, 100, 300, 1e8);
	diffusion_2d_plot(1, 25, 400, 1e7);

	percolation_2d_dimension();

}

