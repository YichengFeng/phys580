#include <iostream>
#include <vector>
#include <cmath>
#include "runge_kutta.h"
#include "planet_spin_orbit.h"

#include <TGraph.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH2.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TF1.h>

using namespace std;

double G_ALG = 4;

//-------------------------------------------------------------------------//

void theta_to_pipi(vector<double> &theta) {

	for(int i=0; i<theta.size(); i++) {
		while(theta[i]>=M_PI) {
			theta[i] -= 2*M_PI;
		}
		while(theta[i]<-M_PI) {
			theta[i] += 2*M_PI;
		}
	}
}

//-------------------------------------------------------------------------//

void theta_to_twopi(vector<double> &theta) {

	for(int i=0; i<theta.size(); i++) {
		while(theta[i]>=2*M_PI) {
			theta[i] -= 2*M_PI;
		}
		while(theta[i]<0) {
			theta[i] += 2*M_PI;
		}
	}
}

//-------------------------------------------------------------------------//

TF1* fit_Lyapunov_exp(const vector<double> &t, const vector<double> &x, double t_start, double t_end) {

	vector<double> tlm;
	vector<double> xlm;
	int n = 0;
	for(int i=1; i<t.size()-1; i++) {
		if(t[i]<t_start || t[i]>t_end) continue;
		if(x[i]<x[i-1] || x[i]<x[i+1]) continue;
		if(x[i]<1e-4) continue;
		//cout << n << "   " << t[i] << "   " << x[i] << endl;
		tlm.push_back(t[i]);
		xlm.push_back(log(x[i]));
		n ++;
	}

	TGraph *g = new TGraph(n, &tlm[0], &xlm[0]);
	//g->SetMarkerColor(kRed);
	//g->Draw("same P");

	gStyle->SetOptFit();
	TF1 *line = new TF1("line", "[0] + [1]*x", t_start, t_end);
	g->Fit(line);

	TF1 *expo = new TF1("expo", "exp([0] + [1]*x)", t_start, t_end);
	expo->SetLineColor(kBlue);
	expo->SetLineStyle(7);
	expo->SetParameters(line->GetParameters());

	TLegend *l_fit_para = new TLegend(0.6,0.9,0.9,0.95);
	l_fit_para->SetHeader(Form("ln(y) = %.4fx %+.4f", line->GetParameter(1), line->GetParameter(0)), "");
	l_fit_para->Draw("same");

	return expo;
}

//-------------------------------------------------------------------------//

void planet_spin_orbit_plot(double e, double theta_start, double omega_start) {

	TString str_e = Form("e%.3d", int(e*100));
	TString str_theta0 = Form("theta%.3d", int(theta_start*100));
	TString str_omega0 = Form("omega%.3d", int(omega_start*100));

	double a = 1/(1-0.123);
	double r_min = a*(1-e);
	double r_max = a*(1+e);
	double f = a*e;
	double b = sqrt(a*a - f*f);

	PlanetSpinOrbit pso(a, e, 0);
	pso.set_alg(G_ALG);
	pso.set_t_end(10);
	pso.set_theta_start(theta_start);
	pso.set_omega_start(omega_start);
	pso.cal();

	vector<double> t = pso.get_t();
	vector< vector<double> > x = pso.get_x();
	int n_stps = pso.get_n_stps();

	theta_to_pipi(x[2]);

	TGraph *g_x_y = new TGraph(n_stps, &x[0][0], &x[1][0]);
	TGraph *g_t_theta = new TGraph(n_stps, &t[0], &x[2][0]);
	TGraph *g_t_omega = new TGraph(n_stps, &t[0], &x[5][0]);
	TGraph *g_theta_omega = new TGraph(n_stps, &x[2][0], &x[5][0]);

	double x_Saturn[1] = {0};
	double y_Saturn[1] = {0};
	TGraph *g_Saturn = new TGraph(1, x_Saturn, y_Saturn);
	g_Saturn->SetMarkerStyle(20);
	g_Saturn->SetMarkerSize(2);
	g_Saturn->SetMarkerColor(kOrange);

	int canvas_nx = 700;
	int canvas_ny = 500;
	double rm = 0.05;
	double rx, ry;
	if(a/canvas_nx < b/canvas_ny) {
		rx = b/canvas_ny*canvas_nx/a;
		ry = 1;
	} else {
		rx = 1;
		ry = a/canvas_nx*canvas_ny/b;
	}
	TH2F *h2frame = new TH2F("frame", "frame", 1,-(r_max+rm*a)*rx,(r_min+rm*a)*rx, 1,-(b+b*rm)*ry,(b+b*rm)*ry);
	h2frame->GetXaxis()->SetTitle("x (HU)");
	h2frame->GetYaxis()->SetTitle("y (HU)");

	TLegend *l_para = new TLegend(0.5,0.9,0.9,0.96);
	l_para->SetNColumns(3);
	l_para->SetFillStyle(0);
	l_para->SetBorderSize(0);
	l_para->AddEntry((TObject*)0, Form("e = %.3f", e), "");
	l_para->AddEntry((TObject*)0, Form("#theta(0) = %.3f", theta_start), "");
	l_para->AddEntry((TObject*)0, Form("#omega(0) = %.3f", omega_start), "");

	TString str_tmp;
	TString str_adj = str_e + "_" + str_theta0 + "_" + str_omega0;

	str_tmp = "x_y_" + str_adj;
	TCanvas *c_x_y = new TCanvas("c_"+str_tmp, "c_"+str_tmp);
	h2frame->GetXaxis()->SetTitle("x (HU)");
	h2frame->GetYaxis()->SetTitle("y (HU)");
	h2frame->Draw();
	g_x_y->SetLineColor(kBlack);
	g_x_y->Draw("L same");
	g_Saturn->Draw("P same");
	TLegend *l_obj = new TLegend(0.1,0.7,0.3,0.9);
	l_obj->SetFillStyle(0);
	l_obj->SetBorderSize(0);
	l_obj->AddEntry(g_x_y, Form("#splitline{Hyperion}{(a = %.3f HU, e = %.3f)}", a, e), "l");
	l_obj->AddEntry(g_Saturn, "Saturn", "p");
	l_obj->Draw("same");
	l_para->Draw("same");
	c_x_y->SaveAs("./plot/"+str_tmp+".pdf");

	str_tmp = "t_theta_" + str_adj;
	TCanvas *c_t_theta = new TCanvas("c_"+str_tmp, "c_"+str_tmp);
	g_t_theta->SetLineColor(kRed);
	g_t_theta->GetXaxis()->SetTitle("t (Hyr)");
	g_t_theta->GetYaxis()->SetTitle("#theta (rad)");
	g_t_theta->Draw("LA");
	l_para->Draw("same");
	c_t_theta->SaveAs("./plot/"+str_tmp+".pdf");

	str_tmp = "t_omega_" + str_adj;
	TCanvas *c_t_omega = new TCanvas("c_"+str_tmp, "c_"+str_tmp);
	g_t_omega->SetLineColor(kBlue);
	g_t_omega->GetXaxis()->SetTitle("t (Hyr)");
	g_t_omega->GetYaxis()->SetTitle("#omega (rad/Hyr)");
	g_t_omega->Draw("LA");
	l_para->Draw("same");
	c_t_omega->SaveAs("./plot/"+str_tmp+".pdf");

	str_tmp = "theta_omega_" + str_adj;
	TCanvas *c_theta_omega = new TCanvas("c_"+str_tmp, "c_"+str_tmp);
	g_theta_omega->SetLineColor(kGreen);
	g_theta_omega->GetXaxis()->SetTitle("#theta (rad)");
	g_theta_omega->GetYaxis()->SetTitle("#omega (rad/Hyr)");
	g_theta_omega->Draw("LA");
	l_para->Draw("same");
	c_theta_omega->SaveAs("./plot/"+str_tmp+".pdf");
}

//-------------------------------------------------------------------------//

void butterfly_effect_plot(double e, double t_end, double fit_t_start, double fit_t_end) {

	TString str_e = Form("e%.3d", int(e*100));

	double a = 1/(1-0.123);
	double r_min = a*(1-e);
	double r_max = a*(1+e);
	double f = a*e;
	double b = sqrt(a*a - f*f);

	double theta_start0 = 0.2;
	double omega_start0 = 0.2;
	PlanetSpinOrbit pso0(a, e, 0);
	pso0.set_alg(G_ALG);
	pso0.set_t_end(t_end);
	pso0.set_theta_start(theta_start0);
	pso0.set_omega_start(omega_start0);
	pso0.cal();
	vector<double> t0 = pso0.get_t();
	vector< vector<double> > x0 = pso0.get_x();
	int n_stps0 = pso0.get_n_stps();
	//theta_to_pipi(x0[2]);

	double theta_start1 = 0.2001;
	double omega_start1 = 0.2;
	PlanetSpinOrbit pso1(a, e, 0);
	pso1.set_alg(G_ALG);
	pso1.set_t_end(t_end);
	pso1.set_theta_start(theta_start1);
	pso1.set_omega_start(omega_start1);
	pso1.cal();
	vector<double> t1 = pso1.get_t();
	vector< vector<double> > x1 = pso1.get_x();
	int n_stps1 = pso1.get_n_stps();
	//theta_to_pipi(x1[2]);

	double theta_start2 = 0.2;
	double omega_start2 = 0.2001;
	PlanetSpinOrbit pso2(a, e, 0);
	pso2.set_alg(G_ALG);
	pso2.set_t_end(t_end);
	pso2.set_theta_start(theta_start2);
	pso2.set_omega_start(omega_start2);
	pso2.cal();
	vector<double> t2 = pso2.get_t();
	vector< vector<double> > x2 = pso2.get_x();
	int n_stps2 = pso2.get_n_stps();
	//theta_to_pipi(x2[2]);

	TString str0 = Form("#theta_{0}(0) = %.4f, #omega_{0}(0) = %.4f", theta_start0, omega_start0);
	TString str1 = Form("#theta_{1}(0) = %.4f, #omega_{1}(0) = %.4f", theta_start1, omega_start1);
	TString str2 = Form("#theta_{2}(0) = %.4f, #omega_{2}(0) = %.4f", theta_start2, omega_start2);
	TString str01 = "#Delta#theta = |#theta_{1}-#theta_{0}|, #Delta#omega = |#omega_{1}-#omega_{0}|";
	TString str02 = "#Delta#theta = |#theta_{2}-#theta_{0}|, #Delta#omega = |#omega_{2}-#omega_{0}|";
	TString stre = Form("e = %.3f", e);

	TLegend *l01 = new TLegend(0.1,0.7,0.4,0.9);
	l01->SetFillStyle(0);
	l01->SetBorderSize(0);
	l01->AddEntry((TObject*)0, stre, "");
	l01->AddEntry((TObject*)0, str0, "");
	l01->AddEntry((TObject*)0, str1, "");
	l01->AddEntry((TObject*)0, str01, "");

	TLegend *l02 = new TLegend(0.1,0.7,0.4,0.9);
	l02->SetFillStyle(0);
	l02->SetBorderSize(0);
	l02->AddEntry((TObject*)0, stre, "");
	l02->AddEntry((TObject*)0, str0, "");
	l02->AddEntry((TObject*)0, str2, "");
	l02->AddEntry((TObject*)0, str02, "");

	vector<double> dtheta01;
	vector<double> domega01;
	vector<double> dtheta02;
	vector<double> domega02;
	for(int i=0; i<n_stps0; i++) {
		dtheta01.push_back(fabs(x1[2][i]-x0[2][i]));
		domega01.push_back(fabs(x1[5][i]-x0[5][i]));
		dtheta02.push_back(fabs(x2[2][i]-x0[2][i]));
		domega02.push_back(fabs(x2[5][i]-x0[5][i]));
	}
	theta_to_twopi(dtheta01);
	theta_to_twopi(domega01);
	theta_to_twopi(dtheta02);
	theta_to_twopi(domega02);

	TGraph *g_t_dtheta01 = new TGraph(n_stps0, &t0[0], &dtheta01[0]);
	TGraph *g_t_domega01 = new TGraph(n_stps0, &t0[0], &domega01[0]);
	TGraph *g_t_dtheta02 = new TGraph(n_stps0, &t0[0], &dtheta02[0]);
	TGraph *g_t_domega02 = new TGraph(n_stps0, &t0[0], &domega02[0]);

	TString str_adj = str_e;
	TString str_tmp;

	str_tmp = "t_dtheta01_" + str_e;
	TCanvas *c_t_dtheta01 = new TCanvas("c_"+str_tmp, "c_"+str_tmp);
	g_t_dtheta01->GetYaxis()->SetRangeUser(1e-5, 7);
	g_t_dtheta01->GetXaxis()->SetTitle("t (Hyr)");
	g_t_dtheta01->GetYaxis()->SetTitle("#Delta#theta");
	g_t_dtheta01->Draw("LA");
	fit_Lyapunov_exp(t0, dtheta01, fit_t_start, fit_t_end)->Draw("same");
	l01->Draw("same");
	c_t_dtheta01->SetLogy();
	c_t_dtheta01->SaveAs("./plot/"+str_tmp+".pdf");

	str_tmp = "t_domega01_" + str_e;
	TCanvas *c_t_domega01 = new TCanvas("c_"+str_tmp, "c_"+str_tmp);
	g_t_domega01->GetYaxis()->SetRangeUser(1e-5, 7);
	g_t_domega01->GetXaxis()->SetTitle("t (Hyr)");
	g_t_domega01->GetYaxis()->SetTitle("#Delta#omega");
	g_t_domega01->Draw("LA");
	fit_Lyapunov_exp(t0, domega01, fit_t_start, fit_t_end)->Draw("same");
	l01->Draw("same");
	c_t_domega01->SetLogy();
	c_t_domega01->SaveAs("./plot/"+str_tmp+".pdf");

	str_tmp = "t_dtheta02_" + str_e;
	TCanvas *c_t_dtheta02 = new TCanvas("c_"+str_tmp, "c_"+str_tmp);
	g_t_dtheta02->GetYaxis()->SetRangeUser(1e-5, 7);
	g_t_dtheta02->GetXaxis()->SetTitle("t (Hyr)");
	g_t_dtheta02->GetYaxis()->SetTitle("#Delta#theta");
	g_t_dtheta02->Draw("LA");
	fit_Lyapunov_exp(t0, dtheta02, fit_t_start, fit_t_end)->Draw("same");
	l02->Draw("same");
	c_t_dtheta02->SetLogy();
	c_t_dtheta02->SaveAs("./plot/"+str_tmp+".pdf");

	str_tmp = "t_domega02_" + str_e;
	TCanvas *c_t_domega02 = new TCanvas("c_"+str_tmp, "c_"+str_tmp);
	g_t_domega02->GetYaxis()->SetRangeUser(1e-5, 7);
	g_t_domega02->GetXaxis()->SetTitle("t (Hyr)");
	g_t_domega02->GetYaxis()->SetTitle("#Delta#omega");
	g_t_domega02->Draw("LA");
	fit_Lyapunov_exp(t0, domega02, fit_t_start, fit_t_end)->Draw("same");
	l02->Draw("same");
	c_t_domega02->SetLogy();
	c_t_domega02->SaveAs("./plot/"+str_tmp+".pdf");
}

//-------------------------------------------------------------------------//

int main() {

	gStyle->SetTitleOffset(1.4,"y");
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);

	planet_spin_orbit_plot(0.000, 0, 0);
	planet_spin_orbit_plot(0.123, 0, 0);
	planet_spin_orbit_plot(0.300, 0, 0);
	planet_spin_orbit_plot(0.600, 0, 0);

	planet_spin_orbit_plot(0.000, 0.2, 0.2);
	planet_spin_orbit_plot(0.123, 0.2, 0.2);
	planet_spin_orbit_plot(0.300, 0.2, 0.2);
	planet_spin_orbit_plot(0.600, 0.2, 0.2);

	butterfly_effect_plot(0.000, 20, 1, 20);
	butterfly_effect_plot(0.123, 10, 0, 10);
	butterfly_effect_plot(0.300, 7 , 0, 7);
	butterfly_effect_plot(0.600, 5 , 0, 5);
}
