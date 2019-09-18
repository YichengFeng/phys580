#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include "runge_kutta.h"
#include "pendulum.h"

#include <TROOT.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TMath.h>

using namespace std;


//-------------------------------------------------------------------------//

void theta_to_twopi(vector<double> &theta) {

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

double ellint_sin_half_theta(double theta_max) {

	double period;
	period = M_PI/2;
	//period *= 1 + pow(theta_max,2)/16 + pow(theta_max,4)*11/3072;
	period *= 1 + pow(theta_max,2)/16 + pow(theta_max,4)*11/3072 + pow(theta_max,6)*173/737280;
}

//-------------------------------------------------------------------------//

void chao_plot(double f_d, double theta_start, double omega_start, double psphi, double psomega) {

	Pendulum pd;
	pd.set_mode(1); // physics pendulum
	pd.set_alg(4);
	//pd.set_theta_start(0.2);
	//pd.set_omega_start(0);
	pd.set_theta_start(theta_start);
	pd.set_omega_start(omega_start);
	pd.set_q(0.5);
	pd.set_F(f_d);
	pd.set_O(2./3);
	pd.set_dt(0.04);
	pd.set_n_periods(400);

	pd.cal();

	vector<double> t = pd.get_t();
	vector< vector<double> > x = pd.get_x();
	vector<double> energy = pd.get_energy();
	int n_stps = pd.get_n_stps();

	theta_to_twopi(x[0]);

	// fixed poincare section
	vector<double> ps_f_d;
	vector<double> ps_theta;
	vector<double> ps_omega;
	vector<double> ps_energy;
	int ps_n_pts = 0;

	//double psphi = 0;
	//double psomega = 2./3;
	double dt = pd.get_dt();

	for(int i_stps=int(n_stps*0.); i_stps<n_stps; i_stps++) {
		if(fabs(psomega*t[i_stps]-psphi-int((psomega*t[i_stps]-psphi)/2./M_PI)*2*M_PI) > 0.5*psomega*dt*0.999999) {
			continue;
		}
		ps_f_d.push_back(f_d);
		ps_theta.push_back(x[0][i_stps]);
		ps_omega.push_back(x[1][i_stps]);
		ps_energy.push_back(energy[i_stps]);
		ps_n_pts ++;
	}

	TString str_f_d = Form("fd_%.3d", int(100*f_d));
	TString str_theta_start = Form("theta0_%.3d", int(100*fabs(theta_start)));
	TString str_omega_start = Form("omega0_%.3d", int(100*fabs(omega_start)));
	TString str_psphi = Form("psphi_%.3d", int(psphi*100));
	TString str_psomega = Form("psomega_%.3d", int(psomega*100));
	TString str_c_tmp;

	int n_tail = int(0.975*n_stps);
	int n_head = int(0.025*n_stps);

	TString str_adj = str_f_d + "_" + str_theta_start + "_" + str_omega_start + "_" + str_psphi + "_" + str_psomega;

	str_c_tmp = "ps_theta_omega_" + str_adj;
	TCanvas *c_ps_theta_omega = new TCanvas(str_c_tmp, str_c_tmp);
	TGraph *g_ps_theta_omega  = new TGraph(ps_n_pts, &ps_theta[0], &ps_omega[0]);
	g_ps_theta_omega->GetXaxis()->SetTitle("#theta");
	g_ps_theta_omega->GetYaxis()->SetTitle("#omega");
	g_ps_theta_omega->SetMarkerStyle(5);
	g_ps_theta_omega->Draw("PA");
	c_ps_theta_omega->SaveAs("./plot/" + str_c_tmp + ".pdf");

	str_c_tmp = "ps_omega_energy_" + str_adj;
	TCanvas *c_ps_omega_energy = new TCanvas(str_c_tmp, str_c_tmp);
	TGraph *g_ps_omega_energy  = new TGraph(ps_n_pts, &ps_omega[0], &ps_energy[0]);
	g_ps_omega_energy->GetXaxis()->SetTitle("#omega");
	g_ps_omega_energy->GetYaxis()->SetTitle("#energy");
	g_ps_omega_energy->SetMarkerStyle(5);
	g_ps_omega_energy->Draw("PA");
	c_ps_omega_energy->SaveAs("./plot/" + str_c_tmp + ".pdf");

	str_c_tmp = "t_theta_head_" + str_adj;
	TCanvas *c_t_theta_head = new TCanvas(str_c_tmp, str_c_tmp);
	TGraph *g_t_theta_head  = new TGraph(n_head, &t[0], &x[0][0]);
	g_t_theta_head->GetXaxis()->SetTitle("t (s)");
	g_t_theta_head->GetYaxis()->SetTitle("#theta");
	g_t_theta_head->SetLineColor(kRed);
	g_t_theta_head->Draw("LA");
	c_t_theta_head->SaveAs("./plot/" + str_c_tmp + ".pdf");

	str_c_tmp = "t_omega_head_" + str_adj;
	TCanvas *c_t_omega_head = new TCanvas(str_c_tmp, str_c_tmp);
	TGraph *g_t_omega_head  = new TGraph(n_head, &t[0], &x[1][0]);
	g_t_omega_head->GetXaxis()->SetTitle("t (s)");
	g_t_omega_head->GetYaxis()->SetTitle("#omega");
	g_t_omega_head->SetLineColor(kBlue);
	g_t_omega_head->Draw("LA");
	c_t_omega_head->SaveAs("./plot/" + str_c_tmp + ".pdf");

	str_c_tmp = "t_energy_head_" + str_adj;
	TCanvas *c_t_energy_head = new TCanvas(str_c_tmp, str_c_tmp);
	TGraph *g_t_energy_head = new TGraph(n_head, &t[0], &energy[0]);
	g_t_energy_head->GetXaxis()->SetTitle("t (s)");
	g_t_energy_head->GetYaxis()->SetTitle("energy");
	g_t_energy_head->SetLineColor(kMagenta);
	g_t_energy_head->Draw("LA");
	c_t_energy_head->SaveAs("./plot/" + str_c_tmp + ".pdf");

	str_c_tmp = "t_theta_tail_" + str_adj;
	TCanvas *c_t_theta_tail = new TCanvas(str_c_tmp, str_c_tmp);
	TGraph *g_t_theta_tail  = new TGraph(n_stps-n_tail, &t[n_tail], &x[0][n_tail]);
	g_t_theta_tail->GetXaxis()->SetTitle("t (s)");
	g_t_theta_tail->GetYaxis()->SetTitle("#theta");
	g_t_theta_tail->SetLineColor(kRed);
	g_t_theta_tail->Draw("LA");
	c_t_theta_tail->SaveAs("./plot/" + str_c_tmp + ".pdf");

	str_c_tmp = "t_omega_tail_" + str_adj;
	TCanvas *c_t_omega_tail = new TCanvas(str_c_tmp, str_c_tmp);
	TGraph *g_t_omega_tail  = new TGraph(n_stps-n_tail, &t[n_tail], &x[1][n_tail]);
	g_t_omega_tail->GetXaxis()->SetTitle("t (s)");
	g_t_omega_tail->GetYaxis()->SetTitle("#omega");
	g_t_omega_tail->SetLineColor(kBlue);
	g_t_omega_tail->Draw("LA");
	c_t_omega_tail->SaveAs("./plot/" + str_c_tmp + ".pdf");

	str_c_tmp = "t_energy_tail_" + str_adj;
	TCanvas *c_t_energy_tail = new TCanvas(str_c_tmp, str_c_tmp);
	TGraph *g_t_energy_tail = new TGraph(n_stps-n_tail, &t[n_tail], &energy[n_tail]);
	g_t_energy_tail->GetXaxis()->SetTitle("t (s)");
	g_t_energy_tail->GetYaxis()->SetTitle("energy");
	g_t_energy_tail->SetLineColor(kMagenta);
	g_t_energy_tail->Draw("LA");
	c_t_energy_tail->SaveAs("./plot/" + str_c_tmp + ".pdf");

	str_c_tmp = "theta_omega_" + str_adj;
	TCanvas *c_theta_omega = new TCanvas(str_c_tmp, str_c_tmp);
	TGraph *g_theta_omega   = new TGraph(n_stps, &x[0][0], &x[1][0]);
	g_theta_omega->GetXaxis()->SetTitle("#theta");
	g_theta_omega->GetYaxis()->SetTitle("#omega");
	g_theta_omega->SetLineColor(kGreen);
	g_theta_omega->Draw("LA");
	g_ps_theta_omega->Draw("P same");
	c_theta_omega->SaveAs("./plot/" + str_c_tmp + ".pdf");

	str_c_tmp = "theta_energy_" + str_adj;
	TCanvas *c_theta_energy = new TCanvas(str_c_tmp, str_c_tmp);
	TGraph *g_theta_energy  = new TGraph(n_stps, &x[0][0], &energy[0]);
	g_theta_energy->GetXaxis()->SetTitle("#theta");
	g_theta_energy->GetYaxis()->SetTitle("energy");
	g_theta_energy->SetLineColor(kOrange);
	g_theta_energy->Draw("LA");
	c_theta_energy->SaveAs("./plot/" + str_c_tmp + ".pdf");

	str_c_tmp = "omega_energy_" + str_adj;
	TCanvas *c_omega_energy = new TCanvas(str_c_tmp, str_c_tmp);
	TGraph *g_omega_energy  = new TGraph(n_stps, &x[1][0], &energy[0]);
	g_omega_energy->GetXaxis()->SetTitle("#omega");
	g_omega_energy->GetYaxis()->SetTitle("energy");
	g_omega_energy->SetLineColor(kPink);
	g_omega_energy->Draw("LA");
	g_ps_omega_energy->Draw("P same");
	c_omega_energy->SaveAs("./plot/" + str_c_tmp + ".pdf");

	str_c_tmp = "theta_omega_tail_" + str_adj;
	TCanvas *c_theta_omega_tail = new TCanvas(str_c_tmp, str_c_tmp);
	TGraph *g_theta_omega_tail   = new TGraph(n_stps-n_tail, &x[0][n_tail], &x[1][n_tail]);
	g_theta_omega_tail->GetXaxis()->SetTitle("#theta");
	g_theta_omega_tail->GetYaxis()->SetTitle("#omega");
	g_theta_omega_tail->SetLineColor(kGreen);
	g_theta_omega_tail->Draw("LA");
	g_ps_theta_omega->Draw("P same");
	c_theta_omega_tail->SaveAs("./plot/" + str_c_tmp + ".pdf");

	str_c_tmp = "theta_energy_tail_" + str_adj;
	TCanvas *c_theta_energy_tail = new TCanvas(str_c_tmp, str_c_tmp);
	TGraph *g_theta_energy_tail  = new TGraph(n_stps-n_tail, &x[0][n_tail], &energy[n_tail]);
	g_theta_energy_tail->GetXaxis()->SetTitle("#theta");
	g_theta_energy_tail->GetYaxis()->SetTitle("energy");
	g_theta_energy_tail->SetLineColor(kOrange);
	g_theta_energy_tail->Draw("LA");
	//g_ps_theta_energy->Draw("P same");
	c_theta_energy_tail->SaveAs("./plot/" + str_c_tmp + ".pdf");

	str_c_tmp = "omega_energy_tail_" + str_adj;
	TCanvas *c_omega_energy_tail = new TCanvas(str_c_tmp, str_c_tmp);
	TGraph *g_omega_energy_tail  = new TGraph(n_stps-n_tail, &x[1][n_tail], &energy[n_tail]);
	g_omega_energy_tail->GetXaxis()->SetTitle("#omega");
	g_omega_energy_tail->GetYaxis()->SetTitle("energy");
	g_omega_energy_tail->SetLineColor(kPink);
	g_omega_energy_tail->Draw("LA");
	g_ps_omega_energy->Draw("P same");
	c_omega_energy_tail->SaveAs("./plot/" + str_c_tmp + ".pdf");
}

//-------------------------------------------------------------------------//

void poincare_section(double psphi, double psomega) {

	TString str_psphi = Form("psphi_%.3d", int(psphi*100));
	TString str_psomega = Form("psomega_%.3d", int(psomega*100));

	double omega_d = 2./3;
	double f_d = 1.30;
	double dt = 0.04;
	double T_d = 2*M_PI/omega_d;

	vector<double> ps_f_d;
	vector<double> ps_theta;
	vector<double> ps_omega;
	vector<double> ps_energy;
	int ps_n_pts = 0;

	double f_d_min = 1.35;
	double f_d_max = 1.49;
	int n_f_d = 200;

	for(int i=0; i<n_f_d; i++) {

		f_d = f_d_min+(f_d_max-f_d_min)/n_f_d*i;

		Pendulum pd;
		pd.set_mode(1); // physics pendulum
		pd.set_alg(4);
		pd.set_theta_start(0.2);
		pd.set_omega_start(0);
		pd.set_q(0.5);
		pd.set_F(f_d);
		pd.set_O(omega_d);
		pd.set_dt(dt);
		pd.set_n_periods(400/omega_d);

		pd.cal();

		vector<double> t = pd.get_t();
		vector< vector<double> > x = pd.get_x();
		vector<double> energy = pd.get_energy();
		int n_stps = pd.get_n_stps();

		theta_to_twopi(x[0]);

		for(int i_stps=300/omega_d/dt; i_stps<n_stps; i_stps++) {
			if(fabs(psomega*t[i_stps]-psphi-int((psomega*t[i_stps]-psphi)/2./M_PI)*2*M_PI) > 0.5*psomega*dt*0.999999) {
				continue;
			}
			ps_f_d.push_back(f_d);
			ps_theta.push_back(x[0][i_stps]);
			ps_omega.push_back(x[1][i_stps]);
			ps_energy.push_back(energy[i_stps]);
			ps_n_pts ++;
		}
	}

	TString str_c_tmp;

	str_c_tmp = "ps_f_d_theta_" + str_psphi + "_" + str_psomega;
	TCanvas *c_f_d_theta = new TCanvas(str_c_tmp, str_c_tmp);
	TGraph *g_f_d_theta = new TGraph(ps_n_pts, &ps_f_d[0], &ps_theta[0]);
	g_f_d_theta->GetXaxis()->SetTitle("f_{d}");
	g_f_d_theta->GetYaxis()->SetTitle("#theta");
	g_f_d_theta->Draw("PA");
	c_f_d_theta->SaveAs("./plot/" + str_c_tmp + ".pdf");

	str_c_tmp = "ps_f_d_omega_" + str_psphi + "_" + str_psomega;
	TCanvas *c_f_d_omega = new TCanvas(str_c_tmp, str_c_tmp);
	TGraph *g_f_d_omega = new TGraph(ps_n_pts, &ps_f_d[0], &ps_omega[0]);
	g_f_d_omega->GetXaxis()->SetTitle("f_{d}");
	g_f_d_omega->GetYaxis()->SetTitle("#omega");
	g_f_d_omega->Draw("PA");
	c_f_d_omega->SaveAs("./plot/" + str_c_tmp + ".pdf");

	str_c_tmp = "ps_f_d_energy_" + str_psphi + "_" + str_psomega;
	TCanvas *c_f_d_energy = new TCanvas(str_c_tmp, str_c_tmp);
	TGraph *g_f_d_energy = new TGraph(ps_n_pts, &ps_f_d[0], &ps_energy[0]);
	g_f_d_energy->GetXaxis()->SetTitle("f_{d}");
	g_f_d_energy->GetYaxis()->SetTitle("scaled energy #frac{2E}{L^{2}m}");
	g_f_d_energy->Draw("PA");
	c_f_d_energy->SaveAs("./plot/" + str_c_tmp + ".pdf");

	str_c_tmp = "ps_theta_omega_" + str_psphi + "_" + str_psomega;
	TCanvas *c_theta_omega = new TCanvas(str_c_tmp, str_c_tmp);
	TGraph *g_theta_omega = new TGraph(ps_n_pts, &ps_theta[0], &ps_omega[0]);
	g_theta_omega->GetXaxis()->SetTitle("#theta");
	g_theta_omega->GetYaxis()->SetTitle("#omega");
	g_theta_omega->Draw("PA");
	c_theta_omega->SaveAs("./plot/" + str_c_tmp + ".pdf");

	str_c_tmp = "ps_theta_energy_" + str_psphi + "_" + str_psomega;
	TCanvas *c_theta_energy = new TCanvas(str_c_tmp, str_c_tmp);
	TGraph *g_theta_energy = new TGraph(ps_n_pts, &ps_theta[0], &ps_energy[0]);
	g_theta_energy->GetXaxis()->SetTitle("#theta");
	g_theta_energy->GetYaxis()->SetTitle("scaled energy #frac{2E}{L^{2}m}");
	g_theta_energy->Draw("PA");
	c_theta_energy->SaveAs("./plot/" + str_c_tmp + ".pdf");

	str_c_tmp = "ps_omega_energy_" + str_psphi + "_" + str_psomega;
	TCanvas *c_omega_energy = new TCanvas(str_c_tmp, str_c_tmp);
	TGraph *g_omega_energy = new TGraph(ps_n_pts, &ps_omega[0], &ps_energy[0]);
	g_omega_energy->GetXaxis()->SetTitle("#omega");
	g_omega_energy->GetYaxis()->SetTitle("scaled energy #frac{2E}{L^{2}m}");
	g_omega_energy->Draw("PA");
	c_omega_energy->SaveAs("./plot/" + str_c_tmp + ".pdf");
}

//-------------------------------------------------------------------------//

int main() {

	gStyle->SetTitleOffset(1.4,"y");
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);

	// (1)
	poincare_section(0, 2./3);

	chao_plot(1.40, 0.2, 0, 0, 2./3);
	chao_plot(1.44, 0.2, 0, 0, 2./3);
	chao_plot(1.20, 0.2, 0, 0, 2./3);

	cout << "(1) done" << endl;

	// (2)
	chao_plot(1.40, 1.2, 0, 0, 2./3);
	chao_plot(1.44, 1.2, 0, 0, 2./3);
	chao_plot(1.20, 1.2, 0, 0, 2./3);

	chao_plot(1.40, 0, 1.2, 0, 2./3);
	chao_plot(1.44, 0, 1.2, 0, 2./3);
	chao_plot(1.20, 0, 1.2, 0, 2./3);

	cout << "(2) done" << endl;

	// (3)
	chao_plot(1.40, 0.2, 0, M_PI/2, 2./3);
	chao_plot(1.44, 0.2, 0, M_PI/2, 2./3);
	chao_plot(1.20, 0.2, 0, M_PI/2, 2./3);

	chao_plot(1.40, 0.2, 0, M_PI, 2./3);
	chao_plot(1.44, 0.2, 0, M_PI, 2./3);
	chao_plot(1.20, 0.2, 0, M_PI, 2./3);

	chao_plot(1.40, 0.2, 0, M_PI*37/180, 2./3);
	chao_plot(1.44, 0.2, 0, M_PI*37/180, 2./3);
	chao_plot(1.20, 0.2, 0, M_PI*37/180, 2./3);

	cout << "(3) done" << endl;

	// (4)
	double omega_1 = M_PI/2./ellint_sin_half_theta(0.2);
	cout << M_PI/2./ellint_sin_half_theta(0.2) << endl;
	chao_plot(1.40, 0.2, 0, 0, 1);
	chao_plot(1.44, 0.2, 0, 0, 1);
	chao_plot(1.20, 0.2, 0, 0, 1);

	chao_plot(1.40, 0.2, 0, 0, omega_1);
	chao_plot(1.44, 0.2, 0, 0, omega_1);
	chao_plot(1.20, 0.2, 0, 0, omega_1);

	cout << "(4) done" << endl;

	return 0;
}

//-------------------------------------------------------------------------//
