#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include "runge_kutta.h"
#include "pendulum.h"
#include "anharmonics.h"

#include <TROOT.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TMath.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <TLine.h>

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

void chao_plot(double f_d, double q, double theta_start, double omega_start, double psphi, double psomega) {

	Pendulum pd;
	pd.set_mode(1); // physics pendulum
	pd.set_alg(4);
	//pd.set_theta_start(0.2);
	//pd.set_omega_start(0);
	pd.set_theta_start(theta_start);
	pd.set_omega_start(omega_start);
	pd.set_q(q);
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
	TString str_q = Form("q_%.3d", int(100*q));
	TString str_theta_start = Form("theta0_%.3d", int(100*fabs(theta_start)));
	TString str_omega_start = Form("omega0_%.3d", int(100*fabs(omega_start)));
	TString str_psphi = Form("psphi_%.3d", int(psphi*100));
	TString str_psomega = Form("psomega_%.3d", int(psomega*100));
	TString str_c_tmp;

	int n_tail = int(0.975*n_stps);
	int n_head = int(0.025*n_stps);

	TString str_adj = str_f_d + "_" + str_q + "_"  + str_theta_start + "_" + str_omega_start + "_" + str_psphi + "_" + str_psomega;

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

void poincare_section_scan_q(double f_d) {

	TString str_f_d = Form("fd_%.3d", int(f_d*100));

	double omega_d = 2./3;
	//double f_d = 1.30;
	double dt = 0.04;
	double T_d = 2*M_PI/omega_d;
	double psphi = 0;
	double psomega = 2./3;
	double q = 0.5;

	vector<double> ps_q;
	vector<double> ps_theta;
	vector<double> ps_omega;
	vector<double> ps_energy;
	int ps_n_pts = 0;

	double q_min = 0;
	double q_max = 1;
	int n_q = 400;

	for(int i=0; i<n_q; i++) {

		q = q_min+(q_max-q_min)/n_q*i;

		Pendulum pd;
		pd.set_mode(1); // physics pendulum
		pd.set_alg(4);
		pd.set_theta_start(0.2);
		pd.set_omega_start(0);
		pd.set_q(q);
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
			ps_q.push_back(q);
			ps_theta.push_back(x[0][i_stps]);
			ps_omega.push_back(x[1][i_stps]);
			ps_energy.push_back(energy[i_stps]);
			ps_n_pts ++;
		}
	}

	TString str_adj = str_f_d;
	TString str_c_tmp;

	str_c_tmp = "ps_q_theta_" + str_adj;
	TCanvas *c_q_theta = new TCanvas(str_c_tmp, str_c_tmp);
	TGraph *g_q_theta = new TGraph(ps_n_pts, &ps_q[0], &ps_theta[0]);
	g_q_theta->GetXaxis()->SetTitle("q");
	g_q_theta->GetYaxis()->SetTitle("#theta");
	g_q_theta->Draw("PA");
	c_q_theta->SaveAs("./plot/" + str_c_tmp + ".pdf");

	str_c_tmp = "ps_q_omega_" + str_adj;
	TCanvas *c_q_omega = new TCanvas(str_c_tmp, str_c_tmp);
	TGraph *g_q_omega = new TGraph(ps_n_pts, &ps_q[0], &ps_omega[0]);
	g_q_omega->GetXaxis()->SetTitle("q");
	g_q_omega->GetYaxis()->SetTitle("#omega");
	g_q_omega->Draw("PA");
	c_q_omega->SaveAs("./plot/" + str_c_tmp + ".pdf");

	str_c_tmp = "ps_q_energy_" + str_adj;
	TCanvas *c_q_energy = new TCanvas(str_c_tmp, str_c_tmp);
	TGraph *g_q_energy = new TGraph(ps_n_pts, &ps_q[0], &ps_energy[0]);
	g_q_energy->GetXaxis()->SetTitle("q");
	g_q_energy->GetYaxis()->SetTitle("scaled energy #frac{2E}{L^{2}m}");
	g_q_energy->Draw("PA");
	c_q_energy->SaveAs("./plot/" + str_c_tmp + ".pdf");
/*
	str_c_tmp = "ps_theta_omega_" + str_adj;
	TCanvas *c_theta_omega = new TCanvas(str_c_tmp, str_c_tmp);
	TGraph *g_theta_omega = new TGraph(ps_n_pts, &ps_theta[0], &ps_omega[0]);
	g_theta_omega->GetXaxis()->SetTitle("#theta");
	g_theta_omega->GetYaxis()->SetTitle("#omega");
	g_theta_omega->Draw("PA");
	c_theta_omega->SaveAs("./plot/" + str_c_tmp + ".pdf");

	str_c_tmp = "ps_theta_energy_" + str_adj;
	TCanvas *c_theta_energy = new TCanvas(str_c_tmp, str_c_tmp);
	TGraph *g_theta_energy = new TGraph(ps_n_pts, &ps_theta[0], &ps_energy[0]);
	g_theta_energy->GetXaxis()->SetTitle("#theta");
	g_theta_energy->GetYaxis()->SetTitle("scaled energy #frac{2E}{L^{2}m}");
	g_theta_energy->Draw("PA");
	c_theta_energy->SaveAs("./plot/" + str_c_tmp + ".pdf");

	str_c_tmp = "ps_omega_energy_" + str_adj;
	TCanvas *c_omega_energy = new TCanvas(str_c_tmp, str_c_tmp);
	TGraph *g_omega_energy = new TGraph(ps_n_pts, &ps_omega[0], &ps_energy[0]);
	g_omega_energy->GetXaxis()->SetTitle("#omega");
	g_omega_energy->GetYaxis()->SetTitle("scaled energy #frac{2E}{L^{2}m}");
	g_omega_energy->Draw("PA");
	c_omega_energy->SaveAs("./plot/" + str_c_tmp + ".pdf");
*/
}

//-------------------------------------------------------------------------//

void anharmonics_amplitude_scan(double alpha, double k=1) {

	const int n_amp = 3;
	double theta_start[n_amp] = {0.2, 0.6, 1.0};
	double omega_start[n_amp] = {0, 0, 0};
	double period[n_amp];
	double dt = 0.02;
	double t_end = 2*M_PI*9;

	TGraph *g_t_theta[n_amp];
	TGraph *g_t_omega[n_amp];

	for(int i_amp=0; i_amp<n_amp; i_amp++) {
		Anharmonics ah;
		ah.set_alpha(alpha);
		ah.set_k(k);
		ah.set_theta_start(theta_start[i_amp]);
		ah.set_omega_start(omega_start[i_amp]);
		ah.set_alg(0); // Euler-Cromer
		ah.set_dt(dt);
		ah.set_t_end(t_end);

		ah.cal();
		period[i_amp] = ah.cal_period();
		cout << setprecision(10) << period[i_amp]/M_PI << " pi" << endl;

		vector<double> t = ah.get_t();
		vector< vector<double> > x = ah.get_x();
		int n_stps = ah.get_n_stps();

		g_t_theta[i_amp] = new TGraph(n_stps, &t[0], &x[0][0]);
		g_t_omega[i_amp] = new TGraph(n_stps, &t[0], &x[1][0]);
	}

	Int_t color[n_amp] = {kBlue, kRed, kGreen};
	TMultiGraph *mg_t_theta = new TMultiGraph();
	TMultiGraph *mg_t_omega = new TMultiGraph();
	TLegend *l_t_theta = new TLegend(0.65,0.7,0.9,0.9);
	l_t_theta->SetFillStyle(0);
	l_t_theta->SetBorderSize(0);
	TLegend *l_t_omega = new TLegend(0.65,0.7,0.9,0.9);
	l_t_omega->SetFillStyle(0);
	l_t_omega->SetBorderSize(0);

	TLegend *l_alpha = new TLegend(0.1,0.8,0.3,0.9);
	l_alpha->SetFillStyle(0);
	l_alpha->SetBorderSize(0);
	l_alpha->AddEntry((TObject*)0, Form("#alpha=%.1f", alpha), "");

	for(int i_amp=0; i_amp<n_amp; i_amp++) {
		TString str_tmp = Form("x_{m}=%.2f: T=%.2f#pi", theta_start[i_amp], period[i_amp]/M_PI);
		g_t_theta[i_amp]->SetLineColor(color[i_amp]);
		g_t_theta[i_amp]->SetMarkerColor(color[i_amp]);
		l_t_theta->AddEntry(g_t_theta[i_amp], str_tmp, "l");
		mg_t_theta->Add(g_t_theta[i_amp]);

		g_t_omega[i_amp]->SetLineColor(color[i_amp]);
		g_t_omega[i_amp]->SetMarkerColor(color[i_amp]);
		l_t_omega->AddEntry(g_t_omega[i_amp], str_tmp, "l");
		mg_t_omega->Add(g_t_omega[i_amp]);
	}

	TString str_alpha = Form("alpha%.1d", int(alpha));
	TString str_k = Form("k%.1d", int(k));
	TString str_adj = str_alpha + "_" + str_k;

	TLine *zero = new TLine(0,0,t_end,0);
	zero->SetLineStyle(7);

	TCanvas *c_t_theta = new TCanvas("t_theta_" + str_adj, "t_theta_" + str_adj);
	mg_t_theta->GetXaxis()->SetTitle("t (s)");
	mg_t_theta->GetYaxis()->SetTitle("x (m)");
	mg_t_theta->GetYaxis()->SetRangeUser(-1.2,1.5);
	mg_t_theta->Draw("LA");
	l_t_theta->Draw("same");
	l_alpha->Draw("same");
	zero->Draw("same");
	c_t_theta->SaveAs("./plot/t_theta_" + str_adj + ".pdf");

	TCanvas *c_t_omega = new TCanvas("t_omega_" + str_adj, "t_omega_" + str_adj);
	mg_t_omega->GetXaxis()->SetTitle("t (s)");
	mg_t_omega->GetYaxis()->SetTitle("v (m/s)");
	mg_t_omega->GetYaxis()->SetRangeUser(-1.2,1.5);
	mg_t_omega->Draw("LA");
	l_t_omega->Draw("same");
	l_alpha->Draw("same");
	zero->Draw("same");
	c_t_omega->SaveAs("./plot/t_omega_" + str_adj + ".pdf");

}

//-------------------------------------------------------------------------//

void chao_delta_start() {

	double dtheta0 = 0.001;
	double f_d = 1.2;
	double theta_start1 = 0.2;
	double theta_start2 = theta_start1 + dtheta0;
	double omega_start1 = 0;
	double omega_start2 = 0;

	Pendulum pd1;
	pd1.set_mode(1); // physics pendulum
	pd1.set_alg(4);
	pd1.set_theta_start(theta_start1);
	pd1.set_omega_start(omega_start1);
	pd1.set_q(0.5);
	pd1.set_F(f_d);
	pd1.set_O(2./3);
	pd1.set_dt(0.02);
	pd1.set_n_periods(20);
	pd1.cal();

	vector<double> t1 = pd1.get_t();
	vector< vector<double> > x1 = pd1.get_x();
	int n_stps1 = pd1.get_n_stps();

	Pendulum pd2;
	pd2.set_mode(1); // physics pendulum
	pd2.set_alg(4);
	pd2.set_theta_start(theta_start2);
	pd2.set_omega_start(omega_start2);
	pd2.set_q(0.5);
	pd2.set_F(f_d);
	pd2.set_O(2./3);
	pd2.set_dt(0.02);
	pd2.set_n_periods(20);
	pd2.cal();

	vector<double> t2 = pd2.get_t();
	vector< vector<double> > x2 = pd2.get_x();
	int n_stps2 = pd2.get_n_stps();

	vector<double> dtheta;
	for(int i=0; i<t1.size(); i++) {
		double tmp_dtheta = fabs(x2[0][i]-x1[0][i]);
		while(tmp_dtheta > 2*M_PI) tmp_dtheta -= 2*M_PI;
		dtheta.push_back(tmp_dtheta);
	}

	theta_to_twopi(x1[0]);
	theta_to_twopi(x2[0]);

	TGraph *g_t_theta_pd1 = new TGraph(n_stps1, &t1[0], &x1[0][0]);
	g_t_theta_pd1->SetLineColor(kBlue);
	g_t_theta_pd1->SetMarkerColor(kBlue);
	g_t_theta_pd1->GetXaxis()->SetTitle("t (s)");
	g_t_theta_pd1->GetYaxis()->SetTitle("#theta (rad)");
	TGraph *g_t_theta_pd2 = new TGraph(n_stps2, &t2[0], &x2[0][0]);
	g_t_theta_pd2->SetLineColor(kRed);
	g_t_theta_pd2->SetMarkerColor(kRed);
	g_t_theta_pd2->GetXaxis()->SetTitle("t (s)");
	g_t_theta_pd2->GetYaxis()->SetTitle("#theta (rad)");
	TGraph *g_t_dtheta = new TGraph(n_stps1, &t1[0], &dtheta[0]);
	g_t_dtheta->GetXaxis()->SetTitle("t (s)");
	g_t_dtheta->GetYaxis()->SetTitle("#Delta#theta (rad)");

	TString str_f_d = Form("fd%.3d", int(100*f_d));
	TString str_adj = str_f_d;

	TCanvas *c_t_theta = new TCanvas("t_theta_" + str_adj, "t_theta_" + str_adj, 980,500);
	g_t_theta_pd1->Draw("LA");
	g_t_theta_pd2->Draw("L same");
	TLegend *l_t_theta = new TLegend(0.4,0.8,0.9,0.9);
	l_t_theta->SetNColumns(2);
	l_t_theta->SetFillStyle(0);
	l_t_theta->SetBorderSize(0);
	l_t_theta->AddEntry(g_t_theta_pd1, Form("#theta_{0}=%.4f", theta_start1), "l");
	l_t_theta->AddEntry(g_t_theta_pd2, Form("#theta_{0}=%.4f", theta_start2), "l");
	l_t_theta->Draw("same");
	c_t_theta->SaveAs("./plot/t_theta_" + str_adj + ".pdf");

	TF1 *func = new TF1("func", "1e-4*exp(0.16997258*x)", 0,75);
	func->SetLineStyle(7);

	TCanvas *c_t_dtheta = new TCanvas("t_dtheta_" + str_adj, "t_dtheta_" + str_adj);
	g_t_dtheta->Draw("LA");
	func->Draw("same");
	c_t_dtheta->SetLogy();
	c_t_dtheta->SaveAs("./plot/t_dtheta_" + str_adj + ".pdf");
}

//-------------------------------------------------------------------------//

void poincare_section_scan_fd(double psphi, double psomega) {

	TString str_psphi = Form("psphi_%.3d", int(psphi*100));
	TString str_psomega = Form("psomega_%.3d", int(psomega*100));
	TString str_q = Form("q_050");

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
	int n_f_d = 400;

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

		for(int i_stps=350/omega_d/dt; i_stps<n_stps; i_stps++) {
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
	TString str_adj = str_q;

	str_c_tmp = "ps_f_d_theta_" + str_adj;
	TCanvas *c_f_d_theta = new TCanvas(str_c_tmp, str_c_tmp);
	TGraph *g_f_d_theta = new TGraph(ps_n_pts, &ps_f_d[0], &ps_theta[0]);
	g_f_d_theta->GetXaxis()->SetTitle("f_{d}");
	g_f_d_theta->GetYaxis()->SetTitle("#theta");
	g_f_d_theta->Draw("PA");
	c_f_d_theta->SaveAs("./plot/" + str_c_tmp + ".pdf");

	str_c_tmp = "ps_f_d_omega_" + str_adj;
	TCanvas *c_f_d_omega = new TCanvas(str_c_tmp, str_c_tmp);
	TGraph *g_f_d_omega = new TGraph(ps_n_pts, &ps_f_d[0], &ps_omega[0]);
	g_f_d_omega->GetXaxis()->SetTitle("f_{d}");
	g_f_d_omega->GetYaxis()->SetTitle("#omega");
	g_f_d_omega->Draw("PA");
	c_f_d_omega->SaveAs("./plot/" + str_c_tmp + ".pdf");

	str_c_tmp = "ps_f_d_energy_" + str_adj;
	TCanvas *c_f_d_energy = new TCanvas(str_c_tmp, str_c_tmp);
	TGraph *g_f_d_energy = new TGraph(ps_n_pts, &ps_f_d[0], &ps_energy[0]);
	g_f_d_energy->GetXaxis()->SetTitle("f_{d}");
	g_f_d_energy->GetYaxis()->SetTitle("scaled energy #frac{2E}{L^{2}m}");
	g_f_d_energy->Draw("PA");
	c_f_d_energy->SaveAs("./plot/" + str_c_tmp + ".pdf");
}

//-------------------------------------------------------------------------//

int main() {

	gStyle->SetTitleOffset(1.4,"y");
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);

	// (1) 
	//poincare_section_scan_q(0.5);
	//poincare_section_scan_q(1.2);
	chao_plot(0.5, 0.05, 0.2, 0, 0, 2./3);
	chao_plot(0.5, 0.12, 0.2, 0, 0, 2./3);
	chao_plot(0.5, 0.50, 0.2, 0, 0, 2./3);
	chao_plot(1.2, 0.10, 0.2, 0, 0, 2./3);
	chao_plot(1.2, 0.37, 0.2, 0, 0, 2./3);
	chao_plot(1.2, 0.70, 0.2, 0, 0, 2./3);

	// (2)
	anharmonics_amplitude_scan(1, 1);
	anharmonics_amplitude_scan(3, 1);

	// (3)

	// (4)
	chao_delta_start();
	poincare_section_scan_fd(0, 2./3);

	// (5)


	return 0;
}

//-------------------------------------------------------------------------//
