#include <iostream>
#include <vector>
#include <cmath>
#include "pendulum.h"

#include <TCanvas.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TMultiGraph.h>

using namespace std;


//-------------------------------------------------------------------------//

void pendulum_simplified(int alg, TString str) {

	gStyle->SetTitleOffset(1.4,"y");
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);

	Pendulum pd;
	pd.set_theta_start(0.08);
	pd.set_omega_start(0);
	pd.set_n_periods(4.5);
	pd.set_alg(alg);

	pd.cal();

	vector<double> t = pd.get_t();
	vector< vector<double> > x = pd.get_x();
	vector< vector<double> > x_exact = pd.get_x_simplified_exact();
	vector< vector<double> > x_error = pd.get_x_simplified_error();
	vector<double> energy = pd.get_energy();
	vector<double> energy_exact = pd.get_energy_exact();
	vector<double> energy_error = pd.get_energy_error();
	int n_stps = pd.get_n_stps();

	TGraph *g_theta = new TGraph(n_stps, &t[0], &x[0][0]);
	TGraph *g_omega = new TGraph(n_stps, &t[0], &x[1][0]);
	TGraph *g_energy = new TGraph(n_stps, &t[0], &energy[0]);
	TGraph *g_theta_exact = new TGraph(n_stps, &t[0], &x_exact[0][0]);
	TGraph *g_omega_exact = new TGraph(n_stps, &t[0], &x_exact[1][0]);
	TGraph *g_energy_exact = new TGraph(n_stps, &t[0], &energy_exact[0]);
	TGraph *g_theta_error = new TGraph(n_stps, &t[0], &x_error[0][0]);
	TGraph *g_omega_error = new TGraph(n_stps, &t[0], &x_error[1][0]);
	TGraph *g_energy_error = new TGraph(n_stps, &t[0], &energy_error[0]);

	TCanvas *c_theta = new TCanvas("theta"+str, "theta"+str);
	g_theta->GetXaxis()->SetTitle("t (s)");
	g_theta->GetYaxis()->SetTitle("#theta (rad)");
	g_theta->SetLineColor(kRed);
	g_theta->Draw("LA");
	g_theta_exact->SetLineColor(kRed);
	g_theta_exact->SetLineWidth(2);
	g_theta_exact->SetLineStyle(7);
	g_theta_exact->Draw("L same");
	TLegend *l_theta = new TLegend(0.5,0.8,0.9,0.9);
	l_theta->SetNColumns(2);
	l_theta->SetFillStyle(0);
	l_theta->AddEntry(g_theta, "approximation", "l");
	l_theta->AddEntry(g_theta_exact, "exact", "l");
	l_theta->Draw("same");
	c_theta->SaveAs("./plot/theta_"+str+".pdf");
	// error
	TCanvas *c_theta_error = new TCanvas("theta_error"+str, "theta_error"+str);
	g_theta_error->GetXaxis()->SetTitle("t (s)");
	g_theta_error->GetYaxis()->SetTitle("#theta error (rad)");
	g_theta_error->SetLineColor(kRed);
	g_theta_error->Draw("LA");
	c_theta_error->SaveAs("./plot/theta_error_"+str+".pdf");

	TCanvas *c_omega = new TCanvas("omega"+str, "omega"+str);
	g_omega->GetXaxis()->SetTitle("t (s)");
	g_omega->GetYaxis()->SetTitle("#omega (rad/s)");
	g_omega->SetLineColor(kBlue);
	g_omega->Draw("LA");
	g_omega_exact->SetLineColor(kBlue);
	g_omega_exact->SetLineWidth(2);
	g_omega_exact->SetLineStyle(7);
	g_omega_exact->Draw("L same");
	TLegend *l_omega = new TLegend(0.5,0.8,0.9,0.9);
	l_omega->SetNColumns(2);
	l_omega->SetFillStyle(0);
	l_omega->AddEntry(g_omega, "approximation", "l");
	l_omega->AddEntry(g_omega_exact, "exact", "l");
	l_omega->Draw("same");
	c_omega->SaveAs("./plot/omega_"+str+".pdf");
	// error
	TCanvas *c_omega_error = new TCanvas("omega_error"+str, "omega_error"+str);
	g_omega_error->GetXaxis()->SetTitle("t (s)");
	g_omega_error->GetYaxis()->SetTitle("#omega error (rad/s)");
	g_omega_error->SetLineColor(kBlue);
	g_omega_error->Draw("LA");
	c_omega_error->SaveAs("./plot/omega_error_"+str+".pdf");

	TCanvas *c_energy = new TCanvas("energy"+str, "energy"+str);
	g_energy->GetXaxis()->SetTitle("t (s)");
	g_energy->GetYaxis()->SetTitle("#frac{2E}{L^{2}m}");
	g_energy->SetLineColor(kGreen);
	g_energy->Draw("LA");
	g_energy_exact->SetLineColor(kGreen);
	g_energy_exact->SetLineWidth(2);
	g_energy_exact->SetLineStyle(7);
	g_energy_exact->Draw("L same");
	TLegend *l_energy = new TLegend(0.5,0.8,0.9,0.9);
	l_energy->SetNColumns(2);
	l_energy->SetFillStyle(0);
	l_energy->AddEntry(g_energy, "approximation", "l");
	l_energy->AddEntry(g_energy_exact, "exact", "l");
	l_energy->Draw("same");
	c_energy->SaveAs("./plot/energy_"+str+".pdf");
	// error
	TCanvas *c_energy_error = new TCanvas("energy_error"+str, "energy_error"+str);
	g_energy_error->GetXaxis()->SetTitle("t (s)");
	g_energy_error->GetYaxis()->SetTitle("#frac{2E}{L^{2}m}");
	g_energy_error->SetLineColor(kGreen);
	g_energy_error->Draw("LA");
	c_energy_error->SaveAs("./plot/energy_error_"+str+".pdf");

	TCanvas *c_theta_omega_energy = new TCanvas("theta_omega_energy"+str, "theta_omega_enrgy"+str, 1900,600);
	c_theta_omega_energy->Divide(3,1,1e-11,1e-11);
	c_theta_omega_energy->cd(1);
	g_theta->Draw("LA");
	g_theta_exact->Draw("L same");
	l_theta->Draw("same");
	c_theta_omega_energy->cd(2);
	g_omega->Draw("LA");
	g_omega_exact->Draw("L same");
	l_omega->Draw("same");
	c_theta_omega_energy->cd(3);
	g_energy->Draw("LA");
	g_energy_exact->Draw("L same");
	l_energy->Draw("same");
	c_theta_omega_energy->SaveAs("./plot/theta_omega_energy_"+str+".pdf");
	
}

//-------------------------------------------------------------------------//

void pendulum_simplified_damped(double q, TString str) {

	Pendulum pd;
	pd.set_theta_start(0.08);
	pd.set_omega_start(0);
	pd.set_n_periods(3.5);
	pd.set_alg(2);
	pd.set_q(q);

	pd.cal();

	vector<double> t = pd.get_t();
	vector< vector<double> > x = pd.get_x();
	vector<double> energy = pd.get_energy();
	int n_stps = pd.get_n_stps();

	TGraph *g_theta = new TGraph(n_stps, &t[0], &x[0][0]);
	TGraph *g_omega = new TGraph(n_stps, &t[0], &x[1][0]);
	TGraph *g_energy = new TGraph(n_stps, &t[0], &energy[0]);

	TLegend *l_damp = new TLegend(0.6,0.6,0.9,0.85);
	l_damp->SetFillStyle(0);
	l_damp->AddEntry((TObject*)0, Form("q = %.2f#sqrt{g/L}",q), "");
	l_damp->AddEntry((TObject*)0, Form("g = 9.8 m/s^{2}"), "");
	l_damp->AddEntry((TObject*)0, Form("L = 9.8 m"), "");
	l_damp->AddEntry((TObject*)0, Form("#Omega_{0} = 1 rad/s"), "");
	l_damp->SetBorderSize(0);

	TCanvas *c_theta = new TCanvas("theta"+str, "theta"+str);
	g_theta->GetXaxis()->SetTitle("t (s)");
	g_theta->GetYaxis()->SetTitle("#theta (rad)");
	g_theta->SetLineColor(kRed);
	g_theta->Draw("LA");
	l_damp->Draw("same");
	c_theta->SaveAs("./plot/theta_"+str+".pdf");

	TCanvas *c_omega = new TCanvas("omega"+str, "omega"+str);
	g_omega->GetXaxis()->SetTitle("t (s)");
	g_omega->GetYaxis()->SetTitle("#omega (rad/s)");
	g_omega->SetLineColor(kBlue);
	g_omega->Draw("LA");
	l_damp->Draw("same");
	c_omega->SaveAs("./plot/omega_"+str+".pdf");

	TCanvas *c_energy = new TCanvas("energy"+str, "energy"+str);
	g_energy->GetXaxis()->SetTitle("t (s)");
	g_energy->GetYaxis()->SetTitle("#frac{2E}{L^{2}m}");
	g_energy->SetLineColor(kGreen);
	g_energy->Draw("LA");
	l_damp->Draw("same");
	c_energy->SaveAs("./plot/energy_"+str+".pdf");
}

//-------------------------------------------------------------------------//

void pendulum_simplified_driven(double od, double fd, TString strod) {

	Pendulum pd;
	pd.set_theta_start(0.08);
	pd.set_omega_start(0);
	pd.set_n_periods(6);
	pd.set_alg(2);
	pd.set_F(fd);
	pd.set_O(od);

	pd.cal();

	vector<double> t = pd.get_t();
	vector< vector<double> > x = pd.get_x();
	vector<double> energy = pd.get_energy();
	int n_stps = pd.get_n_stps();

	TGraph *g_theta = new TGraph(n_stps, &t[0], &x[0][0]);
	TGraph *g_omega = new TGraph(n_stps, &t[0], &x[1][0]);
	TGraph *g_energy = new TGraph(n_stps, &t[0], &energy[0]);

	TLegend *l_damp = new TLegend(0.6,0.6,0.9,0.85);
	l_damp->SetFillStyle(0);
	l_damp->AddEntry((TObject*)0, Form("f_{D} = %.2f", fd), "");
	l_damp->AddEntry((TObject*)0, Form("#Omega_{D} = %.1f rad/s", od), "");
	l_damp->AddEntry((TObject*)0, Form("g = 9.8 m/s^{2}"), "");
	l_damp->AddEntry((TObject*)0, Form("L = 9.8 m"), "");
	l_damp->AddEntry((TObject*)0, Form("#Omega_{0} = 1 rad/s"), "");
	l_damp->SetBorderSize(0);

	TString strfd = Form("%d", int(100*fd));

	TCanvas *c_theta = new TCanvas("theta"+strod+strfd, "theta"+strod+strfd);
	g_theta->GetXaxis()->SetTitle("t (s)");
	g_theta->GetYaxis()->SetTitle("#theta (rad)");
	g_theta->SetLineColor(kRed);
	g_theta->Draw("LA");
	l_damp->Draw("same");
	c_theta->SaveAs("./plot/theta_"+strod+strfd+".pdf");

	TCanvas *c_omega = new TCanvas("omega"+strod+strfd, "omega"+strod+strfd);
	g_omega->GetXaxis()->SetTitle("t (s)");
	g_omega->GetYaxis()->SetTitle("#omega (rad/s)");
	g_omega->SetLineColor(kBlue);
	g_omega->Draw("LA");
	l_damp->Draw("same");
	c_omega->SaveAs("./plot/omega_"+strod+strfd+".pdf");

	TCanvas *c_energy = new TCanvas("energy"+strod+strfd, "energy"+strod+strfd);
	g_energy->GetXaxis()->SetTitle("t (s)");
	g_energy->GetYaxis()->SetTitle("#frac{2E}{L^{2}m}");
	g_energy->SetLineColor(kGreen);
	g_energy->Draw("LA");
	l_damp->Draw("same");
	c_energy->SaveAs("./plot/energy_"+strod+strfd+".pdf");
}

//-------------------------------------------------------------------------//

void pendulum_physics() {

	TString str = "physics";

	Pendulum pd;
	pd.set_theta_start(0.8*M_PI);
	pd.set_omega_start(0);
	pd.set_n_periods(5.5);
	pd.set_mode(1);
	pd.set_alg(2);

	pd.cal();

	vector<double> t = pd.get_t();
	vector< vector<double> > x = pd.get_x();
	vector<double> energy = pd.get_energy();
	int n_stps = pd.get_n_stps();

	TGraph *g_theta = new TGraph(n_stps, &t[0], &x[0][0]);
	TGraph *g_omega = new TGraph(n_stps, &t[0], &x[1][0]);
	TGraph *g_energy = new TGraph(n_stps, &t[0], &energy[0]);

	TCanvas *c_theta = new TCanvas("theta"+str, "theta"+str);
	g_theta->GetXaxis()->SetTitle("t (s)");
	g_theta->GetYaxis()->SetTitle("#theta (rad)");
	g_theta->SetLineColor(kRed);
	g_theta->Draw("LA");
	c_theta->SaveAs("./plot/theta_"+str+".pdf");

	TCanvas *c_omega = new TCanvas("omega"+str, "omega"+str);
	g_omega->GetXaxis()->SetTitle("t (s)");
	g_omega->GetYaxis()->SetTitle("#omega (rad/s)");
	g_omega->SetLineColor(kBlue);
	g_omega->Draw("LA");
	c_omega->SaveAs("./plot/omega_"+str+".pdf");

	TCanvas *c_energy = new TCanvas("energy"+str, "energy"+str);
	g_energy->GetXaxis()->SetTitle("t (s)");
	g_energy->GetYaxis()->SetTitle("#frac{2E}{L^{2}m}");
	g_energy->SetLineColor(kGreen);
	g_energy->Draw("LA");
	c_energy->SaveAs("./plot/energy_"+str+".pdf");
}

//-------------------------------------------------------------------------//

void pendulum_physics_period() {

	int n_stps = 100;
	double theta_max = 0.9*M_PI;
	double dtheta = theta_max/(n_stps+1);

	vector<double> theta;
	vector<double> period;
	vector<double> period_analytical;

	for(int i=1; i<n_stps+1; i++) {
		theta.push_back(i*dtheta);
		Pendulum pd;
		pd.set_theta_start(i*dtheta);
		pd.set_omega_start(0);
		pd.set_n_periods(10);
		pd.set_mode(1);
		pd.set_alg(2);

		pd.cal();
		period.push_back(pd.cal_period());
		period_analytical.push_back(pd.cal_analytical_period());
	}

	TGraph *g_period = new TGraph(n_stps, &theta[0], &period[0]);
	TGraph *g_period_analytical = new TGraph(n_stps, &theta[0], &period_analytical[0]);

	TCanvas *c_period = new TCanvas("physics_period", "physics_period");
	g_period->SetLineColor(kBlue);
	g_period_analytical->SetLineColor(kBlack);
	g_period_analytical->SetLineWidth(2);
	g_period_analytical->SetLineStyle(7);
	TMultiGraph *mg = new TMultiGraph();
	mg->Add(g_period);
	mg->Add(g_period_analytical);
	mg->GetXaxis()->SetTitle("#theta_{m} (rad)");
	mg->GetYaxis()->SetTitle("T (s)");
	mg->Draw("LA");

	TLegend *l_period = new TLegend(0.1,0.7,0.4,0.9);
	l_period->SetFillStyle(0);
	l_period->SetBorderSize(0);
	l_period->AddEntry(g_period, "approximate period", "l");
	l_period->AddEntry(g_period_analytical, "analytical period", "l");
	l_period->Draw("same");

	c_period->SaveAs("./plot/physics_period.pdf");
}

//-------------------------------------------------------------------------//

int main() {

	const int nalg = 3;
	int algs[nalg] = {0, 1, 2};
	TString str_alg[nalg] = {"EulerCromer", "Euler", "RK2"};
	for(int ialg=0; ialg<nalg; ialg++) {
		pendulum_simplified(algs[ialg], str_alg[ialg]);
	}

	const int ndmp = 3;
	double qs[ndmp] = {0.3,2,4};
	TString str_dmp[ndmp] = {"underdamped", "critical", "overdamped"};
	for(int idmp=0; idmp<ndmp; idmp++) {
		pendulum_simplified_damped(qs[idmp], str_dmp[idmp]);
	}

	const int nod = 3;
	double ods[nod] = {0.3,1,1.4};
	const int nfd = 3;
	double fds[nfd] = {0.2,1,1.5};
	TString str_od[nod] = {"slow_drive", "resonance_drive", "quick_drive"};
	for(int iod=0; iod<nod; iod++) {
		for(int ifd=0; ifd<nfd; ifd++) {
			pendulum_simplified_driven(ods[iod], fds[ifd], str_od[iod]);
		}
	}

	pendulum_physics();

	pendulum_physics_period();

	return 0;
}
