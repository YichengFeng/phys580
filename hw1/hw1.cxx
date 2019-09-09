#include <iostream>
#include <vector>
#include "sequential_decay.h"
#include "bicycling.h"
#include "projectile.h"
#include "projectile_3d.h"

#include <TROOT.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TF1.h>
#include <TStyle.h>
#include <TMultiGraph.h>
#include <TH3.h>

using namespace std;


//-------------------------------------------------------------------------//

void sequential_decay_tau_A_B(double tau_A_B, double Na0, double Nb0) {

	SequentialDecay sd;
	SequentialDecay::set_tau_A_B(tau_A_B);

	//sd.set_dt(0.1);
	sd.set_N_A_start(Na0);
	sd.set_N_B_start(Nb0);
	sd.cal();
	vector<double> t = sd.get_t();
	vector< vector<double> > x = sd.get_x();
	int n_stps = sd.get_n_stps();

	sd.cal_exact();
	vector<double> t_exact = sd.get_t_exact();
	vector< vector<double> > x_exact = sd.get_x_exact();
	int n_stps_exact = sd.get_n_stps_exact();
	double N_A_start = sd.get_N_A_start();
	double N_B_start = sd.get_N_B_start();

	TGraph *g_A = new TGraph(n_stps, &t[0], &x[0][0]);
	g_A->SetLineColor(kBlue);
	g_A->SetMarkerColor(kBlue);
	TGraph *g_B = new TGraph(n_stps, &t[0], &x[1][0]);
	g_B->SetLineColor(kRed);
	g_B->SetMarkerColor(kRed);

	TGraph *g_A_exact = new TGraph(n_stps_exact, &t_exact[0], &x_exact[0][0]);
	g_A_exact->SetLineColor(kBlue);
	g_A_exact->SetLineWidth(2);
	g_A_exact->SetLineStyle(7);
	g_A_exact->SetMarkerColor(kBlue);
	TGraph *g_B_exact = new TGraph(n_stps_exact, &t_exact[0], &x_exact[1][0]);
	g_B_exact->SetLineColor(kRed);
	g_B_exact->SetLineWidth(2);
	g_B_exact->SetLineStyle(7);
	g_B_exact->SetMarkerColor(kRed);

	TMultiGraph *mg = new TMultiGraph();
	mg->Add(g_A_exact);
	mg->Add(g_B_exact);
	mg->Add(g_A);
	mg->Add(g_B);

	TString s_c = Form("AB_%.2d_%.0f_%.0f", (int)(tau_A_B*10), N_A_start, N_B_start); 
	TCanvas *c_AB = new TCanvas(s_c, s_c);
	mg->GetXaxis()->SetTitle("scaled time: t/#tau_{B}");
	mg->GetYaxis()->SetTitle("number of nuclei: N");
	mg->Draw("LA");

	TLegend *l_AB = new TLegend(0.65,0.70,0.9,0.9);
	l_AB->SetFillStyle(0);
	l_AB->SetBorderSize(0);
	//l_AB->SetHeader(Form("#tau_{A}/#tau_{B} = %.3f", tau_A_B));
	l_AB->AddEntry(g_A_exact, "exact N_{A}", "l");
	l_AB->AddEntry(g_B_exact, "exact N_{B}", "l");
	l_AB->AddEntry(g_A, "approximate N_{A}", "l");
	l_AB->AddEntry(g_B, "approximate N_{B}", "l");
	l_AB->Draw("same");

	TString s_para1 = Form("#tau_{A}/#tau_{B} = %.3f", tau_A_B);
	TString s_para2 = Form("N_{A}(0) = %.0f, N_{B}(0) = %.0f", N_A_start, N_B_start);
	TString s_para3;
	if(N_A_start/tau_A_B>N_B_start) {
		s_para3 = Form("N_{A}(0)/#tau_{A} > N_{B}(0)/#tau_{B}");
	} else if(N_A_start/tau_A_B<N_B_start) {
		s_para3 = Form("N_{A}(0)/#tau_{A} < N_{B}(0)/#tau_{B}");
	} else {
		s_para3 = Form("N_{A}(0)/#tau_{A} = N_{B}(0)/#tau_{B}");
	}
	TLegend *l_para = new TLegend(0.1,0.75,0.55,0.9);
	l_para->SetFillStyle(0);
	l_para->SetBorderSize(0);
	l_para->AddEntry((TObject*)0, s_para1, "");
	l_para->AddEntry((TObject*)0, s_para2, "");
	l_para->AddEntry((TObject*)0, s_para3, "");
	l_para->Draw("same");

	c_AB->SaveAs("./plot/"+s_c+".pdf");
}

//-------------------------------------------------------------------------//

void sequential_decay_dt() {
	const int ndt = 3;
	double dts[ndt] = {0.02, 0.1, 0.5};
	double style[ndt] = {1, 2, 7};

	double N_A_start = 1000;
	double N_B_start = 300;
	double tau_A_B = 1;

	SequentialDecay::set_tau_A_B(tau_A_B);

	TCanvas *c_dt = new TCanvas("dt", "dt");
	TMultiGraph *mg = new TMultiGraph();
	TGraph *g_A[ndt];
	TGraph *g_B[ndt];
	TLegend *l_error = new TLegend(0.55, 0.7, 0.9, 0.9);
	l_error->SetFillStyle(0);
	l_error->SetNColumns(2);

	for(int idt=0; idt<ndt; idt++) {
		SequentialDecay sd;
		sd.set_N_A_start(1000);
		sd.set_N_B_start(300);
		sd.set_dt(dts[idt]);
		sd.cal();
		sd.cal_exact();
		sd.cal_error();
		int n_stps = sd.get_n_stps();

		vector<double> t = sd.get_t();
		vector< vector<double> > x_error = sd.get_x_error();

		g_A[idt] = new TGraph(n_stps, &t[0], &x_error[0][0]);
		g_A[idt]->SetLineColor(kBlue);
		g_A[idt]->SetMarkerColor(kBlue);
		g_A[idt]->SetLineStyle(style[idt]);
		l_error->AddEntry(g_A[idt], Form("N_{A}, #Deltat = %.2f#tau_{B}",dts[idt]), "l");
		mg->Add(g_A[idt]);
		g_B[idt] = new TGraph(n_stps, &t[0], &x_error[1][0]);
		g_B[idt]->SetLineColor(kRed);
		g_B[idt]->SetMarkerColor(kRed);
		g_B[idt]->SetLineStyle(style[idt]);
		l_error->AddEntry(g_B[idt], Form("N_{B}, #Deltat = %.2f#tau_{B}",dts[idt]), "l");
		mg->Add(g_B[idt]);
	}
	mg->GetXaxis()->SetTitle("scaled time: t/#tau_{B}");
	mg->GetYaxis()->SetTitle("relative error: #DeltaN/N");
	mg->Draw("LA");
	l_error->Draw("same");

	TString s_para1 = Form("#tau_{A}/#tau_{B} = %.3f", tau_A_B);
	TString s_para2 = Form("N_{A}(0) = %.0f, N_{B}(0) = %.0f", N_A_start, N_B_start);
	TString s_para3;
	if(N_A_start/tau_A_B>N_B_start) {
		s_para3 = Form("N_{A}(0)/#tau_{A} > N_{B}(0)/#tau_{B}");
	} else if(N_A_start/tau_A_B<N_B_start) {
		s_para3 = Form("N_{A}(0)/#tau_{A} < N_{B}(0)/#tau_{B}");
	} else {
		s_para3 = Form("N_{A}(0)/#tau_{A} = N_{B}(0)/#tau_{B}");
	}
	TLegend *l_para = new TLegend(0.1,0.75,0.55,0.9);
	l_para->SetFillStyle(0);
	l_para->SetBorderSize(0);
	l_para->AddEntry((TObject*)0, s_para1, "");
	l_para->AddEntry((TObject*)0, s_para2, "");
	l_para->AddEntry((TObject*)0, s_para3, "");
	l_para->Draw("same");

	c_dt->SaveAs("./plot/dt_10_1000_300.pdf");
}

//-------------------------------------------------------------------------//

void bicycling() {
	const int nVar = 3;
	double A[nVar] = {0.33, 0.33, 0.33*0.7};
	double P[nVar] = {400, 234, 234};
	Int_t color[nVar] = {kBlue, kRed, kGreen};

	TCanvas *c = new TCanvas("c_bicycling", "c_bicycling");
	TGraph *g[nVar];
	TMultiGraph *mg = new TMultiGraph();
	TLegend *l = new TLegend(0.7, 0.7, 0.9, 0.9);
	l->SetFillStyle(0);
	l->SetBorderSize(0);

	for(int iVar=0; iVar<nVar; iVar++) {
		Bicycling bc;
		bc.set_A(A[iVar]);
		bc.set_P(P[iVar]);
		bc.cal();
		vector<double> t = bc.get_t();
		vector< vector<double> > x = bc.get_x();
		int n_stps = bc.get_n_stps();

		g[iVar] = new TGraph(n_stps, &t[0], &x[0][0]);
		g[iVar]->SetLineColor(color[iVar]);
		mg->Add(g[iVar]);
		l->AddEntry(g[iVar], Form("A = %.2f m^{2}; P = %.0f W", A[iVar], P[iVar]), "l");
	}
	mg->GetXaxis()->SetTitle("t (s)");
	mg->GetYaxis()->SetTitle("v (m/s)");
	mg->GetYaxis()->SetRangeUser(2, 20);
	mg->Draw("LA");
	l->Draw("same");

	c->SaveAs("./plot/bicycling.pdf");
}

//-------------------------------------------------------------------------//

void projectile(int isOpt) {
	const int nVar = 3; //
	double T0[nVar] = {300, 300, 270};
	int Mod[nVar] = {2, 3, 3};
	Int_t color[nVar] = {kBlue, kRed, kGreen};
	TString strMod[nVar] = {"isothermal model", "adiabatic model", "adiabatic model"};

	TCanvas *c = new TCanvas("c_projectile", "c_projectile");
	TGraph *g[nVar];
	double theta[nVar];
	double range[nVar];
	TMultiGraph *mg = new TMultiGraph();
	TLegend *l = new TLegend(0.55,0.65,0.9,0.9);
	l->SetFillStyle(0);
	l->SetBorderSize(0);
	
	for(int iVar=0; iVar<nVar; iVar++) {
	// cout << "modes:" << endl;
	// cout << "0: no drag; 1: constant drag; 2: isothermal drag; 3: adiabatic drag." << endl;
		Projectile pj;
		pj.set_mode(Mod[iVar]);
		pj.set_T0(T0[iVar]);
		pj.set_v_theta_start(700,45);
		pj.set_rk_order(1); // RK1 aka Euler

		if(isOpt) {
			theta[iVar] = pj.search_theta_for_max_range(700);
		} else {
			theta[iVar] = 45;
			pj.cal();
			pj.cal_range();
		}
		range[iVar] = pj.get_range();
		vector<double> t = pj.get_t();
		vector< vector<double> > x = pj.get_x();
		int n_stps = pj.get_n_stps();
		g[iVar] = new TGraph(n_stps, &x[0][0], &x[1][0]);
		g[iVar]->SetLineColor(color[iVar]);
		mg->Add(g[iVar]);

		TString str = "#splitline{"+strMod[iVar]+Form(" with T0 = %.0f K}{theta = %.1f^{#circ}, range = %.0f m}", T0[iVar], theta[iVar], range[iVar]);
		l->AddEntry(g[iVar], str, "l");
	}

	mg->GetXaxis()->SetTitle("x (m)");
	mg->GetYaxis()->SetTitle("y (m)");
	mg->GetYaxis()->SetRangeUser(0,14000);
	mg->GetXaxis()->SetLimits(0,30000);
	mg->Draw("LA");
	l->Draw("same");
	if(isOpt) {
		c->SaveAs("./plot/projectile_opt.pdf");
	} else {
		c->SaveAs("./plot/projectile_raw.pdf");
	}

}

//-------------------------------------------------------------------------//

void projectile_3d(int isOpt) {
	const int nVar = 4; //
	int Mod[nVar] = {0, 1, 2, 3};
	Int_t color[nVar] = {kBlack, kBlue, kRed, kGreen};
	TString strMod[nVar] = {"no drag", "constant drag", "isothermal model", "adiabatic model"};

	TCanvas *cxy = new TCanvas("cxy_projectile_3d", "cxy_projectile_3d");
	TCanvas *cxz = new TCanvas("cxz_projectile_3d", "cxz_projectile_3d");
	TCanvas *cxyz = new TCanvas("cxyz_projectile_3d", "cxyz_projectile_3d");
	TGraph *gxy[nVar];
	TGraph *gxz[nVar];
	TGraph2D *gxyz[nVar];
	double theta[nVar];
	double range[nVar];
	TMultiGraph *mgxy = new TMultiGraph();
	TMultiGraph *mgxz = new TMultiGraph();
	TLegend *lxy = new TLegend(0.55,0.6,0.9,0.9);
	lxy->SetFillStyle(0);
	lxy->SetBorderSize(0);
	TLegend *lxz = new TLegend(0.1,0.6,0.45,0.9);
	lxz->SetFillStyle(0);
	lxz->SetBorderSize(0);
	
	for(int iVar=0; iVar<nVar; iVar++) {
	// cout << "modes:" << endl;
	// cout << "0: no drag; 1: constant drag; 2: isothermal drag; 3: adiabatic drag." << endl;
		Projectile3D pj;
		//pj.set_dt(0.0005);
		pj.set_mode(Mod[iVar]);
		pj.set_v_theta_start(700,45);
		pj.set_rk_order(4); // RK1 aka Euler
		pj.set_omega_polar(2*M_PI/24/3600, 40+25.0/60, -45);

		if(isOpt) {
			theta[iVar] = pj.search_theta_for_max_range(700);
		} else {
			theta[iVar] = 45;
			pj.cal();
			pj.cal_range();
		}
		range[iVar] = pj.get_range();
		vector<double> t = pj.get_t();
		vector< vector<double> > x = pj.get_x();
		int n_stps = pj.get_n_stps();
		gxy[iVar] = new TGraph(n_stps, &x[0][0], &x[1][0]);
		gxy[iVar]->SetLineColor(color[iVar]);
		gxz[iVar] = new TGraph(n_stps, &x[0][0], &x[2][0]);
		gxz[iVar]->SetLineColor(color[iVar]);
		mgxy->Add(gxy[iVar]);
		mgxz->Add(gxz[iVar]);
		gxyz[iVar] = new TGraph2D(n_stps, &x[2][0], &x[0][0], &x[1][0]);
		gxyz[iVar]->SetLineColor(color[iVar]);
		gxyz[iVar]->SetMarkerColor(color[iVar]);

		TString str = "#splitline{"+strMod[iVar]+Form("}{theta = %.1f^{#circ}, range = %.0f m}", theta[iVar], range[iVar]);
		lxy->AddEntry(gxy[iVar], str, "l");
		lxz->AddEntry(gxz[iVar], str, "l");
	}

	cxy->cd();
	mgxy->GetXaxis()->SetTitle("x (m)");
	mgxy->GetYaxis()->SetTitle("y (m)");
	mgxy->GetYaxis()->SetRangeUser(0,18000);
	mgxy->GetXaxis()->SetLimits(0,55000);
	mgxy->Draw("LA");
	lxy->Draw("same");
	cxz->cd();
	mgxz->GetXaxis()->SetTitle("x (m)");
	mgxz->GetYaxis()->SetTitle("z (m)");
	//mgxz->GetYaxis()->SetRangeUser(0,14000);
	mgxz->GetXaxis()->SetLimits(0,55000);
	mgxz->Draw("LA");
	lxz->Draw("same");
	cxyz->cd();
	TH3F *frame = new TH3F("frame", "frame", 1,0,350, 1,0,55000, 1,0,14000);
	frame->Draw();
	frame->GetXaxis()->SetTitle("z (m)");
	frame->GetYaxis()->SetTitle("x (m)");
	frame->GetZaxis()->SetTitle("y (m)");
	frame->GetXaxis()->SetTitleOffset(1.6);
	frame->GetYaxis()->SetTitleOffset(1.6);
	frame->GetZaxis()->SetTitleOffset(1.6);
	for(int iVar=0; iVar<nVar; iVar++) {
		gxyz[iVar]->Draw("P same");
	}
	lxy->Draw("same");
	if(isOpt) {
		cxy->SaveAs("./plot/projectile_3d_xy_opt.pdf");
		cxz->SaveAs("./plot/projectile_3d_xz_opt.pdf");
		cxyz->SaveAs("./plot/projectile_3d_xyz_opt.pdf");
	} else {
		cxy->SaveAs("./plot/projectile_3d_xy_raw.pdf");
		cxz->SaveAs("./plot/projectile_3d_xz_raw.pdf");
		cxyz->SaveAs("./plot/projectile_3d_xyz_raw.pdf");
	}

}
//-------------------------------------------------------------------------//

int main() {

	gStyle->SetTitleOffset(1.4,"y");
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);

	// (1) sequential decay
	sequential_decay_tau_A_B(1.0/3, 1000, 300);
	sequential_decay_tau_A_B(1.0  , 1000, 300);
	sequential_decay_tau_A_B(3.0  , 1000, 300);
	sequential_decay_tau_A_B(1.0/3, 300, 1000);
	sequential_decay_tau_A_B(1.0  , 300, 1000);
	sequential_decay_tau_A_B(3.0  , 300, 1000);

	sequential_decay_dt();

	// (2) bicycling
	bicycling();
	
	// (3) projectile
	projectile(0);
	projectile(1);

	// (4) projectile_3d
	projectile_3d(0);

	return 0;
}
