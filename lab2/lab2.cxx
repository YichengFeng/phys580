#include <iostream>
#include <vector>
#include <cmath>
#include "projectile.h"

#include <TROOT.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TLine.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TString.h>
#include <TPaveLabel.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>

using namespace std;


void projectile_plot_rk(int rk_order, TString str_method) {

	Projectile::set_g(9.8);
	Projectile::set_B2_m(4e-5);
	Projectile::set_Y0(1e4);
	Projectile::set_a(6.5e-3);
	Projectile::set_T0(293);
	Projectile::set_gamma(1.4);

	Projectile proj;
	proj.set_t_start(0);
	proj.set_v_theta_start(700, 45);
	proj.set_rk_order(rk_order);
	if(rk_order==1) {
		proj.set_dt(0.0005);
	}

	cout << "dt = " << proj.get_dt() << " s" << endl;

	TString str_para1 = Form("g = %.1f m/s^{2}, B2/m = %.1e m^{-1}, Y0 = %.0e m", proj.get_g(), proj.get_B2_m(), proj.get_Y0());
	TString str_para2 = Form("a = %.1e K/m, T0 = %.0f K, #gamma = %.1f, v0 = 700 m/s", proj.get_a(), proj.get_T0(), proj.get_gamma());
	TString str_para3 = Form("#Deltat = %.4f s", proj.get_dt());

	const int nMod = 4;
	TString strMod[nMod] = {"no drag", "constant air density", "isothermal model", "adiabatic model"};
	Int_t colorMod[nMod] = {kBlack, kBlue, kRed, kGreen};

	vector<double> t[nMod];
	vector< vector<double> > x[nMod];
	TGraph *gxy[nMod];
	TMultiGraph *mgxy = new TMultiGraph();
	TCanvas *cxy = new TCanvas(str_method+"xy", str_method+"xy");
	TLegend *lxy = new TLegend(0.64, 0.6, 0.9, 0.9);
	cout << "modes:" << endl;
	cout << "0: no drag; 1: constant drag; 2: isothermal drag; 3: adiabatic drag." << endl;
	for(int iMod=0; iMod<nMod; iMod++) {
		cout << "----- mode " << iMod << " -----" << endl;
		proj.set_mode(iMod);
		proj.cal();
		proj.cal_range();
		cout << "range: " << proj.get_range() << endl;

		t[iMod] = proj.get_t();
		x[iMod] = proj.get_x();
		int npts = proj.get_n_stps();

		gxy[iMod] = new TGraph(npts, &x[iMod][0][0], &x[iMod][1][0]);
		gxy[iMod]->SetLineColor(colorMod[iMod]);
		gxy[iMod]->SetMarkerColor(colorMod[iMod]);
		TString str_tmp = Form("theta: 45.0^{#circ}, range: %.1f m", proj.get_range());
		lxy->AddEntry(gxy[iMod], "#splitline{"+strMod[iMod]+"}{"+str_tmp+"}", "lp");
		mgxy->Add(gxy[iMod]);

	}
	mgxy->GetXaxis()->SetTitle("x (m)");
	mgxy->GetYaxis()->SetTitle("y (m)");
	mgxy->Draw("PLA");
	mgxy->GetXaxis()->SetLimits(0,51000);
	mgxy->GetYaxis()->SetRangeUser(0,17500);
	lxy->Draw("same");
	TLine *zeroline = new TLine(0, 0, 50000, 0);
	zeroline->SetLineStyle(7);
	zeroline->Draw("same");
	TLegend *l_para = new TLegend(0.0,0.7,0.58,0.9);
	l_para->SetFillStyle(0);
	l_para->SetBorderSize(0);
	l_para->AddEntry((TObject*)0, str_method+" method", "");
	l_para->AddEntry((TObject*)0, str_para1, "");
	l_para->AddEntry((TObject*)0, str_para2, "");
	l_para->AddEntry((TObject*)0, str_para3, "");
	l_para->Draw();
	cxy->Modified();
	cxy->SaveAs("./plot/"+str_method+"_xy.pdf");

	vector<double> t_opt[nMod];
	vector< vector<double> > x_opt[nMod];
	TGraph *gxy_opt[nMod];
	TMultiGraph *mgxy_opt = new TMultiGraph();
	TCanvas *cxy_opt = new TCanvas(str_method+"xy_opt", str_method+"xy_opt");
	TLegend *lxy_opt = new TLegend(0.64, 0.6, 0.9, 0.9);
	cout << endl;
	cout << "to maximize the range" << endl;
	for(int iMod=0; iMod<nMod; iMod++) {
		cout << "----- mode " << iMod << " -----" << endl;
		proj.set_mode(iMod);
		double theta_opt = proj.search_theta_for_max_range(700);
		cout << "theta: " << theta_opt << ", max range: " << proj.get_range() << endl;

		t_opt[iMod] = proj.get_t();
		x_opt[iMod] = proj.get_x();
		int npts = proj.get_n_stps();

		gxy_opt[iMod] = new TGraph(npts, &x_opt[iMod][0][0], &x_opt[iMod][1][0]);
		gxy_opt[iMod]->SetLineColor(colorMod[iMod]);
		gxy_opt[iMod]->SetMarkerColor(colorMod[iMod]);
		//lxy_opt->AddEntry(gxy_opt[iMod], strMod[iMod]+Form(" (range: %d m)", (int)proj.get_range()), "lp");
		TString str_tmp = Form("theta: %.1f^{#circ}, range: %.1f m", theta_opt, proj.get_range());
		lxy_opt->AddEntry(gxy_opt[iMod], "#splitline{"+strMod[iMod]+"}{"+str_tmp+"}", "lp");
		mgxy_opt->Add(gxy_opt[iMod]);

	}
	mgxy_opt->GetXaxis()->SetTitle("x (m)");
	mgxy_opt->GetYaxis()->SetTitle("y (m)");
	mgxy_opt->Draw("PLA");
	mgxy_opt->GetXaxis()->SetLimits(0,51000);
	mgxy_opt->GetYaxis()->SetRangeUser(0,17500);
	lxy_opt->Draw("same");
	zeroline->Draw("same");
	l_para->Draw();
	cxy_opt->Modified();
	cxy_opt->SaveAs("./plot/"+str_method+"_xy_opt.pdf");

	cout << endl;

}


int main() {

	gStyle->SetTitleOffset(1.4,"y");
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);

	map<int, TString> str_methods;
	str_methods[1] = "Euler";
	str_methods[2] = "RK2";
	str_methods[4] = "RK4";

	cout << str_methods[1] << endl;
	projectile_plot_rk(1, str_methods[1]);
	cout << str_methods[4] << endl;
	projectile_plot_rk(4, str_methods[4]);

	TLegend *l_drag = new TLegend(0.7,0.7,0.9,0.9);
	TF1 *f_no_drag = new TF1("f_no_drag", "0", 0,15000);
	f_no_drag->SetLineColor(kBlack);
	l_drag->AddEntry(f_no_drag, "no drag", "l");
	TF1 *f_constant_drag = new TF1("f_constant_drag", "1", 0,15000);
	f_constant_drag->SetLineColor(kBlue);
	l_drag->AddEntry(f_constant_drag, "constant air density", "l");
	TF1 *f_isothermal_drag = new TF1("f_isothermal_drag", "exp(-x/10000)", 0,15000);
	f_isothermal_drag->SetLineColor(kRed);
	l_drag->AddEntry(f_isothermal_drag, "isothermal model", "l");
	TF1 *f_adiabatic_drag = new TF1("f_adiabatic_drag", "(1-6.5e-3*x/273)^2.5", 0,15000);
	f_adiabatic_drag->SetLineColor(kGreen);
	l_drag->AddEntry(f_adiabatic_drag, "aiabatic model", "l");

	TCanvas *c_drag = new TCanvas("drag", "drag");
	TH2F *frame = new TH2F("frame", "frame", 1,0,15000, 1,-0.1,1.5);
	frame->GetXaxis()->SetTitle("y (m)");
	frame->GetYaxis()->SetTitle("F_{drag}/(B_{2}(#rho_{0}) v^{2})");
	frame->Draw();
	f_no_drag->Draw("same");
	f_constant_drag->Draw("same");
	f_isothermal_drag->Draw("same");
	f_adiabatic_drag->Draw("same");
	l_drag->Draw("same");
	c_drag->SaveAs("./plot/drag.pdf");

	return 0;
}
