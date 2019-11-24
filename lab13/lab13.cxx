#include <iostream>
#include <vector>
#include <cmath>
#include "square_well_1d.h"
#include "lennard_jones_1d.h"

#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TMultiGraph.h>
#include <TH2.h>

using namespace std;


//-------------------------------------------------------------------------//

void square_well_1d_plot(int parity, double a, double V2, double E) {

	double dx = 0.001;
	double Vmax = 1e5;
	double b = 2;
	double dE = 0.1;

	TString str_par = parity==1?"even":"odd";
	TString str_adj = str_par + Form("_a%.3d_Vm%.3d_E%.4d", int(a*100), int(V2), int(E*100));
	TString str_tmp = "psi_" + str_adj;

	SquareWell1D sw1d(parity, dx, Vmax, a, V2, b, E, dE);
	sw1d.cal_once();
	if(a != 0) sw1d.set_iter(1000);
	sw1d.adjust_until(1e-10);

	sw1d.cal_potential_table();

	vector<double> x = sw1d.get_x();
	vector<double> psi = sw1d.get_psi();
	vector<double> potential = sw1d.get_potential_table();

	TCanvas *c = new TCanvas("c_"+str_tmp, "c_"+str_tmp);
	TH2D *h2frame = new TH2D("frame"+str_tmp, "frame"+str_tmp, 1,-1.2,1.2, 1,-1.2*b,1.2*b);
	h2frame->Draw();
	TGraph *g_psi = new TGraph(x.size(), &x[0], &psi[0]);
	TGraph *g_potential = new TGraph(x.size(), &x[0], &potential[0]);
	g_psi->SetLineColor(kBlue);
	g_potential->SetLineColor(kRed);
	g_psi->Draw("L same");
	g_potential->Draw("L same");
	TLegend *l_g = new TLegend(0.3,0.8,0.7,0.9);
	l_g->SetNColumns(2);
	l_g->SetFillStyle(0);
	l_g->SetBorderSize(0);
	l_g->AddEntry(g_psi, "#psi", "l");
	l_g->AddEntry(g_potential, "V", "l");
	l_g->Draw("same");
	TLegend *l_para = new TLegend(0.1,0.9,0.5,0.96);
	l_para->SetFillStyle(0);
	l_para->SetBorderSize(0);
	l_para->AddEntry((TObject*)0, Form("E = %.4f, #DeltaE = %.2e", sw1d.get_E(), sw1d.get_dE()), "");
	l_para->Draw("same");
	c->SaveAs("./plot/"+str_tmp+".pdf");
}

//-------------------------------------------------------------------------//

void lennard_jones_1d_match_plot(double E, double epsilon = 10) {

	double sigma = 1;
	double dx = 0.01;
	double dE = 0.1;

	TString str_adj = Form("lj_E%.3d", int(fabs(100*E)));
	TString str_tmp;
	if(epsilon > 10) str_tmp = "deep_";
	str_tmp += str_adj;

	LennardJones1D lj1d(epsilon, sigma, dx, E, dE);
	lj1d.adjust_match_until(1e-4);
	lj1d.cal_potential_table();

	vector<double> x = lj1d.get_x();
	vector<double> potential = lj1d.get_potential_table();
	vector<double> xL = lj1d.get_xL();
	vector<double> xR = lj1d.get_xR();
	vector<double> psiL = lj1d.get_psiL();
	vector<double> psiR = lj1d.get_psiR();

	TGraph *g_psiL = new TGraph(xL.size(), &xL[0], &psiL[0]);
	TGraph *g_psiR = new TGraph(xR.size(), &xR[0], &psiR[0]);
	TGraph *g_potential = new TGraph(x.size()-10, &x[10], &potential[10]);

	TCanvas *c = new TCanvas("c_"+str_tmp, "c_"+str_tmp);
	double Lmax = TMath::MaxElement(g_psiL->GetN(), g_psiL->GetY());
	double Lmin = TMath::MinElement(g_psiL->GetN(), g_psiL->GetY());
	double Rmax = TMath::MaxElement(g_psiR->GetN(), g_psiR->GetY());
	double Rmin = TMath::MinElement(g_psiR->GetN(), g_psiR->GetY());
	double ymax = max(max(Lmax, Rmax), 1.5);
	double ymin = min(min(Lmin, Rmin),-1.5);
	TH2D *h2frame = new TH2D("frame"+str_tmp, "frame"+str_tmp, 1,0,5.5, 1,ymin,ymax);
	h2frame->Draw();
	g_psiL->SetLineColor(kBlue);
	g_psiR->SetLineColor(kGreen);
	g_potential->SetLineColor(kRed);
	g_psiL->Draw("L same");
	g_psiR->Draw("L same");
	g_potential->Draw("L same");
	TLegend *l_g = new TLegend(0.7,0.7,0.9,0.9);
	l_g->SetFillStyle(0);
	l_g->SetBorderSize(0);
	l_g->AddEntry(g_psiL, "#psi_{L}", "l");
	l_g->AddEntry(g_psiR, "#psi_{R}", "l");
	l_g->AddEntry(g_potential, "V", "l");
	l_g->Draw("same");
	TLegend *l_para = new TLegend(0.1,0.9,0.5,0.96);
	l_para->SetFillStyle(0);
	l_para->SetBorderSize(0);
	l_para->AddEntry((TObject*)0, Form("E = %.4f, #DeltaE = %.2e", lj1d.get_E(), lj1d.get_dE()), "");
	l_para->Draw("same");
	c->SaveAs("./plot/"+str_tmp+".pdf");
}	

//-------------------------------------------------------------------------//

void lennard_jones_1d_var_plot(double T) {

	double epsilon = 10;
	double sigma = 1;
	double dx = 0.01;
	double dE = 0.1;
	double E = -2;

	TString str_adj = Form("var_lj_T%.4d", int(fabs(1000*T)));
	TString str_tmp = str_adj;

	LennardJones1D lj1d(epsilon, sigma, dx, E, dE);
	lj1d.set_T(T);
	lj1d.set_fr(0.1);
	lj1d.cal_var_attempts(10000);
	lj1d.set_fr(0.01);
	lj1d.cal_var_attempts(5000);
	lj1d.cal_potential_table();

	cout << lj1d.get_varE() << endl;

	vector<double> x = lj1d.get_x();
	vector<double> potential = lj1d.get_potential_table();
	vector<double> psi = lj1d.get_psi();

	TGraph *g_psi = new TGraph(x.size(), &x[0], &psi[0]);
	TGraph *g_potential = new TGraph(x.size()-10, &x[10], &potential[10]);

	TCanvas *c = new TCanvas("c_"+str_tmp, "c_"+str_tmp);
	double Fmax = TMath::MaxElement(g_psi->GetN(), g_psi->GetY());
	double Fmin = TMath::MinElement(g_psi->GetN(), g_psi->GetY());
	double ymax = max(Fmax, 1.5);
	double ymin = min(Fmin,-1.5);
	TH2D *h2frame = new TH2D("frame"+str_tmp, "frame"+str_tmp, 1,0,5.5, 1,ymin,ymax);
	h2frame->Draw();
	g_psi->SetLineColor(kBlue);
	g_potential->SetLineColor(kRed);
	g_psi->Draw("L same");
	g_potential->Draw("L same");
	TLegend *l_g = new TLegend(0.7,0.7,0.9,0.9);
	l_g->SetFillStyle(0);
	l_g->SetBorderSize(0);
	l_g->AddEntry(g_psi, "#psi", "l");
	l_g->AddEntry(g_potential, "V", "l");
	l_g->Draw("same");
	TLegend *l_para = new TLegend(0.1,0.9,0.5,0.96);
	l_para->SetFillStyle(0);
	l_para->SetBorderSize(0);
	l_para->AddEntry((TObject*)0, Form("E = %.4f", lj1d.get_varE()), "");
	l_para->Draw("same");
	c->SaveAs("./plot/"+str_tmp+".pdf");
}	

//-------------------------------------------------------------------------//

int main() {

	gStyle->SetTitleOffset(1.4,"y");
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);
	gStyle->SetOptFit();

	square_well_1d_plot( 1, 0, 0, 1.1);
	square_well_1d_plot(-1, 0, 0, 4.9);
	square_well_1d_plot( 1, 0, 0, 11.);
	square_well_1d_plot(-1, 0, 0, 20.);
	square_well_1d_plot( 1, 0.1, 100, 1.1);
	square_well_1d_plot(-1, 0.1, 100, 4.9);
	square_well_1d_plot( 1, 0.1, 100, 11.);
	square_well_1d_plot(-1, 0.1, 100, 20.);

	lennard_jones_1d_match_plot(-2);
	lennard_jones_1d_match_plot(0.3);
	lennard_jones_1d_match_plot(1.6);

	lennard_jones_1d_match_plot(-60, 30);
	lennard_jones_1d_match_plot(-10, 30);
	lennard_jones_1d_match_plot(-2, 30);
	lennard_jones_1d_match_plot(12, 30);

	lennard_jones_1d_var_plot(0);

	return 0;
}

//-------------------------------------------------------------------------//

