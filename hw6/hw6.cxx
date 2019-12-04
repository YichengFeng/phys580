#include <iostream>
#include <vector>
#include <cmath>
#include "square_well_1d.h"
#include "match_square_well_1d.h"
#include "qm_pow4.h"
#include "lennard_jones_1d.h"

#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TMultiGraph.h>
#include <TH2.h>
#include <TF1.h>

using namespace std;


//-------------------------------------------------------------------------//

void square_well_1d_plot(int parity, double Vmax, double boundary, double a, double V2, double E) {

	double dx = 0.001;
	//double Vmax = 1e5;
	double b = 2;
	double dE = 0.1;
	//double boundary = 4.0;

	TString str_par = parity==1?"even":"odd";
	TString str_adj = str_par + Form("_Vmax%d_a%.3d_Vm%.3d_E%.4d", int(Vmax), int(a*100), int(V2), int(E*100));
	TString str_tmp = "psi_" + str_adj;

	SquareWell1D sw1d(parity, dx, Vmax, a, V2, b, E, dE);
	sw1d.set_boundary(boundary);
	sw1d.cal_once();
	if(a != 0) sw1d.set_iter(1000);
	sw1d.adjust_until(1e-10);

	sw1d.cal_potential_table();

	vector<double> x = sw1d.get_x();
	vector<double> psi = sw1d.get_psi();
	vector<double> potential = sw1d.get_potential_table();

	TCanvas *c = new TCanvas("c_"+str_tmp, "c_"+str_tmp);
	TGraph *g_psi = new TGraph(x.size(), &x[0], &psi[0]);
	TGraph *g_potential = new TGraph(x.size(), &x[0], &potential[0]);
	TF1 *func_tail = new TF1("func_tail", "[0]*exp(-[1]*x)", 1,3);
	g_psi->SetLineColor(kBlue);
	g_potential->SetLineColor(kRed);
	func_tail->SetLineColor(kGreen);
	if(sw1d.get_E() < Vmax) g_psi->Fit(func_tail, "R");
	double ymin = min(TMath::MinElement(g_psi->GetN(), g_psi->GetY()), -b);
	double ymax = max(TMath::MaxElement(g_psi->GetN(), g_psi->GetY()),  b) + 1;
	//TH2D *h2frame = new TH2D("frame"+str_tmp, "frame"+str_tmp, 1,-1.1*boundary,1.1*boundary, 1,1.2*ymin,1.2*ymax);
	TH2D *h2frame = new TH2D("frame"+str_tmp, "frame"+str_tmp, 1,-3.3,3.3, 1,1.2*ymin,1.2*ymax);
	h2frame->Draw();
	g_psi->Draw("L same");
	g_potential->Draw("L same");
	TLegend *l_g = new TLegend(0.3,0.8,0.7,0.9);
	l_g->SetNColumns(2);
	l_g->SetFillStyle(0);
	l_g->SetBorderSize(0);
	l_g->AddEntry(g_psi, "#psi", "l");
	l_g->AddEntry(g_potential, "V", "l");
	l_g->Draw("same");
	TLegend *l_para = new TLegend(0.05,0.9,0.5,0.96);
	l_para->SetFillStyle(0);
	l_para->SetBorderSize(0);
	l_para->AddEntry((TObject*)0, Form("E = %.4f, #DeltaE = %.2e, V_{0} = %.0f", sw1d.get_E(), sw1d.get_dE(), Vmax), "");
	l_para->Draw("same");
	c->SaveAs("./plot/"+str_tmp+".pdf");
}

//-------------------------------------------------------------------------//

void qm_pow4_plot(void) {

	double dx = 0.01;
	double dE = 0.1;
	double E  = 0;

	QMPow4 qmp4(dx, E, dE);
	qmp4.set_T(0);
	qmp4.set_fr(0.1);
	qmp4.cal_var_attempts(10000);
	qmp4.cal_potential_table();

	cout << qmp4.get_varE() << endl;

	vector<double> x = qmp4.get_x();
	vector<double> potential = qmp4.get_potential_table();
	vector<double> psi = qmp4.get_psi();

	TGraph *g_psi = new TGraph(x.size(), &x[0], &psi[0]);
	TGraph *g_potential = new TGraph(x.size(), &x[0], &potential[0]);
	g_psi->SetLineColor(kBlue);
	g_potential->SetLineColor(kRed);

	TCanvas *c_test = new TCanvas("qmpow4", "qmpow4");
	g_psi->Draw("LA");
	g_potential->Draw("L same");
	TLegend *l_test = new TLegend(0.1,0.6,0.3,0.9);
	l_test->SetFillStyle(0);
	l_test->SetBorderSize(0);
	l_test->AddEntry(g_psi, "#psi", "l");
	l_test->AddEntry(g_potential, "V", "l");
	l_test->AddEntry((TObject*)0, Form("E = %.4f", qmp4.get_varE()), "");
	l_test->Draw("same");
	c_test->SaveAs("./plot/qmpow4.pdf");
}

//-------------------------------------------------------------------------//

void match_square_well_1d_plot(double E, int Nmax) {

	//MatchSquareWell1D msw1d(1e5, 0.1, 100, 5e-4, M_PI*M_PI/2, 0.1);
	MatchSquareWell1D msw1d(1e5, 0.1, 100, 5e-4, E, 0.1);
	msw1d.set_Nmax(Nmax);
	msw1d.adjust_match_until(1e-4);
	msw1d.cal_potential_table();

	vector<double> x = msw1d.get_x();
	vector<double> xL = msw1d.get_xL();
	vector<double> xR = msw1d.get_xR();
	vector<double> potential = msw1d.get_potential_table();
	vector<double> psiL = msw1d.get_psiL();
	vector<double> psiR = msw1d.get_psiR();

	TGraph *g_psiL = new TGraph(xL.size(), &xL[0], &psiL[0]);
	g_psiL->SetLineColor(kBlue);
	TGraph *g_psiR = new TGraph(xR.size(), &xR[0], &psiR[0]);
	g_psiR->SetLineColor(kGreen);
	TGraph *g_potential = new TGraph(potential.size(), &x[0], &potential[0]);
	g_potential->SetLineColor(kRed);

	TString str_adj = Form("E%.4d_Nmax%.3d", int(100*E), abs(Nmax));
	TString str_tmp = "msw1d_" + str_adj;

	TCanvas *c_test = new TCanvas("c_"+str_tmp, "c_"+str_tmp);
	TH2D *h2frame = new TH2D("frame_"+str_tmp, "frame_"+str_tmp, 1,-1.3,1.3, 1,-3,3);
	h2frame->Draw();
	g_psiL->Draw("L same");
	g_psiR->Draw("L same");
	g_potential->Draw("L same");
	TLegend *l_psi = new TLegend(0.25,0.8,0.75,0.9);
	l_psi->SetFillStyle(0);
	l_psi->SetBorderSize(0);
	l_psi->SetNColumns(3);
	l_psi->AddEntry(g_psiL, "#psi_{L}", "l");
	l_psi->AddEntry(g_psiR, "#psi_{R}", "l");
	l_psi->AddEntry(g_potential, "V", "l");
	TLegend *l_para = new TLegend(0.1,0.9,0.9,0.96);
	l_para->SetFillStyle(0);
	l_para->SetBorderSize(0);
	l_para->AddEntry((TObject*)0, Form("E = %.4f", msw1d.get_E()), "");
	l_psi->Draw("same");
	l_para->Draw("same");
	c_test->SaveAs("./plot/"+str_tmp+".pdf");
}

//-------------------------------------------------------------------------//

void lennard_jones_1d_plot() {

	LennardJones1D lj1d(10, 1, 0.01, -2, 0.1);
	lj1d.set_T(0);
	lj1d.cal_var_attempts(4000);
	lj1d.cal_potential_table();
	lj1d.var_psi_normalize();

	vector<double> x = lj1d.get_x();
	vector<double> psi = lj1d.get_psi();
	vector<double> potential = lj1d.get_potential_table();

	TGraph *g_psi = new TGraph(x.size(), &x[0], &psi[0]);
	TGraph *g_potential = new TGraph(x.size(), &x[0], &potential[0]);
	g_potential->SetLineColor(kRed);

	LennardJones1D lj1dmatch(10, 1, 0.01, -2, 0.1);
	lj1dmatch.adjust_match_until(1e-4);
	vector<double> xL = lj1dmatch.get_xL();
	vector<double> xR = lj1dmatch.get_xR();
	vector<double> psiL = lj1dmatch.get_psiL();
	vector<double> psiR = lj1dmatch.get_psiR();
	TGraph *g_psiL = new TGraph(xL.size(), &xL[0], &psiL[0]);
	TGraph *g_psiR = new TGraph(xR.size(), &xR[0], &psiR[0]);
	g_psiL->SetLineColor(kBlue);
	g_psiR->SetLineColor(kGreen);

	TCanvas *c = new TCanvas("lj1d", "lj1d");
	TH2D *h2frame = new TH2D("frame_lj1d", "frame_lj1d", 1,0,5.5, 1,-1.5,1.5);
	h2frame->Draw();
	g_psi->Draw("L same");
	g_potential->Draw("L same");
	g_psiL->Draw("L same");
	g_psiR->Draw("L same");

	TLegend *l_var = new TLegend(0.5,0.6,0.7,0.9);
	l_var->SetFillStyle(0);
	//l_var->SetBorderSize(0);
	l_var->SetHeader("#splitline{deterministic}{variational}");
	l_var->AddEntry(g_psi, "#psi", "l");
	l_var->AddEntry(g_potential, "V", "l");
	l_var->AddEntry((TObject*)0, Form("E = %.4f", lj1d.get_varE()), "");
	l_var->AddEntry((TObject*)0, Form("#alpha = %.4f", lj1d.get_varA()), "");
	l_var->Draw("same");
	TLegend *l_match = new TLegend(0.7,0.6,0.9,0.9);
	l_match->SetFillStyle(0);
	//l_match->SetBorderSize(0);
	l_match->SetHeader("#splitline{match}{method}");
	l_match->AddEntry(g_psiL, "#psi_{L}", "l");
	l_match->AddEntry(g_psiR, "#psi_{R}", "l");
	l_match->AddEntry((TObject*)0, Form("E = %.4f", lj1dmatch.get_E()), "");
	l_match->Draw("same");
	c->SaveAs("./plot/lj1d_var.pdf");
}

//-------------------------------------------------------------------------//

int main() {

	gStyle->SetTitleOffset(1.4,"y");
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);
	gStyle->SetOptFit();

	square_well_1d_plot( 1, 10, 3.0, 0, 0, 1.1);
	square_well_1d_plot(-1, 10, 3.0, 0, 0, 4.9);
	square_well_1d_plot( 1, 10, 3.0, 0, 0, 11.);
	square_well_1d_plot(-1, 10, 3.0, 0, 0, 20.);
	square_well_1d_plot( 1, 50, 3.0, 0, 0, 1.1);
	square_well_1d_plot(-1, 50, 3.0, 0, 0, 4.9);
	square_well_1d_plot( 1, 50, 3.0, 0, 0, 11.);
	square_well_1d_plot(-1, 50, 3.0, 0, 0, 15);
	square_well_1d_plot( 1, 50, 3.0, 0, 0, 15);
	square_well_1d_plot(-1, 50, 3.0, 0, 0, 20.);
	square_well_1d_plot( 1, 50, 3.0, 0, 0, 30.);
	square_well_1d_plot(-1, 50, 3.0, 0, 0, 44.);
	square_well_1d_plot( 1, 50, 3.0, 0, 0, 60.);
	square_well_1d_plot( 1, 100, 2.0, 0, 0, 1.1);
	square_well_1d_plot(-1, 100, 2.0, 0, 0, 4.9);
	square_well_1d_plot( 1, 100, 2.0, 0, 0, 11.);
	square_well_1d_plot(-1, 100, 2.0, 0, 0, 20.);
	square_well_1d_plot( 1, 100, 2.0, 0, 0, 30.);
	square_well_1d_plot(-1, 100, 2.0, 0, 0, 44.);
	square_well_1d_plot( 1, 100, 2.0, 0, 0, 60.);
	square_well_1d_plot(-1, 100, 2.0, 0, 0, 60.);
	square_well_1d_plot(-1, 100, 2.0, 0, 0, 80.);

	qm_pow4_plot();
	match_square_well_1d_plot(M_PI*M_PI/2, 1);
	match_square_well_1d_plot(M_PI*M_PI, -999);
	match_square_well_1d_plot(M_PI*M_PI/2, -999);
	match_square_well_1d_plot(M_PI*M_PI*2, 1);
	match_square_well_1d_plot(M_PI*M_PI*2, -999);
	match_square_well_1d_plot(M_PI*M_PI*3, -999);

	lennard_jones_1d_plot();

	return 0;
}

//-------------------------------------------------------------------------//

