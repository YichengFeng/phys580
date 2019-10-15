#include <iostream>
#include <vector>
#include <cmath>
#include "vector_field_2d.h"
#include "loop_magnetic_field.h"

#include <TCanvas.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TF1.h>
#include <TLegend.h>

using namespace std;


void discrete_loop(double I1, TString str_type) {

	LoopMagneticField lmf(30, 0.3, 2, 0.8, 0.1);
	lmf.set_I_loop(1,I1);
	lmf.cal();
	vector< vector<double> > bx = lmf.get_bx();
	vector< vector<double> > by = lmf.get_by();
	vector< vector<double> > bz = lmf.get_bz();

	TString str_adj = "discrete_loop_"+str_type;
	TString str_tmp;

	str_tmp = "bxz_" + str_adj;
	TVectorField2D *vf2_bxz = new TVectorField2D("vf2_"+str_tmp, -1,1, -1,1, bx, bz);
	vf2_bxz->GetXaxis()->SetTitle("x");
	vf2_bxz->GetYaxis()->SetTitle("z");
	TCanvas *c_bxz = new TCanvas("c_"+str_tmp, "c_bxz_"+str_tmp, 600,600);
	vf2_bxz->Draw();
	c_bxz->SaveAs("./plot/"+str_tmp+".pdf");

	const int nz = bz.size();
	const int nx = bz[0].size();
	double z00[nz];
	double bz00[nz];
	for(int i=0; i<nz; i++) {
		z00[i]  = -1+2.0/(nz-1)*i;
		bz00[i] = bz[i][(nx-1)/2];
	}
	TGraph *g_bz00 = new TGraph(nz, z00, bz00);
	g_bz00->GetXaxis()->SetTitle("z (at x=0, y=0)");
	g_bz00->GetYaxis()->SetTitle("B_{z}");

	str_tmp = "bz00_" + str_adj;
	TCanvas *c_bz00 = new TCanvas("c_"+str_tmp, "c_"+str_tmp);
	g_bz00->Draw("LA");
	c_bz00->SaveAs("./plot/"+str_tmp+".pdf");
	
}

void helical_loop(int n_loop) {

	TString str_n_loop = Form("loop%d", n_loop);

	LoopMagneticField lmf(300, 0.3, n_loop, 1.0/n_loop, 0.1);
	lmf.set_mod(1);
	lmf.cal();
	vector< vector<double> > bx = lmf.get_bx();
	vector< vector<double> > by = lmf.get_by();
	vector< vector<double> > bz = lmf.get_bz();

	TString str_adj = "helical_loop_"+str_n_loop;
	TString str_tmp;

	str_tmp = "bxz_" + str_adj;
	TVectorField2D *vf2_bxz = new TVectorField2D("vf2_"+str_tmp, -1,1, -1,1, bx, bz);
	vf2_bxz->GetXaxis()->SetTitle("x");
	vf2_bxz->GetYaxis()->SetTitle("z");
	TCanvas *c_bxz = new TCanvas("c_"+str_tmp, "c_"+str_tmp, 600,600);
	vf2_bxz->Draw();
	c_bxz->SaveAs("./plot/"+str_tmp+".pdf");

	const int nz = bz.size();
	const int nx = bz[0].size();
	double z00[nz];
	double bz00[nz];
	for(int i=0; i<nz; i++) {
		z00[i]  = -1+2.0/(nz-1)*i;
		bz00[i] = bz[i][(nx-1)/2];
	}
	TGraph *g_bz00 = new TGraph(nz, z00, bz00);
	g_bz00->GetXaxis()->SetTitle("z (at x=0, y=0)");
	g_bz00->GetYaxis()->SetTitle("B_{z}");

	str_tmp = "bz00_" + str_adj;
	TCanvas *c_bz00 = new TCanvas("c_"+str_tmp, "c_"+str_tmp);
	g_bz00->Draw("LA");
	TF1 *f_bz00 = new TF1("f_bz00"+str_tmp, Form("2*TMath::Pi()*%d*((0.5+x)/sqrt((0.5+x)*(0.5+x)+0.3*0.3) + (0.5-x)/sqrt((0.5-x)*(0.5-x)+0.3*0.3))",n_loop), -1,1);
	f_bz00->SetLineColor(kRed);
	f_bz00->SetLineStyle(7);
	f_bz00->SetLineWidth(3);
	f_bz00->Draw("same");
	g_bz00->Draw("L same");
	TLegend *l_bz00 = new TLegend(0.1,0.7,0.3,0.9);
	l_bz00->AddEntry(g_bz00, "numerical result", "l");
	l_bz00->AddEntry(f_bz00, "analytical result", "l");
	l_bz00->Draw("same");
	c_bz00->SaveAs("./plot/"+str_tmp+".pdf");
}

int main() {

	gStyle->SetTitleOffset(1.4,"y");
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);

	discrete_loop( 1.0, "same");
	discrete_loop(-1.0, "oppo");
	helical_loop(3);
	helical_loop(30);

	return 0;
}
