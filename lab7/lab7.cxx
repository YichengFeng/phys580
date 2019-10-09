#include <iostream>
#include <vector>
#include <cmath>
#include "capacitor_2d.h"

#include <TH2.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>

using namespace std;


//-------------------------------------------------------------------------//

void capacitor_2d_plot(int alg, double d, double acc) {

	TString str_alg[2] = {"Jacobi", "GaussSeidel"};

	//double d = 0.05;
	//double acc = 1e-6;

	TString str_d = Form("d%.4d", int(d*1000));
	TString str_acc = Form("acc%.1d", int(fabs(log10(acc))));

	Capacitor2D c2d(0.5, 0.5, d, acc);
	c2d.set_alg(alg);
	c2d.cal();
	vector< vector<double> > V = c2d.get_V();
	int nx = V[0].size();
	int ny = V.size();
	double tmp_acc = c2d.get_tmp_acc();
	int n_iter = c2d.get_n_iter();

	cout << tmp_acc << " " << nx << " " << ny << " " << n_iter << endl;

	TH2F *h2_V = new TH2F("capacitor2d", "capacitor2d", nx,-0.5,0.5, ny,-0.5,0.5);
	h2_V->GetXaxis()->SetTitle("x/L");
	h2_V->GetYaxis()->SetTitle("y/L");
	TH2F *h2_minus_V = new TH2F("capacitor2d_minus", "capacitor2d_minus", nx,-0.5,0.5, ny,-0.5,0.5);
	h2_minus_V->GetXaxis()->SetTitle("x/L");
	h2_minus_V->GetYaxis()->SetTitle("y/L");

	for(int j=0; j<ny; j++) {
		for(int i=0; i<nx; i++) {
			h2_V->SetBinContent(i+1, j+1, V[j][i]);
			h2_minus_V->SetBinContent(i+1, j+1, -V[j][i]);
		}
	}

	TString str_adj = str_alg[alg] + "_" + str_d + "_" + str_acc;
	TString str_tmp;

	TLegend *l_iter = new TLegend(0.1,0.9,0.7,0.98);
	l_iter->SetFillStyle(0);
	l_iter->SetBorderSize(0);
	l_iter->AddEntry((TObject*)0, Form("number of iteration: %d",n_iter), "");
	l_iter->AddEntry((TObject*)0, Form("d = %.3fL, accuracy: %.2e < %.2e",d,acc,tmp_acc), "");

	str_tmp = "V_mesh_" + str_adj;
	TCanvas *c_V_mesh = new TCanvas("c_"+str_tmp, "c_"+str_tmp, 600,600);
	h2_V->Draw("SURF2");
	l_iter->Draw("same");
	c_V_mesh->SaveAs("./plot/"+str_tmp+".pdf");

	str_tmp = "V_contour_" + str_adj;
	TCanvas *c_V_contour = new TCanvas("c_"+str_tmp, "c_"+str_tmp, 600,600);
	h2_V->Draw("CONT3");
	l_iter->Draw("same");
	c_V_contour->SaveAs("./plot/"+str_tmp+".pdf");

	str_tmp = "E_field_" + str_adj;
	TCanvas *c_E_field = new TCanvas("c_"+str_tmp, "c_"+str_tmp, 600,600);
	//h2_minus_V->Draw("ARR");
	h2_V->Draw("ARR");
	l_iter->Draw("same");
	c_E_field->SaveAs("./plot/"+str_tmp+".pdf");

	delete h2_V;
	delete h2_minus_V;
}

//-------------------------------------------------------------------------//

void capacitor_2d_rho_plot(int alg) {

	TString str_alg[2] = {"Jacobi", "GaussSeidel"};

	double d = 0.01;
	double acc = 1e-6;

	Capacitor2D c2d(0.5, 0.5, d, acc);
	c2d.set_alg(alg);
	c2d.set_rho(-0.05,0.05,-0.05,0.05,1.0);
	c2d.cal();
	vector< vector<double> > V = c2d.get_V();
	int nx = V[0].size();
	int ny = V.size();
	double tmp_acc = c2d.get_tmp_acc();
	int n_iter = c2d.get_n_iter();

	cout << tmp_acc << " " << nx << " " << ny << " " << n_iter << endl;

	TH2F *h2_V = new TH2F("capacitor2d", "capacitor2d", nx,-0.5,0.5, ny,-0.5,0.5);
	h2_V->GetXaxis()->SetTitle("x/L");
	h2_V->GetYaxis()->SetTitle("y/L");
	TH2F *h2_minus_V = new TH2F("capacitor2d_minus", "capacitor2d_minus", nx,-0.5,0.5, ny,-0.5,0.5);
	h2_minus_V->GetXaxis()->SetTitle("x/L");
	h2_minus_V->GetYaxis()->SetTitle("y/L");

	for(int j=0; j<ny; j++) {
		for(int i=0; i<nx; i++) {
			h2_V->SetBinContent(i+1, j+1, V[j][i]);
			h2_minus_V->SetBinContent(i+1, j+1, -V[j][i]);
		}
	}

	TString str_adj = str_alg[alg];
	TString str_tmp;

	TLegend *l_iter = new TLegend(0.1,0.9,0.7,0.98);
	l_iter->SetFillStyle(0);
	l_iter->SetBorderSize(0);
	l_iter->AddEntry((TObject*)0, Form("number of iteration: %d",n_iter), "");
	l_iter->AddEntry((TObject*)0, Form("d = %.3fL, accuracy: %.2e < %.2e",d,tmp_acc,acc), "");

	str_tmp = "rho_V_mesh_" + str_adj;
	TCanvas *c_V_mesh = new TCanvas("c_"+str_tmp, "c_"+str_tmp, 600,600);
	h2_V->Draw("SURF2");
	l_iter->Draw("same");
	c_V_mesh->SaveAs("./plot/"+str_tmp+".pdf");

	str_tmp = "rho_V_contour_" + str_adj;
	TCanvas *c_V_contour = new TCanvas("c_"+str_tmp, "c_"+str_tmp, 600,600);
	h2_V->Draw("CONT3");
	l_iter->Draw("same");
	c_V_contour->SaveAs("./plot/"+str_tmp+".pdf");

	str_tmp = "rho_E_field_" + str_adj;
	TCanvas *c_E_field = new TCanvas("c_"+str_tmp, "c_"+str_tmp, 600,600);
	//h2_minus_V->Draw("ARR");
	//h2_V->RebinX(2);
	//h2_V->RebinY(2);
	h2_V->Draw("ARR");
	l_iter->Draw("same");
	c_E_field->SaveAs("./plot/"+str_tmp+".pdf");

	delete h2_V;
	delete h2_minus_V;
}

//-------------------------------------------------------------------------//

void capacitor_2d_rho_periodic_plot(int alg) {

	TString str_alg[2] = {"Jacobi", "GaussSeidel"};

	double d = 0.01;
	double acc = 1e-6;

	Capacitor2D c2d(1.0-2*d, 0.5, d, acc);
	c2d.set_alg(alg);
	c2d.set_periodic(1);
	c2d.set_rho(-0.05,0.05,-0.05,0.05,1.0);
	c2d.cal();
	vector< vector<double> > V = c2d.get_V();
	int nx = V[0].size();
	int ny = V.size();
	double tmp_acc = c2d.get_tmp_acc();
	int n_iter = c2d.get_n_iter();

	cout << tmp_acc << " " << nx << " " << ny << " " << n_iter << endl;

	TH2F *h2_V = new TH2F("capacitor2d", "capacitor2d", nx,-0.5,0.5, ny,-0.5,0.5);
	h2_V->GetXaxis()->SetTitle("x/L");
	h2_V->GetYaxis()->SetTitle("y/L");
	TH2F *h2_minus_V = new TH2F("capacitor2d_minus", "capacitor2d_minus", nx,-0.5,0.5, ny,-0.5,0.5);
	h2_minus_V->GetXaxis()->SetTitle("x/L");
	h2_minus_V->GetYaxis()->SetTitle("y/L");

	for(int j=0; j<ny; j++) {
		for(int i=0; i<nx; i++) {
			h2_V->SetBinContent(i+1, j+1, V[j][i]);
			h2_minus_V->SetBinContent(i+1, j+1, -V[j][i]);
		}
	}

	TString str_adj = str_alg[alg];
	TString str_tmp;

	TLegend *l_iter = new TLegend(0.1,0.9,0.7,0.98);
	l_iter->SetFillStyle(0);
	l_iter->SetBorderSize(0);
	l_iter->AddEntry((TObject*)0, Form("number of iteration: %d",n_iter), "");
	l_iter->AddEntry((TObject*)0, Form("d = %.3fL, accuracy: %.2e < %.2e",d,tmp_acc,acc), "");

	str_tmp = "rho_periodic_V_mesh_" + str_adj;
	TCanvas *c_V_mesh = new TCanvas("c_"+str_tmp, "c_"+str_tmp, 600,600);
	h2_V->Draw("SURF2");
	l_iter->Draw("same");
	c_V_mesh->SaveAs("./plot/"+str_tmp+".pdf");

	str_tmp = "rho_periodic_V_contour_" + str_adj;
	TCanvas *c_V_contour = new TCanvas("c_"+str_tmp, "c_"+str_tmp, 600,600);
	h2_V->Draw("CONT3");
	l_iter->Draw("same");
	c_V_contour->SaveAs("./plot/"+str_tmp+".pdf");

	str_tmp = "rho_periodic_E_field_" + str_adj;
	TCanvas *c_E_field = new TCanvas("c_"+str_tmp, "c_"+str_tmp, 600,600);
	//h2_minus_V->Draw("ARR");
	//h2_V->RebinX(2);
	//h2_V->RebinY(2);
	h2_V->Draw("ARR");
	l_iter->Draw("same");
	c_E_field->SaveAs("./plot/"+str_tmp+".pdf");

	delete h2_V;
	delete h2_minus_V;
}

//-------------------------------------------------------------------------//

int main() {

	gStyle->SetTitleOffset(1.4,"y");
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);

	// (1)
	capacitor_2d_plot(0, 0.05, 1e-4);
	capacitor_2d_plot(0, 0.05, 1e-5);
	capacitor_2d_plot(0, 0.05, 1e-6);
	capacitor_2d_plot(0, 0.05, 1e-7);
	capacitor_2d_plot(0, 0.05, 1e-8);
	capacitor_2d_plot(0, 0.05, 1e-9);
	
	capacitor_2d_plot(0,0.005, 1e-6);
	capacitor_2d_plot(0, 0.01, 1e-6);
	capacitor_2d_plot(0, 0.02, 1e-6);
	capacitor_2d_plot(0, 0.04, 1e-6);

	// (2)
	capacitor_2d_plot(1, 0.05, 1e-4);
	capacitor_2d_plot(1, 0.05, 1e-5);
	capacitor_2d_plot(1, 0.05, 1e-6);
	capacitor_2d_plot(1, 0.02, 1e-6);
	capacitor_2d_plot(1, 0.01, 1e-6);

	// (3)
	capacitor_2d_rho_plot(0);
	capacitor_2d_rho_plot(1);

	capacitor_2d_rho_periodic_plot(0);

}
