#include <iostream>
#include <cmath>
#include <vector>
#include <map>
#include "runge_kutta.h"
#include "relativity_precession.h"
#include "two_planet_orbit.h"
#include "capacitor_2d.h"
#include "charge_dist_3d.h"
#include "nintegrate_1d.h"

#include <TGraph.h>
#include <TCanvas.h>
#include <TH2.h>
#include <TStyle.h>
#include <TF1.h>
#include <TLegend.h>

using namespace std;


//-------------------------------------------------------------------------//

double G_MX = 0.1;

//-------------------------------------------------------------------------//

void arc_to_deg(vector<double> &theta) {
	for(int i=0; i<theta.size(); i++) {
		theta[i] = theta[i]/M_PI*180;
	}
}

//-------------------------------------------------------------------------//

double relativity_precession_plot(double alpha) {

	double r_min = 0.307;
	double e = 0.206;
	double a = r_min/(1-e);
	double r_max = a*(1+e);
	double f = a*e;
	double b = sqrt(a*a-f*f);
	double t_end = 2.5;

	RelativityPrecession rp(a, e, 3.30e23/1.99e30);
	rp.set_alpha(alpha);
	rp.set_t_end(t_end);
	rp.set_dt(0.0001);
	rp.set_alg(4);
	rp.cal();
	rp.cal_perihelion();

	vector<double> t = rp.get_t();
	vector< vector<double> > x = rp.get_x();
	int n_stps = rp.get_n_stps();

	vector<double> perihelion_t = rp.get_perihelion_t();
	vector<double> perihelion_x = rp.get_perihelion_x();
	vector<double> perihelion_y = rp.get_perihelion_y();
	vector<double> perihelion_theta = rp.get_perihelion_theta();
	arc_to_deg(perihelion_theta);
	int n_perihelion = perihelion_t.size();

	TGraph *g_x_y = new TGraph(n_stps, &x[0][0], &x[1][0]);
	TGraph *g_perihelion_x_y = new TGraph(n_perihelion, &perihelion_x[0], &perihelion_y[0]);
	TGraph *g_perihelion_t_theta = new TGraph(n_perihelion, &perihelion_t[0], &perihelion_theta[0]);
	double Sun_x[1] = {0};
	double Sun_y[1] = {0};
	TGraph *g_Sun = new TGraph(1, Sun_x, Sun_y);
	g_Sun->SetMarkerStyle(20);
	g_Sun->SetMarkerSize(2);
	g_Sun->SetMarkerColor(kOrange);

	int canvas_nx = 700;
	int canvas_ny = 500;
	double rm = 0.1;
	double rx, ry;
	if(a/canvas_nx < b/canvas_ny) {
		rx = b/canvas_ny*canvas_nx/a;
		ry = 1;
	} else {
		rx = 1;
		ry = a/canvas_nx*canvas_ny/b;
	}
	TH2F *h2frame = new TH2F("frame", "frame", 1,-(r_max+rm*a)*rx,(r_min+rm*a)*rx, 1,-(b+b*rm)*ry,(b+b*rm)*ry);
	h2frame->GetXaxis()->SetTitle("x (AU)");
	h2frame->GetYaxis()->SetTitle("y (AU)");

	TString str_alpha = Form("alpha%.2d", int(10000*alpha));
	TString str_adj = str_alpha;
	TString str_tmp;

	str_tmp = "x_y_" + str_adj;
	TCanvas *c_x_y = new TCanvas("c_"+str_tmp, "c_"+str_tmp);
	h2frame->Draw();
	g_x_y->Draw("L same");
	g_perihelion_x_y->SetMarkerStyle(20);
	g_perihelion_x_y->SetMarkerColor(kBlue);
	g_perihelion_x_y->Draw("P same");
	g_Sun->Draw("P same");
	TLegend *l_obj = new TLegend(0.12,0.7,0.35,0.9);
	l_obj->SetFillStyle(0);
	l_obj->SetBorderSize(0);
	l_obj->AddEntry(g_x_y, "Mercury trajectory", "l");
	l_obj->AddEntry(g_perihelion_x_y, "Mercury perihelions", "p");
	l_obj->AddEntry(g_Sun, "Sun", "p");
	l_obj->Draw();
	c_x_y->SaveAs("./plot/"+str_tmp+".pdf");

	str_tmp = "perihelion_t_theta_" + str_adj;
	TCanvas *c_perihelion_t_theta = new TCanvas("c_"+str_tmp, "c_"+str_tmp);
	g_perihelion_t_theta->SetMarkerStyle(20);
	g_perihelion_t_theta->SetMarkerColor(kBlue);
	g_perihelion_t_theta->GetXaxis()->SetTitle("t (yr)");
	g_perihelion_t_theta->GetYaxis()->SetTitle("#theta (degree)");
	g_perihelion_t_theta->Draw("PA");
	TF1 *func = new TF1("func", "[0]*x", 0, t_end);
	func->SetLineColor(kMagenta);
	g_perihelion_t_theta->Fit(func);
	TLegend *l_fit = new TLegend(0.1,0.7,0.6,0.9);
	l_fit->SetFillStyle(0);
	l_fit->SetBorderSize(0);
	l_fit->AddEntry((TObject*)0, Form("#alpha = %f AU^{2}", alpha), "");
	l_fit->AddEntry((TObject*)0, Form("slope = %f degree/yr", func->GetParameter(0)), "");
	l_fit->Draw("same");
	c_perihelion_t_theta->SaveAs("./plot/"+str_tmp+".pdf");

	return func->GetParameter(0);
}

//-------------------------------------------------------------------------//

int relativity_precession_cal() {

	const int n_alpha = 7;
	double alphas[n_alpha] = {0, 0.0002, 0.0005, 0.0008, 0.0015, 0.0025, 0.0035};
	double dtheta_dt[n_alpha];
	for(int i_alpha=0; i_alpha<n_alpha; i_alpha++){
		dtheta_dt[i_alpha] = relativity_precession_plot(alphas[i_alpha]);
	}
	TGraph *g_alpha_dtheta_dt = new TGraph(n_alpha, alphas, dtheta_dt);
	g_alpha_dtheta_dt->SetMarkerStyle(21);
	g_alpha_dtheta_dt->GetXaxis()->SetTitle("#alpha");
	g_alpha_dtheta_dt->GetYaxis()->SetTitle("d#theta/dt (degree/yr)");
	TCanvas *c_alpha_dtheta_dt = new TCanvas("alpha_dtheta_dt", "alpha_dtheta_dt");
	g_alpha_dtheta_dt->Draw("PA");
	TF1 *func = new TF1("func", "[0]*x", 0, 0.0037);
	g_alpha_dtheta_dt->Fit(func);
	TLegend *l_fit = new TLegend(0.1,0.8,0.45,0.9);
	l_fit->SetFillStyle(0);
	l_fit->SetBorderSize(0);
	l_fit->AddEntry((TObject*)0, Form("slope = %f degree/yr/AU^{2}", func->GetParameter(0)), "");
	l_fit->Draw("same");
	c_alpha_dtheta_dt->SaveAs("./plot/alpha_dtheta_dt.pdf");

	double alpha0 = 1.1e-8;
	cout << endl;
	cout << "Mercury precession angular velocity: " << endl;
	cout << alpha0*func->GetParameter(0)*3600*100 << " arcsecond/century" << endl;
}

//-------------------------------------------------------------------------//

void two_planet_orbit_plot(int mod) {

	TString str_mods[3] = {"real", "binary", "elliptical"};
	double m2s[3] = {3.69397e-8, 3.0027e-6, 3.69397e-8};
	double rdf = 0.02;
	double x1s[3] = {1, sqrt(1-rdf*rdf), 1};
	double x2s[3] = {1+2.570e-3, sqrt(1-rdf*rdf), 1+2.570e-3};
	double y1s[3] = {0, rdf, 0};
	double y2s[3] = {0, -rdf, 0};
	double v1xs[3] = {0, 2*M_PI*sqrt(0.25/rdf*3.0027e-6), 0};
	double v2xs[3] = {0, -2*M_PI*sqrt(0.25/rdf*3.0027e-6), 0};
	double v1ys[3] = {2*M_PI, 2*M_PI, 2*M_PI};
	double v2ys[3] = {2*M_PI*(1+2.570e-3)+2*M_PI/29.53*365.2422*2.570e-3, 2*M_PI, 2*M_PI*(1+2.570e-3)+2*M_PI/29.53*365.2422*2.570e-3*1.2};
	double t_ends[3] = {1.02, 1.2, 1.2};

	TString str_mod = str_mods[mod];
	double m2 = m2s[mod];
	double x1 = x1s[mod];
	double x2 = x2s[mod];
	double y1 = y1s[mod];
	double y2 = y2s[mod];
	double v1x = v1xs[mod];
	double v2x = v2xs[mod];
	double v1y = v1ys[mod];
	double v2y = v2ys[mod];
	double t_end = t_ends[mod];

	TwoPlanetOrbit tpo;
	tpo.set_m2(m2);
	tpo.set_x1_start(x1);
	tpo.set_x2_start(x2);
	tpo.set_y1_start(y1);
	tpo.set_y2_start(y2);
	tpo.set_vx1_start(v1x);
	tpo.set_vx2_start(v2x);
	tpo.set_vy1_start(v1y);
	tpo.set_vy2_start(v2y);
	tpo.set_t_end(t_end);
	tpo.cal();

	vector<double> t = tpo.get_t();
	vector< vector<double> > x = tpo.get_x();
	int n_stps = tpo.get_n_stps();

	vector<double> relative_x;
	vector<double> relative_y;
	vector<double> zoomin_x;
	vector<double> zoomin_y;
	double zoomin = 10;
	for(int i=0; i<n_stps; i++) {
		relative_x.push_back(x[2][i]-x[0][i]);
		relative_y.push_back(x[3][i]-x[1][i]);

		zoomin_x.push_back(x[0][i] + zoomin*(x[2][i]-x[0][i]));
		zoomin_y.push_back(x[1][i] + zoomin*(x[3][i]-x[1][i]));
	}

	double a = 1;
	double b = 1;
	double r_min = 1;
	double r_max = 1;
	int canvas_nx = 700;
	int canvas_ny = 500;
	double rm = 0.1;
	double rx, ry;
	if(a/canvas_nx < b/canvas_ny) {
		rx = b/canvas_ny*canvas_nx/a;
		ry = 1;
	} else {
		rx = 1;
		ry = a/canvas_nx*canvas_ny/b;
	}
	TH2F *h2frame = new TH2F("frame", "frame", 1,-(r_max+rm*a)*rx,(r_min+rm*a)*rx, 1,-(b+b*rm)*ry,(b+b*rm)*ry);
	h2frame->GetXaxis()->SetTitle("x (AU)");
	h2frame->GetYaxis()->SetTitle("y (AU)");

	TGraph *g_Earth_x_y = new TGraph(n_stps, &x[0][0], &x[1][0]);
	TGraph *g_Moon_x_y = new TGraph(n_stps, &x[2][0], &x[3][0]);
	TGraph *g_Moon_zoomin_x_y = new TGraph(n_stps, &zoomin_x[0], &zoomin_y[0]);
	TGraph *g_relative_x_y = new TGraph(n_stps, &relative_x[0], &relative_y[0]);
	double Sun_x[1] = {0};
	double Sun_y[1] = {0};
	TGraph *g_Sun_x_y = new TGraph(1, Sun_x, Sun_y);
	g_Sun_x_y->SetMarkerStyle(20);
	g_Sun_x_y->SetMarkerSize(2);
	g_Sun_x_y->SetMarkerColor(kOrange);

	g_Earth_x_y->SetLineColor(kBlue);
	TCanvas *c_x_y = new TCanvas("x_y", "x_y");
	h2frame->Draw();
	g_Earth_x_y->Draw("L same");
	g_Moon_x_y->Draw("L same");
	g_Sun_x_y->Draw("P same");
	TLegend *l_x_y = new TLegend(0.1,0.7,0.3,0.9);
	l_x_y->SetFillStyle(0);
	l_x_y->SetBorderSize(0);
	l_x_y->SetHeader(str_mod + " case");
	l_x_y->AddEntry(g_Earth_x_y, "Earth", "l");
	l_x_y->AddEntry(g_Moon_x_y, "Moon", "l");
	l_x_y->AddEntry(g_Sun_x_y, "Sun", "p");
	l_x_y->Draw("same");
	c_x_y->SaveAs("./plot/Earth_Moon_"+str_mod+"_x_y.pdf");

	TCanvas *c_zoomin_x_y = new TCanvas("zoomin_x_y", "zoomin_x_y");
	h2frame->Draw();
	g_Earth_x_y->Draw("L same");
	g_Moon_zoomin_x_y->Draw("L same");
	g_Sun_x_y->Draw("P same");
	TLegend *l_zoomin_x_y = new TLegend(0.1,0.7,0.3,0.9);
	l_zoomin_x_y->SetHeader(str_mod + " case");
	l_zoomin_x_y->SetFillStyle(0);
	l_zoomin_x_y->SetBorderSize(0);
	l_zoomin_x_y->AddEntry(g_Earth_x_y, "Earth", "l");
	l_zoomin_x_y->AddEntry(g_Moon_zoomin_x_y, Form("Moon (R_{EM} #times %d)", int(zoomin)), "l");
	l_zoomin_x_y->AddEntry(g_Sun_x_y, "Sun", "p");
	l_zoomin_x_y->Draw("same");
	c_zoomin_x_y->SaveAs("./plot/zoomin_"+str_mod+"_x_y.pdf");

	TCanvas *c_relative_x_y = new TCanvas("relative_x_y", "relative_x_y");
	g_relative_x_y->Draw("LA");
	c_relative_x_y->SaveAs("./plot/relative_"+str_mod+"_x_y.pdf");
}

//-------------------------------------------------------------------------//

void capacitor_2d_plot(int alg, double d, double acc) {

	TString str_alg[3] = {"Jacobi", "GaussSeidel", "SOR"};

	TString str_d = Form("d%.4d", int(d*1000));
	TString str_acc = Form("acc%.1d", int(fabs(log10(acc))));

	Capacitor2D c2d(0.2, 0.3, d, acc);
	c2d.set_alg(alg);
	c2d.cal();
	vector< vector<double> > V = c2d.get_V();
	int nx = V[0].size();
	int ny = V.size();
	double tmp_acc = c2d.get_tmp_acc();
	int n_iter = c2d.get_n_iter();
	double alpha = c2d.get_alpha();

	cout << tmp_acc << " " << nx << " " << ny << " " << n_iter << " " << alpha << endl;

	TH2F *h2_V = new TH2F("capacitor2d", "capacitor2d", nx,-0.5,0.5, ny,-0.5,0.5);
	h2_V->GetXaxis()->SetTitle("x/l");
	h2_V->GetYaxis()->SetTitle("y/l");
	TH2F *h2_minus_V = new TH2F("capacitor2d_minus", "capacitor2d_minus", nx,-0.5,0.5, ny,-0.5,0.5);
	h2_minus_V->GetXaxis()->SetTitle("x/l");
	h2_minus_V->GetYaxis()->SetTitle("y/l");

	for(int j=0; j<ny; j++) {
		for(int i=0; i<nx; i++) {
			h2_V->SetBinContent(i+1, j+1, V[j][i]);
			h2_minus_V->SetBinContent(i+1, j+1, -V[j][i]);
			//cout << V[j][i] << flush;
		}
		//cout << endl;
	}

	TString str_adj = str_alg[alg] + "_" + str_d + "_" + str_acc;
	TString str_tmp;

	TLegend *l_iter = new TLegend(0.1,0.9,0.7,0.98);
	l_iter->SetFillStyle(0);
	l_iter->SetBorderSize(0);
	l_iter->AddEntry((TObject*)0, Form("number of iteration: %d",n_iter), "");
	l_iter->AddEntry((TObject*)0, Form("d = %.3fl, accuracy: %.2e < %.2e",d,acc,tmp_acc), "");

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
	c_E_field->SaveAs("./plot/"+str_tmp+".pdf");

	delete h2_V;
	delete h2_minus_V;
}

//-------------------------------------------------------------------------//

void charge_dist_3d_plot() {

	ChargeDist3D cd3d(0.05, 1e-6);
	cd3d.set_rho(-0.4,-0.4,0,0,0,0, 1);
	cd3d.cal();
	double tmp_acc = cd3d.get_tmp_acc();
	int n_iter = cd3d.get_n_iter();
	double alpha = cd3d.get_alpha();

	vector< vector< vector<double> > > V = cd3d.get_V();
	int nx = V[0][0].size();
	int ny = V[0].size();
	int nz = V.size();

	TH2D *h2_V = new TH2D("charge_dist_3d", "charge_dist_3d", nx,-0.5,0.5, ny,-0.5,0.5);
	h2_V->GetXaxis()->SetTitle("x/L");
	h2_V->GetYaxis()->SetTitle("y/L");
	TH2D *h2_minus_V = new TH2D("charge_dist_3d_minus", "charge_dist_3d_minus", nx,-0.5,0.5, ny,-0.5,0.5);
	h2_minus_V->GetXaxis()->SetTitle("x/L");
	h2_minus_V->GetYaxis()->SetTitle("y/L");
	for(int j=0; j<ny; j++) {
		for(int i=0; i<nx; i++) {
			h2_V->SetBinContent(i+1, j+1, V[(nz+1)/2][j][i]);
			h2_minus_V->SetBinContent(i+1, j+1, -V[(nz+1)/2][j][i]);
		}
	}

	//TString str_adj = str_alg[alg] + "_" + str_d + "_" + str_acc;
	TString str_adj = "tmp3d";
	TString str_tmp;

	TLegend *l_iter = new TLegend(0.1,0.9,0.7,0.98);
	l_iter->SetFillStyle(0);
	l_iter->SetBorderSize(0);
	l_iter->AddEntry((TObject*)0, Form("number of iteration: %d",n_iter), "");
	//l_iter->AddEntry((TObject*)0, Form("d = %.3fL, accuracy: %.2e < %.2e",d,acc,tmp_acc), "");

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
	h2_minus_V->Draw("ARR");
	l_iter->Draw("same");
	c_E_field->SaveAs("./plot/"+str_tmp+".pdf");
	c_E_field->SaveAs("./plot/"+str_tmp+".pdf");

	delete h2_V;
}

//-------------------------------------------------------------------------//

double real_oscillator_period(double x) {

	return sqrt(8.0)/sqrt(cos(x)-cos(G_MX));
}

//-------------------------------------------------------------------------//

void nintegrate_1d_cal() {

	map<double,double> exact;
	exact[0.05] = 6.28417;
	exact[0.1]  = 6.28711;
	exact[0.5]  = 6.38279;
	exact[1.0]  = 6.69998;
	exact[2.0]  = 8.34975;
	exact[3.0]  = 16.1555;
	double T = exact[G_MX];

	NIntegrate1D ni1d(1e-8,G_MX-1e-8, real_oscillator_period);
	ni1d.set_n(1000000);
	double output;

	ni1d.set_alg(0);
	output = ni1d.cal();
	cout << "trapezoidal" << endl;
	cout << "xm = " << G_MX << ", T = " << output << ", err = " << fabs(output-T)/T << endl;;

	ni1d.set_alg(1);
	output = ni1d.cal();
	cout << "Simpson" << endl;
	cout << "xm = " << G_MX << ", T = " << output << ", err = " << fabs(output-T)/T << endl;;
	
}

//-------------------------------------------------------------------------//

int main() {

	gStyle->SetTitleOffset(1.4,"y");
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);

	// (1)
	relativity_precession_cal();

	// (2)
	two_planet_orbit_plot(0);
	two_planet_orbit_plot(1);
	two_planet_orbit_plot(2);

	// (3)
	//capacitor_2d_plot(0, 0.005, 1e-6);
	//capacitor_2d_plot(0, 0.01, 1e-6);
	//capacitor_2d_plot(0, 0.02, 1e-6);
	//capacitor_2d_plot(0, 0.04, 1e-6);
	//capacitor_2d_plot(0, 0.05, 1e-6);
	//
	//capacitor_2d_plot(1, 0.005, 1e-6);
	//capacitor_2d_plot(1, 0.01, 1e-6);
	//capacitor_2d_plot(1, 0.02, 1e-6);
	//capacitor_2d_plot(1, 0.04, 1e-6);
	//capacitor_2d_plot(1, 0.05, 1e-6);

	//capacitor_2d_plot(2, 0.005, 1e-6);
	//capacitor_2d_plot(2, 0.01, 1e-6);
	//capacitor_2d_plot(2, 0.02, 1e-6);
	capacitor_2d_plot(2, 0.04, 1e-6);
	capacitor_2d_plot(2, 0.05, 1e-6);
	capacitor_2d_plot(2, 0.08, 1e-6);
	capacitor_2d_plot(2, 0.1, 1e-6);

	// (4)
	charge_dist_3d_plot();

	// (5)
	G_MX = 0.05;
	nintegrate_1d_cal();
	G_MX = 0.1;
	nintegrate_1d_cal();
	G_MX = 0.5;
	nintegrate_1d_cal();
	G_MX = 1;
	nintegrate_1d_cal();
	G_MX = 2;
	nintegrate_1d_cal();
	G_MX = 3;
	nintegrate_1d_cal();

	return 0;
}

//-------------------------------------------------------------------------//

