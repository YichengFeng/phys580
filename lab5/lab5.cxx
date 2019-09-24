#include <iostream>
#include <vector>
#include <cmath>
#include "runge_kutta.h"
#include "planet_orbit.h"
#include "two_planet_orbit.h"
#include "three_body_2d.h"

#include <TString.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TH2.h>
#include <TMultiGraph.h>

using namespace std;


//-------------------------------------------------------------------------//

void one_planet_orbit_plot(double r_min, double e, double mp, double t_end, TString str_name, bool is_full) {

	TString str_adj = str_name;
	TString str_tmp;	

	// sun
	double x_sun[1] = {0};
	double y_sun[1] = {0};
	TGraph *g_sun = new TGraph(1, x_sun, y_sun);
	g_sun->SetMarkerStyle(20);
	g_sun->SetMarkerSize(2);
	g_sun->SetMarkerColor(kOrange);

	// planet
	double a = r_min/(1-e);
	double r_max = a*(1+e);
	double f = (r_max-r_min)/2;
	double b = sqrt(a*a-f*f);

	PlanetOrbit po_ellipse(a, e, mp);
	//po_ellipse.set_alg(4);
	po_ellipse.set_t_end(t_end);
	po_ellipse.cal();
	vector<double> t_ellipse = po_ellipse.get_t();
	vector< vector<double> > x_ellipse = po_ellipse.get_x();
	int n_stps_ellipse = po_ellipse.get_n_stps();

	double period = po_ellipse.cal_period();

	int canvas_nx = 700;
	int canvas_ny = 500;
	double rm = 0.05;
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

	TGraph *g_ellipse = new TGraph(n_stps_ellipse, &x_ellipse[0][0], &x_ellipse[1][0]);
	g_ellipse->GetXaxis()->SetTitle("x (AU)");
	g_ellipse->GetYaxis()->SetTitle("y (AU)");

	TLegend *l_ellipse = new TLegend(0.15,0.7,0.35,0.9);
	l_ellipse->SetFillStyle(0);
	l_ellipse->SetBorderSize(0);
	l_ellipse->AddEntry(g_ellipse, str_name + " (real)", "l");
	l_ellipse->AddEntry(g_sun, "Sun", "p");

	TLegend *l_period = new TLegend(0.65,0.8,0.9,0.9);
	l_period->SetFillStyle(0);
	l_period->SetBorderSize(0);
	l_period->AddEntry((TObject*)0, Form("T = %.6f yr", period), "");

	str_tmp = "ellipse_" + str_name;
	TCanvas *c_ellipse = new TCanvas("c_" + str_tmp, "c_" + str_tmp, canvas_nx,canvas_ny);
	h2frame->Draw();
	g_sun->Draw("P same");
	g_ellipse->Draw("L same");
	l_ellipse->Draw("same");
	l_period->Draw("same");
	c_ellipse->SaveAs("./plot/" + str_tmp + ".pdf");

	if(!is_full) return;

	PlanetOrbit po_circular(r_min, 0, mp);
	//po_circular.increase_v_start(1.0/sqrt(2.0));
	po_circular.set_t_end(0.5);
	po_circular.cal();
	vector<double> t_circular = po_circular.get_t();
	vector< vector<double> > x_circular = po_circular.get_x();
	int n_stps_circular = po_circular.get_n_stps();

	TGraph *g_circular = new TGraph(n_stps_circular, &x_circular[0][0], &x_circular[1][0]);
	g_circular->GetXaxis()->SetTitle("x (AU)");
	g_circular->GetYaxis()->SetTitle("y (AU)");

	TLegend *l_circular = new TLegend(0.15,0.7,0.35,0.9);
	l_circular->SetFillStyle(0);
	l_circular->SetBorderSize(0);
	l_circular->AddEntry(g_circular, str_name + " (circular)", "l");
	l_circular->AddEntry(g_sun, "Sun", "p");

	str_tmp = "circular_" + str_name;
	TCanvas *c_circular = new TCanvas("c_" + str_tmp, "c_" + str_tmp);
	h2frame->Draw();
	g_sun->Draw("P same");
	g_circular->Draw("L same");
	l_circular->Draw("same");
	c_circular->SaveAs("./plot/" + str_tmp + ".pdf");

	PlanetOrbit po_parabola(a, e, mp);
	po_parabola.increase_v_start(1.0);
	po_parabola.set_t_end(0.2);
	po_parabola.cal();
	vector<double> t_parabola = po_parabola.get_t();
	vector< vector<double> > x_parabola = po_parabola.get_x();
	int n_stps_parabola = po_parabola.get_n_stps();

	TGraph *g_parabola = new TGraph(n_stps_parabola, &x_parabola[0][0], &x_parabola[1][0]);
	g_parabola->GetXaxis()->SetTitle("x (AU)");
	g_parabola->GetYaxis()->SetTitle("y (AU)");

	TLegend *l_parabola = new TLegend(0.15,0.7,0.35,0.9);
	l_parabola->SetFillStyle(0);
	l_parabola->SetBorderSize(0);
	l_parabola->AddEntry(g_parabola, str_name + " (parabola)", "l");
	l_parabola->AddEntry(g_sun, "Sun", "p");

	str_tmp = "parabola_" + str_name;
	TCanvas *c_parabola = new TCanvas("c_" + str_tmp, "c_" + str_tmp);
	g_parabola->Draw("LA");
	g_sun->Draw("P same");
	l_parabola->Draw("same");
	c_parabola->SaveAs("./plot/" + str_tmp + ".pdf");

	PlanetOrbit po_hyperbola(a, e, mp);
	po_hyperbola.increase_v_start(1.1);
	po_hyperbola.set_t_end(0.2);
	po_hyperbola.cal();
	vector<double> t_hyperbola = po_hyperbola.get_t();
	vector< vector<double> > x_hyperbola = po_hyperbola.get_x();
	int n_stps_hyperbola = po_hyperbola.get_n_stps();

	TGraph *g_hyperbola = new TGraph(n_stps_hyperbola, &x_hyperbola[0][0], &x_hyperbola[1][0]);
	g_hyperbola->GetXaxis()->SetTitle("x (AU)");
	g_hyperbola->GetYaxis()->SetTitle("y (AU)");

	TLegend *l_hyperbola = new TLegend(0.15,0.7,0.35,0.9);
	l_hyperbola->SetFillStyle(0);
	l_hyperbola->SetBorderSize(0);
	l_hyperbola->AddEntry(g_hyperbola, str_name + " (hyperbola)", "l");
	l_hyperbola->AddEntry(g_sun, "Sun", "p");

	str_tmp = "hyperbola_" + str_name;
	TCanvas *c_hyperbola = new TCanvas("c_" + str_tmp, "c_" + str_tmp);
	g_hyperbola->Draw("LA");
	g_sun->Draw("P same");
	l_hyperbola->Draw("same");
	c_hyperbola->SaveAs("./plot/" + str_tmp + ".pdf");

	delete h2frame;
}

//-------------------------------------------------------------------------//

void two_planet_orbit_plot(double nMJ) {

	TString str_nmj = Form("nMJ%.4d", int(nMJ));
	TString str_adj = str_nmj;
	TString str_tmp;

	double a1 = 1;
	double e1 = 0.017;
	double f1 = a1*e1;
	double b1 = sqrt(a1*a1 - f1*f1);
	double m1 = 5.97e24/1.99e30;
	double a2 = 5.203;
	double e2 = 0.048;
	double f2 = a2*e2;
	double b2 = sqrt(a2*a2 - f2*f2);
	double m2 = nMJ*1.90e27/1.99e30;

	int canvas_nx = 700;
	int canvas_ny = 500;
	double rm = 0.05;
	double rx, ry;
	if(a2/canvas_nx < b2/canvas_ny) {
		rx = b2/canvas_ny*canvas_nx/a2;
		ry = 1;
	} else {
		rx = 1;
		ry = a2/canvas_nx*canvas_ny/b2;
	}
	TH2F *h2frame = new TH2F("frame", "frame", 1,-(a2+f2+rm*a2)*rx,(a2-f2+rm*a2)*rx, 1,-(b2+b2*rm)*ry,(b2+b2*rm)*ry);
	h2frame->GetXaxis()->SetTitle("x (AU)");
	h2frame->GetYaxis()->SetTitle("y (AU)");

	// sun
	double x_sun[1] = {0};
	double y_sun[1] = {0};
	TGraph *g_sun = new TGraph(1, x_sun, y_sun);
	g_sun->SetMarkerStyle(20);
	g_sun->SetMarkerSize(2);
	g_sun->SetMarkerColor(kOrange);

	TwoPlanetOrbit tpo(a1, e1, m1, a2, e2, m2);
	tpo.set_alg(0);
	tpo.set_t_end(20);
	tpo.set_dt(0.001);
	tpo.cal();
	vector<double> t = tpo.get_t();
	vector< vector<double> > x = tpo.get_x();
	int n_stps = tpo.get_n_stps();

	TGraph *g_e = new TGraph(n_stps, &x[0][0], &x[1][0]);
	g_e->SetLineColor(kBlue);
	g_e->SetMarkerColor(kBlue);

	TGraph *g_j = new TGraph(n_stps, &x[2][0], &x[3][0]);
	g_j->SetLineColor(kRed);
	g_j->SetMarkerColor(kRed);

	TMultiGraph *mg = new TMultiGraph();
	mg->Add(g_e);
	mg->Add(g_j);
	mg->GetXaxis()->SetTitle("x (AU)");
	mg->GetYaxis()->SetTitle("y (AU)");

	TLegend *l_body = new TLegend(0.1,0.7,0.3,0.9);
	l_body->AddEntry(g_sun, "Sun", "l");
	l_body->AddEntry(g_j, "Jupiter", "l");
	l_body->AddEntry(g_e, "Earth", "l");

	TLegend *l_mj = new TLegend(0.6,0.8,0.9,0.9);
	l_mj->SetFillStyle(0);
	l_mj->SetBorderSize(0);
	l_mj->AddEntry((TObject*)0, Form("Jupiter mass #times %d", int(nMJ)), "");

	str_tmp = "Earth_Jupiter_" + str_adj;
	TCanvas *c_ej = new TCanvas("c_" + str_tmp, "c_" + str_tmp);
	h2frame->Draw();
	g_sun->Draw("P same");
	g_e->Draw("L same");
	g_j->Draw("L same");
	l_body->Draw("same");
	l_mj->Draw("same");
	//mg->Draw("L same");
	c_ej->SaveAs("./plot/" + str_tmp + ".pdf");

	delete h2frame;
}

//-------------------------------------------------------------------------//

void three_body_2d_plot(double nMJ) {

	TString str_nmj = Form("nMJ%.4d", int(nMJ));
	TString str_adj = str_nmj;
	TString str_tmp;

	double a1 = 1;
	double e1 = 0.017;
	double f1 = a1*e1;
	double b1 = sqrt(a1*a1 - f1*f1);
	double m1 = 5.97e24/1.99e30;
	double a2 = 5.203;
	double e2 = 0.048;
	double f2 = a2*e2;
	double b2 = sqrt(a2*a2 - f2*f2);
	double m2 = nMJ*1.90e27/1.99e30;
	double a3 = 0;
	double e3 = 0;
	double f3 = a3*e3;
	double b3 = sqrt(a3*a3 - f3*f3);
	double m3 = 1;

	int canvas_nx = 700;
	int canvas_ny = 500;
	double rm = 0.05;
	double rx, ry;
	if(a2/canvas_nx < b2/canvas_ny) {
		rx = b2/canvas_ny*canvas_nx/a2;
		ry = 1;
	} else {
		rx = 1;
		ry = a2/canvas_nx*canvas_ny/b2;
	}
	TH2F *h2frame = new TH2F("frame", "frame", 1,-(a2+f2+rm*a2)*rx,(a2-f2+rm*a2)*rx, 1,-(b2+b2*rm)*ry,(b2+b2*rm)*ry);
	h2frame->GetXaxis()->SetTitle("x (AU)");
	h2frame->GetYaxis()->SetTitle("y (AU)");

	ThreeBody2D tb2d(a1,e1,m1, a2,e2,m2, a3,e3,m3);
	tb2d.set_alg(0);
	tb2d.set_t_end(20);
	tb2d.set_dt(0.001);
	tb2d.cal();
	vector<double> t = tb2d.get_t();
	vector< vector<double> > x = tb2d.get_x();
	int n_stps = tb2d.get_n_stps();

	TGraph *g_e = new TGraph(n_stps, &x[0][0], &x[1][0]);
	g_e->SetLineColor(kBlue);
	g_e->SetMarkerColor(kBlue);

	TGraph *g_j = new TGraph(n_stps, &x[2][0], &x[3][0]);
	g_j->SetLineColor(kRed);
	g_j->SetMarkerColor(kRed);

	TGraph *g_s = new TGraph(n_stps, &x[4][0], &x[5][0]);
	g_s->SetLineColor(kOrange);
	g_s->SetMarkerColor(kOrange);

	TMultiGraph *mg = new TMultiGraph();
	mg->Add(g_e);
	mg->Add(g_j);
	mg->Add(g_s);
	mg->GetXaxis()->SetTitle("x (AU)");
	mg->GetYaxis()->SetTitle("y (AU)");

	TLegend *l_body = new TLegend(0.1,0.7,0.3,0.9);
	l_body->AddEntry(g_s, "Sun", "l");
	l_body->AddEntry(g_j, "Jupiter", "l");
	l_body->AddEntry(g_e, "Earth", "l");

	TLegend *l_mj = new TLegend(0.6,0.8,0.9,0.9);
	l_mj->SetFillStyle(0);
	l_mj->SetBorderSize(0);
	l_mj->AddEntry((TObject*)0, Form("Jupiter mass #times %d", int(nMJ)), "");

	str_tmp = "Earth_Jupiter_Sun_" + str_adj;
	TCanvas *c_ej = new TCanvas("c_" + str_tmp, "c_" + str_tmp);
	h2frame->Draw();
	g_s->Draw("P same");
	g_e->Draw("L same");
	g_j->Draw("L same");
	l_body->Draw("same");
	l_mj->Draw("same");
	//mg->Draw("L same");
	c_ej->SaveAs("./plot/" + str_tmp + ".pdf");

	delete h2frame;
}

//-------------------------------------------------------------------------//

int main() {

	gStyle->SetTitleOffset(1.4,"y");
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);

	// (1)
	one_planet_orbit_plot(0.983, 0.017, 5.97e24/1.99e30, 5, "Earth", true);
	one_planet_orbit_plot(0.307, 0.206, 3.30e23/1.99e30, 5, "Mercury", true);
	one_planet_orbit_plot(29.63, 0.249, 1.27e22/1.99e30, 5, "Pluto", true);

	// (2)
	one_planet_orbit_plot(1.933, 0.572, 0, 20, "Shoemaker_Levy_2", false);
	one_planet_orbit_plot(0.589, 0.967, 0, 120, "Halley", false);

	// (3)
	two_planet_orbit_plot(1);
	two_planet_orbit_plot(10);
	two_planet_orbit_plot(100);
	two_planet_orbit_plot(500);
	two_planet_orbit_plot(1000);

	three_body_2d_plot(1);
	three_body_2d_plot(10);
	three_body_2d_plot(100);
	three_body_2d_plot(500);
	three_body_2d_plot(1000);

	return 0;
}
