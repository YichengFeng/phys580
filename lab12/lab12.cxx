#include <iostream>
#include <vector>
#include <cmath>
#include "molecular_dynamics_2d.h"

#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TH2.h>
#include <TMultiGraph.h>

using namespace std;


//-------------------------------------------------------------------------//

void molecular_dynamics_2d_plot(double factor, double t_change) {

	int N =25;
	double L = 5;

	MolecularDynamics2D md2d(N, L, 0.005, 0.05, 0);

	const int n_ends = 8;
	double t_ends[n_ends] = {0.0,   0.1,  0.2,   4,    13,    16,   22,    24};
	bool  to_draw[n_ends] = {false, true, false, true, false, true, false, true}; 

	TMultiGraph *mg_E = new TMultiGraph();
	TMultiGraph *mg_T = new TMultiGraph();
	TMultiGraph *mg_tag_displace = new TMultiGraph();
	TMultiGraph *mg_tag_distance = new TMultiGraph();

	TString str_vch = Form("Vch%.2d", int(factor*10));
	TString str_adj;
	TString str_tmp;

	for(int i_end=1; i_end<n_ends; i_end++) {

		if(t_ends[i_end] == t_change) md2d.change_velocity(factor);
		md2d.cal_until(t_ends[i_end]);

		vector<double> t = md2d.get_t();
		vector<double> E = md2d.get_E();
		vector<double> T = md2d.get_T();
		vector<double> tag_displace = md2d.get_tag_displace();
		vector<double> tag_distance = md2d.get_tag_distance();
		vector< vector<double> > x = md2d.get_x();
		vector< vector<double> > y = md2d.get_y();
		
		int nt = t.size();

		TGraph *g_E = new TGraph(nt, &t[0], &E[0]);
		mg_E->Add(g_E);
		TGraph *g_T = new TGraph(nt, &t[0], &T[0]);
		mg_T->Add(g_T);
		TGraph *g_tag_displace = new TGraph(nt, &t[0], &tag_displace[0]);
		mg_tag_displace->Add(g_tag_displace);
		TGraph *g_tag_distance = new TGraph(nt, &t[0], &tag_distance[0]);
		mg_tag_distance->Add(g_tag_distance);
		
		if(!to_draw[i_end]) continue;

		str_adj = Form("t%.3d%.3d", int(t_ends[i_end-1]*10), int(t_ends[i_end]*10));
		str_adj = str_adj + "_" + str_vch;
		str_tmp = str_adj + "_xy";
		TCanvas *c_xy = new TCanvas("c_"+str_tmp, "c_"+str_tmp, 600,600);
		TH2D *hframe = new TH2D("frame_"+str_tmp, "frame_"+str_tmp, 1,0,L, 1,0,L);
		hframe->GetXaxis()->SetTitle("x/#sigma");
		hframe->GetYaxis()->SetTitle("y/#sigma");
		hframe->Draw();
		for(int i=0; i<nt; i++) {
			int color = i==nt-1?2:14;
			int style = i==nt-1?20:24;
			double size  = i==nt-1?0.5:0.3;
			TGraph *g_xy = new TGraph(N, &x[i][0], &y[i][0]);
			g_xy->SetMarkerStyle(style);
			g_xy->SetMarkerSize(size);
			g_xy->SetMarkerColor(color);
			g_xy->Draw("P same");
		}
		TLegend *l_time = new TLegend(0.1, 0.9, 0.9, 0.96);
		l_time->SetFillStyle(0);
		l_time->SetBorderSize(0);
		TString str_par = Form("t = %.1f ~ %.1f. (at end, E = %.2f, T = %.2f)", t_ends[i_end-1], t_ends[i_end], E[nt-1], T[nt-1]);
		l_time->AddEntry((TObject*)0, str_par, "");
		l_time->Draw("same");
		c_xy->SaveAs("./plot/"+str_tmp+".pdf");
	}

	str_adj = Form("t%.3d%.3d", int(t_ends[0]*10), int(t_ends[n_ends-1]*10));
	str_adj = str_adj + "_" + str_vch;

	str_tmp = str_adj + "_E";
	TCanvas *c_E = new TCanvas("c_"+str_tmp, "c_"+str_tmp);
	mg_E->GetXaxis()->SetTitle("t");
	mg_E->GetYaxis()->SetTitle("E");
	mg_E->Draw("PLA");
	c_E->SaveAs("./plot/"+str_tmp+".pdf");

	str_tmp = str_adj + "_T";
	TCanvas *c_T = new TCanvas("c_"+str_tmp, "c_"+str_tmp);
	mg_T->GetXaxis()->SetTitle("t");
	mg_T->GetYaxis()->SetTitle("T");
	mg_T->Draw("PLA");
	c_T->SaveAs("./plot/"+str_tmp+".pdf");

	str_tmp = str_adj + "_tag_displace";
	TCanvas *c_tag_displace = new TCanvas("c_"+str_tmp, "c_"+str_tmp);
	mg_tag_displace->GetXaxis()->SetTitle("t");
	mg_tag_displace->GetYaxis()->SetTitle("tag_displace");
	mg_tag_displace->Draw("PLA");
	c_tag_displace->SaveAs("./plot/"+str_tmp+".pdf");

	str_tmp = str_adj + "_tag_distance";
	TCanvas *c_tag_distance = new TCanvas("c_"+str_tmp, "c_"+str_tmp);
	mg_tag_distance->GetXaxis()->SetTitle("t");
	mg_tag_distance->GetYaxis()->SetTitle("tag_distance");
	mg_tag_distance->Draw("PLA");
	c_tag_distance->SaveAs("./plot/"+str_tmp+".pdf");
}

//-------------------------------------------------------------------------//

void molecular_dynamics_2d_set(double L, double vmax) {

	int N =25;
	//double L = 5;

	MolecularDynamics2D md2d(N, L, 0.005, 0, vmax);

	const int n_ends = 10;
	double t_ends[n_ends] = {0.0,   0.1,  0.2,   4,    13,    16,   22,    24,   40,    45};
	bool  to_draw[n_ends] = {false, true, false, true, false, true, false, true, false, true}; 

	TMultiGraph *mg_E = new TMultiGraph();
	TMultiGraph *mg_T = new TMultiGraph();
	TMultiGraph *mg_tag_displace = new TMultiGraph();
	TMultiGraph *mg_tag_distance = new TMultiGraph();

	TString str_adj;
	TString str_tmp;
	TString str_set;

	str_set = Form("L%.2dvmax%.3d", int(L), int(vmax*100));

	for(int i_end=1; i_end<n_ends; i_end++) {

		md2d.cal_until(t_ends[i_end]);

		vector<double> t = md2d.get_t();
		vector<double> E = md2d.get_E();
		vector<double> T = md2d.get_T();
		vector<double> tag_displace = md2d.get_tag_displace();
		vector<double> tag_distance = md2d.get_tag_distance();
		vector< vector<double> > x = md2d.get_x();
		vector< vector<double> > y = md2d.get_y();
		
		int nt = t.size();

		TGraph *g_E = new TGraph(nt, &t[0], &E[0]);
		mg_E->Add(g_E);
		TGraph *g_T = new TGraph(nt, &t[0], &T[0]);
		mg_T->Add(g_T);
		TGraph *g_tag_displace = new TGraph(nt, &t[0], &tag_displace[0]);
		mg_tag_displace->Add(g_tag_displace);
		TGraph *g_tag_distance = new TGraph(nt, &t[0], &tag_distance[0]);
		mg_tag_distance->Add(g_tag_distance);
		
		if(!to_draw[i_end]) continue;

		str_adj = Form("t%.3d%.3d", int(t_ends[i_end-1]*10), int(t_ends[i_end]*10));
		str_adj = str_adj + "_" + str_set;
		str_tmp = str_adj + "_xy";
		TCanvas *c_xy = new TCanvas("c_"+str_tmp, "c_"+str_tmp, 600,600);
		TH2D *hframe = new TH2D("frame_"+str_tmp, "frame_"+str_tmp, 1,0,L, 1,0,L);
		hframe->GetXaxis()->SetTitle("x/#sigma");
		hframe->GetYaxis()->SetTitle("y/#sigma");
		hframe->Draw();
		for(int i=0; i<nt; i++) {
			int color = i==nt-1?2:14;
			int style = i==nt-1?20:24;
			double size  = i==nt-1?0.5:0.3;
			TGraph *g_xy = new TGraph(N, &x[i][0], &y[i][0]);
			g_xy->SetMarkerStyle(style);
			g_xy->SetMarkerSize(size);
			g_xy->SetMarkerColor(color);
			g_xy->Draw("P same");
		}
		TLegend *l_time = new TLegend(0.1, 0.9, 0.9, 0.96);
		l_time->SetFillStyle(0);
		l_time->SetBorderSize(0);
		TString str_par = Form("t = %.1f ~ %.1f. (at end, E = %.2f, T = %.2f)", t_ends[i_end-1], t_ends[i_end], E[nt-1], T[nt-1]);
		l_time->AddEntry((TObject*)0, str_par, "");
		l_time->Draw("same");
		c_xy->SaveAs("./plot/"+str_tmp+".pdf");
	}

	str_adj = Form("t%.3d%.3d", int(t_ends[0]*10), int(t_ends[n_ends-1]*10));
	str_adj = str_adj + "_" + str_set;

	str_tmp = str_adj + "_E";
	TCanvas *c_E = new TCanvas("c_"+str_tmp, "c_"+str_tmp);
	mg_E->GetXaxis()->SetTitle("t");
	mg_E->GetYaxis()->SetTitle("E");
	mg_E->Draw("PLA");
	c_E->SaveAs("./plot/"+str_tmp+".pdf");

	str_tmp = str_adj + "_T";
	TCanvas *c_T = new TCanvas("c_"+str_tmp, "c_"+str_tmp);
	mg_T->GetXaxis()->SetTitle("t");
	mg_T->GetYaxis()->SetTitle("T");
	mg_T->Draw("PLA");
	c_T->SaveAs("./plot/"+str_tmp+".pdf");

	str_tmp = str_adj + "_tag_displace";
	TCanvas *c_tag_displace = new TCanvas("c_"+str_tmp, "c_"+str_tmp);
	mg_tag_displace->GetXaxis()->SetTitle("t");
	mg_tag_displace->GetYaxis()->SetTitle("tag_displace");
	mg_tag_displace->Draw("PLA");
	c_tag_displace->SaveAs("./plot/"+str_tmp+".pdf");

	str_tmp = str_adj + "_tag_distance";
	TCanvas *c_tag_distance = new TCanvas("c_"+str_tmp, "c_"+str_tmp);
	mg_tag_distance->GetXaxis()->SetTitle("t");
	mg_tag_distance->GetYaxis()->SetTitle("tag_distance");
	mg_tag_distance->Draw("PLA");
	c_tag_distance->SaveAs("./plot/"+str_tmp+".pdf");
}

//-------------------------------------------------------------------------//

int main() {

	gStyle->SetTitleOffset(1.4,"y");
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);
	gStyle->SetOptFit();

	//molecular_dynamics_2d_plot(1.0, 13);
	//molecular_dynamics_2d_plot(0.5, 13);
	//molecular_dynamics_2d_plot(1.5, 4);
	//molecular_dynamics_2d_plot(5.0, 16);
	//
	//molecular_dynamics_2d_set(5, 0.02);
	//molecular_dynamics_2d_set(5, 0.05);
	//molecular_dynamics_2d_set(5, 0.10);
	//molecular_dynamics_2d_set(5, 0.20);
	molecular_dynamics_2d_set(10, 2);
	molecular_dynamics_2d_set(10, 5);
	molecular_dynamics_2d_set(10, 10);
	//molecular_dynamics_2d_set(10, 0.20);

	return 0;
}
