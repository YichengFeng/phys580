#include <iostream>
#include "nuclei_decay.h"

#include <TROOT.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TAxis.h>

using namespace std;


int main()
{
	// get data
	const int ndt = 5;
	double dt[ndt] = {0.05, 0.2, 0.8, 1, 1.5};
	NucleiDecay* nd[ndt];
	TGraph* g_euler[ndt];
	TGraph* g_exact[ndt];
	TGraph* g_error[ndt];
	TGraph* g_relative_error[ndt];
	TGraph* g_rk2[ndt];
	TGraph* g_rk2_error[ndt];
	TGraph* g_rk2_relative_error[ndt];
	TGraph* g_rk4[ndt];
	TGraph* g_rk4_error[ndt];
	TGraph* g_rk4_relative_error[ndt];

	const int nal = 3; // number of algorithms
	TString algor_name[nal] = { "Euler", "RK2", "RK4" };
	TGraph* gall_algor[nal][ndt];
	TGraph* gall_error[nal][ndt];
	TGraph* gall_relative_error[nal][ndt];

	for(int i=0; i<ndt; i++) {
		nd[i] = new NucleiDecay(dt[i]);
		nd[i]->cal_g();

		g_euler[i] = nd[i]->get_g_euler();
		g_exact[i] = nd[i]->get_g_exact();
		g_error[i] = nd[i]->get_g_error();
		g_relative_error[i] = nd[i]->get_g_relative_error();
		g_rk2[i] = nd[i]->get_g_rk2();
		g_rk2_error[i] = nd[i]->get_g_rk2_error();
		g_rk2_relative_error[i] = nd[i]->get_g_rk2_relative_error();
		g_rk4[i] = nd[i]->get_g_rk4();
		g_rk4_error[i] = nd[i]->get_g_rk4_error();
		g_rk4_relative_error[i] = nd[i]->get_g_rk4_relative_error();

		gall_algor[0][i] = g_euler[i];
		gall_algor[1][i] = g_rk2[i];
		gall_algor[2][i] = g_rk4[i];

		gall_error[0][i] = g_error[i];
		gall_error[1][i] = g_rk2_error[i];
		gall_error[2][i] = g_rk4_error[i];

		gall_relative_error[0][i] = g_relative_error[i];
		gall_relative_error[1][i] = g_rk2_relative_error[i];
		gall_relative_error[2][i] = g_rk4_relative_error[i];
	}

	// plot
	TLegend *l[nal][ndt];
	TMultiGraph *mg[nal][ndt];
	TCanvas *c[nal][ndt];

	for(int ial=0; ial<nal; ial++) {
	for(int idt=0; idt<ndt; idt++) {
		l[ial][idt] = new TLegend(0.7, 0.7, 0.9, 0.9);
		l[ial][idt]->SetFillStyle(0);
		l[ial][idt]->SetHeader(Form("#Deltat = %.2f#tau", dt[idt]));

		mg[ial][idt] = new TMultiGraph();
		
		g_exact[0]->SetLineColor(kBlack);
		g_exact[0]->SetLineWidth(2);
		g_exact[0]->SetLineStyle(7);
		g_exact[0]->SetMarkerColor(kBlack);
		g_exact[0]->GetXaxis()->SetRangeUser(0, 5);
		l[ial][idt]->AddEntry(g_exact[0], "exact", "l");
		mg[ial][idt]->Add(g_exact[0]);

		gall_algor[ial][idt]->SetLineColor(kBlue);
		gall_algor[ial][idt]->SetMarkerColor(kBlue);
		gall_algor[ial][idt]->SetMarkerStyle(5);
		gall_algor[ial][idt]->GetXaxis()->SetRangeUser(0, 5);
		l[ial][idt]->AddEntry(gall_algor[ial][idt], algor_name[ial], "lp");
		mg[ial][idt]->Add(gall_algor[ial][idt]);
		
		gall_error[ial][idt]->SetLineColor(kGreen);
		gall_error[ial][idt]->SetMarkerColor(kGreen);
		gall_error[ial][idt]->GetXaxis()->SetRangeUser(0, 5);
		l[ial][idt]->AddEntry(gall_error[ial][idt], "error", "lp");
		mg[ial][idt]->Add(gall_error[ial][idt]);
		
		c[ial][idt] = new TCanvas(algor_name[ial]+Form("%.2f", dt[idt]), algor_name[ial]+Form("%.2f", dt[idt]));
		mg[ial][idt]->GetXaxis()->SetLimits(0, 5);
		mg[ial][idt]->Draw("PLA");
		mg[ial][idt]->GetXaxis()->SetLimits(0, 5);
		mg[ial][idt]->GetXaxis()->SetTitle("time scaled by mean life: t/#tau");
		mg[ial][idt]->GetYaxis()->SetTitle("number of nuclei: N");
		l[ial][idt]->Draw("same");

		c[ial][idt]->Modified();
		c[ial][idt]->SaveAs("./plot/"+algor_name[ial]+Form("_dt%.2f.pdf",dt[idt]));
	}
	}

	TCanvas *c_algor[nal];
	TCanvas *c_error[nal];
	TCanvas *c_relative_error[nal];
	TLegend *l_algor[nal];
	TLegend *l_error[nal];
	TLegend *l_relative_error[nal];
	TMultiGraph *mg_algor[nal];
	TMultiGraph *mg_error[nal];
	TMultiGraph *mg_relative_error[nal];
	Int_t color_dt[ndt] = {kBlue, 6, kRed, kOrange, kGreen};
	for(int ial=0; ial<nal; ial++) {
		c_algor[ial] = new TCanvas("algor_"+algor_name[ial], "algor_"+algor_name[ial]);
		l_algor[ial] = new TLegend(0.7, 0.5, 0.9, 0.9);
		l_algor[ial]->SetBorderSize(0);
		l_algor[ial]->SetFillStyle(0);
		l_algor[ial]->SetHeader(algor_name[ial]);
		l_algor[ial]->AddEntry(g_exact[0], "exact", "l");
		mg_algor[ial] = new TMultiGraph();
		mg_algor[ial]->Add(g_exact[0]);
		for(int idt=0; idt<ndt; idt++) {
			gall_algor[ial][idt]->SetLineColor(color_dt[idt]);
			gall_algor[ial][idt]->SetMarkerColor(color_dt[idt]);
			mg_algor[ial]->Add(gall_algor[ial][idt], "l"); 
			l_algor[ial]->AddEntry(gall_algor[ial][idt], Form("#Deltat = %.2f#tau",dt[idt]), "l"); 
		}
		mg_algor[ial]->GetXaxis()->SetTitle("time scaled by mean life: t/#tau");
		mg_algor[ial]->GetYaxis()->SetTitle("number of nuclei: N");
		mg_algor[ial]->Draw("LA");
		l_algor[ial]->Draw("same");
		c_algor[ial]->SaveAs("./plot/algor_"+algor_name[ial]+".pdf");

		c_error[ial] = new TCanvas("error_"+algor_name[ial], "error_"+algor_name[ial]);
		l_error[ial] = new TLegend(0.7, 0.5, 0.9, 0.9);
		l_error[ial]->SetBorderSize(0);
		l_error[ial]->SetFillStyle(0);
		//l_error[ial]->SetHeader(algor_name[ial]+": deviation form the exact");
		l_error[ial]->SetHeader(algor_name[ial]);
		mg_error[ial] = new TMultiGraph();
		for(int idt=0; idt<ndt; idt++) {
			gall_error[ial][idt]->SetLineColor(color_dt[idt]);
			gall_error[ial][idt]->SetMarkerColor(color_dt[idt]);
			mg_error[ial]->Add(gall_error[ial][idt], "l"); 
			l_error[ial]->AddEntry(gall_error[ial][idt], Form("#Deltat = %.2f#tau",dt[idt]), "l"); 
		}
		mg_error[ial]->GetXaxis()->SetTitle("time scaled by mean life: t/#tau");
		mg_error[ial]->GetYaxis()->SetTitle("deviation from the exact: #DeltaN");
		mg_error[ial]->Draw("LA");
		l_error[ial]->Draw("same");
		c_error[ial]->SaveAs("./plot/error_"+algor_name[ial]+".pdf");

		c_relative_error[ial] = new TCanvas("relative_error_"+algor_name[ial], "relative_error_"+algor_name[ial]);
		l_relative_error[ial] = new TLegend(0.15, 0.5, 0.35, 0.9);
		l_relative_error[ial]->SetBorderSize(0);
		l_relative_error[ial]->SetFillStyle(0);
		//l_relative_error[ial]->SetHeader(algor_name[ial]+": relative deviation from the exact");
		l_relative_error[ial]->SetHeader(algor_name[ial]);
		mg_relative_error[ial] = new TMultiGraph();
		for(int idt=0; idt<ndt; idt++) {
			gall_relative_error[ial][idt]->SetLineColor(color_dt[idt]);
			gall_relative_error[ial][idt]->SetMarkerColor(color_dt[idt]);
			mg_relative_error[ial]->Add(gall_relative_error[ial][idt], "l"); 
			l_relative_error[ial]->AddEntry(gall_relative_error[ial][idt], Form("#Deltat = %.2f#tau",dt[idt]), "l"); 
		}
		mg_relative_error[ial]->GetXaxis()->SetTitle("time scaled by mean life: t/#tau");
		mg_relative_error[ial]->GetYaxis()->SetTitle("relative deviation from the exact: #DeltaN/N");
		mg_relative_error[ial]->Draw("LA");
		l_relative_error[ial]->Draw("same");
		c_relative_error[ial]->SaveAs("./plot/relative_error_"+algor_name[ial]+".pdf");
	}

	double fixed_t = 4;
	//----------------------------------------------//
	TCanvas *c_dt[nal];
	TLegend *l_dt[nal];
	TGraph *g_dt[nal];
	double x_dt[nal][ndt];
	double y_dt[nal][ndt];
	TCanvas *c_dt_log = new TCanvas("dt_dependence_log", "dt_dependence_log");
	TLegend *l_dt_log = new TLegend(0.7,0.1,0.9,0.3);
	TMultiGraph *mg_dt_log = new TMultiGraph();
	//----------------------------------------------//
	for(int ial=0; ial<nal; ial++) {
		c_dt[ial] = new TCanvas("dt_"+algor_name[ial], "dt_"+algor_name[ial]);
		l_dt[ial] = new TLegend(0.1,0.8,0.2,0.9);
		l_dt[ial]->SetHeader(algor_name[ial]);

		for(int idt=0; idt<ndt; idt++) {
			int ipt = (int)(fixed_t/dt[idt]);
			double *xpt = gall_error[ial][idt]->GetX();
			double *ypt = gall_error[ial][idt]->GetY();

			x_dt[ial][idt] = dt[idt];
			y_dt[ial][idt] = fabs(ypt[ipt]);
		}

		g_dt[ial] = new TGraph(ndt, x_dt[ial], y_dt[ial]);
		g_dt[ial]->GetXaxis()->SetTitle("step size: #Deltat/#tau");
		g_dt[ial]->GetYaxis()->SetTitle(Form("absolute deviation from the exact: #DeltaN|_{t=%.1f#tau}",fixed_t));
		g_dt[ial]->SetLineColor(kBlue);
		g_dt[ial]->SetMarkerColor(kBlue);
		g_dt[ial]->SetTitle("");
		g_dt[ial]->Draw("PLA");
		l_dt[ial]->Draw("same");

		c_dt[ial]->SaveAs("./plot/error_dt_dependence_"+algor_name[ial]+".pdf");

		g_dt[ial]->SetLineColor(color_dt[ial]);
		g_dt[ial]->SetMarkerColor(color_dt[ial]);
		mg_dt_log->Add(g_dt[ial]);
		l_dt_log->AddEntry(g_dt[ial], algor_name[ial], "lp");
	}

	c_dt_log->cd();
	mg_dt_log->GetXaxis()->SetTitle("step size: #Deltat/#tau");
	mg_dt_log->GetYaxis()->SetTitle(Form("deviation from the exact: |#DeltaN|_{t=%.1f#tau}",fixed_t));
	mg_dt_log->Draw("PLA");
	l_dt_log->Draw("same");
	c_dt_log->SetLogy();
	c_dt_log->SaveAs("./plot/error_dt_dependence_log.pdf");

	//----------------------------------------------//
	TCanvas *c_relative_dt[nal];
	TLegend *l_relative_dt[nal];
	TGraph *g_relative_dt[nal];
	double x_relative_dt[nal][ndt];
	double y_relative_dt[nal][ndt];
	TCanvas *c_relative_dt_log = new TCanvas("relative_dt_dependence_log", "relative_dt_dependence_log");
	TLegend *l_relative_dt_log = new TLegend(0.7,0.1,0.9,0.3);
	TMultiGraph *mg_relative_dt_log = new TMultiGraph();
	//----------------------------------------------//
	for(int ial=0; ial<nal; ial++) {
		c_relative_dt[ial] = new TCanvas("relative_dt_"+algor_name[ial], "relative_dt_"+algor_name[ial]);
		l_relative_dt[ial] = new TLegend(0.1,0.8,0.2,0.9);
		l_relative_dt[ial]->SetHeader(algor_name[ial]);

		for(int idt=0; idt<ndt; idt++) {
			int ipt = (int)(fixed_t/dt[idt]);
			double *xpt = gall_relative_error[ial][idt]->GetX();
			double *ypt = gall_relative_error[ial][idt]->GetY();

			x_relative_dt[ial][idt] = dt[idt];
			y_relative_dt[ial][idt] = fabs(ypt[ipt]);
		}

		g_relative_dt[ial] = new TGraph(ndt, x_relative_dt[ial], y_relative_dt[ial]);
		g_relative_dt[ial]->GetXaxis()->SetTitle("step size: #Deltat/#tau");
		g_relative_dt[ial]->GetYaxis()->SetTitle(Form("relative deviation from the exact: |#DeltaN/N|_{t=%.1f#tau}",fixed_t));
		g_relative_dt[ial]->SetLineColor(kBlue);
		g_relative_dt[ial]->SetMarkerColor(kBlue);
		g_relative_dt[ial]->SetTitle("");
		g_relative_dt[ial]->Draw("PLA");
		l_relative_dt[ial]->Draw("same");

		c_relative_dt[ial]->SaveAs("./plot/relative_error_dt_dependence_"+algor_name[ial]+".pdf");

		g_relative_dt[ial]->SetLineColor(color_dt[ial]);
		g_relative_dt[ial]->SetMarkerColor(color_dt[ial]);
		mg_relative_dt_log->Add(g_relative_dt[ial]);
		l_relative_dt_log->AddEntry(g_relative_dt[ial], algor_name[ial], "lp");
	}

	c_relative_dt_log->cd();
	mg_relative_dt_log->GetXaxis()->SetTitle("step size: #Deltat/#tau");
	mg_relative_dt_log->GetYaxis()->SetTitle(Form("relative deviation from the exact: |#DeltaN/N|_{t=%.1f#tau}",fixed_t));
	mg_relative_dt_log->Draw("PLA");
	l_relative_dt_log->Draw("same");
	c_relative_dt_log->SetLogy();
	c_relative_dt_log->SaveAs("./plot/error_relative_dt_dependence_log.pdf");

	return 0;
}
