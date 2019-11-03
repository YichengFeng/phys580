#include <iostream>
#include <vector>
#include <cmath>
#include "ising_model_2d.h"

#include <TGraph.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TH2.h>
#include <TLegend.h>
#include <TF1.h>
#include <TGraphErrors.h>

using namespace std;


//-------------------------------------------------------------------------//

vector<double> ising_model_2d_plot(double T, double H, int t_max, int t_stb, bool isplot) {

	TString str_adj = Form("T%.3d_H%.3d", int(100*T), int(100*H));
	TString str_tmp;
	int N = 50;

	IsingModel2D im2d(N, T, H, 1, 1);

	vector<double> v_t;
	vector<double> v_m;
	vector<double> v_E;

	for(int t=0; t<=t_max; t++) {

		im2d.cal_until(t);
		v_t.push_back(1.0*t);
		v_m.push_back(im2d.get_m()/N/N);
		v_E.push_back(im2d.get_E()/N/N);

		if(t%10==0 && t<=50 && isplot) {
			vector< vector<int> > spin = im2d.get_spin();

			vector<double> x1;
			vector<double> y1;
			vector<double> x2;
			vector<double> y2;
			for(int i=0; i<N; i++) {
				for(int j=0; j<N; j++) {
					if(spin[i][j]>0) {
						x1.push_back(0.5+i);
						y1.push_back(0.5+j);
					} else {
						x2.push_back(0.5+i);
						y2.push_back(0.5+j);
					}
				}
			}

			str_tmp = "spin_" + str_adj + Form("_t%.3d", t);
			TCanvas *c_spin = new TCanvas("c_"+str_tmp, "c_"+str_tmp, 600,600);
			TH2D *hframe = new TH2D("frame", "frame", 1,0,N, 1,0,N);
			hframe->Draw();
			TGraph *g_up = new TGraph(x1.size(), &x1[0], &y1[0]);
			g_up->SetMarkerStyle(22);
			g_up->SetMarkerColor(kRed);
			TGraph *g_down = new TGraph(x2.size(), &x2[0], &y2[0]);
			g_down->SetMarkerStyle(23);
			g_down->SetMarkerColor(kBlue);
			g_up->Draw("P same");
			g_down->Draw("P same");
			TLegend *l_spin = new TLegend(0.1,0.9,0.5,0.98);
			l_spin->AddEntry(g_up, "spin up", "p");
			l_spin->AddEntry(g_down, "spin down", "p");
			l_spin->SetFillStyle(0);
			l_spin->SetBorderSize(0);
			l_spin->Draw("same");
			TLegend *l_para = new TLegend(0.5,0.9,0.9,0.98);
			l_para->SetNColumns(2);
			l_para->AddEntry((TObject*)0, Form("t = %d", t), "");
			l_para->AddEntry((TObject*)0, Form("H = %.2f", H), "");
			l_para->AddEntry((TObject*)0, Form("T = %.2f", T), "");
			l_para->SetFillStyle(0);
			l_para->SetBorderSize(0);
			l_para->Draw("same");
			c_spin->SaveAs("./plot/" + str_tmp + ".pdf");

			delete hframe;
			delete g_up;
			delete g_down;
			delete c_spin;
		}
	}

	str_tmp = "t_m_" + str_adj;
	TCanvas *c_t_m = new TCanvas("c_"+str_tmp, "c_"+str_tmp);
	TGraph *g_t_m = new TGraph(v_t.size(), &v_t[0], &v_m[0]);
	g_t_m->GetXaxis()->SetTitle("t");
	g_t_m->GetYaxis()->SetTitle("m");
	g_t_m->Draw("LA");
	TF1 *func_t_m = new TF1("func_"+str_tmp, "[0]", 50,1.0*t_max);
	g_t_m->Fit(func_t_m, "R");
	if(isplot) c_t_m->SaveAs("./plot/" + str_tmp + ".pdf");

	str_tmp = "t_E_" + str_adj;
	TCanvas *c_t_E = new TCanvas("c_"+str_tmp, "c_"+str_tmp);
	TGraph *g_t_E = new TGraph(v_t.size(), &v_t[0], &v_E[0]);
	g_t_E->GetXaxis()->SetTitle("t");
	g_t_E->GetYaxis()->SetTitle("E");
	g_t_E->Draw("LA");
	TF1 *func_t_E = new TF1("func_"+str_tmp, "[0]", 50,1.0*t_max);
	g_t_E->Fit(func_t_E, "R");
	if(isplot) c_t_E->SaveAs("./plot/" + str_tmp + ".pdf");

	str_tmp = "m_E_" + str_adj;
	TCanvas *c_m_E = new TCanvas("c_"+str_tmp, "c_"+str_tmp);
	TGraph *g_m_E = new TGraph(v_m.size(), &v_m[0], &v_E[0]);
	g_m_E->GetXaxis()->SetTitle("m");
	g_m_E->GetYaxis()->SetTitle("E");
	g_m_E->Draw("PLA");
	if(isplot) c_m_E->SaveAs("./plot/" + str_tmp + ".pdf");

	vector<double> v_fit;
	v_fit.push_back(func_t_m->GetParameter(0));
	v_fit.push_back(func_t_m->GetParError(0));
	v_fit.push_back(func_t_E->GetParameter(0));
	v_fit.push_back(func_t_E->GetParError(0));

	return v_fit;
}

//-------------------------------------------------------------------------//

void ising_model_2d_T_scan() {

	const int n = 40;
	double Tmin = 1.5;
	double Tmax = 3.5;
	double T[n];
	double Te[n] = {0};
	double m[n];
	double me[n];
	double ms[n];
	double E[n];
	double Ee[n];
	double Es[n];
	for(int i=0; i<n; i++) {
		T[i] = Tmin + (Tmax-Tmin)*i/n;
		vector<double> tmp = ising_model_2d_plot(T[i], 0.01, 3000, 1000, false);
		m[i]  = tmp[0];
		me[i] = tmp[1];
		ms[i] = tmp[1]*sqrt(2000.0);
		E[i]  = tmp[2];
		Ee[i] = tmp[3];
		Es[i] = tmp[3]*sqrt(2000.0);
	}

	TCanvas *c_T_m = new TCanvas("c_T_m", "c_T_m");
	TGraphErrors *g_T_m = new TGraphErrors(n, T, m, Te, me);
	g_T_m->GetXaxis()->SetTitle("T");
	g_T_m->GetYaxis()->SetTitle("m");
	g_T_m->Draw("PLA");
	c_T_m->SaveAs("./plot/T_m.pdf");

	TCanvas *c_T_E = new TCanvas("c_T_E", "c_T_E");
	TGraphErrors *g_T_E = new TGraphErrors(n, T, E, Te, Ee);
	g_T_E->GetXaxis()->SetTitle("T");
	g_T_E->GetYaxis()->SetTitle("E");
	g_T_E->Draw("PLA");
	c_T_E->SaveAs("./plot/T_E.pdf");

	TCanvas *c_T_ms = new TCanvas("c_T_ms", "c_T_ms");
	TGraph *g_T_ms = new TGraph(n, T, ms);
	g_T_ms->GetXaxis()->SetTitle("T");
	g_T_ms->GetYaxis()->SetTitle("#sigma[m]");
	g_T_ms->Draw("PLA");
	c_T_ms->SaveAs("./plot/T_ms.pdf");

	TCanvas *c_T_Es = new TCanvas("c_T_Es", "c_T_Es");
	TGraph *g_T_Es = new TGraph(n, T, Es);
	g_T_Es->GetXaxis()->SetTitle("T");
	g_T_Es->GetYaxis()->SetTitle("#sigma[E]");
	g_T_Es->Draw("PLA");
	c_T_Es->SaveAs("./plot/T_Es.pdf");
}

//-------------------------------------------------------------------------//

int main() {

	gStyle->SetTitleOffset(1.4,"y");
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);
	gStyle->SetOptFit();

	//ising_model_2d_plot(2.00,  0.0, 100, 50, true);
	//ising_model_2d_plot(2.00, -0.2, 100, 50, true);
	//ising_model_2d_plot(2.00,  0.2, 100, 50, true);
	//ising_model_2d_plot(3.00,    0, 100, 50, true);
	//ising_model_2d_plot(3.00, -0.2, 100, 50, true);
	//ising_model_2d_plot(3.00,  0.2, 100, 50, true);

	ising_model_2d_T_scan();
}

