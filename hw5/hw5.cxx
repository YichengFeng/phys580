#include <iostream>
#include <vector>
#include <cmath>
#include "ising_model_2d.h"
#include "molecular_dynamics_2d.h"

#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TStyle.h>
#include <TPaveStats.h>
#include <TMultiGraph.h>
#include <TH1.h>
#include <TH2.h>
#include <THistPainter.h>

using namespace std;


//-------------------------------------------------------------------------//

vector<double> ising_model_2d_plot(int N, double T, double H, int t_max, int t_stb, bool isplot) {

	TString str_adj = Form("T%.3d_H%.3d", int(round(1000*T)), int(100*H));
	TString str_tmp;

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
		}
	}

	str_tmp = "t_m_" + str_adj;
	TCanvas *c_t_m = new TCanvas("c_"+str_tmp, "c_"+str_tmp);
	TGraph *g_t_m = new TGraph(v_t.size(), &v_t[0], &v_m[0]);
	g_t_m->GetXaxis()->SetTitle("t");
	g_t_m->GetYaxis()->SetTitle("m");
	g_t_m->Draw("LA");
	TF1 *func_t_m = new TF1("func_"+str_tmp, "[0]", 50,1.0*t_max);
	//g_t_m->Fit(func_t_m, "R");
	if(isplot) c_t_m->SaveAs("./plot/" + str_tmp + ".pdf");

	str_tmp = "t_E_" + str_adj;
	TCanvas *c_t_E = new TCanvas("c_"+str_tmp, "c_"+str_tmp);
	TGraph *g_t_E = new TGraph(v_t.size(), &v_t[0], &v_E[0]);
	g_t_E->GetXaxis()->SetTitle("t");
	g_t_E->GetYaxis()->SetTitle("E");
	g_t_E->Draw("LA");
	TF1 *func_t_E = new TF1("func_"+str_tmp, "[0]", 50,1.0*t_max);
	//g_t_E->Fit(func_t_E, "R");
	if(isplot) c_t_E->SaveAs("./plot/" + str_tmp + ".pdf");

	str_tmp = "m_E_" + str_adj;
	TCanvas *c_m_E = new TCanvas("c_"+str_tmp, "c_"+str_tmp);
	TGraph *g_m_E = new TGraph(v_m.size(), &v_m[0], &v_E[0]);
	g_m_E->GetXaxis()->SetTitle("m");
	g_m_E->GetYaxis()->SetTitle("E");
	g_m_E->Draw("PLA");
	if(isplot) c_m_E->SaveAs("./plot/" + str_tmp + ".pdf");

	double Em = 0, E2m = 0;
	double Mm = 0, M2m = 0;
	for(int t=t_stb; t<t_max; t++) {
		Em  += v_E[t];
		E2m += v_E[t]*v_E[t];
		Mm  += v_m[t];
		M2m += v_m[t]*v_m[t];
	}
	Em  /= (t_max-t_stb);
	E2m /= (t_max-t_stb);
	Mm  /= (t_max-t_stb);
	M2m /= (t_max-t_stb);

	vector<double> v_fit;
	//v_fit.push_back(func_t_m->GetParameter(0));
	//v_fit.push_back(func_t_m->GetParError(0));
	//v_fit.push_back(func_t_E->GetParameter(0));
	//v_fit.push_back(func_t_E->GetParError(0));
	v_fit.push_back(Mm);
	v_fit.push_back(sqrt(M2m-Mm*Mm)/sqrt(t_max-t_stb));
	v_fit.push_back(Em);
	v_fit.push_back(sqrt(E2m-Em*Em)/sqrt(t_max-t_stb));

	return v_fit;
}

//-------------------------------------------------------------------------//

void ising_model_2d_beta() {

	const int n = 80;
	const double Tmax = 2.8;
	const double Tmin = 2.0;

	double T[n];
	double Te[n] = {0};
	double M[n];
	double Me[n];
	double Ms[n];
	double E[n];
	double Ee[n];
	double Es[n];
	double C[n];
	double Chi[n];

	int t_max = 4000;
	int t_stb = 1000;
	int N = 50;

	for(int i=0; i<n; i++) {
		double Ttmp = Tmin + (Tmax - Tmin)/n*i;
		T[i] = Ttmp;
		vector<double> v_tmp = ising_model_2d_plot(N, Ttmp, 0.01, t_max, t_stb, false);
		M[i]  = v_tmp[0];
		Me[i] = v_tmp[1];
		Ms[i] = v_tmp[1]*sqrt(t_max-t_stb);
		E[i]  = v_tmp[2];
		Ee[i] = v_tmp[3];
		Es[i] = v_tmp[3]*sqrt(t_max-t_stb);
		C[i]  = N*N*Es[i]*Es[i]/Ttmp/Ttmp;
		Chi[i] = N*N*Ms[i]*Ms[i]/Ttmp;
	}

	TCanvas *c_T_M = new TCanvas("c_T_M", "c_T_M");
	TGraphErrors *g_T_M = new TGraphErrors(n, T, M, Te, Me);
	g_T_M->GetXaxis()->SetTitle("T");
	g_T_M->GetYaxis()->SetTitle("M");
	g_T_M->Draw("PLA");
	c_T_M->SaveAs("./plot/T_M.pdf");

	TCanvas *c_T_E = new TCanvas("c_T_E", "c_T_E");
	TGraphErrors *g_T_E = new TGraphErrors(n, T, E, Te, Ee);
	g_T_E->GetXaxis()->SetTitle("T");
	g_T_E->GetYaxis()->SetTitle("E");
	g_T_E->Draw("PLA");
	c_T_E->SaveAs("./plot/T_E.pdf");

	TCanvas *c_T_Ms = new TCanvas("c_T_Ms", "c_T_Ms");
	TGraph *g_T_Ms = new TGraph(n, T, Ms);
	g_T_Ms->GetXaxis()->SetTitle("T");
	g_T_Ms->GetYaxis()->SetTitle("#sigma[M]");
	g_T_Ms->Draw("PLA");
	c_T_Ms->SaveAs("./plot/T_Ms.pdf");

	TCanvas *c_T_Es = new TCanvas("c_T_Es", "c_T_Es");
	TGraph *g_T_Es = new TGraph(n, T, Es);
	g_T_Es->GetXaxis()->SetTitle("T");
	g_T_Es->GetYaxis()->SetTitle("#sigma[E]");
	g_T_Es->Draw("PLA");
	c_T_Es->SaveAs("./plot/T_Es.pdf");

	TCanvas *c_T_C = new TCanvas("c_T_C", "c_T_C");
	TGraph *g_T_C = new TGraph(n, T, C);
	g_T_C->GetXaxis()->SetTitle("T");
	g_T_C->GetYaxis()->SetTitle("C");
	g_T_C->Draw("PLA");
	c_T_C->SaveAs("./plot/T_C.pdf");

	TCanvas *c_T_Chi = new TCanvas("c_T_Chi", "c_T_Chi");
	TGraph *g_T_Chi = new TGraph(n, T, Chi);
	g_T_Chi->GetXaxis()->SetTitle("T");
	g_T_Chi->GetYaxis()->SetTitle("#chi");
	g_T_Chi->Draw("PLA");
	c_T_Chi->SaveAs("./plot/T_Chi.pdf");

	int imax = 0;
	double Cmax = 0;
	for(int i=0; i<n; i++) {
		if(Cmax >= C[i]) continue;
		Cmax = C[i];
		imax = i;
	}
	double Tc = T[imax];
	cout << "Tc = " << Tc << endl;

	double logTdif[n] = {0};
	double logM[n] = {0};
	double logC[n] = {0};
	double logChi[n] = {0};
	for(int i=0; i<n; i++) {
		if(T[i] == Tc) continue;
		logTdif[i] = log(fabs(T[i] - Tc));
		logM[i] = log(M[i]);
		logC[i] = log(C[i]);
		logChi[i] = log(Chi[i]);
	}

	TCanvas *c_log_T_M = new TCanvas("c_log_T_M", "c_log_T_M");
	TGraph *g_log_T_M = new TGraph(imax, logTdif, logM);
	g_log_T_M->GetXaxis()->SetTitle("ln|T-T_{c}|");
	g_log_T_M->GetYaxis()->SetTitle("ln(M)");
	g_log_T_M->SetMarkerStyle(20);
	g_log_T_M->Draw("PA");
	TF1 *func_log_T_M = new TF1("func_log_T_M", "[0] + [1]*x", -2,0);
	//g_log_T_M->Fit(func_log_T_M, "R");
	g_log_T_M->Fit(func_log_T_M);
	gPad->Update();
	TPaveStats *st_log_T_M = (TPaveStats*)g_log_T_M->FindObject("stats");
	st_log_T_M->SetX1NDC(0.1);
	st_log_T_M->SetX2NDC(0.4);
	st_log_T_M->SetY1NDC(0.75);
	st_log_T_M->SetY2NDC(0.9);
	c_log_T_M->SaveAs("./plot/log_T_M.pdf");

	TCanvas *c_log_T_C = new TCanvas("c_log_T_C", "c_log_T_C");
	TGraph *g_log_T_C = new TGraph(imax, logTdif, logC);
	g_log_T_C->GetXaxis()->SetTitle("ln|T-T_{c}|");
	g_log_T_C->GetYaxis()->SetTitle("ln(C)");
	g_log_T_C->SetMarkerStyle(20);
	g_log_T_C->Draw("PA");
	TF1 *func_log_T_C = new TF1("func_log_T_C", "[0] + [1]*x", -2,0);
	//g_log_T_C->Fit(func_log_T_C, "R");
	g_log_T_C->Fit(func_log_T_C);
	gPad->Update();
	TPaveStats *st_log_T_C = (TPaveStats*)g_log_T_C->FindObject("stats");
	st_log_T_C->SetX1NDC(0.6);
	st_log_T_C->SetX2NDC(0.9);
	st_log_T_C->SetY1NDC(0.75);
	st_log_T_C->SetY2NDC(0.9);
	c_log_T_C->SaveAs("./plot/log_T_C.pdf");

	TCanvas *c_log_T_Chi = new TCanvas("c_log_T_Chi", "c_log_T_Chi");
	TGraph *g_log_T_Chi = new TGraph(imax, logTdif, logChi);
	g_log_T_Chi->GetXaxis()->SetTitle("ln|T-T_{c}|");
	g_log_T_Chi->GetYaxis()->SetTitle("ln(#chi)");
	g_log_T_Chi->SetMarkerStyle(20);
	g_log_T_Chi->Draw("PA");
	TF1 *func_log_T_Chi = new TF1("func_log_T_Chi", "[0] + [1]*x", -2,0);
	//g_log_T_Chi->Fit(func_log_T_Chi, "R");
	g_log_T_Chi->Fit(func_log_T_Chi);
	gPad->Update();
	TPaveStats *st_log_T_Chi = (TPaveStats*)g_log_T_Chi->FindObject("stats");
	st_log_T_Chi->SetX1NDC(0.6);
	st_log_T_Chi->SetX2NDC(0.9);
	st_log_T_Chi->SetY1NDC(0.75);
	st_log_T_Chi->SetY2NDC(0.9);
	c_log_T_Chi->SaveAs("./plot/log_T_Chi.pdf");
}

//-------------------------------------------------------------------------//

void ising_model_2d_specific_heat() {

	const int n = 200;
	const double Tmax = 4;
	const double Tmin = 1;
	const double dT = (Tmax - Tmin)/n;

	double T[n];
	double Te[n] = {0};
	double M[n];
	double Me[n];
	double Ms[n];
	double E[n];
	double Ee[n];
	double Es[n];
	double C_fluctuation[n];
	double C_differentiation[n];

	int t_max = 3000;
	int t_stb = 1000;
	int N = 10;
	double DeltaE = 0;

	for(int i=0; i<n; i++) {
		double Ttmp = Tmin + dT*i;
		T[i] = Ttmp;
		vector<double> v_tmp = ising_model_2d_plot(N, Ttmp, 0.00, t_max, t_stb, false);
		M[i]  = v_tmp[0];
		Me[i] = v_tmp[1];
		Ms[i] = v_tmp[1]*sqrt(t_max-t_stb);
		E[i]  = v_tmp[2];
		Ee[i] = v_tmp[3];
		Es[i] = v_tmp[3]*sqrt(t_max-t_stb);
		C_fluctuation[i]  = N*N*Es[i]*Es[i]/Ttmp/Ttmp;
	}

	for(int i=0; i<n; i++) {
		if(i==0)   C_differentiation[i] = (E[i+1]-E[i])/dT;
		if(i==n-1) C_differentiation[i] = (E[i]-E[i-1])/dT;
		C_differentiation[i] = (E[i+1]-E[i-1])/dT/2;
		if(!isnan(C_fluctuation[i])) {
			DeltaE += C_fluctuation[i];
		} else {
			cout << i << endl;
		}
	}
	DeltaE *= dT;

	cout << "------- DeltaE " << DeltaE << ", " << E[n-1]-E[0] << endl;

	TString str_tmp;
	TString str_adj = "";

	str_tmp = "T_C_fluctuation" + str_adj;
	TCanvas *c_T_C_fluctuation = new TCanvas(str_tmp, str_tmp);
	TGraph *g_T_C_fluctuation = new TGraph(n, &T[0], &C_fluctuation[0]);
	g_T_C_fluctuation->GetXaxis()->SetTitle("T");
	g_T_C_fluctuation->GetYaxis()->SetTitle("C from E fluctuation");
	g_T_C_fluctuation->GetYaxis()->SetRangeUser(-0.5,2.5);
	g_T_C_fluctuation->Draw("PLA");
	c_T_C_fluctuation->SaveAs("./plot/"+str_tmp+".pdf");

	str_tmp = "T_C_differentiation" + str_adj;
	TCanvas *c_T_C_differentiation = new TCanvas(str_tmp, str_tmp);
	TGraph *g_T_C_differentiation = new TGraph(n-2, &T[1], &C_differentiation[1]);
	g_T_C_differentiation->GetXaxis()->SetTitle("T");
	g_T_C_differentiation->GetYaxis()->SetTitle("C from E differentiation");
	g_T_C_differentiation->GetYaxis()->SetRangeUser(-0.5,2.5);
	g_T_C_differentiation->Draw("PLA");
	c_T_C_differentiation->SaveAs("./plot/"+str_tmp+".pdf");
}


//-------------------------------------------------------------------------//

void ising_model_2d_scaling() {

	const int n = 100;
	const double Tmax = 3;
	const double Tmin = 1.5;
	const double dT = (Tmax - Tmin)/n;
	const int nH = 5;
	const double dH = 0.01;
	double H[nH];
	Int_t color[nH] = {kBlack, kBlue, kRed, kGreen, kOrange};

	double T[nH][n];
	double Te[nH][n] = {0};
	double M[nH][n];
	double Me[nH][n];
	double Ms[nH][n];
	double E[nH][n];
	double Ee[nH][n];
	double Es[nH][n];
	double C[nH][n];
	double Chi[nH][n];
	int    iTc[nH] = {0};
	double Tc[nH];

	int t_max = 3000;
	int t_stb = 1000;
	int N = 20;

	for(int iH=0; iH<nH; iH++) {
		H[iH] = 0.01*(iH+1);
		for(int i=0; i<n; i++) {
			double Ttmp = Tmin + dT*i;
			T[iH][i] = Ttmp;
			vector<double> v_tmp = ising_model_2d_plot(N, Ttmp, H[iH], t_max, t_stb, false);
			M[iH][i]  = v_tmp[0];
			Me[iH][i] = v_tmp[1];
			Ms[iH][i] = v_tmp[1]*sqrt(t_max-t_stb);
			E[iH][i]  = v_tmp[2];
			Ee[iH][i] = v_tmp[3];
			Es[iH][i] = v_tmp[3]*sqrt(t_max-t_stb);
			C[iH][i]  = N*N*Es[iH][i]*Es[iH][i]/Ttmp/Ttmp;
			Chi[iH][i] = Ms[iH][i]*Ms[iH][i]*N*N/Ttmp;
		}
		int imax = 0;
		double Cmax = 0;
		for(int i=1; i<n; i++) {
			if(Cmax >= C[iH][i]) continue;
			Cmax = C[iH][i];
			imax = i;
		}
		iTc[iH] = imax;
		Tc[iH] = T[iH][imax];
	}

	vector<double> T_scaled[nH];
	vector<double> M_scaled[nH];
	vector<double> H_scaled[nH];
	vector<double> Chi_scaled[nH];
	for(int iH=0; iH<nH; iH++) {
		for(int i=0; i<n; i++) {
			if(T[iH][i] == Tc[iH]) continue;
			if(M[iH][i]  <0) continue;
			if(Chi[iH][i]<0) continue;
			double Ttmp = fabs(T[iH][i] - Tc[iH])/Tc[iH];
			T_scaled[iH].push_back(Ttmp);
			M_scaled[iH].push_back(M[iH][i]*pow(Ttmp,-0.125));
			H_scaled[iH].push_back(H[iH]*pow(Ttmp, -1.875));
			Chi_scaled[iH].push_back(Chi[iH][i]*pow(Ttmp, 1.75));
		}
	}

	TString str_tmp;
	TString str_adj = "";

	str_tmp = "T_M_variousH" + str_adj;
	TCanvas *c_T_M_variousH = new TCanvas("c_"+str_tmp, "c_"+str_tmp);
	TLegend *l_T_M_variousH = new TLegend(0.7,0.6,0.9,0.9);
	l_T_M_variousH->SetFillStyle(0);
	l_T_M_variousH->SetBorderSize(0);
	TMultiGraph *mg_T_M_variousH = new TMultiGraph();
	TGraph *g_T_M_variousH[nH];
	for(int iH=0; iH<nH; iH++) {
		g_T_M_variousH[iH] = new TGraph(n, T[iH], M[iH]);
		g_T_M_variousH[iH]->SetLineColor(color[iH]);
		g_T_M_variousH[iH]->SetMarkerColor(color[iH]);
		g_T_M_variousH[iH]->SetMarkerStyle(20);
		g_T_M_variousH[iH]->SetMarkerSize(0.8);
		l_T_M_variousH->AddEntry(g_T_M_variousH[iH], Form("H = %.2f", H[iH]), "lp");
		mg_T_M_variousH->Add(g_T_M_variousH[iH]);
	}
	mg_T_M_variousH->GetXaxis()->SetTitle("T");
	mg_T_M_variousH->GetYaxis()->SetTitle("m");
	mg_T_M_variousH->Draw("PA");
	l_T_M_variousH->Draw("same");
	c_T_M_variousH->SaveAs("./plot/"+str_tmp+".pdf");

	str_tmp = "H_M_scaled" + str_adj;
	TCanvas *c_H_M_scaled = new TCanvas("c_"+str_tmp, "c_"+str_tmp);
	TLegend *l_H_M_scaled = new TLegend(0.7,0.1,0.9,0.4);
	l_H_M_scaled->SetFillStyle(0);
	l_H_M_scaled->SetBorderSize(0);
	TMultiGraph *mg_H_M_scaled = new TMultiGraph();
	TGraph *g_H_M_scaled[nH];
	for(int iH=0; iH<nH; iH++) {
		g_H_M_scaled[iH] = new TGraph(H_scaled[iH].size(), &H_scaled[iH][0], &M_scaled[iH][0]);
		g_H_M_scaled[iH]->SetLineColor(color[iH]);
		g_H_M_scaled[iH]->SetMarkerColor(color[iH]);
		g_H_M_scaled[iH]->SetMarkerStyle(20);
		g_H_M_scaled[iH]->SetMarkerSize(0.8);
		l_H_M_scaled->AddEntry(g_H_M_scaled[iH], Form("H = %.2f", H[iH]), "lp");
		mg_H_M_scaled->Add(g_H_M_scaled[iH]);
	}
	mg_H_M_scaled->GetXaxis()->SetTitle("h/|t|^{15/8}");
	mg_H_M_scaled->GetYaxis()->SetTitle("m/|t|^{1/8}");
	mg_H_M_scaled->Draw("PA");
	l_H_M_scaled->Draw("same");
	c_H_M_scaled->SetLogx();
	c_H_M_scaled->SetLogy();
	c_H_M_scaled->SaveAs("./plot/"+str_tmp+".pdf");

	str_tmp = "H_Chi_scaled" + str_adj;
	TCanvas *c_H_Chi_scaled = new TCanvas("c_"+str_tmp, "c_"+str_tmp);
	TMultiGraph *mg_H_Chi_scaled = new TMultiGraph();
	TLegend *l_H_Chi_scaled = new TLegend(0.7,0.6,0.9,0.9);
	l_H_Chi_scaled->SetFillStyle(0);
	l_H_Chi_scaled->SetBorderSize(0);
	TGraph *g_H_Chi_scaled[nH];
	for(int iH=0; iH<nH; iH++) {
		g_H_Chi_scaled[iH] = new TGraph(H_scaled[iH].size(), &H_scaled[iH][0], &Chi_scaled[iH][0]);
		g_H_Chi_scaled[iH]->SetLineColor(color[iH]);
		g_H_Chi_scaled[iH]->SetMarkerColor(color[iH]);
		g_H_Chi_scaled[iH]->SetMarkerStyle(20);
		g_H_Chi_scaled[iH]->SetMarkerSize(0.8);
		l_H_Chi_scaled->AddEntry(g_H_Chi_scaled[iH], Form("H = %.2f", H[iH]), "lp");
		mg_H_Chi_scaled->Add(g_H_Chi_scaled[iH]);
	}
	mg_H_Chi_scaled->GetXaxis()->SetTitle("#chi|t|^{7/4}");
	mg_H_Chi_scaled->GetYaxis()->SetTitle("m/|t|^{1/8}");
	mg_H_Chi_scaled->Draw("PA");
	l_H_Chi_scaled->Draw("same");
	c_H_Chi_scaled->SetLogx();
	c_H_Chi_scaled->SetLogy();
	c_H_Chi_scaled->SaveAs("./plot/"+str_tmp+".pdf");
}

//-------------------------------------------------------------------------//

void molecular_dynamics_2d_plot() {

	TString str_adj;
	TString str_tmp;

	int N =20;
	double L = 10;

	MolecularDynamics2D md2d(N, L, 0.005, 0.05, 0);
	//MolecularDynamics2D md2d(N, L, 0.02, 0.1, 0);
	//md2d.set_rec(1);

	const int n_ends = 4;
	double t_ends[n_ends] = {0.0,   20.0, 40.0, 60.0};
	bool  to_draw[n_ends] = {false, true, true, true}; 

	TMultiGraph *mg_E = new TMultiGraph();
	TMultiGraph *mg_T = new TMultiGraph();
	TMultiGraph *mg_tag_displace = new TMultiGraph();
	TMultiGraph *mg_tag_distance = new TMultiGraph();
	TH1D *hvx[n_ends-1];
	TH1D *hvy[n_ends-1];
	TH1D *hv[n_ends-1];
	str_adj = Form("t%.3d%.3d", int(t_ends[0]*10), int(t_ends[n_ends-1]*10));
	str_tmp = str_adj + "_v";
	TH1D *hsumv = new TH1D("hsum_"+str_tmp, "hsum_"+str_tmp, 40,0,8);

	for(int i_end=1; i_end<n_ends; i_end++) {

		md2d.cal_until(t_ends[i_end]);

		vector<double> t = md2d.get_t();
		vector<double> E = md2d.get_E();
		vector<double> T = md2d.get_T();
		vector<double> tag_displace = md2d.get_tag_displace();
		vector<double> tag_distance = md2d.get_tag_distance();
		vector< vector<double> > x = md2d.get_x();
		vector< vector<double> > y = md2d.get_y();
		vector< vector<double> > vx = md2d.get_vx();
		vector< vector<double> > vy = md2d.get_vy();
		
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

		str_tmp = str_adj + "_vx";
		hvx[i_end-1] = new TH1D("h_"+str_tmp, "h_"+str_tmp, 40,-5,5);
		str_tmp = str_adj + "_vy";
		hvy[i_end-1] = new TH1D("h_"+str_tmp, "h_"+str_tmp, 40,-5,5);
		str_tmp = str_adj + "_v";
		hv[i_end-1]  = new TH1D("h_"+str_tmp, "h_"+str_tmp, 40,0,8);
		for(int i=0; i<vx.size(); i++) {
			for(int j=0; j<vx[i].size(); j++) {
				hvx[i_end-1]->Fill(vx[i][j]);
				hvy[i_end-1]->Fill(vy[i][j]);
				hv[i_end-1] ->Fill(sqrt(vx[i][j]*vx[i][j]+vy[i][j]*vy[i][j]));
				hsumv->Fill(sqrt(vx[i][j]*vx[i][j]+vy[i][j]*vy[i][j]));
			}
		}
		hvx[i_end-1]->Scale(4.0/hvx[i_end-1]->GetEntries());
		hvy[i_end-1]->Scale(4.0/hvx[i_end-1]->GetEntries());
		hv[i_end-1] ->Scale(5.0/hv[i_end-1] ->GetEntries());
	}
	hsumv->Scale(5.0/hsumv->GetEntries());

	str_adj = Form("t%.3d%.3d", int(t_ends[0]*10), int(t_ends[n_ends-1]*10));

	str_tmp = str_adj + "_vx";
	TCanvas *c_vx = new TCanvas("c_"+str_tmp, "c_"+str_tmp);
	TLegend *l_vx = new TLegend(0.1,0.7,0.3,0.9);
	for(int i=0; i<n_ends-1; i++) {
		hvx[i]->SetMarkerStyle(24+i);
		hvx[i]->SetMarkerColor(2+i);
		hvx[i]->SetLineColor(2+i);
		l_vx->AddEntry(hvx[i], Form("t = %.1f ~ %.1f", i*20.0, (i+1)*20.0), "lp");
		hvx[i]->GetXaxis()->SetTitle("v_{x}");
		hvx[i]->GetYaxis()->SetTitle("P(v_{x})");
		hvx[i]->Draw("same L");
	}
	l_vx->Draw("same");
	c_vx->SaveAs("./plot/"+str_tmp+".pdf");

	str_tmp = str_adj + "_vy";
	TCanvas *c_vy = new TCanvas("c_"+str_tmp, "c_"+str_tmp);
	TLegend *l_vy = new TLegend(0.1,0.7,0.3,0.9);
	for(int i=0; i<n_ends-1; i++) {
		hvy[i]->SetMarkerStyle(24+i);
		hvy[i]->SetMarkerColor(2+i);
		hvy[i]->SetLineColor(2+i);
		l_vy->AddEntry(hvy[i], Form("t = %.1f ~ %.1f", i*20.0, (i+1)*20.0), "lp");
		hvy[i]->GetXaxis()->SetTitle("v_{y}");
		hvy[i]->GetYaxis()->SetTitle("P(v_{y})");
		hvy[i]->Draw("L same");
	}
	l_vy->Draw("same");
	c_vy->SaveAs("./plot/"+str_tmp+".pdf");

	str_tmp = str_adj + "_v";
	TCanvas *c_v = new TCanvas("c_"+str_tmp, "c_"+str_tmp);
	TLegend *l_v = new TLegend(0.65,0.6,0.9,0.8);
	l_v->SetFillStyle(0);
	l_v->SetBorderSize(0);
	TF1 *func_v = new TF1(Form("func%d",0), "[0]*x*exp(-[0]*0.5*x*x)", 0,8);
	func_v->SetLineColor(1);
	func_v->SetParameter(0,1);
	hsumv->GetXaxis()->SetTitle("v");
	hsumv->GetYaxis()->SetTitle("P(v)");
	hsumv->Draw("");
	hsumv->Fit(func_v);
	l_v->AddEntry(hsumv, "t = 0.0 ~ 60.0", "lp");
	for(int i=0; i<n_ends-1; i++) {
		hv[i]->SetMarkerStyle(24+i);
		hv[i]->SetMarkerColor(2+i);
		hv[i]->SetLineColor(2+i);
		l_v->AddEntry(hv[i], Form("t = %.1f ~ %.1f", i*20.0, (i+1)*20.0), "lp");
		hv[i]->GetXaxis()->SetTitle("v");
		hv[i]->GetYaxis()->SetTitle("P(v)");
		hv[i]->Draw("same L");
	}
	l_v->Draw("same");
	c_v->SaveAs("./plot/"+str_tmp+".pdf");

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

	//ising_model_2d_beta();
	//ising_model_2d_specific_heat();
	//ising_model_2d_scaling();

	molecular_dynamics_2d_plot();

	return 0;
}

//-------------------------------------------------------------------------//

