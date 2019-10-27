#include <iostream>
#include <vector>
#include <cmath>
#include "percolation_2d.h"

#include <TCanvas.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TH2.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TF1.h>

using namespace std;


//-------------------------------------------------------------------------//

vector<double> percolation_2d_plot(double p, int L, int trial, bool tosave = true) {

	TString str_tmp;
	TString str_adj;
	str_adj = Form("L%d",L);

	Percolation2D pc2d(p, L, trial);
	pc2d.cal_trials();

	cout << endl;
	cout << "P(p) = " <<  pc2d.get_P_average() << endl;
	cout << "S(p) = " <<  pc2d.get_S_average() << endl;
	cout << "cluster number: " << pc2d.get_n_cluster() << endl;
	cout << "occupied sites: " << pc2d.get_n_occupied() << endl;

	vector<double> vPS;
	vPS.push_back(pc2d.get_P_average());
	vPS.push_back(pc2d.get_S_average());

	vector< vector<bool> > occupied = pc2d.get_occupied();
	vector<double> occupied_x;
	vector<double> occupied_y;
	for(int iy=0; iy<occupied.size(); iy++) {
		for(int ix=0; ix<occupied[iy].size(); ix++) {
			if(occupied[iy][ix]){
				occupied_x.push_back(ix+0.5);
				occupied_y.push_back(iy+0.5);
			}
		}
	}

	vector< vector<int> > cluster_x = pc2d.get_cluster_x();
	vector< vector<int> > cluster_y = pc2d.get_cluster_y();
	vector<int> scid = pc2d.get_spanning_cluster_id();
	vector<double> scx;
	vector<double> scy;
	for(int i=0; i<cluster_x[scid[0]].size(); i++) {
		scx.push_back(cluster_x[scid[0]][i]+0.5);
		scy.push_back(cluster_y[scid[0]][i]+0.5);
	}

	str_tmp = "occupied_" + str_adj;
	TCanvas *c_occupied = new TCanvas("c_"+str_tmp, "c_"+str_tmp, 600,600);
	TH2D *hframe = new TH2D("frame", "frame", 1,0,L, 1,0,L);
	hframe->GetXaxis()->SetTickSize(0);
	hframe->GetYaxis()->SetTickSize(0);
	hframe->Draw();
	TGraph *g_occ = new TGraph(occupied_x.size(), &occupied_x[0], &occupied_y[0]);
	g_occ->SetMarkerStyle(21);
	g_occ->SetMarkerSize(600*0.8/L/8);
	g_occ->Draw("P");
	TGraph *g_spn = new TGraph(scx.size(), &scx[0], &scy[0]);
	g_spn->SetMarkerStyle(21);
	g_spn->SetMarkerColor(kOrange);
	g_spn->SetMarkerSize(600*0.8/L/8);
	g_spn->Draw("P");
	if(tosave) c_occupied->SaveAs("./plot/"+str_tmp+".pdf");

	delete hframe;

	return vPS;
}

//-------------------------------------------------------------------------//

void percolation_2d_Lscan_plot() {

	TString str_tmp;
	TString str_adj = "";

	const int nL = 7;
	double Ls[nL] = {50, 100, 200, 400, 600, 800, 1000};
	double Ps[nL];
	double Ss[nL];
	double logLs[nL];
	double logPs[nL];
	double logSs[nL];
	for(int i=0; i<nL; i++) {
		vector<double> tmpv = percolation_2d_plot(0.593, (int)Ls[i], 50);
		Ps[i] = tmpv[0];
		Ss[i] = tmpv[1];
		logLs[i] = log(1.0*Ls[i]);
		logPs[i] = log(1.0*Ps[i]);
		logSs[i] = log(1.0*Ss[i]);
	}

	str_tmp = "L_P" + str_adj;
	TCanvas *c_L_P = new TCanvas(str_tmp, str_tmp);
	TGraph *g_L_P = new TGraph(nL, Ls, Ps);
	g_L_P->GetXaxis()->SetTitle("L");
	g_L_P->GetYaxis()->SetTitle("P(p_{c})");
	g_L_P->SetMarkerStyle(20);
	g_L_P->Draw("PLA");
	c_L_P->SaveAs("./plot/"+str_tmp+".pdf");

	str_tmp = "logL_logP" + str_adj;
	TCanvas *c_logL_logP = new TCanvas(str_tmp, str_tmp);
	TGraph *g_logL_logP = new TGraph(nL, logLs, logPs);
	g_logL_logP->GetXaxis()->SetTitle("ln(L)");
	g_logL_logP->GetYaxis()->SetTitle("ln(P(p_{c}))");
	g_logL_logP->SetMarkerStyle(20);
	g_logL_logP->Draw("PLA");
	TF1 *func_logL_logP = new TF1("func_logL_logP", "[0]+[1]*x");
	func_logL_logP->SetLineStyle(7);
	func_logL_logP->SetLineColor(kRed);
	g_logL_logP->Fit(func_logL_logP);
	TLegend *l_logL_logP = new TLegend(0.1,0.9,0.9,0.96);
	l_logL_logP->SetFillStyle(0);
	l_logL_logP->SetBorderSize(0);
	l_logL_logP->SetNColumns(2);
	l_logL_logP->AddEntry((TObject*)0, Form("slope = %.5f",func_logL_logP->GetParameter(1)), "");
	l_logL_logP->AddEntry((TObject*)0, Form("intercept = %.5f",func_logL_logP->GetParameter(0)), "");
	l_logL_logP->Draw("same");
	c_logL_logP->SaveAs("./plot/"+str_tmp+".pdf");

	str_tmp = "L_S" + str_adj;
	TCanvas *c_L_S = new TCanvas(str_tmp, str_tmp);
	TGraph *g_L_S = new TGraph(nL, Ls, Ss);
	g_L_S->GetXaxis()->SetTitle("L");
	g_L_S->GetYaxis()->SetTitle("S(p_{c})");
	g_L_S->SetMarkerStyle(20);
	g_L_S->Draw("PLA");
	c_L_S->SaveAs("./plot/"+str_tmp+".pdf");

	str_tmp = "logL_logS" + str_adj;
	TCanvas *c_logL_logS = new TCanvas(str_tmp, str_tmp);
	TGraph *g_logL_logS = new TGraph(nL, logLs, logSs);
	g_logL_logS->GetXaxis()->SetTitle("ln(L)");
	g_logL_logS->GetYaxis()->SetTitle("ln(S(p_{c}))");
	g_logL_logS->SetMarkerStyle(20);
	g_logL_logS->Draw("PLA");
	TF1 *func_logL_logS = new TF1("func_logL_logS", "[0]+[1]*x");
	func_logL_logS->SetLineStyle(7);
	func_logL_logS->SetLineColor(kRed);
	g_logL_logS->Fit(func_logL_logS);
	TLegend *l_logL_logS = new TLegend(0.1,0.9,0.9,0.96);
	l_logL_logS->SetNColumns(2);
	l_logL_logS->AddEntry((TObject*)0, Form("slope = %.5f",func_logL_logS->GetParameter(1)), "");
	l_logL_logS->AddEntry((TObject*)0, Form("intercept = %.5f",func_logL_logS->GetParameter(0)), "");
	l_logL_logS->Draw("same");
	l_logL_logS->SetFillStyle(0);
	l_logL_logS->SetBorderSize(0);
	c_logL_logS->SaveAs("./plot/"+str_tmp+".pdf");

}

//-------------------------------------------------------------------------//

void percolation_2d_pscan_plot(int range) {

	TString str_tmp;
	TString str_adj;

	const int np = 10;
	double *ps;
	//double pslow[np] = {0.338, 0.466, 0.53, 0.562, 0.578, 0.586, 0.590, 0.592};
	//double pshigh[np] = {0.6, 0.602, 0.605, 0.610, 0.624, 0.656, 0.72, 0.848};
	double pslow[np];
	double pshigh[np];
	double pl = 0.23;
	double ph = 0.95;
	double pc = 0.593;
	double dph = 0.003;
	double dpl = 0.10;
	for(int i=0; i<np; i++) {
		pslow[i]  = pc - dpl*pow((pc-pl)/dpl, 1.0*i/np);
		pshigh[i] = pc + dph*pow((ph-pc)/dph, 1.0*i/np);
	} 
	if(range == 0) {
		str_adj = "low";
		ps = pslow;
	} else {
		str_adj = "high";
		ps = pshigh;
	}
	double Ps[np];
	double Ss[np];
	double logps[np];
	double logPs[np];
	double logSs[np];
	for(int i=0; i<np; i++) {
		vector<double> tmpv = percolation_2d_plot(ps[i], 100, 200, false);
		Ps[i] = tmpv[0];
		Ss[i] = tmpv[1];
		logps[i] = log(1.0*fabs(ps[i]-0.593));
		logPs[i] = log(1.0*Ps[i]);
		logSs[i] = log(1.0*Ss[i]);
	}

	str_tmp = "p_P_" + str_adj;
	TCanvas *c_p_P = new TCanvas(str_tmp, str_tmp);
	TGraph *g_p_P = new TGraph(np, ps, Ps);
	g_p_P->GetXaxis()->SetTitle("p");
	g_p_P->GetYaxis()->SetTitle("P(p)");
	g_p_P->SetMarkerStyle(20);
	g_p_P->Draw("PLA");
	c_p_P->SaveAs("./plot/"+str_tmp+".pdf");

	str_tmp = "logp_logP_" + str_adj;
	TCanvas *c_logp_logP = new TCanvas(str_tmp, str_tmp);
	TGraph *g_logp_logP = new TGraph(np, logps, logPs);
	g_logp_logP->GetXaxis()->SetTitle("ln|p-p_{c}|");
	g_logp_logP->GetYaxis()->SetTitle("ln(P(p))");
	g_logp_logP->SetMarkerStyle(20);
	g_logp_logP->Draw("PLA");
	TF1 *func_logp_logP = new TF1("func_logp_logP", "[0]+[1]*x");
	func_logp_logP->SetLineStyle(7);
	func_logp_logP->SetLineColor(kRed);
	g_logp_logP->Fit(func_logp_logP);
	TLegend *l_logp_logP = new TLegend(0.1,0.9,0.9,0.96);
	l_logp_logP->SetFillStyle(0);
	l_logp_logP->SetBorderSize(0);
	l_logp_logP->SetNColumns(2);
	l_logp_logP->AddEntry((TObject*)0, Form("slope = %.5f",func_logp_logP->GetParameter(1)), "");
	l_logp_logP->AddEntry((TObject*)0, Form("intercept = %.5f",func_logp_logP->GetParameter(0)), "");
	l_logp_logP->Draw("same");
	c_logp_logP->SaveAs("./plot/"+str_tmp+".pdf");

	str_tmp = "p_S_" + str_adj;
	TCanvas *c_p_S = new TCanvas(str_tmp, str_tmp);
	TGraph *g_p_S = new TGraph(np, ps, Ss);
	g_p_S->GetXaxis()->SetTitle("p");
	g_p_S->GetYaxis()->SetTitle("S(p)");
	g_p_S->SetMarkerStyle(20);
	g_p_S->Draw("PLA");
	c_p_S->SaveAs("./plot/"+str_tmp+".pdf");

	str_tmp = "logp_logS_" + str_adj;
	TCanvas *c_logp_logS = new TCanvas(str_tmp, str_tmp);
	TGraph *g_logp_logS = new TGraph(np, logps, logSs);
	g_logp_logS->GetXaxis()->SetTitle("ln|p-p_{c}|");
	g_logp_logS->GetYaxis()->SetTitle("ln(S(p))");
	g_logp_logS->SetMarkerStyle(20);
	g_logp_logS->Draw("PLA");
	TF1 *func_logp_logS = new TF1("func_logp_logS", "[0]+[1]*x");
	func_logp_logS->SetLineStyle(7);
	func_logp_logS->SetLineColor(kRed);
	g_logp_logS->Fit(func_logp_logS);
	TLegend *l_logp_logS = new TLegend(0.1,0.9,0.9,0.96);
	l_logp_logS->SetNColumns(2);
	l_logp_logS->AddEntry((TObject*)0, Form("slope = %.5f",func_logp_logS->GetParameter(1)), "");
	l_logp_logS->AddEntry((TObject*)0, Form("intercept = %.5f",func_logp_logS->GetParameter(0)), "");
	l_logp_logS->Draw("same");
	l_logp_logS->SetFillStyle(0);
	l_logp_logS->SetBorderSize(0);
	c_logp_logS->SaveAs("./plot/"+str_tmp+".pdf");

}

//-------------------------------------------------------------------------//

int main() {

	gStyle->SetTitleOffset(1.4,"y");
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);

	//percolation_2d_Lscan_plot();
	percolation_2d_pscan_plot(0);
	percolation_2d_pscan_plot(1);

	return 0;
}
