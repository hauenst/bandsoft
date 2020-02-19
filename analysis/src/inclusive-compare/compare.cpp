
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <unistd.h>
#include <ctype.h>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"

#include "constants.h"

TTree* readTree(TFile* inFile);
void label1D(TH1D* data, TH1D* sim, TString xlabel, TString ylabel);
void label1D_ratio(TH1D* data, TH1D* sim, TString xlabel, TString ylabel);
using namespace std;
int main(int argc, char ** argv){

	if (argc<4){
		cerr << "Wrong number of arguments used.\n\tPlease instead use: ./code [SimFile] [DataFile] [OutputFile]\n";
		return -1;
	}

	TFile * inSimFile = new TFile(argv[1]);
	TFile * inDatFile = new TFile(argv[2]);
	TTree * inSimTree = (TTree*)	readTree(inSimFile);
	TTree * inDatTree = (TTree*)	readTree(inDatFile);

	
	// Define histograms for the output
	TFile * outFile = new TFile(argv[3],"RECREATE");
	TH1D * h1dat_p_e	= new TH1D("h1dat_p_e",		"",50,0,10);
	TH1D * h1sim_p_e	= new TH1D("h1sim_p_e",		"",50,0,10);

	TH1D * h1dat_px_e	= new TH1D("h1dat_px_e",		"",50,-5,5);
	TH1D * h1sim_px_e	= new TH1D("h1sim_px_e",		"",50,-5,5);
	TH1D * h1dat_py_e	= new TH1D("h1dat_py_e",		"",50,-5,5);
	TH1D * h1sim_py_e	= new TH1D("h1sim_py_e",		"",50,-5,5);
	TH1D * h1dat_pz_e	= new TH1D("h1dat_pz_e",		"",50,0,10);
	TH1D * h1sim_pz_e	= new TH1D("h1sim_pz_e",		"",50,0,10);

	TH1D * h1dat_th_e	= new TH1D("h1dat_th_e",	"",80,0,40);
	TH1D * h1sim_th_e	= new TH1D("h1sim_th_e",	"",80,0,40);

	TH1D * h1dat_ph_e	= new TH1D("h1dat_ph_e",	"",80,-200,200);
	TH1D * h1sim_ph_e	= new TH1D("h1sim_ph_e",	"",80,-200,200);

	TH1D * h1dat_q		= new TH1D("h1dat_q",		"",50,0,10);
	TH1D * h1sim_q		= new TH1D("h1sim_q",		"",50,0,10);

	TH1D * h1dat_th_q	= new TH1D("h1dat_th_q",	"",50,0,50);
	TH1D * h1sim_th_q	= new TH1D("h1sim_th_q",	"",50,0,50);

	TH1D * h1dat_ph_q	= new TH1D("h1dat_ph_q",	"",50,-200,200);
	TH1D * h1sim_ph_q	= new TH1D("h1sim_ph_q",	"",50,-200,200);

	TH1D * h1dat_nu		= new TH1D("h1dat_nu",		"",50,0,10);
	TH1D * h1sim_nu		= new TH1D("h1sim_nu",		"",50,0,10);

	TH1D * h1dat_Q2		= new TH1D("h1dat_Q2",		"",50,0,10);
	TH1D * h1sim_Q2		= new TH1D("h1sim_Q2",		"",50,0,10);

	TH1D * h1dat_xB		= new TH1D("h1dat_xB",		"",50,0,1);
	TH1D * h1sim_xB		= new TH1D("h1sim_xB",		"",50,0,1);

	TH1D * h1dat_W		= new TH1D("h1dat_W",		"",50,0,5);
	TH1D * h1sim_W		= new TH1D("h1sim_W",		"",50,0,5);


	// Cuts used for the plots
	TCut inclusive = "weight*(lV > 15 && lW > 15 && p_e < Ebeam && E_tot > 0.25 && e_vtz < 10 && e_vtz > -15 && E_tot/p_e < 0.3 && E_tot/p_e > 0.15 && p_e>2)";

	TCanvas * c_p_e = new TCanvas("c_p_e");
	c_p_e->Divide(1,2);
	c_p_e->cd(1);
	inDatTree->Draw("p_e >> h1dat_p_e"								,inclusive);
	inSimTree->Draw("p_e >> h1sim_p_e"								,inclusive);
	label1D(	h1dat_p_e,	h1sim_p_e,	"|p_e| [GeV]",	"Counts"	);
	c_p_e->cd(2);
	label1D_ratio(	h1dat_p_e,      h1sim_p_e,      "|p_e| [GeV]",  "Data/Sim"	);
	c_p_e->Update();

	TCanvas * c_p_comp = new TCanvas("c_p_comp");
	c_p_comp->Divide(2,2);
	c_p_comp->cd(1);
	inDatTree->Draw("p_e*sin(theta_e)*cos(phi_e) >> h1dat_px_e"                                    	,inclusive);
        inSimTree->Draw("p_e*sin(theta_e)*cos(phi_e) >> h1sim_px_e"                                    	,inclusive);
	label1D(        h1dat_px_e,	h1sim_px_e,	"p_e X [GeV]",	"Counts"	);
	c_p_comp->cd(2);
	inDatTree->Draw("p_e*sin(theta_e)*sin(phi_e) >> h1dat_py_e"                                    	,inclusive);
        inSimTree->Draw("p_e*sin(theta_e)*sin(phi_e) >> h1sim_py_e"                                    	,inclusive);
	label1D(        h1dat_py_e,	h1sim_py_e,	"p_e Y [GeV]",	"Counts"	);
	c_p_comp->cd(3);
	inDatTree->Draw("p_e*cos(theta_e) >> h1dat_pz_e"                                               	,inclusive);
        inSimTree->Draw("p_e*cos(theta_e) >> h1sim_pz_e"                                               	,inclusive);
	label1D(        h1dat_pz_e,	h1sim_pz_e,	"p_e Z [GeV]",	"Counts"	);
	c_p_comp->cd(2);
	c_p_comp->Update();
	

	TCanvas * c_th_e = new TCanvas("c_th_e");
	c_th_e->Divide(1,2);
	c_th_e->cd(1);
	inDatTree->Draw("theta_e*180./TMath::Pi() >> h1dat_th_e"					,inclusive);
	inSimTree->Draw("theta_e*180./TMath::Pi() >> h1sim_th_e"					,inclusive);
	label1D(        h1dat_th_e,	h1sim_th_e,	"Theta e [deg.]",	"Counts"	);
	c_th_e->cd(2);
	label1D_ratio(	h1dat_th_e,     h1sim_th_e,     "Theta e [deg.]",       "Data/Sim"	);
	c_th_e->Update();	

	TCanvas * c_ph_e = new TCanvas("c_ph_e");
	c_ph_e->Divide(1,2);
	c_ph_e->cd(1);
	inDatTree->Draw("phi_e*180./TMath::Pi() >> h1dat_ph_e"						,inclusive);
	inSimTree->Draw("phi_e*180./TMath::Pi() >> h1sim_ph_e"						,inclusive);
	label1D(        h1dat_ph_e,	h1sim_ph_e,	"Phi e [deg.]",	"Counts"	);
	c_ph_e->cd(2);
	label1D_ratio(	h1dat_ph_e,     h1sim_ph_e,     "Phi e [deg.]", "Data/Sim"	);
	c_ph_e->Update();

	TCanvas * c_q = new TCanvas("c_q");
	c_q->Divide(1,2);
	c_q->cd(1);
	inDatTree->Draw("q >> h1dat_q"									,inclusive);
	inSimTree->Draw("q >> h1sim_q"									,inclusive);
	label1D(        h1dat_q,	h1sim_q,	"|q| [GeV]",	"Counts"	);
	c_q->cd(2);
	label1D_ratio(	h1dat_q,        h1sim_q,        "|q| [GeV]",    "Data/Sim"	);
	c_q->Update();

	TCanvas * c_th_q = new TCanvas("c_th_q");
	c_th_q->Divide(1,2);
	c_th_q->cd(1);
	inDatTree->Draw("theta_q*180./TMath::Pi() >> h1dat_th_q"					,inclusive);
	inSimTree->Draw("theta_q*180./TMath::Pi() >> h1sim_th_q"					,inclusive);
	label1D(        h1dat_th_q,	h1sim_th_q,	"Theta q [deg.]",	"Counts"	);
	c_th_q->cd(2);
	label1D_ratio(	h1dat_th_q,     h1sim_th_q,     "Theta q [deg.]",       "Data/Sim"	);
	c_th_q->Update();

	TCanvas * c_ph_q = new TCanvas("c_ph_q");
	c_ph_q->Divide(1,2);
	c_ph_q->cd(1);
	inDatTree->Draw("phi_q*180./TMath::Pi() >> h1dat_ph_q"						,inclusive);
	inSimTree->Draw("phi_q*180./TMath::Pi() >> h1sim_ph_q"						,inclusive);
	label1D(        h1dat_ph_q,	h1sim_ph_q,	"Phi q [deg.]",	"Counts"	);
	c_ph_q->cd(2);
	label1D_ratio(	h1dat_ph_q,     h1sim_ph_q,     "Phi q [deg.]", "Data/Sim"	);
	c_ph_q->Update();

	TCanvas * c_nu = new TCanvas("c_nu");
	c_nu->Divide(1,2);
	c_nu->cd(1);
	inDatTree->Draw("nu >> h1dat_nu"								,inclusive);
	inSimTree->Draw("nu >> h1sim_nu"								,inclusive);
	label1D(        h1dat_nu,	h1sim_nu,	"Nu [GeV]",	"Counts"	);
	c_nu->cd(2);
	label1D_ratio(	h1dat_nu,       h1sim_nu,       "Nu [GeV]",     "Data/Sim"	);
	c_nu->Update();
	
	TCanvas * c_Q2 = new TCanvas("c_Q2");
	c_Q2->Divide(1,2);
	c_Q2->cd(1);
	inDatTree->Draw("Q2 >> h1dat_Q2"								,inclusive);
	inSimTree->Draw("Q2 >> h1sim_Q2"								,inclusive);
	label1D(        h1dat_Q2,	h1sim_Q2,	"Q2 [GeV^2]",	"Counts"	);
	c_Q2->cd(2);
	label1D_ratio(	h1dat_Q2,       h1sim_Q2,       "Q2 [GeV^2]",   "Data/Sim"	);
	c_Q2->Update();

	TCanvas * c_xB = new TCanvas("c_xB");
	c_xB->Divide(1,2);
	c_xB->cd(1);
	inDatTree->Draw("xB >> h1dat_xB"								,inclusive);
	inSimTree->Draw("xB >> h1sim_xB"								,inclusive);
	label1D(        h1dat_xB,	h1sim_xB,	"xB",	"Counts"	);
	c_xB->cd(2);
	label1D_ratio(	h1dat_xB,       h1sim_xB,       "xB",   "Data/Sim"	);
	c_xB->Update();

	TCanvas * c_W = new TCanvas("c_W");
	c_W->Divide(1,2);
	c_W->cd(1);
	inDatTree->Draw("sqrt(W2) >> h1dat_W"								,inclusive);
	inSimTree->Draw("sqrt(W2) >> h1sim_W"								,inclusive);
	label1D(        h1dat_W,	h1sim_W,	"W [GeV]",	"Counts"	);
	c_W->cd(2);
	label1D_ratio(	h1dat_W,        h1sim_W,        "W [GeV]",      "Data/Sim"	);
	c_W->Update();


	TCanvas * c0 = new TCanvas("c0","c0",700,900);
	c0 -> Modified();
	c0 -> Update(); 
	c0 -> Print(Form("comparison_test.pdf(" ));
	c_p_e	->Print(Form("comparison_test.pdf(" ));
	c_p_comp->Print(Form("comparison_test.pdf(" ));
	c_th_e	->Print(Form("comparison_test.pdf(" ));
	c_ph_e	->Print(Form("comparison_test.pdf(" ));
	c_q	->Print(Form("comparison_test.pdf(" ));
	c_th_q	->Print(Form("comparison_test.pdf(" ));
	c_ph_q	->Print(Form("comparison_test.pdf(" ));
	c_nu	->Print(Form("comparison_test.pdf(" ));
	c_Q2	->Print(Form("comparison_test.pdf(" ));
	c_xB	->Print(Form("comparison_test.pdf(" ));
	c_W	->Print(Form("comparison_test.pdf(" ));
	c0 -> Print(Form("comparison_test.pdf)" ));
	
	

	outFile->cd();
	c_p_e->Write();
	c_p_comp->Write();
	c_th_e->Write();
	c_ph_e->Write();
	c_q->Write();
	c_th_q->Write();
	c_ph_q->Write();
	c_nu->Write();
	c_Q2->Write();
	c_xB->Write();
	c_W->Write();
	outFile->Close();

	return 1;
}

TTree* readTree(TFile* inFile){
	TTree * inTree = (TTree*)inFile->Get("skim");
	double Ebeam		= 0;
	double gated_charge	= 0;
	double livetime		= 0;
	double starttime	= 0;
	int ePid		= 0;
	int eCharge		= 0;
	int eStatus		= 0;
	double eTime		= 0;
	double eBeta 		= 0;
	double eChi2pid		= 0;
	double E_tot		= 0;
	double E_pcal		= 0;
	double t_e		= 0;
	double dL_e		= 0;
	double lU		= 0;
	double lV		= 0;
	double lW		= 0;
	double e_vtx		= 0;
	double e_vty		= 0;
	double e_vtz		= 0;
	double p_e		= 0;
	double theta_e		= 0;
	double phi_e		= 0;
	double q		= 0;
	double theta_q		= 0;
	double phi_q		= 0;
	double nu		= 0;
	double Q2		= 0;
	double xB		= 0;
	double W2		= 0;
	double weight		= 0;
	inTree->Branch("Ebeam"		,&Ebeam			);
	inTree->Branch("gated_charge"	,&gated_charge		);
	inTree->Branch("livetime"	,&livetime		);
	inTree->Branch("starttime"	,&starttime		);
	inTree->Branch("ePid"		,&ePid			);
	inTree->Branch("eCharge"	,&eCharge		);
	inTree->Branch("eStatus"	,&eStatus		);
	inTree->Branch("eTime"		,&eTime			);
	inTree->Branch("eBeta"		,&eBeta 		);
	inTree->Branch("eChi2pid"	,&eChi2pid		);
	inTree->Branch("E_tot"		,&E_tot			);
	inTree->Branch("E_pcal"		,&E_pcal		);
	inTree->Branch("t_e"		,&t_e			);
	inTree->Branch("dL_e"		,&dL_e			);
	inTree->Branch("lU"		,&lU			);
	inTree->Branch("lV"		,&lV			);
	inTree->Branch("lW"		,&lW			);
	inTree->Branch("e_vtx"		,&e_vtx			);
	inTree->Branch("e_vty"		,&e_vty			);
	inTree->Branch("e_vtz"		,&e_vtz			);
	inTree->Branch("p_e"		,&p_e			);
	inTree->Branch("theta_e"	,&theta_e		);
	inTree->Branch("phi_e"		,&phi_e			);
	inTree->Branch("q"		,&q			);
	inTree->Branch("theta_q"	,&theta_q		);
	inTree->Branch("phi_q"		,&phi_q			);
	inTree->Branch("nu"		,&nu			);
	inTree->Branch("Q2"		,&Q2			);
	inTree->Branch("xB"		,&xB			);
	inTree->Branch("W2"		,&W2			);
	inTree->Branch("weight"		,&weight		);
	
	return inTree;
}

void label1D(TH1D* data, TH1D* sim, TString xlabel, TString ylabel){
	data->SetLineColor(1);
	data->SetLineWidth(3);
	data->SetStats(0);

	sim->SetLineColor(9);
	sim->SetLineWidth(3);
	sim->SetStats(0);
	sim->Scale(data->Integral() / sim->Integral() );

	data->Draw("h,p");
	sim->Draw("hist,same");

	double max1 = data->GetMaximum()*1.1;
	double max2 = sim->GetMaximum()*1.1;
	data->GetYaxis()->SetRangeUser(0,max(max1,max2));
	
	data->GetXaxis()->SetTitle(xlabel);
	data->GetYaxis()->SetTitle(ylabel);

	TLegend * legend = new TLegend(0.7,0.8,0.9,0.9);
	//legend->AddEntry(data,"Radiation On","f");
	//legend->AddEntry(sim,"Radiation Off","f");
	legend->AddEntry(data,"Data","f");
	legend->AddEntry(sim,"Sim","f");
	legend->Draw();

	return;
}
void label1D_ratio(TH1D* data, TH1D* sim, TString xlabel, TString ylabel){
	TH1D * data_copy = (TH1D*) data->Clone();
	TH1D * sim_copy = (TH1D*) sim->Clone();
	
	data_copy->SetLineColor(1);
	data_copy->SetLineWidth(3);
	data_copy->SetStats(0);

	sim_copy->SetLineColor(9);
	sim_copy->SetLineWidth(3);
	sim_copy->SetStats(0);
	sim_copy->Scale(data_copy->Integral() / sim_copy->Integral() );

	data_copy->Divide(sim_copy);
	data_copy->Draw("ep");
	TLine* line = new TLine(data_copy->GetXaxis()->GetBinCenter(1), 1., data_copy->GetXaxis()->GetBinCenter(data_copy->GetXaxis()->GetNbins()), 1.);
	line->SetLineWidth(3);
	line->SetLineColor(2);
	line->Draw("same");

	double max1 = data_copy->GetMaximum()*1.1;
	data_copy->GetYaxis()->SetRangeUser(0,2);
	
	data_copy->GetXaxis()->SetTitle(xlabel);
	data_copy->GetYaxis()->SetTitle(ylabel);

	return;
}
