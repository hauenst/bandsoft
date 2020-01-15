#include <cstdlib>
#include <iostream>

#include "TRint.h"
#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TCanvas.h"

#include "reader.h"
#include "bank.h"

#include "BParticle.h"
#include "BCalorimeter.h"
#include "BScintillator.h"
#include "BBand.h"
#include "BEvent.h"

using namespace std;

// Forward-declaring functions
void PrettyTH1F(TH1F * h1,TString titx,TString tity,int color);
void PrettyTH2F(TH2F * h2,TString titx,TString tity);
double fn_Emiss(double Pmiss, double omega, double M_tar, double Enuc, double Mnuc);

// ========================================================================================================================================
int main(int argc, char** argv) {

#ifdef WITHRINT
	TRint *myapp = new TRint("RootSession",&argc,argv,NULL,0);
#else
	TApplication *myapp = new TApplication("myapp",0,0);
#endif

	// ----------------------------------------------------------------------------------
	// Useful variables
	double mp      = 0.93827; //GeV (proton mass      )
	double mPiC    = 0.13957; //GeV (charged pion mass)
	double mD      = 1.8756;  //GeV (deuteron mass    )
	double rad2deg = 180./3.14159;

	// ----------------------------------------------------------------------------------
	// Getting input arguments
	double Ebeam, mtar;

	if(argc>=3){
		if(atoi(argv[1])==1){
			cout << "Will assume this hipo file corresponds to: Ebeam =  6.4 GeV, target = H (i.e. RGA)" << endl;
			Ebeam = 6.4; //GeV
			mtar  = mp;
		}
		else if(atoi(argv[1])==2){
			cout << "Will assume this hipo file corresponds to: Ebeam = 10.6 GeV, target = D (i.e. RGB)" << endl;
			Ebeam = 10.41; //GeV
			mtar  = mD;
		}
	}
	else {
		cout << "=========================\nRun this code as:\n./code A path/to/input/file\n" << endl;
		cout << "where: A = 1 -> Ebeam =  6.4 GeV, target = H (i.e. RGA)" << endl;
		cout << "         = 2 -> Ebeam = 10.6 GeV, target = D (i.e. RGB)" << endl;
		cout << "=========================" << endl;
		exit(0);
	}

	TVector3 V3_Ebeam(0,0,Ebeam);
	TLorentzVector V4_Ebeam(V3_Ebeam,Ebeam);
	TLorentzVector V4_mtar(0,0,0,mtar);

	// ----------------------------------------------------------------------------------
	// Event selection cuts
	double cut_ep      =     2; //GeV
	double cut_chi2pid =     5;
	double cut_min_vz  =   -15; //cm
	double cut_max_vz  =    10; //cm
	double cut_W       =     0; //GeV
	double cut_uvw     =    15; //cm
	double cut_Epcal   = 0.060; //GeV (60 MeV)
	double cut_tof_e   =    10; //ns

	TH1I *h1_triggerbits = new TH1I("h1_triggerbits","Neg. Hadron Sector 1 and Positive Hadron Sector 4",16,-0.5,15.5);
        h1_triggerbits->GetXaxis()->SetTitle("Trigger bit");

        TH1F *h1_phi_plus = new TH1F("h1_phi_plus","Phi distribution positive hadron",360,-180,180);
        h1_phi_plus->GetXaxis()->SetTitle("#phi [degree]");

        TH1F *h1_phi_minus = new TH1F("h1_phi_minus","Phi distribution negative hadron",360,-180,180);
        h1_phi_minus->GetXaxis()->SetTitle("#phi [degree]");

        TH1F *h1_mom_plus = new TH1F("h1_mom_plus","Momentum distribution positive hadron",200,0,10);
        h1_mom_plus->GetXaxis()->SetTitle("Momentum [GeV]");

        TH1F *h1_mom_minus = new TH1F("h1_mom_minus","Momentum distribution negative hadron",200,0,10);
        h1_mom_minus->GetXaxis()->SetTitle("Momentum [GeV]");

        TH1F *h1_the_plus = new TH1F("h1_the_plus","Theta distribution positive hadron",50,5,50);
        h1_the_plus->GetXaxis()->SetTitle("#Theta [degree]");

        TH1F *h1_the_minus = new TH1F("h1_the_minus","Theta distribution negative hadron",50,5,50);
        h1_the_minus->GetXaxis()->SetTitle("#Theta [degree]");

        TH2F *h2_the_mom_plus = new TH2F("h2_the_mom_plus","Theta-Mom dist positive hadron",40,5,50,40,0,10);
        h2_the_mom_plus->GetXaxis()->SetTitle("#Theta [degree]");
        h2_the_mom_plus->GetYaxis()->SetTitle("Momentum [GeV]");

        TH2F *h2_the_mom_minus = new TH2F("h2_the_mom_minus","Theta-Mom dist negative hadron",40,5,50,40,0,10);
        h2_the_mom_minus->GetXaxis()->SetTitle("#Theta [degree]");
        h2_the_mom_minus->GetYaxis()->SetTitle("Momentum [GeV]");

        TH1F *h1_pcal_energy_pos[6];
	TH1F *h1_pcal_energy_neg[6];

	TH1F *h1_ecal_energy_pos[6];
	TH1F *h1_ecal_energy_neg[6];

	TH1F *h1_ecal_tot_e_pos[6];
	TH1F *h1_ecal_tot_e_neg[6];
	for (int i = 0; i < 6; i++) {
		int pos_sector = i+1; // 1 , 2 , 3 , 4 , 5 , 6 ,
		int neg_sector = i<=2 ? 4+i : i-2;  // 4 , 5 , 6 , 1, 2, 3;
	//	cout << i  << " pos sector " << pos_sector << " neg se " << neg_sector << endl; 
		h1_pcal_energy_pos[i] = new TH1F(Form("h1_pcal_energy_pos_%i",i),Form("PCal energy sector %i for pos. hadron in %i",i+1,i+1),100,0.01,0.2);
                h1_pcal_energy_neg[i] = new TH1F(Form("h1_pcal_energy_neg_%i",i),Form("PCal energy sector %i for neg. hadron in %i",i+1,i+1),100,0.01,0.2);
		h1_pcal_energy_pos[i]->GetXaxis()->SetTitle("PCal energy [GeV]");
                h1_pcal_energy_neg[i]->GetXaxis()->SetTitle("PCal energy [GeV]");
	
		h1_ecal_energy_pos[i] = new TH1F(Form("h1_ecal_energy_pos_%i",i),Form("EC outer energy sector %i for pos. hadron in %i",i+1,i+1),100,0.01,0.2);
                h1_ecal_energy_neg[i] = new TH1F(Form("h1_ecal_energy_neg_%i",i),Form("EC outer energy sector %i for neg. hadron in %i",i+1,i+1),100,0.01,0.2);
                h1_ecal_energy_pos[i]->GetXaxis()->SetTitle("ECal energy [GeV]");
                h1_ecal_energy_neg[i]->GetXaxis()->SetTitle("ECal energy [GeV]");

		h1_ecal_tot_e_pos[i] = new TH1F(Form("h1_ecal_tot_e_pos_%i",i),Form("EC total energy sector %i for pos. hadron in %i",i+1,i+1),100,0.01,0.2);
                h1_ecal_tot_e_neg[i] = new TH1F(Form("h1_ecal_tot_e_neg_%i",i),Form("EC total energy sector %i for neg. hadron in %i",i+1,i+1),100,0.01,0.2);
                h1_ecal_tot_e_pos[i]->GetXaxis()->SetTitle("ECal tot. energy [GeV]");
                h1_ecal_tot_e_neg[i]->GetXaxis()->SetTitle("ECal tot. energy [GeV]");
	}


	TH1F *h1_pcal_energy_sec1 = new TH1F("h1_pcal_energy_sec1","PCal energy sector 1",100,0.01,0.2);
        h1_pcal_energy_sec1 ->GetXaxis()->SetTitle("PCal Energy [GeV]");

        TH1F *h1_pcal_energy_sec4 = new TH1F("h1_pcal_energy_sec4","PCal energy sector 4",100,0.01,0.2);
        h1_pcal_energy_sec4->GetXaxis()->SetTitle("PCal Energy [GeV]");

        TH1F *h1_ecal_energy_sec1 = new TH1F("h1_ecal_energy_sec1","EC outer energy sector 1",100,0.01,0.2);
        h1_ecal_energy_sec1 ->GetXaxis()->SetTitle("ECal Energy [GeV]");

        TH1F *h1_ecal_energy_sec4 = new TH1F("h1_ecal_energy_sec4","EC outer energy sector 4",100,0.01,0.2);
        h1_ecal_energy_sec4 ->GetXaxis()->SetTitle("ECal Energy [GeV]");

	TH1F *h1_track_chi2 = new TH1F("h1_track_chi2","Track Chi2 from REC::Tracks",400,0,500);

	for( int file_i = 2 ; file_i < argc ; file_i++) { 
	   TString inputFile = argv[file_i];
	   // ----------------------------------------------------------------------------------
	   // ----------------------------------------------------------------------------------
	   // Opening input HIPO file
	   hipo::reader reader;
	   reader.open(inputFile);

	   //Read Dictionary of Hipo File  // new hipo4
           hipo::dictionary  factory;      // new hipo4
           reader.readDictionary(factory); // new hipo4
           //factory.show();               // new hipo4

	   BEvent        event       (factory.getSchema("REC::Event"       ));
	   BParticle     particles   (factory.getSchema("REC::Particle"    ));
	   BCalorimeter  calo        (factory.getSchema("REC::Calorimeter" ));
	   BScintillator scintillator(factory.getSchema("REC::Scintillator"));

	   hipo::bank run_config  (factory.getSchema("RUN::config"));
	   hipo::bank rec_track   (factory.getSchema("REC::Track"));
	   //One also needs a hipo::event object which is called from the reader for each event to get
           //the information for each bank
           hipo::event readevent;  // new hipo4

	   int event_counter = 0;
	   int goodevents = 0;

	   // ----------------------------------------------------------------------------------
	   // Loop over events and print them on the screen
	   while(reader.next()==true){

		//Reader has to load information about event in hipo::event class
                reader.read(readevent); // new hipo4

                //Load explicitly all information for each bank for the event
                readevent.getStructure(event       );   // new hipo4
                readevent.getStructure(particles   );   // new hipo4
                readevent.getStructure(calo        );   // new hipo4
                readevent.getStructure(scintillator);   // new hipo4
		readevent.getStructure(rec_track);
		readevent.getStructure(run_config);
		int* trigger_bits1 = new int[64];
      			
		long TriggerWord = run_config.getLong("trigger",0);
		//cout << "Triggerword" << TriggerWord << endl;
		for (int i=63; i>=0; i--) {
           		   if((TriggerWord & (((long) 1) << i)) != 0) trigger_bits1[i] = 1;
              		   else trigger_bits1[i] = 0;
          	}
	 
	
                //Now everything is loaded and can be used as before with HIPO3 files. There is only one difference:
                //The number of hits in each bank is determined by the function "getRows()" and not by "getSize" as before.
		if(event_counter%1000000==0) cout << "event: " << event_counter << endl;
		event_counter++;

		// Particle bank
		int pid0       = particles.getPid    (0);	// electron candidate id assigned by clas
		// -------------------------------------------------------------------------
		// Fill some histograms before cutting on good electrons
		// - (DONE)	p > 2 GeV
		// - (DONE)	p < Ebeam
		// - (DONE)	TOF > 10 ns //(need bank 330) event.json
		// - (DONE)	vz: -15 to 10 cm
		// - (DONE)	W > 0 GeV
		// - 		Sampling fraction +/-3sigma about E/p vs. p mean
		// - (DONE)	~15 cm fiducial cuts on U, V, W to contain full shower (need bank 332 lu, lv, lw)
		// - (DONE)	abs(chisq PID) < 5 (goodness of PID from EB)
		// - (DONE)	PCAL > 60 MeV (to remove min-i) (bank 332 layers 1(PCAL), 4(EC inner), 7(EC outter))
		// -------------------------------------------------------------------------
		// Only keep events for which the first particle is an electron
/*		if(             (pid0!=11              )||
				(chr0!=-1              )||
				(chi2pid>=cut_chi2pid  )||
				(ep<=cut_ep            )||
				(ep>=Ebeam             )||
				(V3_ev.Z()>cut_max_vz  )||
				(V3_ev.Z()<cut_min_vz  )||
				(lU<cut_uvw            )||
				(lV<cut_uvw            )||
				(lW<cut_uvw            )||
				(Epcal<cut_Epcal       )||
				(TMath::Sqrt(W2)<=cut_W)||
				(tof_e<cut_tof_e       )
		  ) continue;

*/		// -------------------------------------------------------------------------
		//I dont want to have an electron
		if (   pid0==11 ) continue;

		//if ( event_counter > 100) continue;
		int nParticles = particles.getRows();	
		int nPlus =  0;
		int nMinus = 0;
		int tmp_fast_pip_idx = -1;
		int tmp_fast_pim_idx = -1;
		int sectorid_plus = -1;
		int sectorid_minus= -1;	

		int nTracks = rec_track.getRows();
		//cout << "NTracks " << nTracks << endl;
		if (nTracks < 2) continue;
		//First track infos
		int sector_tr0 = rec_track.getInt("sector",0);
		int charge_tr0 = rec_track.getInt("q",0);
		float chi_tr0  = rec_track.getFloat("chi2",0);
		int pindex_tr0 = rec_track.getInt("pindex",0);
		if (chi_tr0 > 150 || pindex_tr0 != 0 || sector_tr0 == 0 || charge_tr0 == 0) continue;
		if (charge_tr0 == -1) {
		   nMinus++;
		   sectorid_minus = sector_tr0;
		   tmp_fast_pim_idx = pindex_tr0;
		}
		else if (charge_tr0 == 1) {
  		   nPlus++;
		   sectorid_plus = sector_tr0;
		   tmp_fast_pip_idx = pindex_tr0;
		}
		else { cout << "This should not happen. Charge of first track after cuts is " << charge_tr0 << endl; } 
		for(int track = 1; track < nTracks; track++) {
			int pindex_track = rec_track.getInt("pindex",track);
			float chi2 = rec_track.getFloat("chi2",track);
			int sector = rec_track.getInt("sector",track);
			int charge = rec_track.getInt("q",track); 
			h1_track_chi2->Fill(chi2);
			if (chi2 > 150) continue;

		//	if (sector > 0)			cout << "charge" << charge << " sector " << sector << " pindex " << pindex_track << " chi2 " << chi2 << endl;
	
			if (sector > 0 && abs(sector - sector_tr0) == 3 && charge == -1) {
				if (sectorid_minus < 0) {
				  sectorid_minus = sector;
				}
				nMinus++;
				if (tmp_fast_pim_idx < 0 ) {
                                  tmp_fast_pim_idx = pindex_track;
                                }
			}
			if (sector > 0 && abs(sector - sector_tr0) == 3 && charge == 1 ) {
				if (sectorid_plus < 0) {
				  sectorid_plus = sector;
				}
				nPlus++;
				if (tmp_fast_pip_idx < 0 ) {
                                  tmp_fast_pip_idx = pindex_track;
                                }
			}
		}

		//if ( sectorid_plus>0 && sectorid_minus>0 && abs(sectorid_minus - sectorid_plus) == 3 )  cout << "found sector match " << endl;

		if (nPlus>=1 && nMinus>=1 && sectorid_plus>0 && sectorid_minus>0 && abs(sectorid_minus - sectorid_plus) == 3) { 
			//cout << nPiPlus << " " << nPiMinus << " found one goodevents is " << goodevents << endl; 
			goodevents++; 
		
//			int scintRows = scintillator.getRows();
//			int caloRows  = calo.getRows();

/*			for (int row_scint = 0; row_scint < scintRows; row_scint++) {
				int pindex = scintillator.getPindex( row_scint);
				int sector = scintillator.getSector( row_scint);
				int detector = scintillator.getDetector(row_scint);
				//FTOF detector ID is 12
				if (pindex == tmp_fast_pip_idx && detector == 12) sectorid_plus = sector; 
				if (pindex == tmp_fast_pim_idx && detector == 12) sectorid_minus = sector;				
			}
*/
			//Loop over all sector combinations for plus and minus (each dihadron bit) and fill PCal and ECal outer energy distribution for each combination
			for (int i_sector = 1; i_sector < 7; i_sector++) {
				//i_sector defines sector with negative hadron
				//oppo_sector is sector with positive hadron
				int oppo_sector = i_sector<=3 ? 3+i_sector : i_sector-3; //Sector numbers 4 , 5, 6, 1, 2, 3
				if ( (sectorid_plus == oppo_sector && sectorid_minus == i_sector) ) {
					h1_pcal_energy_neg[i_sector-1]->Fill(calo.getPcalE(tmp_fast_pim_idx));
					h1_pcal_energy_pos[oppo_sector-1]->Fill(calo.getPcalE(tmp_fast_pip_idx));

					h1_ecal_energy_neg[i_sector-1]->Fill(calo.getECoutE(tmp_fast_pim_idx));
                                        h1_ecal_energy_pos[oppo_sector-1]->Fill(calo.getECoutE(tmp_fast_pip_idx));

					h1_ecal_tot_e_neg[i_sector-1]->Fill(calo.getECinE(tmp_fast_pim_idx)+calo.getECoutE(tmp_fast_pim_idx));
                                        h1_ecal_tot_e_pos[oppo_sector-1]->Fill(calo.getECinE(tmp_fast_pip_idx)+calo.getECoutE(tmp_fast_pip_idx));
				}
			}
			if ( (sectorid_plus == 4 && sectorid_minus == 1) ) { 
	//			std::cout << "FOund trigger 7 combination" << std::endl;//should be trigger bit 7 dominated
				//cout << "trigger bits 1 -16 are " ;
				for (int i = 0; i < 15; i++) { 
				//	cout << trigger_bits1[i] << " " ;
					if (trigger_bits1[i]!=0) {			
						h1_triggerbits->Fill(i);
					} 
				}
				h1_phi_plus->Fill(180/TMath::Pi()*(particles.getV3P(tmp_fast_pip_idx).Phi())); 
				h1_phi_minus->Fill(180/TMath::Pi()*(particles.getV3P(tmp_fast_pim_idx).Phi()));
			
				h1_the_plus->Fill(180/TMath::Pi()*(particles.getV3P(tmp_fast_pip_idx).Theta()));
                                h1_the_minus->Fill(180/TMath::Pi()*(particles.getV3P(tmp_fast_pim_idx).Theta()));

				h1_mom_plus->Fill(particles.getV3P(tmp_fast_pip_idx).Mag());
                                h1_mom_minus->Fill(particles.getV3P(tmp_fast_pim_idx).Mag());

				h2_the_mom_plus->Fill(180/TMath::Pi()*(particles.getV3P(tmp_fast_pip_idx).Theta()), particles.getV3P(tmp_fast_pip_idx).Mag());
                                h2_the_mom_minus->Fill(180/TMath::Pi()*(particles.getV3P(tmp_fast_pim_idx).Theta()), particles.getV3P(tmp_fast_pim_idx).Mag());
			
				h1_pcal_energy_sec1->Fill(calo.getPcalE(tmp_fast_pim_idx));
                        	h1_pcal_energy_sec4->Fill(calo.getPcalE(tmp_fast_pip_idx));

				h1_ecal_energy_sec1->Fill(calo.getECoutE(tmp_fast_pim_idx));
				h1_ecal_energy_sec4->Fill(calo.getECoutE(tmp_fast_pip_idx));
			}
		}


	   } // End of while loop over events

	cout << goodevents << " from total " << event_counter <<  endl;
	} //end of loop over file_i
	// -------------------------------------------------------------------------------------------

	TCanvas *c1 = new TCanvas("c1","");
	c1->cd();
	h1_triggerbits->Draw("BAR");
//	h1_track_chi2->Draw("");

	c1->SaveAs("triggerbit.pdf");

        TCanvas *c2 = new TCanvas("c2","");
        c2->Divide(2,2);
	c2->cd(1);
        h1_phi_plus->Draw("E");
        c2->cd(2);
        h1_phi_minus->Draw("E");
        c2->cd(3);
        h1_the_plus->Draw("E");
        c2->cd(4);
        h1_the_minus->Draw("E");
	c2->SaveAs("1D_phithe.pdf");
	
	TCanvas *c3 = new TCanvas("c3","");
        c3->Divide(2,2);
        c3->cd(1);
        h1_mom_plus->Draw("E");
        c3->cd(2);
        h1_mom_minus->Draw("E");
        c3->cd(3);
        h2_the_mom_plus->Draw("COLZ");
        c3->cd(4);
        h2_the_mom_minus->Draw("COLZ");
	c3->SaveAs("2D_themom.pdf");
	
	/*TCanvas *c4 = new TCanvas("c4","");
        c4->Divide(2,2);
        c4->cd(1);
        h1_pcal_energy_sec1->Draw("");
        c4->cd(2);
        h1_pcal_energy_sec4->Draw("");
        c4->cd(3);
        h1_ecal_energy_sec1->Draw("");
        c4->cd(4);
        h1_ecal_energy_sec4->Draw("");
	c4->SaveAs("energy_ecal-pcal.pdf");
*/
	//PCAL EC OUT energy plots for neg. hadron in Sector 1 and positive hadron in Sector 4
	TCanvas *ctrig1 = new TCanvas("ctrig1","");
        ctrig1->Divide(2,3);
        ctrig1->cd(1);
        h1_pcal_energy_neg[0]->Draw("");
        ctrig1->cd(2);
        h1_pcal_energy_pos[3]->Draw("");
        ctrig1->cd(3);
        h1_ecal_energy_neg[0]->Draw("");
        ctrig1->cd(4);
        h1_ecal_energy_pos[3]->Draw("");
	ctrig1->cd(5);
	h1_ecal_tot_e_neg[0]->Draw("");
        ctrig1->cd(6);
        h1_ecal_tot_e_pos[3]->Draw("");
        ctrig1->SaveAs("energy_ecal-pcal_trig1.pdf");

	//PCAL EC OUT energy plots for neg. hadron in Sector 2 and positive hadron in Sector 5
	//
        TCanvas *ctrig2 = new TCanvas("ctrig2","");
        ctrig2->Divide(2,3);
        ctrig2->cd(1);
        h1_pcal_energy_neg[1]->Draw("");
        ctrig2->cd(2);
        h1_pcal_energy_pos[4]->Draw("");
        ctrig2->cd(3);
        h1_ecal_energy_neg[1]->Draw("");
        ctrig2->cd(4);
        h1_ecal_energy_pos[4]->Draw("");
 	ctrig2->cd(5);
        h1_ecal_tot_e_neg[1]->Draw("");
        ctrig2->cd(6);
        h1_ecal_tot_e_pos[4]->Draw("");
        ctrig2->SaveAs("energy_ecal-pcal_trig2.pdf");

	//PCAL EC OUT energy plots for neg. hadron in Sector 3 and positive hadron in Sector 6
	TCanvas *ctrig3 = new TCanvas("ctrig3","");
        ctrig3->Divide(2,3);
        ctrig3->cd(1);
        h1_pcal_energy_neg[2]->Draw("");
        ctrig3->cd(2);
        h1_pcal_energy_pos[5]->Draw("");
        ctrig3->cd(3);
        h1_ecal_energy_neg[2]->Draw("");
        ctrig3->cd(4);
        h1_ecal_energy_pos[5]->Draw("");
        ctrig3->cd(5);
        h1_ecal_tot_e_neg[2]->Draw("");
        ctrig3->cd(6);
        h1_ecal_tot_e_pos[5]->Draw("");
        ctrig3->SaveAs("energy_ecal-pcal_trig3.pdf");

        //PCAL EC OUT energy plots for neg. hadron in Sector 4 and positive hadron in Sector 1
        TCanvas *ctrig4 = new TCanvas("ctrig4","");
        ctrig4->Divide(2,3);
        ctrig4->cd(1);
        h1_pcal_energy_neg[3]->Draw("");
        ctrig4->cd(2);
        h1_pcal_energy_pos[0]->Draw("");
        ctrig4->cd(3);
        h1_ecal_energy_neg[3]->Draw("");
        ctrig4->cd(4);
        h1_ecal_energy_pos[0]->Draw("");
        ctrig4->cd(5);
        h1_ecal_tot_e_neg[3]->Draw("");
        ctrig4->cd(6);
        h1_ecal_tot_e_pos[0]->Draw("");
        ctrig4->SaveAs("energy_ecal-pcal_trig4.pdf");

        //PCAL EC OUT energy plots for neg. hadron in Sector 4 and positive hadron in Sector 2
        TCanvas *ctrig5 = new TCanvas("ctrig5","");
        ctrig5->Divide(2,3);
        ctrig5->cd(1);
        h1_pcal_energy_neg[4]->Draw("");
        ctrig5->cd(2);
        h1_pcal_energy_pos[1]->Draw("");
        ctrig5->cd(3);
        h1_ecal_energy_neg[4]->Draw("");
        ctrig5->cd(4);
        h1_ecal_energy_pos[1]->Draw("");
        ctrig5->cd(5);
        h1_ecal_tot_e_neg[4]->Draw("");
        ctrig5->cd(6);
        h1_ecal_tot_e_pos[1]->Draw("");
        ctrig5->SaveAs("energy_ecal-pcal_trig5.pdf");

        //PCAL EC OUT energy plots for neg. hadron in Sector 6 and positive hadron in Sector 3
        TCanvas *ctrig6 = new TCanvas("ctrig6","");
        ctrig6->Divide(2,3);
        ctrig6->cd(1);
        h1_pcal_energy_neg[5]->Draw("");
        ctrig6->cd(2);
        h1_pcal_energy_pos[2]->Draw("");
        ctrig6->cd(3);
        h1_ecal_energy_neg[5]->Draw("");
        ctrig6->cd(4);
        h1_ecal_energy_pos[2]->Draw("");
        ctrig6->cd(5);
        h1_ecal_tot_e_neg[5]->Draw("");
        ctrig6->cd(6);
        h1_ecal_tot_e_pos[2]->Draw("");
	ctrig6->SaveAs("energy_ecal-pcal_trig6.pdf");


	myapp -> Run();
	return 0;
}
// ========================================================================================================================================
void PrettyTH1F(TH1F * h1,TString titx,TString tity,int color) {
	h1 -> GetXaxis() -> SetTitle(titx);
	h1 -> GetYaxis() -> SetTitle(tity);
	h1 -> SetLineColor(color);
	h1 -> SetLineWidth(2);
}
// ========================================================================================================================================
void PrettyTH2F(TH2F * h2,TString titx,TString tity) {
	h2 -> GetXaxis() -> SetTitle(titx);
	h2 -> GetYaxis() -> SetTitle(tity);
}
// ========================================================================================================================================
double fn_Emiss(double Pmiss, double omega, double M_tar, double Enuc, double Mnuc){
	// Calculates missing energy
	// Takes as input: missing momentum, transfer energy, struck nucleon energy, and struck nucleon mass.
	double Tb   = omega + M_tar - Enuc - TMath::Sqrt(pow(omega + M_tar - Enuc,2)- Pmiss*Pmiss );	// Kinetic energy of A-1 system
	double Tnuc = Enuc - Mnuc;                                                             			// Kinetic energy of struck neutron
	return omega - Tnuc - Tb;
}
