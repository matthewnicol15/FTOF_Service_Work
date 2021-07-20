#include "DC_Fiducial_Cuts_CLAS12.cxx"
#include <iostream>
#include <cstdlib>
#include <chrono>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <vector>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TPaletteAxis.h>
#include <TBenchmark.h>
#include "clas12reader.h"
#include <stdlib.h>
#include "Riostream.h"
#include "TLine.h"
#include "TVirtualPad.h"
#include "TClass.h"
#include "TVirtualX.h"
#include "TMath.h"
#include "TStyle.h"

using namespace clas12;

void SetLorentzVector(TLorentzVector &p4,clas12::region_part_ptr rp){
  p4.SetXYZM(rp->par()->getPx(),rp->par()->getPy(),
  rp->par()->getPz(),p4.M());

}

void FTOF_2pi_Missing_Pim(){

  // Data files to process
<<<<<<< HEAD
  TString inputFile1("/home/matthewn/links/RGB_Spring_2019_Inbending_dst/*.hipo");
=======
  TString inputFile1("/lustre19/expphy/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass1/v0/dst/train/skim4/skim4_0050*.hipo");
>>>>>>> 746002c24cbd1e4ca1e3da8f0b430e15d2e783da
  // TString inputFile1("/lustre19/expphy/volatile/clas12/rg-a/production/recon/fall2018/torus+1/pass1/v1/dst/train/skim4/*.hipo");
  // TString inputFile1("/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass1/v0/dst/train/skim4/skim4_005117.hipo");
  // TString inputFile2("/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass1/v0/dst/train/skim4/skim4_005124.hipo");
  // TString inputFile3("/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass1/v0/dst/train/skim4/skim4_005125.hipo");
  // TString inputFile4("/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass1/v0/dst/train/skim4/skim4_005126.hipo");
  // TString inputFile5("/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass1/v0/dst/train/skim4/skim4_005128.hipo");
  // TString inputFile6("/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass1/v0/dst/train/skim4/skim4_005130.hipo");

  // Creating a TChain of all the input files
  TChain fake("hipo");
  // Adding the different input files to the TChain
  fake.Add(inputFile1.Data());
  // fake.Add(inputFile2.Data());
  // fake.Add(inputFile3.Data());
  // fake.Add(inputFile4.Data());
  // fake.Add(inputFile5.Data());
  // fake.Add(inputFile6.Data());


  // Access information from PDG e.g. particle masses
  auto db=TDatabasePDG::Instance();
  TLorentzVector beam(0,0,10.6,10.6); // beam Px,Py,Pz,E
<<<<<<< HEAD
  TLorentzVector target(0,0,0,1.8756); // target Px,Py,Pz,E
  // TLorentzVector target(0,0,0,db->GetParticle(2212)->Mass()); // target Px,Py,Pz,E
=======
  TLorentzVector target(0,0,0,db->GetParticle(2212)->Mass()); // target Px,Py,Pz,E
>>>>>>> 746002c24cbd1e4ca1e3da8f0b430e15d2e783da
  TLorentzVector el(0,0,0,db->GetParticle(11)->Mass()); // scattered e^- Px,Py,Pz,E
  TLorentzVector pr(0,0,0,db->GetParticle(2212)->Mass()); // proton Px,Py,Pz,E
  TLorentzVector pip(0,0,0,db->GetParticle(211)->Mass()); // pi^+ Px,Py,Pz,E

  Int_t part_pid;

  // x,y,z positions for 3 DC layers
  Double_t part_DC_c1x,part_DC_c1y,part_DC_c1z;
  Double_t part_DC_c2x,part_DC_c2y,part_DC_c2z;
  Double_t part_DC_c3x,part_DC_c3y,part_DC_c3z;

  // Retrieving list of files
  auto files=fake.GetListOfFiles();
  // Gets total events in all files for run dependence binning
  Int_t Bins = files->GetEntries();
  // Output file location and name
<<<<<<< HEAD
  TFile fileOutput1("/volatile/clas12/matthewn/FTOF/FTOF_Efficiency_RGB_SPING2019_dst_Inbending_2pi_misspim_binning_20072021_01.root","recreate");

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Create histograms here

=======
  TFile fileOutput1("/volatile/clas12/matthewn/FTOF/FTOF_Efficiency_RGA_FALL2018_skim4_Inbending_50_2pi_misspim_FTOF2_13072021_01.root","recreate");

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Create histograms here
>>>>>>> 746002c24cbd1e4ca1e3da8f0b430e15d2e783da
  // 2 pi event histograms
  auto* hmass=new TH1F("pimmass","Missing Mass e' p #pi^{+};MM(e'p#pi^{+}) [GeV];Counts",200,-1,1);
  auto* hmass_cuts=new TH1F("hmass_cuts","Missing Mass e' p #pi^{+};MM(e'p#pi^{+}) [GeV];Counts",200,-1,1);
  auto* hdeltaP=new TH1F("DeltaMomentum","Momentum difference of #pi^{-} detected and reconstructed;#Delta P [GeV];Counts",400,-2,2);
  auto* hdeltaTheta=new TH1F("DeltaTheta","#theta difference of #pi^{-} detected and reconstructed;#Delta #theta [deg];Counts",360,-180,180);
  auto* hdeltaPhi=new TH1F("DeltaPhi","#phi difference of #pi^{-} detected and reconstructed;#Delta #phi [deg];Counts",360,-180,180);

  // Checking the counter vs L
  auto* hlvscounter1A=new TH2F("hlvscounter1A","L vs counter FTOF1A;Counter;L [cm]",24,0,24,500,0,500);
  auto* hlvscounter1B=new TH2F("hlvscounter1B","L vs counter FTOF1B;Counter;L [cm]",63,0,63,500,0,500);
  auto* hlvscounter2=new TH2F("hlvscounter2","L vs counter FTOF2;Counter;L [cm]",6,0,6,300,600,900);

  // Look at momentum distribution of particles
  auto* hMomentum_1A_p=new TH1F("hMomentum_1A_p", "Momentum of positive tracks in FTOF1A; Momentum [GeV/c]; Counts;", 200,0,8);
  auto* hMomentum_1B_p=new TH1F("hMomentum_1B_p", "Momentum of positive tracks in FTOF1A; Momentum [GeV/c]; Counts;", 200,0,8);
  auto* hMomentum_2_p=new TH1F("hMomentum_2_p", "Momentum of positive tracks in FTOF1A; Momentum [GeV/c]; Counts;", 200,0,8);
  auto* hMomentum_1A_n=new TH1F("hMomentum_1A_n", "Momentum of negative tracks in FTOF1A; Momentum [GeV/c]; Counts;", 200,0,8);
  auto* hMomentum_1B_n=new TH1F("hMomentum_1B_n", "Momentum of negative tracks in FTOF1A; Momentum [GeV/c]; Counts;", 200,0,8);
  auto* hMomentum_2_n=new TH1F("hMomentum_2_n", "Momentum of negative tracks in FTOF1A; Momentum [GeV/c]; Counts;", 200,0,8);
  auto* hEnergy_1a=new TH1F("hEnergy_1a","Energy 1A",510,-10,500);
  auto* hEnergy_1b=new TH1F("hEnergy_1b","Energy 1A",510,-10,500);
  auto *h_z_vertex=new TH1D("h_z_vertex","z vertex",100,-15,15);
  auto *h_Pid=new TH1D("h_Pid","PID after cuts",10000,-5000,5000);

  // Histograms for timing
  auto* h_PCAL_time=new TH1F("h_PCAL_time","PCAL time;time [ns];Counts",400,0,200);
  auto* h_PCAL_beta=new TH2D("h_PCAL_beta","PCAL beta;P [GeV];beta",1100,0,11,400,-2,2);
  auto* h_TOF=new TH1F("h_TOF","PCAL TOF;time [ns];Counts",400,0,200);
  auto* h_beta=new TH2D("h_beta","beta;P [GeV];beta",1100,0,11,400,-2,2);
  auto* h_PCAL_path=new TH1F("h_PCAL_path","PCAL path;path [cm];Counts",1000,0,1000);
  auto* h_FTOF_time_1A=new TH1F("h_FTOF_time_1A","FTOF1A time;time [ns];Counts",400,0,200);
  auto* h_FTOF_time_1B=new TH1F("h_FTOF_time_1B","FTOF1A time;time [ns];Counts",400,0,200);
  auto* h_start_time=new TH1F("h_start_time","Start time;Start time [ns];Counts",400,0,200);

  // Arrays of histograms [layer][charge][sector]
  TH3F *h_Trajectories[3][2][6]; // Trajectories from DC
  TH3F *h_Tracks[3][2][6]; // Trajectories from DC with energy deposited in FTOF
  TH3F *h_Efficiency[3][2][6]; // Tracks divided by Trajectories

  // Looping over the FTOF layers
  for(Int_t i_detector=0;i_detector<3;i_detector++){
    // Looping over negative and positive particles
    for(Int_t i_charge=0;i_charge<2;i_charge++){
      // Looping over the sectors
      for(Int_t i_sector=0;i_sector<6;i_sector++){

        //create a string which we can append integers to, which allows us to define a number of histograms in a for loop
        // Histogram names
        ostringstream Traj_name_stream;
        ostringstream Tracks_name_stream;

        // Histogram Titles
        ostringstream Traj_title_stream;
        ostringstream Tracks_title_stream;

        // Defining the histogram name strings
        Traj_name_stream<<"h_Traj_Det_"<<i_detector<<"_Charge_"<<i_charge<<"_Sec_"<<i_sector;
        Tracks_name_stream<<"h_Tracks_Det_"<<i_detector<<"_Charge_"<<i_charge<<"_Sec_"<<i_sector;

        // Defining the histogram title strings
        if (i_detector==0){
          Traj_title_stream<<"Trajectories FTOF1A Charge "<<2*i_charge-1<<" Sec "<<i_sector+1<<"; Run no.; L [cm]";

          //convert the stringstream to a string and define our histograms in each element of the array
          h_Trajectories[i_detector][i_charge][i_sector] = new TH3F(Traj_name_stream.str().c_str(),"", Bins,0,Bins, 25, 54.51, 434.01, 50, -250, 250);
          h_Tracks[i_detector][i_charge][i_sector] = new TH3F(Tracks_name_stream.str().c_str(),"", Bins,0,Bins, 25, 54.51, 434.01, 50, -250 , 250);

        }
        else if (i_detector==1){
          Traj_title_stream<<"Trajectories FTOF2 Charge "<<2*i_charge-1<<" Sec "<<i_sector+1<<"; Run no.; L [cm]";

          //convert the stringstream to a string and define our histograms in each element of the array
          h_Trajectories[i_detector][i_charge][i_sector] = new TH3F(Traj_name_stream.str().c_str(),"", Bins,0,Bins, 65, 42.7, 440.5, 50, -250, 250);
          h_Tracks[i_detector][i_charge][i_sector] = new TH3F(Tracks_name_stream.str().c_str(),"", Bins,0,Bins, 65, 42.7, 440.5, 50, -250 , 250);
        }
        else if (i_detector==2){
          Traj_title_stream<<"Trajectories FTOF2 Charge "<<2*i_charge-1<<" Sec "<<i_sector+1<<"; Run no.; L [cm]";

          //convert the stringstream to a string and define our histograms in each element of the array
          h_Trajectories[i_detector][i_charge][i_sector] = new TH3F(Traj_name_stream.str().c_str(),"", Bins,0,Bins, 7, 693, 847, 50, -250, 250);
          h_Tracks[i_detector][i_charge][i_sector] = new TH3F(Tracks_name_stream.str().c_str(),"", Bins,0,Bins, 7, 693, 847, 50, -250 , 250);
        }


        // Setting the title for each histogram as it did not work when put in the lines above
        h_Trajectories[i_detector][i_charge][i_sector]->SetTitle(Traj_title_stream.str().c_str());
        h_Tracks[i_detector][i_charge][i_sector]->SetTitle(Tracks_title_stream.str().c_str());
      }
    }
  }

  // Distance between trajectory point and scintillator hit
  auto* h_radia_residual_1A = new TH1D("h_radia_residual_1A","Radius Residual FTOF1A",300,0,30);
  auto* h_radia_residual_1B = new TH1D("h_radia_residual_1B","Radius Residual FTOF1A",300,0,30);
  auto* h_radia_residual_2 = new TH1D("h_radia_residual_2","Radius Residual FTOF1A",300,0,30);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Creating variables for different layers of FTOF

  // Positions and angles
  Double_t x_1a, x_1b, x_2, y_1a, y_1b, y_2, z_1a, z_1b, z_2; // DC trajectory x,y,z positions
  Double_t d_1a, d_1b, d_2; // Distance to xy position (used for calculating perpendicular distance)
  Double_t x_FTOF, y_FTOF, z_FTOF; // Scintillator x,y,z positions
  Double_t L_1a, L_1b, L_2; // Perpendicular distance to sector in lab frame
  Double_t L_det_1a, L_det_1b, L_det_2; // Perpendicular distance to sector in detector frame

  Double_t alpha_1a, alpha_1b, alpha_2; // angle to hit (used to determine sector)
  Double_t L_Perp_1a, L_Perp_1b, L_Perp_2; // Distance along counter
  Double_t L_Theta; // Angle at middle of sector (check if hit is left or right of the middle of the sector)
  Double_t radia_residual; // Distance between trajectory bank and scintillator hit

  Double_t Component; // Getting specific counter hit in scintillator

  // Energy Values
  Double_t FTOF_1A_Energy, FTOF_1B_Energy, FTOF_1_Energy, FTOF_2_Energy; // FTOF energy depositions
  Double_t PCAL_Energy, ECIN_Energy, ECOUT_Energy, ECAL_Energy; // ECAL energy depositions

  // Timing and beta
  Double_t PCAL_Time = 0;
  Double_t PCAL_path = 0;
  Double_t PCAL_beta = 0;
  Double_t beta = 0;
  Double_t start_time = 0;
  Double_t TOF = 0;
  Double_t Momentum = 0;

  // Missing pim
  TLorentzVector misspim;
  // Negative particle set to pi^-
  TLorentzVector pim;

  // Run Number
  Int_t runno; // runno as a integer
  vector<TString> v_Runno; // runno as a vector of string
  Int_t Binno=0; // Count the number of files this corresponds to the number of x bins
  // Status
  Int_t Status, Calorimeter_Hits; // Status is used to find out if there is a ECAL hit

  // Particle numbers for 2pi events
  Int_t negative, positive, nonelectron;
  // Creating variables for comparing detected and missing pion
  Double_t DeltaP, DeltaTheta, DeltaPhi;
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // 1st Loop over data to find 2pi events

  //Loop over files
  for(Int_t i=0;i<files->GetEntries();i++){

    Binno++; // Count the number of files, therefore the number of x bins

    //create the event reader
    clas12reader c12(files->At(i)->GetTitle());

    while(c12.next()==true){
      // Used to access the particle bank for each event
      auto particles = c12.getDetParticles();

      // Setting variables to zero at start  of each event
      negative = 0;
      positive = 0;
      nonelectron = 0;
      DeltaP = 0;
      DeltaTheta = 0;
      DeltaPhi = 0;


      auto electrons=c12.getByID(11);
      auto protons=c12.getByID(2212);
      auto pips=c12.getByID(211);

      // Looping over all particles in this event
      for(auto& p : particles){

        // Looking at negative particles
        if(p->par()->getCharge() < 0){
          // negative particles not electron are set to pi^-
          if(p->par()->getPid() != 11){
            nonelectron++;
            pim.SetXYZM(p->par()->getPx(),p->par()->getPy(),p->par()->getPz(),db->GetParticle(-211)->Mass());
          }
        }
        // Looking at positive particles
        else if(p->par()->getCharge() > 0) positive++;
      }

      // Getting 2pi events
      if(nonelectron != 1 || electrons.size() != 1 || protons.size() != 1 || pips.size() != 1 || positive == 2)continue;
      {
        // cout<<pim.Rho()<<endl;
        // Apply cut on pion momentum to look at momentum dependence
        // if(pim.Rho()<2.5 /*|| pim.Rho()>2.5*/)continue;
        // cout<<"test "<<pim.Rho()<<endl;



        SetLorentzVector(el,electrons[0]);
        SetLorentzVector(pr,protons[0]);
        SetLorentzVector(pip,pips[0]);

        if(fabs(electrons[0]->par()->getChi2Pid())>3 || fabs(protons[0]->par()->getChi2Pid())>3 || fabs(pips[0]->par()->getChi2Pid())>3)continue;


        misspim = beam + target - el - pr - pip;
        hmass->Fill(misspim.M2());

        DeltaP = misspim.Rho() - pim.Rho();
        DeltaTheta = TMath::RadToDeg()* (misspim.Theta() - pim.Theta());
        DeltaPhi = TMath::RadToDeg()* (misspim.Phi() - pim.Phi());


        // Checking the background after cuts
        if(fabs(DeltaP) < 0.3 && fabs(DeltaTheta) < 10 && fabs(DeltaPhi) < 10 ){

          hmass_cuts->Fill(misspim.M2());
        }

        // Cut on missing mass of the pi^-
        if(misspim.M2() > - 0.1 && misspim.M2() < 0.2){


          // Plotting pi^- variables
          hdeltaP->Fill(DeltaP);
          hdeltaTheta->Fill(DeltaTheta);
          hdeltaPhi->Fill(DeltaPhi);

          // Other cuts to be applied and tuned
          if(fabs(DeltaP) < 0.3 && fabs(DeltaTheta) < 10 && fabs(DeltaPhi) < 10 ){

            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // 2nd Loop over data to get detector information

            //Set the particle index to 0 and loop through the particles
            int pindex=0;
            for(auto& p : particles){

              //get information from the different detectors
              switch(p->getRegion()) {
                // Only look at particles in FD as FTOF is only in FD
                case FD :

                // Increase the particle index with each loop of the particles
                pindex++;

                //Ignore the first particle (trigger), any neutrals and electron
                if (pindex==1 || p->par()->getCharge()==0 || p->par()->getPid()==11) continue;

                runno = c12.runconfig()->getRun(); // Getting the run number
                Status = p->par()->getStatus(); // Getting the status

                // Checking the z vertex before applying a cut
                h_z_vertex->Fill(p->par()->getVz());
                if(p->par()->getVz() > 2 || p->par()->getVz() < -9)continue;

                // Setting all energy values to zero
                FTOF_1A_Energy = 0;
                FTOF_1B_Energy = 0;
                FTOF_2_Energy = 0;
                PCAL_Energy = 0;
                ECIN_Energy = 0;
                ECOUT_Energy = 0;

                // Setting the energy values to a double
                if(p->sci(FTOF1A)->getEnergy()) FTOF_1A_Energy = p->sci(FTOF1A)->getEnergy();
                if(p->sci(FTOF1B)->getEnergy()) FTOF_1B_Energy = p->sci(FTOF1B)->getEnergy();
                if(p->sci(FTOF2)->getEnergy()) FTOF_2_Energy = p->sci(FTOF2)->getEnergy();
                if(p->cal(PCAL)->getEnergy()) PCAL_Energy = p->cal(PCAL)->getEnergy();
                if(p->cal(ECIN)->getEnergy()) ECIN_Energy = p->cal(ECIN)->getEnergy();
                if(p->cal(ECOUT)->getEnergy()) ECOUT_Energy = p->cal(ECOUT)->getEnergy();

                ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                // DC fiducial cuts

                // if(( p->traj(DC,6)->getLayer() != 6 || p->traj(DC,18)->getLayer() != 18 || p->traj(DC,36)->getLayer() != 36))continue;

                // if(p->traj(DC,6)->getDetector() == 6 && p->traj(DC,6)->getLayer() == 6){
                //   part_DC_c1x = p->traj(DC,6)->getX();
                //   part_DC_c1y = p->traj(DC,6)->getY();
                //   part_DC_c1z = p->traj(DC,6)->getZ();
                // }
                // if(p->traj(DC,18)->getDetector() == 6 && p->traj(DC,18)->getLayer() == 18){
                //   part_DC_c2x = p->traj(DC,18)->getX();
                //   part_DC_c2y = p->traj(DC,18)->getY();
                //   part_DC_c2z = p->traj(DC,18)->getZ();
                // }
                // if(p->traj(DC,36)->getDetector() == 6 && p->traj(DC,36)->getLayer() == 36){
                //   part_DC_c3x = p->traj(DC,36)->getX();
                //   part_DC_c3y = p->traj(DC,36)->getY();
                //   part_DC_c3z = p->traj(DC,36)->getZ();
                // }
                //
                // part_pid = 0;
                // // Positive particles asssigned to pi^+, negative to pi^-
                //
                // if(p->par()->getPid())
                // part_pid = p->par()->getPid();
                //
                // else if(p->par()->getCharge()>0) part_pid = 2;
                // else if(p->par()->getCharge()<0) part_pid = 3;
                //
                // int part_DC_sector = determineSector(part_DC_c2x, part_DC_c2y,part_DC_c2z);

                // Use this cut if looking at inbending data
                // if(!DC_fiducial_cut_theta_phi(part_DC_c1x,part_DC_c1y,part_DC_c1z,part_DC_c2x, part_DC_c2y, part_DC_c2z,part_DC_c3x,part_DC_c3y, part_DC_c3z,part_pid,part_DC_sector,1) ||
                // !DC_fiducial_cut_theta_phi(part_DC_c1x,part_DC_c1y,part_DC_c1z,part_DC_c2x, part_DC_c2y, part_DC_c2z,part_DC_c3x,part_DC_c3y, part_DC_c3z,part_pid,part_DC_sector,2) ||
                // !DC_fiducial_cut_theta_phi(part_DC_c1x,part_DC_c1y,part_DC_c1z,part_DC_c2x, part_DC_c2y, part_DC_c2z,part_DC_c3x, part_DC_c3y, part_DC_c3z,part_pid,part_DC_sector,3))continue;

                // // Use this cut if looking at outbending data
                // if(!DC_fiducial_cut_XY(part_DC_c1x,part_DC_c1y,part_DC_c1z,part_DC_c2x, part_DC_c2y, part_DC_c2z,part_DC_c3x, part_DC_c3y, part_DC_c3z,part_pid,part_DC_sector,1) ||
                // !DC_fiducial_cut_XY(part_DC_c1x,part_DC_c1y,part_DC_c1z,part_DC_c2x, part_DC_c2y, part_DC_c2z,part_DC_c3x, part_DC_c3y, part_DC_c3z, part_pid,part_DC_sector,2) ||
                // !DC_fiducial_cut_XY(part_DC_c1x,part_DC_c1y,part_DC_c1z,part_DC_c2x, part_DC_c2y, part_DC_c2z,part_DC_c3x, part_DC_c3y, part_DC_c3z, part_pid,part_DC_sector,3)) continue;


                ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


                // Skips any particle with no hits in the calorimeter (only applied to 1A and 2 due to angular coverage)
                // if(PCAL_Energy > 0 && ECIN_Energy > 0){

                FTOF_1_Energy = FTOF_1A_Energy + FTOF_1B_Energy;
                ECAL_Energy = PCAL_Energy + ECIN_Energy + ECOUT_Energy;

                // Defining momentum and beta of particles
                Momentum = p->par()->getP();
                beta = p->par()->getBeta();
                start_time = c12.event()->getStartTime();
                h_beta->Fill(Momentum,p->par()->getBeta());

                // Removing particles with unphysical beta
                // if(beta < 0.4 || beta >1.1)continue;

                // PCAL and start time information
                if(p->cal(PCAL)->getEnergy()){
                  PCAL_Time = p->cal(PCAL)->getTime();
                  h_PCAL_time->Fill(PCAL_Time);
                  TOF = PCAL_Time - start_time;
                  h_TOF->Fill(TOF);
                  PCAL_path = p->cal(PCAL)->getPath();
                  PCAL_beta = (PCAL_path / TOF) / 30.0;
                  h_PCAL_beta->Fill(Momentum,PCAL_beta);
                  h_PCAL_path->Fill(PCAL_path);
                }

                h_start_time->Fill(c12.event()->getStartTime());

                h_Pid->Fill(p->par()->getPid());

                ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                //FTOF1A

                // Requiring a hit in FTOF1B
                // if(p->sci(FTOF1B)->getEnergy()>0){

                if(p->traj(FTOF,FTOF1A)->getDetector()==12 && p->traj(FTOF,FTOF1A)->getLayer()==1){
                  if(p->par()->getCharge()>0) hMomentum_1A_p->Fill(p->par()->getP());
                  else if(p->par()->getCharge()<0) hMomentum_1A_n->Fill(p->par()->getP());

                  // Getting the x-, y- and z- co-ordinates from DC
                  x_1a =  p->traj(FTOF, FTOF1A)->getX();
                  y_1a =  p->traj(FTOF, FTOF1A)->getY();
                  z_1a = p->traj(FTOF, FTOF1A)->getZ();

                  // Getting x-, y- and z- co-ordinates from FTOF hit
                  if(p->sci(FTOF1A)->getEnergy()){
                    x_FTOF =  p->sci(FTOF1A)->getX();
                    y_FTOF =  p->sci(FTOF1A)->getY();
                    z_FTOF =  p->sci(FTOF1A)->getZ();

                    // FTOF time
                    h_FTOF_time_1A->Fill(p->sci(FTOF1A)->getTime()-start_time);

                    // if(p->sci(FTOF1A)->getTime()-start_time < 20 || p->sci(FTOF1A)->getTime()-start_time > 40)continue;

                    // Distance between DC and FTOF co-ordinates
                    radia_residual = sqrt(pow(x_1a-x_FTOF,2) + pow(y_1a-y_FTOF,2) + pow(z_1a-z_FTOF,2));
                    h_radia_residual_1A->Fill(radia_residual);
                  }

                  // Calculating d, distance to hit from (0,0) to (x,y)
                  d_1a = sqrt(pow(x_1a,2) + pow(y_1a,2));

                  // Calculating alpha, angle from (0,0) to hit
                  alpha_1a = TMath::RadToDeg()*atan(y_1a/x_1a);

                  // Calculating length from (0,0) along centre of sector
                  // For sectors 1 and 4, L is just the x value as their centre is the x-axis
                  if(fabs(alpha_1a) < 30) L_1a = fabs(x_1a);
                  else L_1a = (d_1a * cos(TMath::DegToRad()*(fabs(fabs(alpha_1a)-60))));

                  // Getting L in detector plane by accounting for tilt of FTOF1A (25 deg)
                  L_det_1a = L_1a / cos(TMath::DegToRad()*25);
                  // Calculating the distance along the width of the counter
                  L_Perp_1a = pow(pow(d_1a,2) - pow(L_1a,2), 0.5);

                  // Resetting the component to unphysical value
                  Component = -10;

                  // Checking if there is energy deposition in FTOF1A
                  if(p->sci(FTOF1A)->getEnergy()){
                    // Getting the component number for the counter hit
                    Component = p->sci(FTOF1A)->getComponent();
                    // Plotting L vs the counter number
                    hlvscounter1A->Fill(Component,L_det_1a);
                  }

                  // Determine which sector the hit is in using alpha
                  // Look at positive x (sectors 1,2 and 6)
                  if(x_1a > 0){
                    // Sector 2
                    if(alpha_1a > 30){
                      // Checking to see which side of the middle the hit is
                      L_Theta = 60;
                      if(L_Theta - alpha_1a < 0)L_Perp_1a = - L_Perp_1a;

                      // Positive particles
                      if(p->par()->getCharge()>0){
                        h_Trajectories[0][1][1]->Fill(i,L_det_1a, L_Perp_1a);
                      }

                      // Negative particles
                      else if(p->par()->getCharge()<0){
                        h_Trajectories[0][0][1]->Fill(i,L_det_1a, L_Perp_1a);
                      }
                      if(p->sci(FTOF1A)->getEnergy())hEnergy_1a->Fill(p->sci(FTOF1A)->getEnergy());

                      // Check if there is energy deposited on the scintillator
                      if(p->sci(FTOF1A)->getEnergy()){
                        // Positive particles
                        if(p->par()->getCharge()>0){
                          h_Tracks[0][1][1]->Fill(i,L_det_1a, L_Perp_1a);
                        }

                        // Negative particles
                        else if(p->par()->getCharge()<0){
                          h_Tracks[0][0][1]->Fill(i,L_det_1a, L_Perp_1a);
                        }
                      }
                    }

                    // Sector 1
                    else if(fabs(alpha_1a) < 30) {

                      if(alpha_1a < 0)L_Perp_1a = - L_Perp_1a;

                      // Positive particles
                      if(p->par()->getCharge()>0){
                        h_Trajectories[0][1][0]->Fill(i,L_det_1a, L_Perp_1a);
                      }

                      // Negative particles
                      else if(p->par()->getCharge()<0){
                        h_Trajectories[0][0][0]->Fill(i,L_det_1a, L_Perp_1a);
                      }

                      // Check if there is energy deposited on the scintillator
                      if(p->sci(FTOF1A)->getEnergy()){
                        if(p->par()->getCharge()>0){
                          h_Tracks[0][1][0]->Fill(i,L_det_1a, L_Perp_1a);
                        }

                        else if(p->par()->getCharge()<0){
                          h_Tracks[0][0][0]->Fill(i,L_det_1a, L_Perp_1a);
                        }
                      }
                    }

                    // Sector 6

                    else if(alpha_1a < -30) {

                      L_Theta = 60;
                      if(L_Theta + alpha_1a < 0)L_Perp_1a = - L_Perp_1a;

                      // Positive particles
                      if(p->par()->getCharge()>0){
                        h_Trajectories[0][1][5]->Fill(i,L_det_1a, L_Perp_1a);
                      }

                      // Negative particles
                      else if(p->par()->getCharge()<0){
                        h_Trajectories[0][0][5]->Fill(i,L_det_1a, L_Perp_1a);
                      }

                      // Check if there is energy deposited on the scintillator
                      if(p->sci(FTOF1A)->getEnergy()){
                        if(p->par()->getCharge()>0){
                          h_Tracks[0][1][5]->Fill(i,L_det_1a, L_Perp_1a);
                        }

                        else if(p->par()->getCharge()<0){
                          h_Tracks[0][0][5]->Fill(i,L_det_1a, L_Perp_1a);
                        }
                      }
                    }
                  }

                  // Look at negative x (sectors 3,4 and 5)
                  else if(x_1a < 0){

                    // Sector 5

                    if(alpha_1a > 30){

                      L_Theta = 60;
                      if(L_Theta - alpha_1a < 0)L_Perp_1a = - L_Perp_1a;

                      // Positive particles
                      if(p->par()->getCharge()>0){
                        h_Trajectories[0][1][4]->Fill(i,L_det_1a, L_Perp_1a);
                      }

                      // Negative particles
                      else if(p->par()->getCharge()<0){
                        h_Trajectories[0][0][4]->Fill(i,L_det_1a, L_Perp_1a);
                      }

                      // Check if there is energy deposited on the scintillator
                      if(p->sci(FTOF1A)->getEnergy()){
                        if(p->par()->getCharge()>0){
                          h_Tracks[0][1][4]->Fill(i,L_det_1a, L_Perp_1a);
                        }

                        else if(p->par()->getCharge()<0){
                          h_Tracks[0][0][4]->Fill(i,L_det_1a, L_Perp_1a);
                        }
                      }
                    }

                    // Sector 4
                    else if(fabs(alpha_1a) < 30) {

                      if(alpha_1a < 0)L_Perp_1a = - L_Perp_1a;

                      // Positive particles
                      if(p->par()->getCharge()>0){
                        h_Trajectories[0][1][3]->Fill(i,L_det_1a, L_Perp_1a);
                      }

                      // Negative particles
                      else if(p->par()->getCharge()<0){
                        h_Trajectories[0][0][3]->Fill(i,L_det_1a, L_Perp_1a);
                      }

                      // Check if there is energy deposited on the scintillator
                      if(p->sci(FTOF1A)->getEnergy()){
                        if(p->par()->getCharge()>0){
                          h_Tracks[0][1][3]->Fill(i,L_det_1a, L_Perp_1a);
                        }

                        else if(p->par()->getCharge()<0){
                          h_Tracks[0][0][3]->Fill(i,L_det_1a, L_Perp_1a);
                        }
                      }
                    }

                    // Sector 3
                    else if(alpha_1a < -30) {

                      L_Theta = 60;
                      if(L_Theta + alpha_1a < 0)L_Perp_1a = - L_Perp_1a;

                      // Positive particles
                      if(p->par()->getCharge()>0){
                        h_Trajectories[0][1][2]->Fill(i,L_det_1a, L_Perp_1a);
                      }

                      // Negative particles
                      else if(p->par()->getCharge()<0){
                        h_Trajectories[0][0][2]->Fill(i,L_det_1a, L_Perp_1a);
                      }

                      // Check if there is energy deposited on the scintillator
                      if(p->sci(FTOF1A)->getEnergy()){
                        if(p->par()->getCharge()>0){
                          h_Tracks[0][1][2]->Fill(i,L_det_1a, L_Perp_1a);
                        }

                        else if(p->par()->getCharge()<0){
                          h_Tracks[0][0][2]->Fill(i,L_det_1a, L_Perp_1a);
                        }
                      }
                    }
                  }
                }
                // }
                ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                //FTOF1B

                // Requiring a hit in FTOF1A
                // if(p->sci(FTOF1A)->getEnergy()>0){

                if(p->traj(FTOF,FTOF1B)->getDetector()==12 && p->traj(FTOF,FTOF1B)->getLayer()==2){
                  if(p->par()->getCharge()>0) hMomentum_1B_p->Fill(p->par()->getP());
                  else if(p->par()->getCharge()<0) hMomentum_1B_n->Fill(p->par()->getP());

                  // Getting the x-, y- and z- co-ordinates from DC
                  x_1b =  p->traj(FTOF, FTOF1B)->getX();
                  y_1b =  p->traj(FTOF, FTOF1B)->getY();
                  z_1b = p->traj(FTOF, FTOF1B)->getZ();

                  // Getting x-, y- and z- co-ordinates from FTOF hit
                  if(p->sci(FTOF1B)->getEnergy()){
                    x_FTOF =  p->sci(FTOF1B)->getX();
                    y_FTOF =  p->sci(FTOF1B)->getY();
                    z_FTOF =  p->sci(FTOF1B)->getZ();

                    // FTOF time
                    h_FTOF_time_1B->Fill(p->sci(FTOF1B)->getTime() - start_time);

                    // if(p->sci(FTOF1B)->getTime()-start_time < 20 || p->sci(FTOF1B)->getTime()-start_time > 40)continue;

                    // Distance between DC and FTOF co-ordinates
                    radia_residual = sqrt(pow(x_1b-x_FTOF,2) + pow(y_1b-y_FTOF,2) + pow(z_1b-z_FTOF,2));
                    h_radia_residual_1B->Fill(radia_residual);
                  }

                  // Calculating d, distance to hit from (0,0) to (x,y)
                  d_1b = sqrt(pow(x_1b,2) + pow(y_1b,2));

                  // Calculating alpha, angle from (0,0) to hit
                  alpha_1b = TMath::RadToDeg()*atan(y_1b/x_1b);

                  // Calculating length from (0,0) along centre of sector
                  // For sectors 1 and 4, L is just the x value as their centre is the x-axis
                  if(fabs(alpha_1b) < 30) L_1b = fabs(x_1b);
                  else L_1b = (d_1b * cos(TMath::DegToRad()*(fabs(fabs(alpha_1b)-60))));

                  // Getting L in detector plane by accounting for tilt of FTOF1B (25 deg)
                  L_det_1b = L_1b / cos(TMath::DegToRad()*25);
                  // Calculating the distance along the width of the counter
                  L_Perp_1b = pow(pow(d_1b,2) - pow(L_1b,2), 0.5);

                  // Resetting the component to unphysical value
                  Component = -10;

                  // Checking if there is energy deposition in FTOF1B
                  if(p->sci(FTOF1B)->getEnergy()){
                    // Getting the component number for the counter hit
                    Component = p->sci(FTOF1B)->getComponent();
                    // Plotting L vs the counter number
                    hlvscounter1B->Fill(Component,L_det_1b);
                  }

                  // Determine which sector the hit is in using alpha
                  // Look at positive x (sectors 1,2 and 6)
                  if(x_1b > 0){
                    // Sector 2
                    if(alpha_1b > 30){
                      // Checking to see which side of the middle the hit is
                      L_Theta = 60;
                      if(L_Theta - alpha_1b < 0)L_Perp_1b = - L_Perp_1b;

                      // Positive particles
                      if(p->par()->getCharge()>0){
                        h_Trajectories[1][1][1]->Fill(i,L_det_1b, L_Perp_1b);
                      }

                      // Negative particles
                      else if(p->par()->getCharge()<0){
                        h_Trajectories[1][0][1]->Fill(i,L_det_1b, L_Perp_1b);
                      }
                      // Check if there is energy deposited on the scintillator
                      if(p->sci(FTOF1B)->getEnergy()){
                        // Positive particles
                        if(p->par()->getCharge()>0){
                          h_Tracks[1][1][1]->Fill(i,L_det_1b, L_Perp_1b);
                        }

                        // Negative particles
                        else if(p->par()->getCharge()<0){
                          h_Tracks[1][0][1]->Fill(i,L_det_1b, L_Perp_1b);
                        }
                      }
                    }

                    // Sector 1
                    else if(fabs(alpha_1b) < 30) {

                      if(alpha_1b < 0)L_Perp_1b = - L_Perp_1b;

                      // Positive particles
                      if(p->par()->getCharge()>0){
                        h_Trajectories[1][1][0]->Fill(i,L_det_1b, L_Perp_1b);
                      }

                      // Negative particles
                      else if(p->par()->getCharge()<0){
                        h_Trajectories[1][0][0]->Fill(i,L_det_1b, L_Perp_1b);
                      }

                      // Check if there is energy deposited on the scintillator
                      if(p->sci(FTOF1B)->getEnergy()){
                        if(p->par()->getCharge()>0){
                          h_Tracks[1][1][0]->Fill(i,L_det_1b, L_Perp_1b);
                        }

                        else if(p->par()->getCharge()<0){
                          h_Tracks[1][0][0]->Fill(i,L_det_1b, L_Perp_1b);
                        }
                      }
                    }

                    // Sector 6

                    else if(alpha_1b < -30) {

                      L_Theta = 60;
                      if(L_Theta + alpha_1b < 0)L_Perp_1b = - L_Perp_1b;

                      // Positive particles
                      if(p->par()->getCharge()>0){
                        h_Trajectories[1][1][5]->Fill(i,L_det_1b, L_Perp_1b);
                      }

                      // Negative particles
                      else if(p->par()->getCharge()<0){
                        h_Trajectories[1][0][5]->Fill(i,L_det_1b, L_Perp_1b);
                      }

                      // Check if there is energy deposited on the scintillator
                      if(p->sci(FTOF1B)->getEnergy()){
                        if(p->par()->getCharge()>0){
                          h_Tracks[1][1][5]->Fill(i,L_det_1b, L_Perp_1b);
                        }

                        else if(p->par()->getCharge()<0){
                          h_Tracks[1][0][5]->Fill(i,L_det_1b, L_Perp_1b);
                        }
                      }
                    }
                  }

                  // Look at negative x (sectors 3,4 and 5)
                  else if(x_1b < 0){

                    // Sector 5

                    if(alpha_1b > 30){

                      L_Theta = 60;
                      if(L_Theta - alpha_1b < 0)L_Perp_1b = - L_Perp_1b;

                      // Positive particles
                      if(p->par()->getCharge()>0){
                        h_Trajectories[1][1][4]->Fill(i,L_det_1b, L_Perp_1b);
                      }

                      // Negative particles
                      else if(p->par()->getCharge()<0){
                        h_Trajectories[1][0][4]->Fill(i,L_det_1b, L_Perp_1b);
                      }

                      // Check if there is energy deposited on the scintillator
                      if(p->sci(FTOF1B)->getEnergy()){
                        if(p->par()->getCharge()>0){
                          h_Tracks[1][1][4]->Fill(i,L_det_1b, L_Perp_1b);
                        }

                        else if(p->par()->getCharge()<0){
                          h_Tracks[1][0][4]->Fill(i,L_det_1b, L_Perp_1b);
                        }
                      }
                    }

                    // Sector 4
                    else if(fabs(alpha_1b) < 30) {

                      if(alpha_1b < 0)L_Perp_1b = - L_Perp_1b;

                      // Positive particles
                      if(p->par()->getCharge()>0){
                        h_Trajectories[1][1][3]->Fill(i,L_det_1b, L_Perp_1b);
                      }

                      // Negative particles
                      else if(p->par()->getCharge()<0){
                        h_Trajectories[1][0][3]->Fill(i,L_det_1b, L_Perp_1b);
                      }

                      // Check if there is energy deposited on the scintillator
                      if(p->sci(FTOF1B)->getEnergy()){
                        if(p->par()->getCharge()>0){
                          h_Tracks[1][1][3]->Fill(i,L_det_1b, L_Perp_1b);
                        }

                        else if(p->par()->getCharge()<0){
                          h_Tracks[1][0][3]->Fill(i,L_det_1b, L_Perp_1b);
                        }
                      }
                    }

                    // Sector 3
                    else if(alpha_1b < -30) {

                      L_Theta = 60;
                      if(L_Theta + alpha_1b < 0)L_Perp_1b = - L_Perp_1b;

                      // Positive particles
                      if(p->par()->getCharge()>0){
                        h_Trajectories[1][1][2]->Fill(i,L_det_1b, L_Perp_1b);
                      }

                      // Negative particles
                      else if(p->par()->getCharge()<0){
                        h_Trajectories[1][0][2]->Fill(i,L_det_1b, L_Perp_1b);
                      }

                      // Check if there is energy deposited on the scintillator
                      if(p->sci(FTOF1B)->getEnergy()){
                        if(p->par()->getCharge()>0){
                          h_Tracks[1][1][2]->Fill(i,L_det_1b, L_Perp_1b);
                        }

                        else if(p->par()->getCharge()<0){
                          h_Tracks[1][0][2]->Fill(i,L_det_1b, L_Perp_1b);
                        }
                      }
                    }
                  }
                }
                // }
                // }
                ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                //FTOF2
                if(p->traj(FTOF,FTOF2)->getDetector()==12 && p->traj(FTOF,FTOF2)->getLayer()==3){
                  if(p->par()->getCharge()>0) hMomentum_2_p->Fill(p->par()->getP());
                  else if(p->par()->getCharge()<0) hMomentum_2_n->Fill(p->par()->getP());


                  // Getting the x-, y- and z- co-ordinates from DC
                  x_2 =  p->traj(FTOF, FTOF2)->getX();
                  y_2 =  p->traj(FTOF, FTOF2)->getY();
                  z_2 = p->traj(FTOF, FTOF2)->getZ();

                  // Getting x-, y- and z- co-ordinates from FTOF hit
                  if(p->sci(FTOF2)->getEnergy()){
                    x_FTOF =  p->sci(FTOF2)->getX();
                    y_FTOF =  p->sci(FTOF2)->getY();
                    z_FTOF =  p->sci(FTOF2)->getZ();

                    // Distance between DC and FTOF co-ordinates
                    radia_residual = sqrt(pow(x_2-x_FTOF,2) + pow(y_2-y_FTOF,2) + pow(z_2-z_FTOF,2));
                    h_radia_residual_2->Fill(radia_residual);
                  }

                  // Calculating d, distance to hit from (0,0) to (x,y)
                  d_2 = sqrt(pow(x_2,2) + pow(y_2,2));

                  // Calculating alpha, angle from (0,0) to hit
                  alpha_2 = TMath::RadToDeg()*atan(y_2/x_2);

                  // Calculating length from (0,0) along centre of sector
                  // For sectors 1 and 4, L is just the x value as their centre is the x-axis
                  if(fabs(alpha_2) < 30) L_2 = fabs(x_2);
                  else L_2 = (d_2 * cos(TMath::DegToRad()*(fabs(fabs(alpha_2)-60))));

                  // Getting L in detector plane by accounting for tilt of FTOF2 (58.11 deg)
                  L_det_2 = L_2 / cos(TMath::DegToRad()*58.11);
                  // Calculating the distance along the width of the counter
                  L_Perp_2 = pow(pow(d_2,2) - pow(L_2,2), 0.5);


                  // Resetting the component to unphysical value
                  Component = -10;

                  // Checking if there is energy deposition in FTOF2
                  if(p->sci(FTOF2)->getEnergy()){
                    // Getting the component number for the counter hit
                    Component = p->sci(FTOF2)->getComponent();
                    // Plotting L vs the counter number
                    hlvscounter2->Fill(Component,L_det_2);
                  }


                  // Determine which sector the hit is in using alpha
                  // Look at positive x (sectors 1,2 and 6)
                  if(x_2 > 0){
                    // Sector 2
                    if(alpha_2 > 30){
                      // Checking to see which side of the middle the hit is
                      L_Theta = 60;
                      if(L_Theta - alpha_2 < 0)L_Perp_2 = - L_Perp_2;

                      // Positive particles
                      if(p->par()->getCharge()>0){
                        h_Trajectories[2][1][1]->Fill(i,L_det_2, L_Perp_2);
                      }

                      // Negative particles
                      else if(p->par()->getCharge()<0){
                        h_Trajectories[2][0][1]->Fill(i,L_det_2, L_Perp_2);
                      }
                      // Check if there is energy deposited on the scintillator
                      if(p->sci(FTOF2)->getEnergy()){
                        // Positive particles
                        if(p->par()->getCharge()>0){
                          h_Tracks[2][1][1]->Fill(i,L_det_2, L_Perp_2);
                        }

                        // Negative particles
                        else if(p->par()->getCharge()<0){
                          h_Tracks[2][0][1]->Fill(i,L_det_2, L_Perp_2);
                        }
                      }
                    }

                    // Sector 1
                    else if(fabs(alpha_2) < 30) {

                      if(alpha_2 < 0)L_Perp_2 = - L_Perp_2;

                      // Positive particles
                      if(p->par()->getCharge()>0){
                        h_Trajectories[2][1][0]->Fill(i,L_det_2, L_Perp_2);
                      }

                      // Negative particles
                      else if(p->par()->getCharge()<0){
                        h_Trajectories[2][0][0]->Fill(i,L_det_2, L_Perp_2);
                      }

                      // Check if there is energy deposited on the scintillator
                      if(p->sci(FTOF2)->getEnergy()){
                        if(p->par()->getCharge()>0){
                          h_Tracks[2][1][0]->Fill(i,L_det_2, L_Perp_2);
                        }

                        else if(p->par()->getCharge()<0){
                          h_Tracks[2][0][0]->Fill(i,L_det_2, L_Perp_2);
                        }
                      }
                    }

                    // Sector 6

                    else if(alpha_2 < -30) {

                      L_Theta = 60;
                      if(L_Theta + alpha_2 < 0)L_Perp_2 = - L_Perp_2;

                      // Positive particles
                      if(p->par()->getCharge()>0){
                        h_Trajectories[2][1][5]->Fill(i,L_det_2, L_Perp_2);
                      }

                      // Negative particles
                      else if(p->par()->getCharge()<0){
                        h_Trajectories[2][0][5]->Fill(i,L_det_2, L_Perp_2);
                      }

                      // Check if there is energy deposited on the scintillator
                      if(p->sci(FTOF2)->getEnergy()){
                        if(p->par()->getCharge()>0){
                          h_Tracks[2][1][5]->Fill(i,L_det_2, L_Perp_2);
                        }

                        else if(p->par()->getCharge()<0){
                          h_Tracks[2][0][5]->Fill(i,L_det_2, L_Perp_2);
                        }
                      }
                    }
                  }

                  // Look at negative x (sectors 3,4 and 5)
                  else if(x_2 < 0){

                    // Sector 5

                    if(alpha_2 > 30){

                      L_Theta = 60;
                      if(L_Theta - alpha_2 < 0)L_Perp_2 = - L_Perp_2;

                      // Positive particles
                      if(p->par()->getCharge()>0){
                        h_Trajectories[2][1][4]->Fill(i,L_det_2, L_Perp_2);
                      }

                      // Negative particles
                      else if(p->par()->getCharge()<0){
                        h_Trajectories[2][0][4]->Fill(i,L_det_2, L_Perp_2);
                      }

                      // Check if there is energy deposited on the scintillator
                      if(p->sci(FTOF2)->getEnergy()){
                        if(p->par()->getCharge()>0){
                          h_Tracks[2][1][4]->Fill(i,L_det_2, L_Perp_2);
                        }

                        else if(p->par()->getCharge()<0){
                          h_Tracks[2][0][4]->Fill(i,L_det_2, L_Perp_2);
                        }
                      }
                    }

                    // Sector 4
                    else if(fabs(alpha_2) < 30) {

                      if(alpha_2 < 0)L_Perp_2 = - L_Perp_2;

                      // Positive particles
                      if(p->par()->getCharge()>0){
                        h_Trajectories[2][1][3]->Fill(i,L_det_2, L_Perp_2);
                      }

                      // Negative particles
                      else if(p->par()->getCharge()<0){
                        h_Trajectories[2][0][3]->Fill(i,L_det_2, L_Perp_2);
                      }

                      // Check if there is energy deposited on the scintillator
                      if(p->sci(FTOF2)->getEnergy()){
                        if(p->par()->getCharge()>0){
                          h_Tracks[2][1][3]->Fill(i,L_det_2, L_Perp_2);
                        }

                        else if(p->par()->getCharge()<0){
                          h_Tracks[2][0][3]->Fill(i,L_det_2, L_Perp_2);
                        }
                      }
                    }

                    // Sector 3
                    else if(alpha_2 < -30) {

                      L_Theta = 60;
                      if(L_Theta + alpha_2 < 0)L_Perp_2 = - L_Perp_2;

                      // Positive particles
                      if(p->par()->getCharge()>0){
                        h_Trajectories[2][1][2]->Fill(i,L_det_2, L_Perp_2);
                      }

                      // Negative particles
                      else if(p->par()->getCharge()<0){
                        h_Trajectories[2][0][2]->Fill(i,L_det_2, L_Perp_2);
                      }

                      // Check if there is energy deposited on the scintillator
                      if(p->sci(FTOF2)->getEnergy()){
                        if(p->par()->getCharge()>0){
                          h_Tracks[2][1][2]->Fill(i,L_det_2, L_Perp_2);
                        }

                        else if(p->par()->getCharge()<0){
                          h_Tracks[2][0][2]->Fill(i,L_det_2, L_Perp_2);
                        }
                      }
                    }
                  }
                }
                // }
              }
            }
          }
        }
      }
    }

    v_Runno.push_back(to_string(runno)); // Converting runno integer to a string
  }

  //saving the file
  fileOutput1.Write();

}
