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

void DC_Fiducial_Cuts_Test(){

  // Data files to process
  TString inputFile("/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass1/v0/dst/train/skim4/skim4_005032.hipo");
  // TFile fileOutput1("/u/home/matthewn/Documents/Service_Work/Sectors/L_Dependence/FTOF_Efficiency_RGA_RGB_FULL_RUN_CALORIMETER_HITS_ON_231120_02.root","recreate");

  bool inbending = true;
  bool outbending = false;

  Int_t part_pid;

  auto* h1xy=new TH2D("h1xy","",1000,-500,500,1000,-500,500);
  auto* h2xy=new TH2D("h2xy","",1000,-500,500,1000,-500,500);
  auto* h3xy=new TH2D("h3xy","",1000,-500,500,1000,-500,500);
  auto* h1xy_after=new TH2D("h1xy_after","",1000,-500,500,1000,-500,500);
  auto* h2xy_after=new TH2D("h2xy_after","",1000,-500,500,1000,-500,500);
  auto* h3xy_after=new TH2D("h3xy_after","",1000,-500,500,1000,-500,500);

  TChain fake("hipo");
  fake.Add(inputFile.Data());

  Double_t part_DC_c1x,part_DC_c1y,part_DC_c1z;
  Double_t part_DC_c2x,part_DC_c2y,part_DC_c2z;
  Double_t part_DC_c3x,part_DC_c3y,part_DC_c3z;

  auto files=fake.GetListOfFiles();

  //Loop over files
  for(Int_t i=0;i<files->GetEntries();i++){
    //create the event reader
    clas12reader c12(files->At(i)->GetTitle());

    while(c12.next()==true){
      auto particles = c12.getDetParticles();




      int pindex=0;
      for(auto& p : particles){


        cout<<p->traj(CTOF,0)->getX()<<" "<<p->traj(CTOF,1)->getX()<<" "<<p->traj(CTOF,2)->getX()<<endl;
        // cout<<p->traj(CTOF,0)->getX()<<endl;
        // cout<<p->traj(CTOF,2)->getX()<<endl;
        // cout<<p->traj(CTOF,CTOF0)->getX()<<endl;
        // cout<<p->traj(CTOF,CTOF1)->getX()<<endl;
        // cout<<p->traj(CTOF,CTOF2)->getX()<<endl;

        if(p->traj(DC,6)->getDetector() == 6 && p->traj(DC,6)->getLayer() == 6){
          part_DC_c1x = p->traj(DC,6)->getX();
          part_DC_c1y = p->traj(DC,6)->getY();
          part_DC_c1z = p->traj(DC,6)->getZ();
        }
        if(p->traj(DC,18)->getDetector() == 6 && p->traj(DC,18)->getLayer() == 18){
          part_DC_c2x = p->traj(DC,18)->getX();
          part_DC_c2y = p->traj(DC,18)->getY();
          part_DC_c2z = p->traj(DC,18)->getZ();
        }
        if(p->traj(DC,36)->getDetector() == 6 && p->traj(DC,36)->getLayer() == 36){
          part_DC_c3x = p->traj(DC,36)->getX();
          part_DC_c3y = p->traj(DC,36)->getY();
          part_DC_c3z = p->traj(DC,36)->getZ();
        }

        part_pid = 0;
        // Positive particles asssigned to pi^+, negative to pi^-

        if(p->par()->getCharge()==0)continue;
        if(p->par()->getPid())
        part_pid = p->par()->getPid();

        else if(p->par()->getCharge()>0) part_pid = 2;
        else if(p->par()->getCharge()<0) part_pid = 3;

        int part_DC_sector = determineSector(part_DC_c2x, part_DC_c2y,part_DC_c2z);


        h1xy->Fill(part_DC_c1x,part_DC_c1y);
        h2xy->Fill(part_DC_c2x,part_DC_c2y);
        h3xy->Fill(part_DC_c3x,part_DC_c3y);

        // Use this cut if looking at inbending data
        if(!DC_fiducial_cut_theta_phi(part_DC_c1x,part_DC_c1y,part_DC_c1z,part_DC_c2x, part_DC_c2y, part_DC_c2z,part_DC_c3x,part_DC_c3y, part_DC_c3z,part_pid,part_DC_sector,1) ||
        !DC_fiducial_cut_theta_phi(part_DC_c1x,part_DC_c1y,part_DC_c1z,part_DC_c2x, part_DC_c2y, part_DC_c2z,part_DC_c3x,part_DC_c3y, part_DC_c3z,part_pid,part_DC_sector,2) ||
        !DC_fiducial_cut_theta_phi(part_DC_c1x,part_DC_c1y,part_DC_c1z,part_DC_c2x, part_DC_c2y, part_DC_c2z,part_DC_c3x, part_DC_c3y, part_DC_c3z,part_pid,part_DC_sector,3))continue;



        // // Use this cut if looking at outbending data
        // if(!DC_fiducial_cut_XY(part_DC_c1x,part_DC_c1y,part_DC_c1z,part_DC_c2x, part_DC_c2y, part_DC_c2z,part_DC_c3x, part_DC_c3y, part_DC_c3z, part_pid,part_DC_sector,1) ||
        //!DC_fiducial_cut_XY(part_DC_c1x,part_DC_c1y,part_DC_c1z,part_DC_c2x, part_DC_c2y, part_DC_c2z,part_DC_c3x, part_DC_c3y, part_DC_c3z, part_pid,part_DC_sector,2) ||
        //!DC_fiducial_cut_XY(part_DC_c1x,part_DC_c1y,part_DC_c1z,part_DC_c2x, part_DC_c2y, part_DC_c2z,part_DC_c3x, part_DC_c3y, part_DC_c3z, part_pid,part_DC_sector,3)) continue;

        h1xy_after->Fill(part_DC_c1x,part_DC_c1y);
        h2xy_after->Fill(part_DC_c2x,part_DC_c2y);
        h3xy_after->Fill(part_DC_c3x,part_DC_c3y);
      }
    }
  }
  //saving the file
  // fileOutput1.Write();
}
