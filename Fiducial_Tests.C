{


  TFile *f1=new TFile("/mnt/f/PhD/FTOF/2pi/miss_pim/FTOF_Efficiency_RGA_FALL2018_skim4_Inbending_2pi_misspim_withoutcuts_ftoffiducial_tight_21092021_01.root");
  // TFile *f1=new TFile("/media/mn688/Elements1/PhD/FTOF/2pi/miss_pim/FTOF_Efficiency_RGA_FALL2018_skim4_Inbending_2pi_misspim_withoutcuts_ftoffiducial_loose_21092021_01.root");

  TH3F *h_Traj[3][2];
  TH3F *h_Tracks[3][2];


  h_Traj[0][0] = (TH3F*) f1->Get("h_Traj_Det_0_Charge_0_Sec_2");
  h_Tracks[0][0] = (TH3F*) f1->Get("h_Tracks_Det_0_Charge_0_Sec_2");


  TH2D *Traj_yz = (TH2D*) h_Traj[0][0]->Project3D("yz"); // denominator
  TH2D *Tracks_yz = (TH2D*) h_Tracks[0][0]->Project3D("yz"); // numerator


  Tracks_yz->Divide(Traj_yz);
  Tracks_yz->SetMaximum(1.0);
  Tracks_yz->SetMinimum(0.0);
  Tracks_yz->Draw("colz");

///////////////////////////////////////////////////////////////////////////////////////////////
  // FTOF fiducial cuts

  // FTOF1A

  // Loose fiducial cut lines
  TLine *l1=new TLine(16.15,77,179.68,420.6);
  TLine *l2=new TLine(-16.15,77,-179.68,420.6);
  TLine *l3=new TLine(-16.15,77,16.15,77);
  TLine *l4=new TLine(-179.68,420.6,179.68,420.6);

  // Tight fiducial cut lines
  TLine *l5=new TLine(-10,80,-160,410);
  TLine *l6=new TLine(10,80,160,410);
  TLine *l7=new TLine(-10,80,10,80);
  TLine *l8=new TLine(-160,410,160,410);

  // FTOF1B

  // Loose fiducial cut lines
  TLine *l9=new TLine(8.635,50.61,203.95,434.27);
  TLine *l10=new TLine(-8.635,50.61,-203.95,434.27);
  TLine *l11=new TLine(-8.635,50.61,8.635,50.61);
  TLine *l12=new TLine(-203.95,434.27,203.95,434.27);

  // Tight fiducial cut lines
  TLine *l13=new TLine(8.635,75,150,390);
  TLine *l14=new TLine(-8.635,75,-150,390);
  TLine *l15=new TLine(-8.635,75,8.635,75);
  TLine *l16=new TLine(-150,390,150,390);

  // FTOF2

  // Loose fiducial cut lines
  // TLine *l17=new TLine(16.15,77,179.68,420.6);
  // TLine *l18=new TLine(-16.15,77,-179.68,420.6);
  // TLine *l19=new TLine(-16.15,77,16.15,77);
  // TLine *l20=new TLine(-179.68,420.6,179.68,420.6);
  //
  // // Tight fiducial cut lines
  // TLine *l21=new TLine(-10,80,-160,410);
  // TLine *l22=new TLine(10,80,160,410);
  // TLine *l23=new TLine(-10,80,10,80);
  // TLine *l24=new TLine(-160,410,160,410);


  // l1->Draw("same");
  // l2->Draw("same");
  // l3->Draw("same");
  // l4->Draw("same");
  //
  // l5->Draw("same");
  // l6->Draw("same");
  // l7->Draw("same");
  // l8->Draw("same");

  l9->Draw("same");
  l10->Draw("same");
  l11->Draw("same");
  l12->Draw("same");

  l13->Draw("same");
  l14->Draw("same");
  l15->Draw("same");
  l16->Draw("same");

}
