{
  // f1 for missing pi^-, first polarity
  TFile *f1=new TFile("/media/mn688/Elements1/PhD/FTOF/RGA_Spring19_Inbending_22022022_01/RGA_Spring19_Inbending_Total.root");
  // f2 for missing pi^+ and p, opposite polarity
  TFile *f2=new TFile("/media/mn688/Elements1/PhD/FTOF/RGA_Fall18_Inbending/RGA_Fall18_Inbending_Total.root");
  TFile *f3=new TFile("/media/mn688/Elements1/PhD/FTOF/RGA_Fall18_Outbending/RGA_Fall18_Outbending_Total.root");

  // To look at different detector layers change Det_x to:
  // x = 0 for FTOF1A, x = 1 for FTOF1B, and x = 2 for FTOF2

  // To look at different sectors change Sec_x to:
  // For sector n, use Sec_n-1
  cout<<"1"<<endl;
  auto *pim_Num3 = (TH2F*)f1->Get("h_Numerator_Topo_3_Det_0_Charge_0_Sec_3");
  auto *pim_Num = (TH1F*)pim_Num3->ProjectionX();
  cout<<"2"<<endl;

  auto *pim_Den3 = (TH2F*)f1->Get("h_Denominator_Topo_3_Det_0_Charge_0_Sec_3");
  auto *pim_Den = (TH1F*)pim_Den3->ProjectionX();
  cout<<"3"<<endl;

  pim_Num->SetName("pim_Num");
  pim_Num->SetTitle("pim_Num_title");
  pim_Den->SetName("pim_Den");
  pim_Den->SetTitle("pim_Den_title");
  cout<<"4"<<endl;


  auto *pip_Num3 = (TH2F*)f2->Get("h_Numerator_Topo_3_Det_0_Charge_0_Sec_3");
  auto *pip_Num = (TH1F*)pip_Num3->ProjectionX();

  auto *pip_Den3 = (TH2F*)f2->Get("h_Denominator_Topo_3_Det_0_Charge_0_Sec_3");
  auto *pip_Den = (TH1F*)pip_Den3->ProjectionX();
  cout<<"5"<<endl;


  pip_Num->SetName("pip_Num");
  pip_Num->SetTitle("pip_Num_title");
  pip_Den->SetName("pip_Den");
  pip_Den->SetTitle("pip_Den_title");
  cout<<"6"<<endl;

  auto *proton_Num3 = (TH2F*)f3->Get("h_Numerator_Topo_3_Det_0_Charge_1_Sec_3");
  auto *proton_Num = (TH1F*)proton_Num3->ProjectionX();
  //
  auto *proton_Den3 = (TH2F*)f3->Get("h_Denominator_Topo_3_Det_0_Charge_1_Sec_3");
  auto *proton_Den = (TH1F*)proton_Den3->ProjectionX();

  pim_Num->Sumw2();
  pip_Den->Sumw2();
  proton_Num->Sumw2();

  pim_Num->Divide(pim_Den);
  pip_Num->Divide(pip_Den);
  proton_Num->Divide(proton_Den);

  for(Int_t l = 2; l < pim_Num->GetNbinsX(); l++){

    // Changing the bin labels to match the counter numbers
    ostringstream test;
    test<< to_string(l-1);

    pim_Num->GetXaxis()->SetBinLabel(l,test.str().c_str());
  }


  pim_Num->SetMinimum(0.75);
  pim_Num->SetMaximum(1.1);
  pip_Num->SetLineColor(kRed);
  proton_Num->SetLineColor(kGreen);

  auto *c1=new TCanvas();
  c1->cd();
  pim_Num->Draw();
  pip_Num->Draw("same");
  proton_Num->Draw("same");
}
