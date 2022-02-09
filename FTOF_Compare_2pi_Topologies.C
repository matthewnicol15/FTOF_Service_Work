{
  // f1 for missing pi^-, first polarity
  TFile *f1=new TFile("");
  // f2 for missing pi^+ and p, opposite polarity
  TFile *f2=new TFile("");

  // To look at different detector layers change Det_x to:
  // FTOF1A = Det_0, FTOF1B = Det_1, FTOF2 = Det_2

  // To look at different sectors change Sec_x to:
  // For sector n, use Sec_n-1

  auto *pim_Num3 = (TH3F*)f1->Get("h_Numerator_Topo_0_Det_0_Charge_0_Sec_5");
  auto *pim_Num = (TH1F*)pim_Num3->Project3D("y");

  auto *pim_Den3 = (TH3F*)f1->Get("h_Denominator_Topo_0_Det_0_Charge_0_Sec_5");
  auto *pim_Den = (TH1F*)pim_Den3->Project3D("y");

  pim_Num->SetName("pim_Num");
  cout<<"test1"<<endl;
  pim_Num->SetTitle("pim_Num_title");
  cout<<"test1b"<<endl;
  pim_Den->SetName("pim_Den");
  cout<<"test1c"<<endl;
  pim_Den->SetTitle("pim_Den_title");


  auto *pip_Num3 = (TH3F*)f2->Get("h_Numerator_Topo_1_Det_0_Charge_1_Sec_5");
  auto *pip_Num = (TH1F*)pip_Num3->Project3D("y");

  auto *pip_Den3 = (TH3F*)f2->Get("h_Denominator_Topo_1_Det_0_Charge_1_Sec_5");
  auto *pip_Den = (TH1F*)pip_Den3->Project3D("y");
  cout<<"test2"<<endl;


  pip_Num->SetName("pip_Num");
  pip_Num->SetTitle("pip_Num_title");
  pip_Den->SetName("pip_Den");
  pip_Den->SetTitle("pip_Den_title");
  cout<<"test3"<<endl;

  auto *proton_Num3 = (TH3F*)f2->Get("h_Numerator_Topo_2_Det_0_Charge_1_Sec_5");
  auto *proton_Num = (TH1F*)proton_Num3->Project3D("y");

  auto *proton_Den3 = (TH3F*)f2->Get("h_Denominator_Topo_2_Det_0_Charge_1_Sec_5");
  auto *proton_Den = (TH1F*)proton_Den3->Project3D("y");
  cout<<"test4"<<endl;

  pim_Num->Sumw2();
  pip_Den->Sumw2();
  proton_Num->Sumw2();

  pim_Num->Divide(pim_Den);
  pip_Num->Divide(pip_Den);
  proton_Num->Divide(proton_Den);
  cout<<"test5"<<endl;

  for(Int_t l = 2; l < pim_Num->GetNbinsX(); l++){

    // Changing the bin labels to match the counter numbers
    ostringstream test;
    test<< to_string(l-1);

    pim_Num->GetXaxis()->SetBinLabel(l,test.str().c_str());
  }

  cout<<"test6"<<endl;

  pim_Num->SetMinimum(0.0);
  pim_Num->SetMaximum(1.1);
  pip_Num->SetLineColor(kRed);
  proton_Num->SetLineColor(kGreen);

  auto *c1=new TCanvas();
  c1->cd();
  pim_Num->Draw();
  pip_Num->Draw("same");
  proton_Num->Draw("same");
}
