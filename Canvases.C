{

  // TFile *f1=new TFile("/mnt/f/PhD/Service_Work/Sectors/L_Dependence/FTOF_Efficiency_RGB_ECAL_1_ON_Beta_FULL_RUN_CALORIMETER_HITS_ON.root");
  TFile *f1=new TFile("/mnt/f/PhD/FTOF/2pi/miss_pim/FTOF_Efficiency_RGA_FALL2018_skim4_Good_Runs_Inbending_2pi_misspim_0__1_5Momentum_FTOF2_29062021_01.root");
  TFile *f2=new TFile("/mnt/f/PhD/FTOF/2pi/miss_pim/FTOF_Efficiency_RGA_FALL2018_skim4_Good_Runs_Inbending_2pi_misspim_1_5__2_5Momentum_FTOF2_29062021_01.root");
  TFile *f3=new TFile("/mnt/f/PhD/FTOF/2pi/miss_pim/FTOF_Efficiency_RGA_FALL2018_skim4_Good_Runs_Inbending_2pi_misspim_2_5__Momentum_FTOF2_29062021_01.root");

  // Creating arrays of histograms for the numerator and denominator
  // Low Momentum
  // Denominator
  TH3F *h_Traj_low_0[3][2];
  TH3F *h_Traj_low_1[3][2];
  TH3F *h_Traj_low_2[3][2];
  TH3F *h_Traj_low_3[3][2];
  TH3F *h_Traj_low_4[3][2];
  TH3F *h_Traj_low_5[3][2];
  // Clones of denominator to avoid using the same pointers later on
  TH3F *h_Traj_low_0_clone[3][2];
  TH3F *h_Traj_low_1_clone[3][2];
  TH3F *h_Traj_low_2_clone[3][2];
  TH3F *h_Traj_low_3_clone[3][2];
  TH3F *h_Traj_low_4_clone[3][2];
  TH3F *h_Traj_low_5_clone[3][2];
  // Numerator
  TH3F *h_Tracks_low_0[3][2];
  TH3F *h_Tracks_low_1[3][2];
  TH3F *h_Tracks_low_2[3][2];
  TH3F *h_Tracks_low_3[3][2];
  TH3F *h_Tracks_low_4[3][2];
  TH3F *h_Tracks_low_5[3][2];
  // Clones of numerator to avoid using the same pointers later on
  TH3F *h_Tracks_low_0_clone[3][2];
  TH3F *h_Tracks_low_1_clone[3][2];
  TH3F *h_Tracks_low_2_clone[3][2];
  TH3F *h_Tracks_low_3_clone[3][2];
  TH3F *h_Tracks_low_4_clone[3][2];
  TH3F *h_Tracks_low_5_clone[3][2];

  // Med Momentum
  // Denominator
  TH3F *h_Traj_med_0[3][2];
  TH3F *h_Traj_med_1[3][2];
  TH3F *h_Traj_med_2[3][2];
  TH3F *h_Traj_med_3[3][2];
  TH3F *h_Traj_med_4[3][2];
  TH3F *h_Traj_med_5[3][2];
  // Clones of denominator to avoid using the same pointers later on
  TH3F *h_Traj_med_0_clone[3][2];
  TH3F *h_Traj_med_1_clone[3][2];
  TH3F *h_Traj_med_2_clone[3][2];
  TH3F *h_Traj_med_3_clone[3][2];
  TH3F *h_Traj_med_4_clone[3][2];
  TH3F *h_Traj_med_5_clone[3][2];
  // Numerator
  TH3F *h_Tracks_med_0[3][2];
  TH3F *h_Tracks_med_1[3][2];
  TH3F *h_Tracks_med_2[3][2];
  TH3F *h_Tracks_med_3[3][2];
  TH3F *h_Tracks_med_4[3][2];
  TH3F *h_Tracks_med_5[3][2];
  // Clones of numerator to avoid using the same pointers later on
  TH3F *h_Tracks_med_0_clone[3][2];
  TH3F *h_Tracks_med_1_clone[3][2];
  TH3F *h_Tracks_med_2_clone[3][2];
  TH3F *h_Tracks_med_3_clone[3][2];
  TH3F *h_Tracks_med_4_clone[3][2];
  TH3F *h_Tracks_med_5_clone[3][2];

  // High Momentum
  // Denominator
  TH3F *h_Traj_high_0[3][2];
  TH3F *h_Traj_high_1[3][2];
  TH3F *h_Traj_high_2[3][2];
  TH3F *h_Traj_high_3[3][2];
  TH3F *h_Traj_high_4[3][2];
  TH3F *h_Traj_high_5[3][2];
  // Clones of denominator to avoid using the same pointers later on
  TH3F *h_Traj_high_0_clone[3][2];
  TH3F *h_Traj_high_1_clone[3][2];
  TH3F *h_Traj_high_2_clone[3][2];
  TH3F *h_Traj_high_3_clone[3][2];
  TH3F *h_Traj_high_4_clone[3][2];
  TH3F *h_Traj_high_5_clone[3][2];
  // Numerator
  TH3F *h_Tracks_high_0[3][2];
  TH3F *h_Tracks_high_1[3][2];
  TH3F *h_Tracks_high_2[3][2];
  TH3F *h_Tracks_high_3[3][2];
  TH3F *h_Tracks_high_4[3][2];
  TH3F *h_Tracks_high_5[3][2];
  // Clones of numerator to avoid using the same pointers later on
  TH3F *h_Tracks_high_0_clone[3][2];
  TH3F *h_Tracks_high_1_clone[3][2];
  TH3F *h_Tracks_high_2_clone[3][2];
  TH3F *h_Tracks_high_3_clone[3][2];
  TH3F *h_Tracks_high_4_clone[3][2];
  TH3F *h_Tracks_high_5_clone[3][2];

  // Creating arrays of canvases [FTOF layer] [charge]
  TCanvas* canvas_yz[3][2]; // canvases for L vs LPerp
  TCanvas* canvas_yx[3][2]; // canvases for L vs run no
  TCanvas* canvas_zx[3][2]; // canvases for LPerp vs run no

  gStyle->SetOptStat(0); // Remove statistics from canvases

  // Assigning histograms
  // Looping over the FTOF layers
  for(Int_t i_detector=0;i_detector<2;i_detector++){

    // Looping over negative and positive particles
    for(Int_t i_charge=0;i_charge<2;i_charge++){
      // Creating canvas names and file names
      ostringstream Canvas_yz_name_stream;
      ostringstream Canvas_yz_file_name_stream;
      ostringstream Canvas_yx_name_stream;
      ostringstream Canvas_yx_file_name_stream;
      ostringstream Canvas_zx_name_stream;
      ostringstream Canvas_zx_file_name_stream;

      // Setting the Title for each canvas
      Canvas_yz_name_stream<<"canvas_det_"<<i_detector<<"_Charge_"<<2*i_charge-1;

      // Setting the file name for each canvas
      // FTOF1A
      if(i_detector==0){
        Canvas_yz_file_name_stream<<"/mnt/f/PhD/Service_Work/Sectors/L_Dependence/Canvases/2pi/RGA_FALL2018_Inbending_skim4_momentum_2pi_290621/RGA_FALL2018_skim4_pim_Inbending_1D_det_FTOF1A_L_LPerp_Charge_"<<2*i_charge-1<<"_290621_01.pdf";
        Canvas_yx_file_name_stream<<"/mnt/f/PhD/Service_Work/Sectors/L_Dependence/Canvases/2pi/RGA_FALL2018_Inbending_skim4_momentum_2pi_290621/RGA_FALL2018_skim4_pim_Inbending_1D_det_FTOF1A_L_Run_Charge_"<<2*i_charge-1<<"_290621_01.pdf";
        Canvas_zx_file_name_stream<<"/mnt/f/PhD/Service_Work/Sectors/L_Dependence/Canvases/2pi/RGA_FALL2018_Inbending_skim4_momentum_2pi_290621/RGA_FALL2018_skim4_pim_Inbending_1D_det_FTOF1A_LPerp_Run_Charge_"<<2*i_charge-1<<"_290621_01.pdf";
      }

      // FTOF1B
      if(i_detector==1){
        Canvas_yz_file_name_stream<<"/mnt/f/PhD/Service_Work/Sectors/L_Dependence/Canvases/2pi/RGA_FALL2018_Inbending_skim4_momentum_2pi_290621/RGA_FALL2018_skim4_pim_Inbending_1D_det_FTOF1B_L_LPerp_Charge_"<<2*i_charge-1<<"_290621_01.pdf";
        Canvas_yx_file_name_stream<<"/mnt/f/PhD/Service_Work/Sectors/L_Dependence/Canvases/2pi/RGA_FALL2018_Inbending_skim4_momentum_2pi_290621/RGA_FALL2018_skim4_pim_Inbending_1D_det_FTOF1B_L_Run_Charge_"<<2*i_charge-1<<"_290621_01.pdf";
        Canvas_zx_file_name_stream<<"/mnt/f/PhD/Service_Work/Sectors/L_Dependence/Canvases/2pi/RGA_FALL2018_Inbending_skim4_momentum_2pi_290621/RGA_FALL2018_skim4_pim_Inbending_1D_det_FTOF1B_LPerp_Run_Charge_"<<2*i_charge-1<<"_290621_01.pdf";
      }

      // FTOF2
      if(i_detector==2){
        Canvas_yz_file_name_stream<<"/mnt/f/PhD/Service_Work/Sectors/L_Dependence/Canvases/2pi/RGA_FALL2018_Inbending_skim4_momentum_2pi_290621/RGA_FALL2018_skim4_pim_Inbending_1D_det_FTOF2_L_LPerp_Charge_"<<2*i_charge-1<<"_290621_01.pdf";
        Canvas_yx_file_name_stream<<"/mnt/f/PhD/Service_Work/Sectors/L_Dependence/Canvases/2pi/RGA_FALL2018_Inbending_skim4_momentum_2pi_290621/RGA_FALL2018_skim4_pim_Inbending_1D_det_FTOF2_L_Run_Charge_"<<2*i_charge-1<<"_290621_01.pdf";
        Canvas_zx_file_name_stream<<"/mnt/f/PhD/Service_Work/Sectors/L_Dependence/Canvases/2pi/RGA_FALL2018_Inbending_skim4_momentum_2pi_290621/RGA_FALL2018_skim4_pim_Inbending_1D_det_FTOF2_LPerp_Run_Charge_"<<2*i_charge-1<<"_290621_01.pdf";
      }

      // creating canvas for each FTOF layer for L vs LPerp
      canvas_yz[i_detector][i_charge] = new TCanvas(Canvas_yz_name_stream.str().c_str(),"", 2500,1600);
      canvas_yz[i_detector][i_charge]->Divide(3,2);

      // creating canvas for each FTOF layer for L vs Run
      canvas_yx[i_detector][i_charge] = new TCanvas(Canvas_yx_name_stream.str().c_str(),"", 2500,1600);
      canvas_yx[i_detector][i_charge]->Divide(3,2);

      // creating canvas for each FTOF layer for LPerp vs Run
      canvas_zx[i_detector][i_charge] = new TCanvas(Canvas_zx_name_stream.str().c_str(),"", 2500,1600);
      canvas_zx[i_detector][i_charge]->Divide(3,2);

      // Looping over the sectors
      // for(Int_t i_sector=0;i_sector<2;i_sector++){
      // Creating histogram names
      ostringstream Traj_name_stream_0,Traj_name_stream_1, Traj_name_stream_2, Traj_name_stream_3,
      Traj_name_stream_4, Traj_name_stream_5; // denominator
      ostringstream Tracks_name_stream_0, Tracks_name_stream_1, Tracks_name_stream_2,
      Tracks_name_stream_3, Tracks_name_stream_4, Tracks_name_stream_5; // numerator
      ostringstream Efficiency_name_stream; // efficiency

      // Setting the histogram names
      // Denominator
      Traj_name_stream_0<<"h_Traj_Det_"<<i_detector<<"_Charge_"<<i_charge<<"_Sec_0";  // denominator
      Traj_name_stream_1<<"h_Traj_Det_"<<i_detector<<"_Charge_"<<i_charge<<"_Sec_1";  // denominator
      Traj_name_stream_2<<"h_Traj_Det_"<<i_detector<<"_Charge_"<<i_charge<<"_Sec_2";  // denominator
      Traj_name_stream_3<<"h_Traj_Det_"<<i_detector<<"_Charge_"<<i_charge<<"_Sec_3";  // denominator
      Traj_name_stream_4<<"h_Traj_Det_"<<i_detector<<"_Charge_"<<i_charge<<"_Sec_4";  // denominator
      Traj_name_stream_5<<"h_Traj_Det_"<<i_detector<<"_Charge_"<<i_charge<<"_Sec_5";  // denominator
      // Numerator
      Tracks_name_stream_0<<"h_Tracks_Det_"<<i_detector<<"_Charge_"<<i_charge<<"_Sec_0"; // numerator
      Tracks_name_stream_1<<"h_Tracks_Det_"<<i_detector<<"_Charge_"<<i_charge<<"_Sec_1"; // numerator
      Tracks_name_stream_2<<"h_Tracks_Det_"<<i_detector<<"_Charge_"<<i_charge<<"_Sec_2"; // numerator
      Tracks_name_stream_3<<"h_Tracks_Det_"<<i_detector<<"_Charge_"<<i_charge<<"_Sec_3"; // numerator
      Tracks_name_stream_4<<"h_Tracks_Det_"<<i_detector<<"_Charge_"<<i_charge<<"_Sec_4"; // numerator
      Tracks_name_stream_5<<"h_Tracks_Det_"<<i_detector<<"_Charge_"<<i_charge<<"_Sec_5"; // numerator

      // Histogram title for efficiency plots
      ostringstream Efficiency_yz_title_stream_0, Efficiency_yz_title_stream_1,
      Efficiency_yz_title_stream_2, Efficiency_yz_title_stream_3, Efficiency_yz_title_stream_4,
      Efficiency_yz_title_stream_5;
      // ostringstream Efficiency_yx_title_stream,Efficiency_zx_title_stream;

      //FTOF1A
      if (i_detector==0){
        Efficiency_yz_title_stream_0<<"Efficiency FTOF1A vs L Charge "<<2*i_charge-1<<" Sec 1; L [cm]; #epsilon";
        Efficiency_yz_title_stream_1<<"Efficiency FTOF1A vs L Charge "<<2*i_charge-1<<" Sec 2; L [cm]; #epsilon";
        Efficiency_yz_title_stream_2<<"Efficiency FTOF1A vs L Charge "<<2*i_charge-1<<" Sec 3; L [cm]; #epsilon";
        Efficiency_yz_title_stream_3<<"Efficiency FTOF1A vs L Charge "<<2*i_charge-1<<" Sec 4; L [cm]; #epsilon";
        Efficiency_yz_title_stream_4<<"Efficiency FTOF1A vs L Charge "<<2*i_charge-1<<" Sec 5; L [cm]; #epsilon";
        Efficiency_yz_title_stream_5<<"Efficiency FTOF1A vs L Charge "<<2*i_charge-1<<" Sec 6; L [cm]; #epsilon";
        // Efficiency_yz_title_stream<<"Efficiency FTOF1A L vs L_{Perp} Charge "<<2*i_charge-1<<" Sec "<<i_sector+1<<"; L_{Perp} [cm]; L [cm]";
        // Efficiency_yx_title_stream<<"Efficiency FTOF1A L vs Run Charge "<<2*i_charge-1<<" Sec "<<i_sector+1<<"; Run no.; L [cm]";
        // Efficiency_zx_title_stream<<"Efficiency FTOF1A L_{Perp} vs Run Charge "<<2*i_charge-1<<" Sec "<<i_sector+1<<"; Run no.; L_{Perp} [cm]";
      }

      // FTOF1B
      else if (i_detector==1) {
        Efficiency_yz_title_stream_0<<"Efficiency FTOF1B vs L Charge "<<2*i_charge-1<<" Sec 1; L [cm]; #epsilon";
        Efficiency_yz_title_stream_1<<"Efficiency FTOF1B vs L Charge "<<2*i_charge-1<<" Sec 2; L [cm]; #epsilon";
        Efficiency_yz_title_stream_2<<"Efficiency FTOF1B vs L Charge "<<2*i_charge-1<<" Sec 3; L [cm]; #epsilon";
        Efficiency_yz_title_stream_3<<"Efficiency FTOF1B vs L Charge "<<2*i_charge-1<<" Sec 4; L [cm]; #epsilon";
        Efficiency_yz_title_stream_4<<"Efficiency FTOF1B vs L Charge "<<2*i_charge-1<<" Sec 5; L [cm]; #epsilon";
        Efficiency_yz_title_stream_5<<"Efficiency FTOF1B vs L Charge "<<2*i_charge-1<<" Sec 6; L [cm]; #epsilon";
        // Efficiency_yz_title_stream<<"Efficiency FTOF1B L vs L_{Perp} Charge "<<2*i_charge-1<<" Sec "<<i_sector+1<<"; L_{Perp} [cm]; L [cm]";
        // Efficiency_yx_title_stream<<"Efficiency FTOF1B L vs Run Charge "<<2*i_charge-1<<" Sec "<<i_sector+1<<"; Run no.; L [cm]";
        // Efficiency_zx_title_stream<<"Efficiency FTOF1B L_{Perp} vs Run Charge "<<2*i_charge-1<<" Sec "<<i_sector+1<<"; Run no.; L_{Perp} [cm]";
      }

      // FTOF2
      else if (i_detector==2) {
        Efficiency_yz_title_stream_0<<"Efficiency FTOF2 vs L Charge "<<2*i_charge-1<<" Sec 1; L [cm]; #epsilon";
        Efficiency_yz_title_stream_1<<"Efficiency FTOF2 vs L Charge "<<2*i_charge-1<<" Sec 2; L [cm]; #epsilon";
        Efficiency_yz_title_stream_2<<"Efficiency FTOF2 vs L Charge "<<2*i_charge-1<<" Sec 3; L [cm]; #epsilon";
        Efficiency_yz_title_stream_3<<"Efficiency FTOF2 vs L Charge "<<2*i_charge-1<<" Sec 4; L [cm]; #epsilon";
        Efficiency_yz_title_stream_4<<"Efficiency FTOF2 vs L Charge "<<2*i_charge-1<<" Sec 5; L [cm]; #epsilon";
        Efficiency_yz_title_stream_5<<"Efficiency FTOF2 vs L Charge "<<2*i_charge-1<<" Sec 6; L [cm]; #epsilon";
        // Efficiency_yz_title_stream<<"Efficiency FTOF2 L vs L_{Perp} Charge "<<2*i_charge-1<<" Sec "<<i_sector+1<<";  L_{Perp} [cm]; L [cm]";
        Efficiency_yx_title_stream<<"Efficiency FTOF2 L vs Run Charge "<<2*i_charge-1<<" Sec "<<i_sector+1<<"; Run no.; L [cm]";
        Efficiency_zx_title_stream<<"Efficiency FTOF2 L_{Perp} vs Run Charge "<<2*i_charge-1<<" Sec "<<i_sector+1<<"; Run no.; L_{Perp} [cm]";
      }

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Getting histograms from input files

      // Low momentum
      // Getting the denominator histograms from 1st file
      h_Traj_low_0[i_detector][i_charge] = (TH3F*) f1->Get(Traj_name_stream_0.str().c_str());
      h_Traj_low_1[i_detector][i_charge] = (TH3F*) f1->Get(Traj_name_stream_1.str().c_str());
      h_Traj_low_2[i_detector][i_charge] = (TH3F*) f1->Get(Traj_name_stream_2.str().c_str());
      h_Traj_low_3[i_detector][i_charge] = (TH3F*) f1->Get(Traj_name_stream_3.str().c_str());
      h_Traj_low_4[i_detector][i_charge] = (TH3F*) f1->Get(Traj_name_stream_4.str().c_str());
      h_Traj_low_5[i_detector][i_charge] = (TH3F*) f1->Get(Traj_name_stream_5.str().c_str());
      // Clone the denominator histograms
      h_Traj_low_0_clone[i_detector][i_charge] = (TH3F*) h_Traj_low_0[i_detector][i_charge]->Clone("h_Traj_low_0");
      h_Traj_low_1_clone[i_detector][i_charge] = (TH3F*) h_Traj_low_1[i_detector][i_charge]->Clone("h_Traj_low_1");
      h_Traj_low_2_clone[i_detector][i_charge] = (TH3F*) h_Traj_low_2[i_detector][i_charge]->Clone("h_Traj_low_2");
      h_Traj_low_3_clone[i_detector][i_charge] = (TH3F*) h_Traj_low_3[i_detector][i_charge]->Clone("h_Traj_low_3");
      h_Traj_low_4_clone[i_detector][i_charge] = (TH3F*) h_Traj_low_4[i_detector][i_charge]->Clone("h_Traj_low_4");
      h_Traj_low_5_clone[i_detector][i_charge] = (TH3F*) h_Traj_low_5[i_detector][i_charge]->Clone("h_Traj_low_5");

      // Getting the numerator histograms from 1st file
      h_Tracks_low_0[i_detector][i_charge] = (TH3F*) f1->Get(Tracks_name_stream_0.str().c_str());
      h_Tracks_low_1[i_detector][i_charge] = (TH3F*) f1->Get(Tracks_name_stream_1.str().c_str());
      h_Tracks_low_2[i_detector][i_charge] = (TH3F*) f1->Get(Tracks_name_stream_2.str().c_str());
      h_Tracks_low_3[i_detector][i_charge] = (TH3F*) f1->Get(Tracks_name_stream_3.str().c_str());
      h_Tracks_low_4[i_detector][i_charge] = (TH3F*) f1->Get(Tracks_name_stream_4.str().c_str());
      h_Tracks_low_5[i_detector][i_charge] = (TH3F*) f1->Get(Tracks_name_stream_5.str().c_str());
      // Clone the numerator histograms
      h_Tracks_low_0_clone[i_detector][i_charge] = (TH3F*) h_Tracks_low_0[i_detector][i_charge]->Clone("h_Tracks_low_0");
      h_Tracks_low_1_clone[i_detector][i_charge] = (TH3F*) h_Tracks_low_1[i_detector][i_charge]->Clone("h_Tracks_low_1");
      h_Tracks_low_2_clone[i_detector][i_charge] = (TH3F*) h_Tracks_low_2[i_detector][i_charge]->Clone("h_Tracks_low_2");
      h_Tracks_low_3_clone[i_detector][i_charge] = (TH3F*) h_Tracks_low_3[i_detector][i_charge]->Clone("h_Tracks_low_3");
      h_Tracks_low_4_clone[i_detector][i_charge] = (TH3F*) h_Tracks_low_4[i_detector][i_charge]->Clone("h_Tracks_low_4");
      h_Tracks_low_5_clone[i_detector][i_charge] = (TH3F*) h_Tracks_low_5[i_detector][i_charge]->Clone("h_Tracks_low_5");

      // med momentum
      // Getting the denominator histograms from 1st file
      h_Traj_med_0[i_detector][i_charge] = (TH3F*) f2->Get(Traj_name_stream_0.str().c_str());
      h_Traj_med_1[i_detector][i_charge] = (TH3F*) f2->Get(Traj_name_stream_1.str().c_str());
      h_Traj_med_2[i_detector][i_charge] = (TH3F*) f2->Get(Traj_name_stream_2.str().c_str());
      h_Traj_med_3[i_detector][i_charge] = (TH3F*) f2->Get(Traj_name_stream_3.str().c_str());
      h_Traj_med_4[i_detector][i_charge] = (TH3F*) f2->Get(Traj_name_stream_4.str().c_str());
      h_Traj_med_5[i_detector][i_charge] = (TH3F*) f2->Get(Traj_name_stream_5.str().c_str());
      // Clone the denominator histograms
      h_Traj_med_0_clone[i_detector][i_charge] = (TH3F*) h_Traj_med_0[i_detector][i_charge]->Clone("h_Traj_med_0");
      h_Traj_med_1_clone[i_detector][i_charge] = (TH3F*) h_Traj_med_1[i_detector][i_charge]->Clone("h_Traj_med_1");
      h_Traj_med_2_clone[i_detector][i_charge] = (TH3F*) h_Traj_med_2[i_detector][i_charge]->Clone("h_Traj_med_2");
      h_Traj_med_3_clone[i_detector][i_charge] = (TH3F*) h_Traj_med_3[i_detector][i_charge]->Clone("h_Traj_med_3");
      h_Traj_med_4_clone[i_detector][i_charge] = (TH3F*) h_Traj_med_4[i_detector][i_charge]->Clone("h_Traj_med_4");
      h_Traj_med_5_clone[i_detector][i_charge] = (TH3F*) h_Traj_med_5[i_detector][i_charge]->Clone("h_Traj_med_5");

      // Getting the numerator histograms from 1st file
      h_Tracks_med_0[i_detector][i_charge] = (TH3F*) f2->Get(Tracks_name_stream_0.str().c_str());
      h_Tracks_med_1[i_detector][i_charge] = (TH3F*) f2->Get(Tracks_name_stream_1.str().c_str());
      h_Tracks_med_2[i_detector][i_charge] = (TH3F*) f2->Get(Tracks_name_stream_2.str().c_str());
      h_Tracks_med_3[i_detector][i_charge] = (TH3F*) f2->Get(Tracks_name_stream_3.str().c_str());
      h_Tracks_med_4[i_detector][i_charge] = (TH3F*) f2->Get(Tracks_name_stream_4.str().c_str());
      h_Tracks_med_5[i_detector][i_charge] = (TH3F*) f2->Get(Tracks_name_stream_5.str().c_str());
      // Clone the numerator histograms
      h_Tracks_med_0_clone[i_detector][i_charge] = (TH3F*) h_Tracks_med_0[i_detector][i_charge]->Clone("h_Tracks_med_0");
      h_Tracks_med_1_clone[i_detector][i_charge] = (TH3F*) h_Tracks_med_1[i_detector][i_charge]->Clone("h_Tracks_med_1");
      h_Tracks_med_2_clone[i_detector][i_charge] = (TH3F*) h_Tracks_med_2[i_detector][i_charge]->Clone("h_Tracks_med_2");
      h_Tracks_med_3_clone[i_detector][i_charge] = (TH3F*) h_Tracks_med_3[i_detector][i_charge]->Clone("h_Tracks_med_3");
      h_Tracks_med_4_clone[i_detector][i_charge] = (TH3F*) h_Tracks_med_4[i_detector][i_charge]->Clone("h_Tracks_med_4");
      h_Tracks_med_5_clone[i_detector][i_charge] = (TH3F*) h_Tracks_med_5[i_detector][i_charge]->Clone("h_Tracks_med_5");

      // High momentum
      // Getting the denominator histograms from 2nd file
      h_Traj_high_0[i_detector][i_charge] = (TH3F*) f3->Get(Traj_name_stream_0.str().c_str());
      h_Traj_high_1[i_detector][i_charge] = (TH3F*) f3->Get(Traj_name_stream_1.str().c_str());
      h_Traj_high_2[i_detector][i_charge] = (TH3F*) f3->Get(Traj_name_stream_2.str().c_str());
      h_Traj_high_3[i_detector][i_charge] = (TH3F*) f3->Get(Traj_name_stream_3.str().c_str());
      h_Traj_high_4[i_detector][i_charge] = (TH3F*) f3->Get(Traj_name_stream_4.str().c_str());
      h_Traj_high_5[i_detector][i_charge] = (TH3F*) f3->Get(Traj_name_stream_5.str().c_str());
      // Clone the denominator histograms
      h_Traj_high_0_clone[i_detector][i_charge] = (TH3F*) h_Traj_high_0[i_detector][i_charge]->Clone("h_Traj_high_0");
      h_Traj_high_1_clone[i_detector][i_charge] = (TH3F*) h_Traj_high_1[i_detector][i_charge]->Clone("h_Traj_high_1");
      h_Traj_high_2_clone[i_detector][i_charge] = (TH3F*) h_Traj_high_2[i_detector][i_charge]->Clone("h_Traj_high_2");
      h_Traj_high_3_clone[i_detector][i_charge] = (TH3F*) h_Traj_high_3[i_detector][i_charge]->Clone("h_Traj_high_3");
      h_Traj_high_4_clone[i_detector][i_charge] = (TH3F*) h_Traj_high_4[i_detector][i_charge]->Clone("h_Traj_high_4");
      h_Traj_high_5_clone[i_detector][i_charge] = (TH3F*) h_Traj_high_5[i_detector][i_charge]->Clone("h_Traj_high_5");

      // Getting the numerator histograms from 3rd file
      h_Tracks_high_0[i_detector][i_charge] = (TH3F*) f3->Get(Tracks_name_stream_0.str().c_str());
      h_Tracks_high_1[i_detector][i_charge] = (TH3F*) f3->Get(Tracks_name_stream_1.str().c_str());
      h_Tracks_high_2[i_detector][i_charge] = (TH3F*) f3->Get(Tracks_name_stream_2.str().c_str());
      h_Tracks_high_3[i_detector][i_charge] = (TH3F*) f3->Get(Tracks_name_stream_3.str().c_str());
      h_Tracks_high_4[i_detector][i_charge] = (TH3F*) f3->Get(Tracks_name_stream_4.str().c_str());
      h_Tracks_high_5[i_detector][i_charge] = (TH3F*) f3->Get(Tracks_name_stream_5.str().c_str());
      // Clone the numerator histograms
      h_Tracks_high_0_clone[i_detector][i_charge] = (TH3F*) h_Tracks_high_0[i_detector][i_charge]->Clone("h_Tracks_high_0");
      h_Tracks_high_1_clone[i_detector][i_charge] = (TH3F*) h_Tracks_high_1[i_detector][i_charge]->Clone("h_Tracks_high_1");
      h_Tracks_high_2_clone[i_detector][i_charge] = (TH3F*) h_Tracks_high_2[i_detector][i_charge]->Clone("h_Tracks_high_2");
      h_Tracks_high_3_clone[i_detector][i_charge] = (TH3F*) h_Tracks_high_3[i_detector][i_charge]->Clone("h_Tracks_high_3");
      h_Tracks_high_4_clone[i_detector][i_charge] = (TH3F*) h_Tracks_high_4[i_detector][i_charge]->Clone("h_Tracks_high_4");
      h_Tracks_high_5_clone[i_detector][i_charge] = (TH3F*) h_Tracks_high_5[i_detector][i_charge]->Clone("h_Tracks_high_5");

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Making 1D projections of the numerator and denominator histograms

      // Low momentum
      TH1D *Traj_y_low_0 = (TH1D*) h_Traj_low_0_clone[i_detector][i_charge]->Project3D("y"); // denominator
      TH1D *Traj_y_low_1 = (TH1D*) h_Traj_low_1_clone[i_detector][i_charge]->Project3D("y"); // denominator
      TH1D *Traj_y_low_2 = (TH1D*) h_Traj_low_2_clone[i_detector][i_charge]->Project3D("y"); // denominator
      TH1D *Traj_y_low_3 = (TH1D*) h_Traj_low_3_clone[i_detector][i_charge]->Project3D("y"); // denominator
      TH1D *Traj_y_low_4 = (TH1D*) h_Traj_low_4_clone[i_detector][i_charge]->Project3D("y"); // denominator
      TH1D *Traj_y_low_5 = (TH1D*) h_Traj_low_5_clone[i_detector][i_charge]->Project3D("y"); // denominator
      TH1D *Tracks_y_low_0 = (TH1D*) h_Tracks_low_0_clone[i_detector][i_charge]->Project3D("y"); // numerator
      TH1D *Tracks_y_low_1 = (TH1D*) h_Tracks_low_1_clone[i_detector][i_charge]->Project3D("y"); // numerator
      TH1D *Tracks_y_low_2 = (TH1D*) h_Tracks_low_2_clone[i_detector][i_charge]->Project3D("y"); // numerator
      TH1D *Tracks_y_low_3 = (TH1D*) h_Tracks_low_3_clone[i_detector][i_charge]->Project3D("y"); // numerator
      TH1D *Tracks_y_low_4 = (TH1D*) h_Tracks_low_4_clone[i_detector][i_charge]->Project3D("y"); // numerator
      TH1D *Tracks_y_low_5 = (TH1D*) h_Tracks_low_5_clone[i_detector][i_charge]->Project3D("y"); // numerator

      // Med momentum
      TH1D *Traj_y_med_0 = (TH1D*) h_Traj_med_0_clone[i_detector][i_charge]->Project3D("y"); // denominator
      TH1D *Traj_y_med_1 = (TH1D*) h_Traj_med_1_clone[i_detector][i_charge]->Project3D("y"); // denominator
      TH1D *Traj_y_med_2 = (TH1D*) h_Traj_med_2_clone[i_detector][i_charge]->Project3D("y"); // denominator
      TH1D *Traj_y_med_3 = (TH1D*) h_Traj_med_3_clone[i_detector][i_charge]->Project3D("y"); // denominator
      TH1D *Traj_y_med_4 = (TH1D*) h_Traj_med_4_clone[i_detector][i_charge]->Project3D("y"); // denominator
      TH1D *Traj_y_med_5 = (TH1D*) h_Traj_med_5_clone[i_detector][i_charge]->Project3D("y"); // denominator
      TH1D *Tracks_y_med_0 = (TH1D*) h_Tracks_med_0_clone[i_detector][i_charge]->Project3D("y"); // numerator
      TH1D *Tracks_y_med_1 = (TH1D*) h_Tracks_med_1_clone[i_detector][i_charge]->Project3D("y"); // numerator
      TH1D *Tracks_y_med_2 = (TH1D*) h_Tracks_med_2_clone[i_detector][i_charge]->Project3D("y"); // numerator
      TH1D *Tracks_y_med_3 = (TH1D*) h_Tracks_med_3_clone[i_detector][i_charge]->Project3D("y"); // numerator
      TH1D *Tracks_y_med_4 = (TH1D*) h_Tracks_med_4_clone[i_detector][i_charge]->Project3D("y"); // numerator
      TH1D *Tracks_y_med_5 = (TH1D*) h_Tracks_med_5_clone[i_detector][i_charge]->Project3D("y"); // numerator

      // High momentum
      TH1D *Traj_y_high_0 = (TH1D*) h_Traj_high_0_clone[i_detector][i_charge]->Project3D("y"); // denominator
      TH1D *Traj_y_high_1 = (TH1D*) h_Traj_high_1_clone[i_detector][i_charge]->Project3D("y"); // denominator
      TH1D *Traj_y_high_2 = (TH1D*) h_Traj_high_2_clone[i_detector][i_charge]->Project3D("y"); // denominator
      TH1D *Traj_y_high_3 = (TH1D*) h_Traj_high_3_clone[i_detector][i_charge]->Project3D("y"); // denominator
      TH1D *Traj_y_high_4 = (TH1D*) h_Traj_high_4_clone[i_detector][i_charge]->Project3D("y"); // denominator
      TH1D *Traj_y_high_5 = (TH1D*) h_Traj_high_5_clone[i_detector][i_charge]->Project3D("y"); // denominator
      TH1D *Tracks_y_high_0 = (TH1D*) h_Tracks_high_0_clone[i_detector][i_charge]->Project3D("y"); // numerator
      TH1D *Tracks_y_high_1 = (TH1D*) h_Tracks_high_1_clone[i_detector][i_charge]->Project3D("y"); // numerator
      TH1D *Tracks_y_high_2 = (TH1D*) h_Tracks_high_2_clone[i_detector][i_charge]->Project3D("y"); // numerator
      TH1D *Tracks_y_high_3 = (TH1D*) h_Tracks_high_3_clone[i_detector][i_charge]->Project3D("y"); // numerator
      TH1D *Tracks_y_high_4 = (TH1D*) h_Tracks_high_4_clone[i_detector][i_charge]->Project3D("y"); // numerator
      TH1D *Tracks_y_high_5 = (TH1D*) h_Tracks_high_5_clone[i_detector][i_charge]->Project3D("y"); // numerator

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Calculating errors

      // Low momentum
      Tracks_y_low_0->Sumw2();
      Tracks_y_low_1->Sumw2();
      Tracks_y_low_2->Sumw2();
      Tracks_y_low_3->Sumw2();
      Tracks_y_low_4->Sumw2();
      Tracks_y_low_5->Sumw2();

      // Med momentum
      Tracks_y_med_0->Sumw2();
      Tracks_y_med_1->Sumw2();
      Tracks_y_med_2->Sumw2();
      Tracks_y_med_3->Sumw2();
      Tracks_y_med_4->Sumw2();
      Tracks_y_med_5->Sumw2();

      // High momentum
      Tracks_y_high_0->Sumw2();
      Tracks_y_high_1->Sumw2();
      Tracks_y_high_2->Sumw2();
      Tracks_y_high_3->Sumw2();
      Tracks_y_high_4->Sumw2();
      Tracks_y_high_5->Sumw2();


      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Dividing the numerator histogram by the denominator histogram

      // Low momentum
      Tracks_y_low_0->Divide(Traj_y_low_0);
      Tracks_y_low_1->Divide(Traj_y_low_1);
      Tracks_y_low_2->Divide(Traj_y_low_2);
      Tracks_y_low_3->Divide(Traj_y_low_3);
      Tracks_y_low_4->Divide(Traj_y_low_4);
      Tracks_y_low_5->Divide(Traj_y_low_5);

      // Med momentum
      Tracks_y_med_0->Divide(Traj_y_med_0);
      Tracks_y_med_1->Divide(Traj_y_med_1);
      Tracks_y_med_2->Divide(Traj_y_med_2);
      Tracks_y_med_3->Divide(Traj_y_med_3);
      Tracks_y_med_4->Divide(Traj_y_med_4);
      Tracks_y_med_5->Divide(Traj_y_med_5);

      // High momentum
      Tracks_y_high_0->Divide(Traj_y_high_0);
      Tracks_y_high_1->Divide(Traj_y_high_1);
      Tracks_y_high_2->Divide(Traj_y_high_2);
      Tracks_y_high_3->Divide(Traj_y_high_3);
      Tracks_y_high_4->Divide(Traj_y_high_4);
      Tracks_y_high_5->Divide(Traj_y_high_5);

      // // Making 2D projections of the trajectory and track histograms L vs run
      // TH2D *Traj_yx = (TH2D*) h_Traj[i_detector][i_charge][i_sector]->Project3D("yx");
      // TH2D *Tracks_yx = (TH2D*) h_Tracks[i_detector][i_charge][i_sector]->Project3D("yx");
      //
      // // Dividing the numerator 2D histogram by the denominator 2D histogram
      // Tracks_yx->Divide(Traj_yx);
      //
      // // Making 2D projections of the trajectory and track histograms Lperp vs run
      // TH2D *Traj_zx = (TH2D*) h_Traj[i_detector][i_charge][i_sector]->Project3D("zx");
      // TH2D *Tracks_zx = (TH2D*) h_Tracks[i_detector][i_charge][i_sector]->Project3D("zx");
      //
      // // Dividing the numerator 2D histogram by the denominator 2D histogram
      // Tracks_zx->Divide(Traj_zx);

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Setting titles and aligning them on canvases

      // Low momentum
      Tracks_y_low_0->GetXaxis()->SetTitleSize(0.045);
      Tracks_y_low_1->GetXaxis()->SetTitleSize(0.045);
      Tracks_y_low_2->GetXaxis()->SetTitleSize(0.045);
      Tracks_y_low_3->GetXaxis()->SetTitleSize(0.045);
      Tracks_y_low_4->GetXaxis()->SetTitleSize(0.045);
      Tracks_y_low_5->GetXaxis()->SetTitleSize(0.045);
      Tracks_y_low_0->GetYaxis()->SetTitleSize(0.045);
      Tracks_y_low_1->GetYaxis()->SetTitleSize(0.045);
      Tracks_y_low_2->GetYaxis()->SetTitleSize(0.045);
      Tracks_y_low_3->GetYaxis()->SetTitleSize(0.045);
      Tracks_y_low_4->GetYaxis()->SetTitleSize(0.045);
      Tracks_y_low_5->GetYaxis()->SetTitleSize(0.045);
      Tracks_y_low_0->GetYaxis()->SetTitleOffset(1.1);
      Tracks_y_low_1->GetYaxis()->SetTitleOffset(1.1);
      Tracks_y_low_2->GetYaxis()->SetTitleOffset(1.1);
      Tracks_y_low_3->GetYaxis()->SetTitleOffset(1.1);
      Tracks_y_low_4->GetYaxis()->SetTitleOffset(1.1);
      Tracks_y_low_5->GetYaxis()->SetTitleOffset(1.1);

      // Med momentum
      Tracks_y_med_0->GetXaxis()->SetTitleSize(0.045);
      Tracks_y_med_1->GetXaxis()->SetTitleSize(0.045);
      Tracks_y_med_2->GetXaxis()->SetTitleSize(0.045);
      Tracks_y_med_3->GetXaxis()->SetTitleSize(0.045);
      Tracks_y_med_4->GetXaxis()->SetTitleSize(0.045);
      Tracks_y_med_5->GetXaxis()->SetTitleSize(0.045);
      Tracks_y_med_0->GetYaxis()->SetTitleSize(0.045);
      Tracks_y_med_1->GetYaxis()->SetTitleSize(0.045);
      Tracks_y_med_2->GetYaxis()->SetTitleSize(0.045);
      Tracks_y_med_3->GetYaxis()->SetTitleSize(0.045);
      Tracks_y_med_4->GetYaxis()->SetTitleSize(0.045);
      Tracks_y_med_5->GetYaxis()->SetTitleSize(0.045);
      Tracks_y_med_0->GetYaxis()->SetTitleOffset(1.1);
      Tracks_y_med_1->GetYaxis()->SetTitleOffset(1.1);
      Tracks_y_med_2->GetYaxis()->SetTitleOffset(1.1);
      Tracks_y_med_3->GetYaxis()->SetTitleOffset(1.1);
      Tracks_y_med_4->GetYaxis()->SetTitleOffset(1.1);
      Tracks_y_med_5->GetYaxis()->SetTitleOffset(1.1);

      // High momentum
      Tracks_y_high_0->GetXaxis()->SetTitleSize(0.045);
      Tracks_y_high_1->GetXaxis()->SetTitleSize(0.045);
      Tracks_y_high_2->GetXaxis()->SetTitleSize(0.045);
      Tracks_y_high_3->GetXaxis()->SetTitleSize(0.045);
      Tracks_y_high_4->GetXaxis()->SetTitleSize(0.045);
      Tracks_y_high_5->GetXaxis()->SetTitleSize(0.045);
      Tracks_y_high_0->GetYaxis()->SetTitleSize(0.045);
      Tracks_y_high_1->GetYaxis()->SetTitleSize(0.045);
      Tracks_y_high_2->GetYaxis()->SetTitleSize(0.045);
      Tracks_y_high_3->GetYaxis()->SetTitleSize(0.045);
      Tracks_y_high_4->GetYaxis()->SetTitleSize(0.045);
      Tracks_y_high_5->GetYaxis()->SetTitleSize(0.045);
      Tracks_y_high_0->GetYaxis()->SetTitleOffset(1.1);
      Tracks_y_high_1->GetYaxis()->SetTitleOffset(1.1);
      Tracks_y_high_2->GetYaxis()->SetTitleOffset(1.1);
      Tracks_y_high_3->GetYaxis()->SetTitleOffset(1.1);
      Tracks_y_high_4->GetYaxis()->SetTitleOffset(1.1);
      Tracks_y_high_5->GetYaxis()->SetTitleOffset(1.1);

      // Setting titles
      // Low momentum
      Tracks_y_low_0->SetTitle(Efficiency_yz_title_stream_0.str().c_str());
      Tracks_y_low_1->SetTitle(Efficiency_yz_title_stream_1.str().c_str());
      Tracks_y_low_2->SetTitle(Efficiency_yz_title_stream_2.str().c_str());
      Tracks_y_low_3->SetTitle(Efficiency_yz_title_stream_3.str().c_str());
      Tracks_y_low_4->SetTitle(Efficiency_yz_title_stream_4.str().c_str());
      Tracks_y_low_5->SetTitle(Efficiency_yz_title_stream_5.str().c_str());

      // Med momentum
      Tracks_y_med_0->SetTitle(Efficiency_yz_title_stream_0.str().c_str());
      Tracks_y_med_1->SetTitle(Efficiency_yz_title_stream_1.str().c_str());
      Tracks_y_med_2->SetTitle(Efficiency_yz_title_stream_2.str().c_str());
      Tracks_y_med_3->SetTitle(Efficiency_yz_title_stream_3.str().c_str());
      Tracks_y_med_4->SetTitle(Efficiency_yz_title_stream_4.str().c_str());
      Tracks_y_med_5->SetTitle(Efficiency_yz_title_stream_5.str().c_str());

      // High momentum
      Tracks_y_high_0->SetTitle(Efficiency_yz_title_stream_0.str().c_str());
      Tracks_y_high_1->SetTitle(Efficiency_yz_title_stream_1.str().c_str());
      Tracks_y_high_2->SetTitle(Efficiency_yz_title_stream_2.str().c_str());
      Tracks_y_high_3->SetTitle(Efficiency_yz_title_stream_3.str().c_str());
      Tracks_y_high_4->SetTitle(Efficiency_yz_title_stream_4.str().c_str());
      Tracks_y_high_5->SetTitle(Efficiency_yz_title_stream_5.str().c_str());

      // Tracks_yx->GetXaxis()->SetTitleSize(0.045);
      // Tracks_yx->GetYaxis()->SetTitleSize(0.045);
      // Tracks_yx->GetYaxis()->SetTitleOffset(1.1);
      //
      // Tracks_zx->GetXaxis()->SetTitleSize(0.045);
      // Tracks_zx->GetYaxis()->SetTitleSize(0.045);
      // Tracks_zx->GetYaxis()->SetTitleOffset(1.1);

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Setting min, max and ranges for histograms
      // Low momentum
      // Setting Minimum
      Tracks_y_low_0->SetMinimum(0.5);
      Tracks_y_low_1->SetMinimum(0.5);
      Tracks_y_low_2->SetMinimum(0.5);
      Tracks_y_low_3->SetMinimum(0.5);
      Tracks_y_low_4->SetMinimum(0.5);
      Tracks_y_low_5->SetMinimum(0.5);
      // Setting Maximum
      Tracks_y_low_0->SetMaximum(1.1);
      Tracks_y_low_1->SetMaximum(1.1);
      Tracks_y_low_2->SetMaximum(1.1);
      Tracks_y_low_3->SetMaximum(1.1);
      Tracks_y_low_4->SetMaximum(1.1);
      Tracks_y_low_5->SetMaximum(1.1);
      // Setting range
      Tracks_y_low_0->GetXaxis()->SetRangeUser(0,500);
      Tracks_y_low_1->GetXaxis()->SetRangeUser(0,500);
      Tracks_y_low_2->GetXaxis()->SetRangeUser(0,500);
      Tracks_y_low_3->GetXaxis()->SetRangeUser(0,500);
      Tracks_y_low_4->GetXaxis()->SetRangeUser(0,500);
      Tracks_y_low_5->GetXaxis()->SetRangeUser(0,500);

      // Med momentum
      // Setting Minimum
      Tracks_y_med_0->SetMinimum(0.5);
      Tracks_y_med_1->SetMinimum(0.5);
      Tracks_y_med_2->SetMinimum(0.5);
      Tracks_y_med_3->SetMinimum(0.5);
      Tracks_y_med_4->SetMinimum(0.5);
      Tracks_y_med_5->SetMinimum(0.5);
      // Setting Maximum
      Tracks_y_med_0->SetMaximum(1.1);
      Tracks_y_med_1->SetMaximum(1.1);
      Tracks_y_med_2->SetMaximum(1.1);
      Tracks_y_med_3->SetMaximum(1.1);
      Tracks_y_med_4->SetMaximum(1.1);
      Tracks_y_med_5->SetMaximum(1.1);
      // Setting range
      Tracks_y_med_0->GetXaxis()->SetRangeUser(0,500);
      Tracks_y_med_1->GetXaxis()->SetRangeUser(0,500);
      Tracks_y_med_2->GetXaxis()->SetRangeUser(0,500);
      Tracks_y_med_3->GetXaxis()->SetRangeUser(0,500);
      Tracks_y_med_4->GetXaxis()->SetRangeUser(0,500);
      Tracks_y_med_5->GetXaxis()->SetRangeUser(0,500);

      // High momentum
      // Setting Minimum
      Tracks_y_high_0->SetMinimum(0.5);
      Tracks_y_high_1->SetMinimum(0.5);
      Tracks_y_high_2->SetMinimum(0.5);
      Tracks_y_high_3->SetMinimum(0.5);
      Tracks_y_high_4->SetMinimum(0.5);
      Tracks_y_high_5->SetMinimum(0.5);
      // Setting Maximum
      Tracks_y_high_0->SetMaximum(1.1);
      Tracks_y_high_1->SetMaximum(1.1);
      Tracks_y_high_2->SetMaximum(1.1);
      Tracks_y_high_3->SetMaximum(1.1);
      Tracks_y_high_4->SetMaximum(1.1);
      Tracks_y_high_5->SetMaximum(1.1);
      // Setting range
      Tracks_y_high_0->GetXaxis()->SetRangeUser(0,500);
      Tracks_y_high_1->GetXaxis()->SetRangeUser(0,500);
      Tracks_y_high_2->GetXaxis()->SetRangeUser(0,500);
      Tracks_y_high_3->GetXaxis()->SetRangeUser(0,500);
      Tracks_y_high_4->GetXaxis()->SetRangeUser(0,500);
      Tracks_y_high_5->GetXaxis()->SetRangeUser(0,500);

      // Setting colours
      // Med momentum
      Tracks_y_med_0->SetLineColor(kRed);
      Tracks_y_med_1->SetLineColor(kRed);
      Tracks_y_med_2->SetLineColor(kRed);
      Tracks_y_med_3->SetLineColor(kRed);
      Tracks_y_med_4->SetLineColor(kRed);
      Tracks_y_med_5->SetLineColor(kRed);
      // High momentum
      Tracks_y_high_0->SetLineColor(kGreen);
      Tracks_y_high_1->SetLineColor(kGreen);
      Tracks_y_high_2->SetLineColor(kGreen);
      Tracks_y_high_3->SetLineColor(kGreen);
      Tracks_y_high_4->SetLineColor(kGreen);
      Tracks_y_high_5->SetLineColor(kGreen);

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Plotting L vs LPerp efficieny for each sector on canvas

      // Sector 1
      canvas_yz[i_detector][i_charge]->cd(1);
      Tracks_y_low_0->Draw();
      Tracks_y_med_0->Draw("same");
      Tracks_y_high_0->Draw("same");

      // Sector 2
      canvas_yz[i_detector][i_charge]->cd(2);
      Tracks_y_low_1->Draw();
      Tracks_y_med_1->Draw("same");
      Tracks_y_high_1->Draw("same");

      // Sector 3
      canvas_yz[i_detector][i_charge]->cd(3);
      Tracks_y_low_2->Draw();
      Tracks_y_med_2->Draw("same");
      Tracks_y_high_2->Draw("same");

      // Sector 4
      canvas_yz[i_detector][i_charge]->cd(4);
      Tracks_y_low_3->Draw();
      Tracks_y_med_3->Draw("same");
      Tracks_y_high_3->Draw("same");

      // Sector 5
      canvas_yz[i_detector][i_charge]->cd(5);
      Tracks_y_low_4->Draw();
      Tracks_y_med_4->Draw("same");
      Tracks_y_high_4->Draw("same");

      // Sector 6
      canvas_yz[i_detector][i_charge]->cd(6);
      Tracks_y_low_5->Draw();
      Tracks_y_med_5->Draw("same");
      Tracks_y_high_5->Draw("same");



      // // Plotting L vs run efficieny for each sector on canvas
      // canvas_yx[i_detector][i_charge]->cd(i_sector+1);
      // Tracks_yx->SetTitle(Efficiency_yx_title_stream.str().c_str());
      // Tracks_yx->SetLabelSize(0);
      // Tracks_yx->SetMaximum(1);
      // Tracks_yx->Draw("colz");
      //
      // // Plotting L vs run efficiency for each sector on canvas
      // canvas_zx[i_detector][i_charge]->cd(i_sector+1);
      // Tracks_zx->SetTitle(Efficiency_zx_title_stream.str().c_str());
      // Tracks_zx->SetLabelSize(0);
      // Tracks_zx->SetMaximum(1);
      // Tracks_zx->Draw("colz");

      // }
      // Saving the canvases
      canvas_yz[i_detector][i_charge]->SaveAs(Canvas_yz_file_name_stream.str().c_str());
      // canvas_yx[i_detector][i_charge]->SaveAs(Canvas_yx_file_name_stream.str().c_str());
      // canvas_zx[i_detector][i_charge]->SaveAs(Canvas_zx_file_name_stream.str().c_str());
    }
  }
}
