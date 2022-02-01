// Make sure to set the Data, Date and Version string before running

{
  // Give the input analysis file to determine efficieny for
  TFile *f1=new TFile("/media/mn688/Elements1/PhD/FTOF/FTOF_Efficiency_RGB_Spring_2019_63_dst_unified_17012022_01.root");

  //////////////////////////////////////////////////////////////////////////////
  ////Define variables for naming and limits ///////////
  //////////////////////////////////////////////////////////////////////////////

  // Information for canvas and histogram name
  ostringstream Data;
  ostringstream Date;
  ostringstream Version;
  ostringstream Topology;
  ostringstream Layer;
  ostringstream Charge;
  ostringstream Sector;
  // Setting the strings for canvas name
  Data<<"RGB_Spring2019_63_dst";
  Date<<"01022022";
  Version<<"03";




  Int_t detector_min, detector_max; // Limits for FTOF layers
  Int_t charge_min, charge_max; // Limits for charge

  //////////////////////////////////////////////////////////////////////////////
  ////Creating arrays of histograms for the numerator and denominator///////////
  //////////////////////////////////////////////////////////////////////////////

  // Denominator
  TH1F *Denominator_y[4][3][2][6];
  TH1F *pi_Denominator_y[3][2][6];
  TH1F *Overall_Denominator_y[3][2][6];

  // Numerator
  TH1F *Numerator_y[4][3][2][6];
  TH1F *pi_Numerator_y[3][2][6];
  TH1F *Overall_Numerator_y[3][2][6];


  //////////////////////////////////////////////////////////////////////////////
  ////Creating arrays of canvases for efficieny plots //////////////////////////
  //////////////////////////////////////////////////////////////////////////////

  // Creating arrays of canvases [FTOF layer] [charge]
  TCanvas* canvas_y[4][3][2]; // canvases for efficieny vs L

  gStyle->SetOptStat(0); // Remove statistics from canvases

  //////////////////////////////////////////////////////////////////////////////
  ////Loop over topology, layer, charge, and sector ///////////
  //////////////////////////////////////////////////////////////////////////////

  // Looping over topology
  for(Int_t i_topology=0;i_topology < 5;i_topology++){

    // Defining topology specific strings and parameters

    // Missing pi^-
    if(i_topology == 0) {
      // Limiting the charge and layers to look at based on topology
      charge_min = 0; // Start at negatives
      charge_max = 1; // Only do negatives for pi-
      detector_min = 0; // Start at 1A
      detector_max = 3; // Do all layers (1A, 1B, and 2) for pi-

      // Name the topology
      Topology.str("pim");

    }

    // Missing pi^+ and missing proton
    else if(i_topology == 1 || i_topology == 2) {
      // Limiting the charge and layers to look at based on topology
      charge_min = 1; // Start at positives
      charge_max = 2; // Only do positives for pi+ and p
      detector_min = 0; // Start at 1A
      detector_max = 3; // Do all layers (1A, 1B, and 2) for pi+ and p

      // Name the topology
      if(i_topology == 1) Topology.str("pip");
      if(i_topology == 2) Topology.str("proton");
    }

    // Single track events
    else if(i_topology == 3) {
      // Limiting the charge and layers to look at based on topology
      charge_min = 0; // Start at negatives
      charge_max = 2; // Do positives and negatives for single track
      detector_min = 0; // Start at 1A
      detector_max = 2; // Only do 1A and 1B for single track

      // Name the topology
      Topology.str("single_track");
    }

    // Overall result combining single track and 2 pi method
    else if(i_topology == 4){
      // Limiting the charge and layers to look at based on topology
      charge_min = 0; // Start at negatives
      charge_max = 2; // Do positives and negatives for overall
      detector_min = 0; // Start at 1A
      detector_max = 3; // Do all layers (1A, 1B, and 2) for overall

      // Name the topology
      Topology.str("overall");
    }

    // Looping over the FTOF layers
    for(Int_t i_detector = detector_min; i_detector < detector_max; i_detector++){

      // Name the FTOF layer
      if(i_detector == 0) Layer.str("1A");
      else if(i_detector == 1) Layer.str("1B");
      else if(i_detector == 2) Layer.str("2");


      // Looping over charge
      for(Int_t i_charge = charge_min; i_charge < charge_max; i_charge++){

        // Name the charge
        if(i_charge == 0) Charge.str("negatives");
        else if(i_charge == 1) Charge.str("positives");

        // Creating strings for canvas names and file names
        ostringstream canvas_y_name_stream;
        ostringstream canvas_y_file_name_stream;


        // Setting the name for each canvas
        canvas_y_name_stream<<Topology.str().c_str()<<"_"<<Layer.str().c_str()<<"_"<<Charge.str().c_str()<<"_L";

        // Setting the file name for each canvas
        canvas_y_file_name_stream<<"/media/mn688/Elements1/PhD/FTOF/Efficiency_Canvases/RGB_Spring_2019_Inbending_01022022/"<<Data.str().c_str()<<"_"<<Topology.str().c_str()<<"_"<<Layer.str().c_str()<<"_"<<Charge.str().c_str()<<"_"<<Date.str().c_str()<<"_"<<Version.str().c_str()<<".pdf";


        // creating canvas for each topology, FTOF layer, and charge
        canvas_y[i_topology][i_detector][i_charge] = new TCanvas(canvas_y_name_stream.str().c_str(),"", 2500,1600);
        canvas_y[i_topology][i_detector][i_charge]->Divide(3,2);


        // Looping over the sectors
        for(Int_t i_sector = 0; i_sector < 6; i_sector++){

          // Printing out the loop
          cout<<"topology "<<i_topology<<" detector "<<i_detector<<" charge "<<i_charge<<" sector "<<i_sector<<endl;

          // Creating strings for the histograms and canvases
          ostringstream Denominator_name_stream, Numerator_name_stream;
          ostringstream pi_Denominator_name_stream, pi_Numerator_name_stream;
          ostringstream Efficiency_title_stream;
          ostringstream Efficiency_overall_title_stream;

          //////////////////////////////////////////////////////////////////////////////
          // Setting the histogram names  //////////////////////////////////////////////
          //////////////////////////////////////////////////////////////////////////////

          // For 2pi and single track topologies
          if(i_topology != 4){
            Denominator_name_stream<<"h_Denominator_Topo_"<<i_topology<<"_Det_"<<i_detector<<"_Charge_"<<i_charge<<"_Sec_"<<i_sector;
            Numerator_name_stream<<"h_Numerator_Topo_"<<i_topology<<"_Det_"<<i_detector<<"_Charge_"<<i_charge<<"_Sec_"<<i_sector;
          }

          // For overall result
          else{

            // FTOF 2
            if(i_detector == 2){

              // negatively charged particles
              if(i_charge == 0) {

                // Give histogram name from pim
                Denominator_name_stream<<"h_Denominator_Topo_0_Det_"<<i_detector<<"_Charge_"<<i_charge<<"_Sec_"<<i_sector;
                Numerator_name_stream<<"h_Numerator_Topo_0_Det_"<<i_detector<<"_Charge_"<<i_charge<<"_Sec_"<<i_sector;

              }
              // positively charge particles
              else{
                // Give histogram name from pip
                Denominator_name_stream<<"h_Denominator_Topo_1_Det_"<<i_detector<<"_Charge_"<<i_charge<<"_Sec_"<<i_sector;
                Numerator_name_stream<<"h_Numerator_Topo_1_Det_"<<i_detector<<"_Charge_"<<i_charge<<"_Sec_"<<i_sector;
              }
            }

            // FTOF 1A and 1B
            else{
              // Give histogram name from single track
              Denominator_name_stream<<"h_Denominator_Topo_3_Det_"<<i_detector<<"_Charge_"<<i_charge<<"_Sec_"<<i_sector;
              Numerator_name_stream<<"h_Numerator_Topo_3_Det_"<<i_detector<<"_Charge_"<<i_charge<<"_Sec_"<<i_sector;

              // Negatives
              if(i_charge == 0){
                // Get histogram name from missing pi-
                pi_Denominator_name_stream<<"h_Denominator_Topo_0_Det_"<<i_detector<<"_Charge_"<<i_charge<<"_Sec_"<<i_sector;
                pi_Numerator_name_stream<<"h_Numerator_Topo_0_Det_"<<i_detector<<"_Charge_"<<i_charge<<"_Sec_"<<i_sector;
              }
              // Positives
              else{
                // Get histogram name from missing pi+
                pi_Denominator_name_stream<<"h_Denominator_Topo_1_Det_"<<i_detector<<"_Charge_"<<i_charge<<"_Sec_"<<i_sector;
                pi_Numerator_name_stream<<"h_Numerator_Topo_1_Det_"<<i_detector<<"_Charge_"<<i_charge<<"_Sec_"<<i_sector;
              }
            }
          }

          //////////////////////////////////////////////////////////////////////////////
          // Setting the efficieny titles //////////////////////////////////////////////
          //////////////////////////////////////////////////////////////////////////////

          // For 2 pi topologies
          if(i_topology < 3){
            Efficiency_title_stream<<Topology.str().c_str()<<" Efficiency of FTOF"<<Layer.str().c_str()<<" Sector "<<i_sector + 1<<"; Counter; #epsilon";
          }
          // single track topology
          else if(i_topology == 3){
            Efficiency_title_stream<<Topology.str().c_str()<<" Efficiency of "<<Charge.str().c_str()<<" FTOF"<<Layer.str().c_str()<<" Sector "<<i_sector + 1<<"; Counter; #epsilon";
          }
          // overall result
          else if(i_topology == 4){
            Efficiency_overall_title_stream<<"Overall efficiency of "<<Charge.str().c_str()<<" FTOF"<<Layer.str().c_str()<<" Sector "<<i_sector + 1<<"; Counter; #epsilon";
          }

          //////////////////////////////////////////////////////////////////////////////
          // Getting histograms from input file  ///////////////////////////////////////
          //////////////////////////////////////////////////////////////////////////////

          // Getting the histograms from data file
          TH3F *h_Denominator = (TH3F*) f1->Get(Denominator_name_stream.str().c_str());
          TH3F *h_Numerator = (TH3F*) f1->Get(Numerator_name_stream.str().c_str());

          // Take 1d projection
          Denominator_y[i_topology][i_detector][i_charge][i_sector] = (TH1F*) h_Denominator->Project3D("y"); // denominator
          Numerator_y[i_topology][i_detector][i_charge][i_sector] = (TH1F*) h_Numerator->Project3D("y"); // numerator

          // For overall topology layer 1A and 1B
          if(i_topology == 4 && (i_detector == 0 || i_detector == 1)){

            // Also get the 2pi histogram to use for the parts not covered by
            // calorimeter
            TH3F *pi_Denominator = (TH3F*) f1->Get(pi_Denominator_name_stream.str().c_str());
            TH3F *pi_Numerator = (TH3F*) f1->Get(pi_Numerator_name_stream.str().c_str());

            // Take the 1D projection
            pi_Denominator_y[i_detector][i_charge][i_sector] = (TH1F*) pi_Denominator->Project3D("y");
            pi_Numerator_y[i_detector][i_charge][i_sector] = (TH1F*) pi_Numerator->Project3D("y");

            // Remove the 3D histograms once the 1D projections are taken
            delete pi_Denominator;
            delete pi_Numerator;
          }

          //////////////////////////////////////////////////////////////////////////////
          // Get results from 2pi and single track method for overall efficieny ////////
          //////////////////////////////////////////////////////////////////////////////

          if(i_topology == 4){

            // Creating strings for names and titles of histograms for overall efficieny
            ostringstream histogram_name_num, histogram_title_num;
            ostringstream histogram_name_den, histogram_title_den;

            // Setting the names and titles of histograms for overall efficieny
            histogram_name_num<<"overall_hist_num_"<<i_detector<<"_"<<i_charge<<"_"<<i_sector;
            histogram_title_num<<"overall_hist_num_"<<i_detector<<"_"<<i_charge<<"_"<<i_sector;
            histogram_name_den<<"overall_hist_den_"<<i_detector<<"_"<<i_charge<<"_"<<i_sector;
            histogram_title_den<<"overall_hist_den_"<<i_detector<<"_"<<i_charge<<"_"<<i_sector;

            // For FTOF1A and FTOF1B get single track result where calorimeter
            // covers and 2 pi values where it doesn't
            if(i_detector == 0 || i_detector == 1){

              // Change bin content for those outside of calorimeter acceptance
              // For FTOF 1A change last counter to 2 pi method
              if(i_detector == 0){
                // Define the overall histograms
                Overall_Numerator_y[i_detector][i_charge][i_sector] = new TH1F(histogram_name_num.str().c_str(),histogram_title_num.str().c_str(),25, 61.98, 441.48);
                Overall_Denominator_y[i_detector][i_charge][i_sector] = new TH1F(histogram_name_den.str().c_str(),histogram_title_den.str().c_str(),25, 61.98, 441.48);


                for(Int_t bin_0 = 2; bin_0 < 25; bin_0++){
                  // These bins are covered by the calorimeter, so single track
                  // is used
                  if(bin_0 < 24){
                    // Getting denominator values from single track
                    Overall_Denominator_y[i_detector][i_charge][i_sector]->SetBinContent
                    (bin_0,Denominator_y[3][i_detector][i_charge][i_sector]->GetBinContent(bin_0));

                    // Getting numerator values from single track
                    Overall_Numerator_y[i_detector][i_charge][i_sector]->SetBinContent
                    (bin_0,Numerator_y[3][i_detector][i_charge][i_sector]->GetBinContent(bin_0));
                  }
                  // These bins are not covered by the calorimeter, so 2pi
                  // is used
                  else{
                    // Getting denominator values from 2 pi
                    Overall_Denominator_y[i_detector][i_charge][i_sector]->SetBinContent
                    (bin_0,pi_Denominator_y[i_detector][i_charge][i_sector]->GetBinContent(bin_0));

                    // Getting numerator values from 2pi
                    Overall_Numerator_y[i_detector][i_charge][i_sector]->SetBinContent
                    (bin_0,pi_Numerator_y[i_detector][i_charge][i_sector]->GetBinContent(bin_0));
                  }
                }
              }

              // For FTOF 1B change last 5 counters to 2 pi method
              if(i_detector == 1){
                // Define the overall histograms
                Overall_Numerator_y[i_detector][i_charge][i_sector] = new TH1F(histogram_name_num.str().c_str(),histogram_title_num.str().c_str(),65, 46.04, 443.84);
                Overall_Denominator_y[i_detector][i_charge][i_sector] = new TH1F(histogram_title_den.str().c_str(),histogram_title_den.str().c_str(),65, 46.04, 443.84);

                for(Int_t bin_1 = 2; bin_1 < 65; bin_1++){
                  // These bins are covered by the calorimeter, so single track
                  // is used
                  if(bin_1 < 57){
                    // Getting denominator values from single track
                    Overall_Denominator_y[i_detector][i_charge][i_sector]->SetBinContent
                    (bin_1,Denominator_y[3][i_detector][i_charge][i_sector]->GetBinContent(bin_1));

                    // Getting numerator values from single track
                    Overall_Numerator_y[i_detector][i_charge][i_sector]->SetBinContent
                    (bin_1,Numerator_y[3][i_detector][i_charge][i_sector]->GetBinContent(bin_1));
                  }

                  else{
                    // Getting denominator values from 2pi
                    Overall_Denominator_y[i_detector][i_charge][i_sector]->SetBinContent
                    (bin_1,pi_Denominator_y[i_detector][i_charge][i_sector]->GetBinContent(bin_1));

                    // Getting numerator values from 2pi
                    Overall_Numerator_y[i_detector][i_charge][i_sector]->SetBinContent
                    (bin_1,pi_Numerator_y[i_detector][i_charge][i_sector]->GetBinContent(bin_1));
                  }
                }
              }
            }

            // FTOF 2
            if(i_detector == 2){
              // Define overall histograms
              Overall_Numerator_y[i_detector][i_charge][i_sector] = new TH1F(histogram_name_num.str().c_str(),histogram_title_num.str().c_str(),7, 693, 847);
              Overall_Denominator_y[i_detector][i_charge][i_sector] = new TH1F(histogram_name_den.str().c_str(),histogram_title_den.str().c_str(),7, 693, 847);

              // Loop over counters
              for(Int_t bin_2 = 0; bin_2 < 7; bin_2++){

                // For negative efficieny get the pi- results
                if(i_charge == 0){
                  // Getting denominator from 2pi
                  Overall_Denominator_y[i_detector][i_charge][i_sector]->SetBinContent
                  (bin_2,Denominator_y[0][i_detector][i_charge][i_sector]->GetBinContent(bin_2));

                  // Getting numerator from 2pi
                  Overall_Numerator_y[i_detector][i_charge][i_sector]->SetBinContent
                  (bin_2,Numerator_y[0][i_detector][i_charge][i_sector]->GetBinContent(bin_2));
                }

                // For positive efficieny get the pi+ results
                else{
                  // Getting denominator from 2pi
                  Overall_Denominator_y[i_detector][i_charge][i_sector]->SetBinContent
                  (bin_2,Denominator_y[1][i_detector][i_charge][i_sector]->GetBinContent(bin_2));

                  // Getting numerator from 2pi
                  Overall_Numerator_y[i_detector][i_charge][i_sector]->SetBinContent
                  (bin_2,Numerator_y[1][i_detector][i_charge][i_sector]->GetBinContent(bin_2));
                }
              }
            }

            // Determine errors for overall numerator
            Overall_Numerator_y[i_detector][i_charge][i_sector]->Sumw2();
            // Get ratio between numerator and denominator for efficieny
            Overall_Numerator_y[i_detector][i_charge][i_sector]->Divide(Overall_Denominator_y[i_detector][i_charge][i_sector]);
          }


          //////////////////////////////////////////////////////////////////////////////
          // Determine efficieny and error for 2pi and single track method /////////////
          //////////////////////////////////////////////////////////////////////////////

          // 2pi and single track
          else{
            // Determine errors
            Numerator_y[i_topology][i_detector][i_charge][i_sector]->Sumw2();
          }
          // Dividing the numerator histogram by the denominator histogram
          Numerator_y[i_topology][i_detector][i_charge][i_sector]->Divide(Denominator_y[i_topology][i_detector][i_charge][i_sector]);

          //////////////////////////////////////////////////////////////////////
          // Setting presentation on canvas   //////////////////////////////////
          //////////////////////////////////////////////////////////////////////

          // // Setting titles
          Numerator_y[i_topology][i_detector][i_charge][i_sector]->SetTitle(Efficiency_title_stream.str().c_str());

          // Aligning Title
          Numerator_y[i_topology][i_detector][i_charge][i_sector]->GetXaxis()->SetTitleSize(0.045);
          Numerator_y[i_topology][i_detector][i_charge][i_sector]->GetYaxis()->SetTitleSize(0.045);
          Numerator_y[i_topology][i_detector][i_charge][i_sector]->GetYaxis()->SetTitleOffset(1.1);


          // Setting minimum and maximum
          Numerator_y[i_topology][i_detector][i_charge][i_sector]->SetMinimum(0.75);
          Numerator_y[i_topology][i_detector][i_charge][i_sector]->SetMaximum(1.1);

          // Setting ranges
          // FTOF 1
          if(i_detector<2)  Numerator_y[i_topology][i_detector][i_charge][i_sector]->GetXaxis()->SetRangeUser(0,500);

          // FTOF 2
          else  Numerator_y[i_topology][i_detector][i_charge][i_sector]->GetXaxis()->SetRangeUser(600,900);

          // Setting colour
          Numerator_y[i_topology][i_detector][i_charge][i_sector]->SetLineColor(kBlue);


          // Set presenation on canvas for overall plots
          if(i_topology == 4 ){
            // // Setting titles
            Overall_Numerator_y[i_detector][i_charge][i_sector]->SetTitle(Efficiency_overall_title_stream.str().c_str());

            // Aligning Title
            Overall_Numerator_y[i_detector][i_charge][i_sector]->GetXaxis()->SetTitleSize(0.045);
            Overall_Numerator_y[i_detector][i_charge][i_sector]->GetYaxis()->SetTitleSize(0.045);
            Overall_Numerator_y[i_detector][i_charge][i_sector]->GetYaxis()->SetTitleOffset(1.1);

            // Setting minimum and maximum
            Overall_Numerator_y[i_detector][i_charge][i_sector]->SetMinimum(0.75);
            Overall_Numerator_y[i_detector][i_charge][i_sector]->SetMaximum(1.1);

            // Setting ranges
            // FTOF 1
            if(i_detector < 2)
              Overall_Numerator_y[i_detector][i_charge][i_sector]->GetXaxis()->SetRangeUser(0,500);

            // FTOF 2
            else
              Overall_Numerator_y[i_detector][i_charge][i_sector]->GetXaxis()->SetRangeUser(600,900);

            // Setting colour
            Overall_Numerator_y[i_detector][i_charge][i_sector]->SetLineColor(kBlue);

          }


          //////////////////////////////////////////////////////////////////////
          // Renaming the bin labels to match the counters  ///////////
          //////////////////////////////////////////////////////////////////////

          // Loop over bins in histogram
          for(Int_t l = 2; l < Numerator_y[i_topology][i_detector][i_charge][i_sector]->GetNbinsX(); l++){

            // Changing the bin labels to match the counter numbers
            ostringstream test;
            test<< to_string(l-1);

            // Single track and 2pi
            if(i_topology < 4)
            Numerator_y[i_topology][i_detector][i_charge][i_sector]->GetXaxis()->SetBinLabel(l,test.str().c_str());
            // overall method
            else
            Overall_Numerator_y[i_detector][i_charge][i_sector]->GetXaxis()->SetBinLabel(l,test.str().c_str());
          }

          //////////////////////////////////////////////////////////////////////
          // Plotting efficieny against L for each sector on canvas  ///////////
          //////////////////////////////////////////////////////////////////////

          // Selecting correct canvas and drawing histogram
          canvas_y[i_topology][i_detector][i_charge]->cd(i_sector+1);

          // Single track and 2pi
          if(i_topology!=4){
            Numerator_y[i_topology][i_detector][i_charge][i_sector]->Draw();
          }

          // overall method
          else{
            Overall_Numerator_y[i_detector][i_charge][i_sector]->Draw();
          }

          // Remove the 3D histograms
          delete h_Numerator;
          delete h_Denominator;
        }
        // Saving the canvases
        canvas_y[i_topology][i_detector][i_charge]->SaveAs(canvas_y_file_name_stream.str().c_str());

        // Clearing all ostringstreams
        Data.clear();
        Date.clear();
        Version.clear();
        Topology.clear();
        Layer.clear();
        Charge.clear();

      }
    }
  }
}
