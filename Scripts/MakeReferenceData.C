#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

enum e_Planes {kPlane0,kPlane1,kPlane2};
std::vector<std::string> str_Planes = {"Plane0","Plane1","Plane2"};

const std::vector<std::vector<double>> dEdx_bins = {
  {0.000,0.500,1.000,1.500,2.000,2.500,3.000,3.500,4.000,4.500,5.000,5.500,6.000,6.500,7.000,7.500,8.000,9.000,10.000,12.000,15.000,20.000,25.000,30.000,35.000,40.000,50.000},
  {0.000,0.500,1.000,1.500,2.000,2.500,3.000,3.500,4.000,4.500,5.000,5.500,6.000,6.500,7.000,7.500,8.000,9.000,10.000,12.000,15.000,20.000,25.000,30.000,35.000,40.000,50.000},
  {0.000,0.500,1.000,1.500,2.000,2.500,3.000,3.500,4.000,4.500,5.000,5.500,6.000,6.500,7.000,7.500,8.000,9.000,10.000,12.000,15.000,20.000,25.000,30.000,35.000,40.000,50.000}
};

const std::vector<std::vector<double>> ResRange_bins = {
  {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30},
  {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30},
  {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30}
};

const std::vector<std::vector<double>> Angle_bins = {
  {0.0,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90},
  {0.0,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90},
  {0.0,2.5,5,7.5,10,12.5,15,20,30,40,50,90}
};

const std::vector<int> SupportedPDGs = {3222,3112,321,2212,13,211};
const std::vector<std::string> SupportedPDGs_str = {"SigmaP","SigmaM","Kaon","Proton","Muon","Pion"};

const float TrackAngleTrueRecoCut = 30; // Maximum angle between true and reco directions of tracks to be used to make reference
const float TrackLengthCut = 2; // Minimum length of track to be used to generate reference

void MakeReferenceData(){

  // Setup the histograms
  std::vector<std::vector<TH3D*>> h_dEdx_ResidualRange_Angle_v(3); 

  // 1D histograms of each variable to help guide binning choices
  std::vector<std::vector<TH1D*>> h_dEdx_v(3),h_ResidualRange_v(3),h_Angle_v(3);

  for(int i_pl=0;i_pl<3;i_pl++){  
    std::string plane = str_Planes.at(i_pl);
    const int angle_bins = Angle_bins.at(i_pl).size()-1; 
    const int rr_bins = ResRange_bins.at(i_pl).size()-1; 
    const int dedx_bins = dEdx_bins.at(i_pl).size()-1; 
    std::vector<double> tmp = Angle_bins.at(i_pl); double* angle_binning = &tmp[0]; // apparently this can't be done as a single line  
    std::vector<double> tmp2 = dEdx_bins.at(i_pl); double* dedx_binning = &tmp2[0];   
    std::vector<double> tmp3 = ResRange_bins.at(i_pl); double* rr_binning = &tmp3[0];   

    for(size_t i_pdg=0;i_pdg<SupportedPDGs.size();i_pdg++){
      std::string pdg = SupportedPDGs_str.at(i_pdg);
      h_dEdx_ResidualRange_Angle_v.at(i_pl).push_back(new TH3D(("h_dEdx_ResidualRange_Angle_"+plane+"_"+pdg).c_str(),(pdg+";Residual Range (cm);dE/dx (MeV/cm);Angle (deg)").c_str(),rr_bins,rr_binning,dedx_bins,dedx_binning,angle_bins,angle_binning));
      h_dEdx_v.at(i_pl).push_back(new TH1D(("h_dEdx_"+plane+"_"+pdg).c_str(),";dE/dx (MeV/cm);",dedx_bins,dedx_binning));
      h_ResidualRange_v.at(i_pl).push_back(new TH1D(("h_ResidualRange_"+plane+"_"+pdg).c_str(),";Residual Range (cm);",rr_bins,rr_binning));  
      h_Angle_v.at(i_pl).push_back(new TH1D(("h_Angle_"+plane+"_"+pdg).c_str(),";Angle (deg);",angle_bins,angle_binning));
    } 
  }


  // Load the trees containing the data
  TFile* f_in = TFile::Open("dEdxTrees.root");
  TTree* t_in = static_cast<TTree*>(f_in->Get("ana/OutputTree"));
  vector<int>     *TrackTruePDG=0;
  vector<float>   *TrackLength=0;
  vector<float>   *TrackAngleTrueReco=0;
  vector<vector<vector<float>>*> ResidualRange(3,0);
  vector<vector<vector<float>>*> dEdx(3,0);
  vector<vector<vector<float>>*> Pitch(3,0);
  t_in->SetBranchAddress("TrackTruePDG",&TrackTruePDG);
  t_in->SetBranchAddress("TrackLength",&TrackLength);
  t_in->SetBranchAddress("TrackAngleTrueReco",&TrackAngleTrueReco);
  t_in->SetBranchAddress("ResidualRange_Plane0",&ResidualRange.at(kPlane0));
  t_in->SetBranchAddress("dEdx_Plane0",&dEdx.at(kPlane0));
  t_in->SetBranchAddress("Pitch_Plane0",&Pitch.at(kPlane0));
  t_in->SetBranchAddress("ResidualRange_Plane1",&ResidualRange.at(kPlane1));
  t_in->SetBranchAddress("dEdx_Plane1",&dEdx.at(kPlane1));
  t_in->SetBranchAddress("Pitch_Plane1",&Pitch.at(kPlane1));
  t_in->SetBranchAddress("ResidualRange_Plane2",&ResidualRange.at(kPlane2));
  t_in->SetBranchAddress("dEdx_Plane2",&dEdx.at(kPlane2));
  t_in->SetBranchAddress("Pitch_Plane2",&Pitch.at(kPlane2));

  Long64_t nentries = t_in->GetEntriesFast();

  // Fill the histograms
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    t_in->GetEntry(jentry);

    if(jentry % 10000 == 0) std::cout << "Event " << jentry << "/" << nentries << std::endl;

    for(size_t i_tr=0;i_tr<TrackTruePDG->size();i_tr++){

      // Apply quality cuts
      if(TrackLength->at(i_tr) < TrackLengthCut) continue;
      if(TrackAngleTrueReco->at(i_tr) > TrackAngleTrueRecoCut) continue;

      // Only use tracks in the list of pdgs to make reference for
      int pdg_index = -1;
      for(size_t i_pdg=0;i_pdg<SupportedPDGs.size();i_pdg++) 
        if(SupportedPDGs.at(i_pdg) == TrackTruePDG->at(i_tr))
          pdg_index = i_pdg;
      if(pdg_index == -1) continue;

      // Fill histograms
      for(size_t i_pl=0;i_pl<ResidualRange.size();i_pl++){
        std::vector<float>* rr = &ResidualRange.at(i_pl)->at(i_tr);
        std::vector<float>* pitch = &Pitch.at(i_pl)->at(i_tr);
        std::vector<float>* dedx = &dEdx.at(i_pl)->at(i_tr);
        for(size_t i_p=0;i_p<rr->size();i_p++){
          h_dEdx_ResidualRange_Angle_v.at(i_pl).at(pdg_index)->Fill(rr->at(i_p),dedx->at(i_p),180/3.142*acos(0.3/pitch->at(i_p)));
          h_dEdx_v.at(i_pl).at(pdg_index)->Fill(dedx->at(i_p));
          h_ResidualRange_v.at(i_pl).at(pdg_index)->Fill(rr->at(i_p));
          h_Angle_v.at(i_pl).at(pdg_index)->Fill(180/3.142*acos(0.3/pitch->at(i_p)));
        }
      }

    }

  }

  f_in->Close();

  // Normalise all of the histograms to 1 to create likelihood profile in 3D - want each slice in rr/pitch space to have pm of 1
  for(size_t i_pl=0;i_pl<h_dEdx_ResidualRange_Angle_v.size();i_pl++){
    for(size_t i_pdg=0;i_pdg<SupportedPDGs.size();i_pdg++){
      TH3D* h = h_dEdx_ResidualRange_Angle_v.at(i_pl).at(i_pdg);
      for(int i_rr=1;i_rr<h->GetNbinsX()+1;i_rr++){
        for(int i_p=1;i_p<h->GetNbinsZ()+1;i_p++){
          // Calculate integral of each slice of angle/rr space
          double mass = 0.0;
          for(int i_dd=0;i_dd<h->GetNbinsY()+1;i_dd++) mass += h->GetBinContent(i_rr,i_dd,i_p);
          for(int i_dd=0;i_dd<h->GetNbinsY()+1;i_dd++) h->SetBinContent(i_rr,i_dd,i_p,h->GetBinContent(i_rr,i_dd,i_p)/mass);
        }
      }
    }
  }

  // Draw plots of each slice in pitch space for each particle
  TCanvas* c = new TCanvas("c","c");
  gSystem->Exec("mkdir -p Plots/");

  for(size_t i_pl=0;i_pl<h_dEdx_ResidualRange_Angle_v.size();i_pl++){

    std::string plane = str_Planes.at(i_pl); 

    // Draw this 2D histograms showing the dE/dx vs RR curves for each slice in
    // pitch space
    for(size_t i_pdg=0;i_pdg<SupportedPDGs.size();i_pdg++){

      std::string particle = SupportedPDGs_str.at(i_pdg); 
      TH3D h_3d = *h_dEdx_ResidualRange_Angle_v.at(i_pl).at(i_pdg); 

      for(int i_z=1;i_z<h_3d.GetNbinsZ()+1;i_z++){
        h_3d.GetZaxis()->SetRange(i_z,i_z);
        TH2D* h_proj = static_cast<TH2D*>(h_3d.Project3D("yx"));
        std::string title = particle + " , Angle (" + std::to_string(h_3d.GetZaxis()->GetBinLowEdge(i_z)) + "," + std::to_string(h_3d.GetZaxis()->GetBinLowEdge(i_z+1)) + ")";
        std::string label = "Plots/" + particle + "_" + std::to_string(i_z) + "_" + plane + ".png";
        h_proj->Draw("colz");
        h_proj->SetStats(0);
        h_proj->SetTitle(title.c_str());
        c->Print(label.c_str());
        c->Clear();
      }

      h_dEdx_v.at(i_pl).at(i_pdg)->SetLineColor(1);
      h_dEdx_v.at(i_pl).at(i_pdg)->SetLineWidth(2);
      h_dEdx_v.at(i_pl).at(i_pdg)->Draw("e1");
      h_dEdx_v.at(i_pl).at(i_pdg)->SetStats(0);   
      c->Print(("Plots/dEdx_" + particle + "_" + plane + ".png").c_str());  
      c->Clear(); 

      h_ResidualRange_v.at(i_pl).at(i_pdg)->SetLineColor(1);
      h_ResidualRange_v.at(i_pl).at(i_pdg)->SetLineWidth(2);
      h_ResidualRange_v.at(i_pl).at(i_pdg)->Draw("e1");
      h_ResidualRange_v.at(i_pl).at(i_pdg)->SetStats(0);   
      c->Print(("Plots/ResidualRange_" + particle + "_" + plane + ".png").c_str());  
      c->Clear(); 

      h_Angle_v.at(i_pl).at(i_pdg)->SetLineColor(1);
      h_Angle_v.at(i_pl).at(i_pdg)->SetLineWidth(2);
      h_Angle_v.at(i_pl).at(i_pdg)->Draw("e1");
      h_Angle_v.at(i_pl).at(i_pdg)->SetStats(0);   
      c->Print(("Plots/Angle_" + particle + "_" + plane + ".png").c_str());  
      c->Clear(); 

    }

  }

  // Write the histograms to file to serve as reference tables
  TFile* f_out = TFile::Open("dEdx_Reference.root","RECREATE");  
  for(size_t i_pl=0;i_pl<h_dEdx_ResidualRange_Angle_v.size();i_pl++){
    std::string plane = str_Planes.at(i_pl);
    for(size_t i_pdg=0;i_pdg<SupportedPDGs.size();i_pdg++){
      const string pdg = std::to_string(SupportedPDGs.at(i_pdg)); 
      h_dEdx_ResidualRange_Angle_v.at(i_pl).at(i_pdg)->Write((pdg + "_" + plane).c_str());
    }
  }
  f_out->Close();

}
