#ifndef _PIDCalc_cxx_
#define _PIDCalc_cxx_

#include "ubana/FlexiPID/Alg/PIDCalc.h"

using namespace FlexiPID;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

PIDCalc::PIDCalc(const fhicl::ParameterSet& p) :
  PIDReferenceHists(p.get<std::string>("PIDReferenceHists")),
  SupportedPDGs(p.get<std::vector<int>>("SupportedPDGs")),
  Debug(p.get<bool>("Debug",false))
{

  TFile* f = TFile::Open(PIDReferenceHists.c_str());

  h_dEdx_Reference.resize(kInvalid);

  for(int i_pl=0;i_pl<kInvalid;i_pl++){
    std::string plane = "Plane" + std::to_string(i_pl);
    for(int pdg : SupportedPDGs){
      std::string name = std::to_string(pdg) + "_" + plane;
      h_dEdx_Reference.at(i_pl).push_back(static_cast<TH3D*>(f->Get(name.c_str()))); 
    }
  }

  f->Close();

  if(Debug) std::cout << "PIDCalc: Loaded PID reference hists from " << PIDReferenceHists << std::endl;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

double PIDCalc::GetGenericLLRPID(std::vector<art::Ptr<anab::Calorimetry>> calo_v,std::pair<int,int> hypotheses) const {

  if(PIDReferenceHists == "")
    throw cet::exception("PIDCalc:") << "Generic LLR PID reference hists not loaded" << std::endl;

  int pdg_index_first = -1;
  int pdg_index_second = -1;
  for(size_t i_pdg=0;i_pdg<SupportedPDGs.size();i_pdg++){
    if(SupportedPDGs.at(i_pdg) == hypotheses.first) pdg_index_first = i_pdg;
    if(SupportedPDGs.at(i_pdg) == hypotheses.second) pdg_index_second = i_pdg;
  }

  // TODO: Add a check to catch if hypothesis isn't in reference

  double llr = 0;

  for(auto const &calo : calo_v){

    auto const &plane = calo->PlaneID().Plane;

    if(plane != kPlane0 && plane != kPlane1 && plane != kPlane2) continue;

    TH3D* h_hyp_first = h_dEdx_Reference.at(plane).at(pdg_index_first);
    TH3D* h_hyp_second = h_dEdx_Reference.at(plane).at(pdg_index_second);

    auto const &dedx_values = calo->dEdx();
    auto const &rr = calo->ResidualRange();
    auto const &pitch = calo->TrkPitchVec();
    for(size_t i_p=0;i_p<dedx_values.size();i_p++){ 
      int bin_x = h_hyp_first->GetXaxis()->FindBin(rr.at(i_p));
      int bin_y = h_hyp_first->GetYaxis()->FindBin(dedx_values.at(i_p));
      int bin_z = h_hyp_first->GetZaxis()->FindBin(180/3.142*acos(0.3/pitch.at(i_p)));
      double l_first = h_hyp_first->GetBinContent(bin_x,bin_y,bin_z);
      double l_second = h_hyp_second->GetBinContent(bin_x,bin_y,bin_z);
      //std::cout << l_first <<  "  " << l_second << std::endl;
      if(l_first > 0 && l_second > 0 && !std::isnan(l_first) && !std::isnan(l_second)) llr += log(l_first) - log(l_second);
    } 
  }

  //std::cout << "Score: " << 2/3.1415*atan(llr) << std::endl;

  return 2.0/3.1415*atan(llr);

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif
