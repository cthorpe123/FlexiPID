#ifndef _PIDCalc_h_
#define _PIDCalc_h_

#include <string>
#include <vector>
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "cetlib_except/exception.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include <TH3.h>
#include <TFile.h>

namespace FlexiPID {

  class PIDCalc {

    public: 

      PIDCalc(const fhicl::ParameterSet& p);

      double GetGenericLLRPID(std::vector<art::Ptr<anab::Calorimetry>> calo_v,std::pair<int,int> hypotheses) const;

    private:

      enum Planes {kPlane0,kPlane1,kPlane2,kInvalid};

      const std::string PIDReferenceHists;
      const std::vector<int> SupportedPDGs;

      std::vector<std::vector<TH3D*>> h_dEdx_Reference;

      void LoadGenericLLRPID();

  };

}

#endif
