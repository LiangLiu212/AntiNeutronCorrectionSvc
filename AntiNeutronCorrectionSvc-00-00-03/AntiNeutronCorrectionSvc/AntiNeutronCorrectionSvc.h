//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++        Anti-neutron Correction 						++
//++ 	author Liang Liu  : liangzy@mail.ustc.edu.cn		++
//++	Fri Aug  7 11:31:27 CST 2020						++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#ifndef Anti_Neutron_Correction_Svc
#define Anti_Neutron_Correction_Svc
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>
#include "TGraph2D.h"
#include "TGraph.h"
#include <TRandom.h>
#include <TH2D.h>
#include <vector>
#include "GaudiKernel/IInterface.h"
#include "GaudiKernel/Kernel.h"
#include "GaudiKernel/Service.h"
#include "GaudiKernel/IService.h"
#include "GaudiKernel/IIncidentListener.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "AntiNeutronCorrectionSvc/IAntiNeutronCorrectionSvc.h"
#include "EmcRecEventModel/RecEmcShower.h"
#include "CLHEP/Vector/LorentzVector.h"
using CLHEP::HepLorentzVector;

#include "TString.h"

typedef std::vector<int> Vint;

using namespace std;

class AntiNeutronTrk;

class AntiNeutronCorrectionSvc : public Service, virtual public IAntiNeutronCorrectionSvc
{
		public:
				AntiNeutronCorrectionSvc(const std::string& name, ISvcLocator* svcLoc);
				virtual ~AntiNeutronCorrectionSvc();
				virtual StatusCode queryInterface(const InterfaceID& riid, void** ppvIF);
				virtual StatusCode initialize();
				virtual StatusCode finalize();


				AntiNeutronTrk *setAntiNeutronTrkR(const HepLorentzVector initp4);
				//++++++++++++++++++++++++++++++++++++++++++++++++++++++
				//		initp4: anti-neutron truth momentum
				//		initpost: anti-neutron truth initial position
				//		IP:	interaction point of this event
				//		statistical_uncertainty :  when the value of statistical_uncertainty is 0, using the centre value
				//		of the efficiency surface to estimate the anti-neutron selection efficiency; when this value 
				//		equals 1, efficiency = centre value + gaussian(0, sigma). sigma is the statistical uncertainty 
				//		of the efficiency from Jpsi-> p nbar pi control sample.
				//
				void *setAntiNeutronTrk(const HepLorentzVector initp4, const HepLorentzVector initposi, const HepLorentzVector IP,  int statistical_uncertainty = 0);
				void *setErrorMatrix(RecEmcShower *nbarTrk);
				RecEmcShower *getNbarShower();
				AntiNeutronTrk *setEnergy(double nbar_energy);
				AntiNeutronTrk *setJpsirunNo(const int runNo);
				bool isAntiNeutronCorrectionValid();
				int EMC_module();


		private:
				std::string m_eff0912_path;
				std::string m_error_matrix0912_path;
				std::string m_eff18_path;
				std::string m_error_matrix18_path;
				std::string m_eff_path;
				std::string m_error_matrix_path;
				int m_rdmSeed;
				double m_Nbar_Energy;
				double m_Nbar_SecMom;
				int m_Nbar_Hits;
				AntiNeutronTrk *m_nbarTrk;


};
#endif
