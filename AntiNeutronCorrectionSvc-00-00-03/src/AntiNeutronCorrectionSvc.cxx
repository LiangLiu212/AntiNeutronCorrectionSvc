
#ifndef Anti_Neutron_Correction_Svc_C
#define Anti_Neutron_Correction_Svc_C
#include <algorithm>

#include "GaudiKernel/Kernel.h"
#include "GaudiKernel/IInterface.h"
#include "GaudiKernel/StatusCode.h"
#include "GaudiKernel/SvcFactory.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/PropertyMgr.h"


#include "GaudiKernel/IIncidentSvc.h"
#include "GaudiKernel/Incident.h"
#include "GaudiKernel/IIncidentListener.h"

#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/Bootstrap.h"


#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TGraph.h"
#include "TH2D.h"
#include "TString.h"
#include "TGraph2D.h"
#include "TVector3.h"
#include "TRandom.h"

#include "CLHEP/Vector/LorentzVector.h"
#include "AntiNeutronCorrectionSvc/AntiNeutronCorrectionSvc.h"
#include "AntiNeutronCorrectionSvc/AntiNeutronTrk.h"
#endif
using CLHEP::HepLorentzVector;

//TString Itoa(int idx, int base);

using namespace std;
//double myfunc(double phi1, double theta1, double phi2, double theta2, double angle);


AntiNeutronCorrectionSvc::AntiNeutronCorrectionSvc(const std::string& name, ISvcLocator* svcLoc)
		: Service(name, svcLoc)
{
		declareProperty("eff0912_path", m_eff0912_path = "");
		declareProperty("error_matrix0912_path", m_error_matrix0912_path = "");
		declareProperty("eff18_path", m_eff18_path = "");
		declareProperty("error_matrix18_path", m_error_matrix18_path = "");
		declareProperty("eff_path", m_eff_path = "");
		declareProperty("error_matrix_path", m_error_matrix_path = "");
		declareProperty("Nbar_Energy", m_Nbar_Energy);
		declareProperty("Nbar_SecMom", m_Nbar_SecMom);
		declareProperty("Nbar_Hits", m_Nbar_Hits);
		declareProperty("rdmSeed", m_rdmSeed);
}

AntiNeutronCorrectionSvc::~AntiNeutronCorrectionSvc()
{
}

StatusCode AntiNeutronCorrectionSvc::queryInterface(const InterfaceID& riid, void** ppvIF)
{
  if (IAntiNeutronCorrectionSvc::interfaceID().versionMatch(riid)) {
    *ppvIF = dynamic_cast<IAntiNeutronCorrectionSvc*>(this);
  }
  else {
    return Service::queryInterface(riid, ppvIF);
  }
  return StatusCode::SUCCESS;
}



StatusCode AntiNeutronCorrectionSvc::initialize()
{
  MsgStream log(messageService(), name());
  log << MSG::INFO << "@initialize() " << endreq;
//  cout << "Hellow world from initialize() AntiNeutronCorrectionSvc " << m_eff_path << endl;

  StatusCode sc = Service::initialize();

  if (sc.isSuccess()) cout << "Hellow world from AntiNeutronCorrectionSvc!" << endl;
	const char* eff0912path = m_eff0912_path.c_str();
	const char* err0912path = m_error_matrix0912_path.c_str();
	const char* eff18path = m_eff18_path.c_str();
	const char* err18path = m_error_matrix18_path.c_str();
	const char* effpath = m_eff_path.c_str();
	const char* errpath = m_error_matrix_path.c_str();
	m_nbarTrk = AntiNeutronTrk::instance();
//    m_nbarTrk->openFile(eff0912path, err0912path, eff18path, err18path,  m_rdmSeed);
	m_nbarTrk->openFile(effpath, errpath, m_rdmSeed);
	m_nbarTrk->setEventCut(m_Nbar_Energy, m_Nbar_SecMom, m_Nbar_Hits);
   return sc;
}

StatusCode AntiNeutronCorrectionSvc::finalize()
{
  MsgStream log(messageService(), name());
  log << MSG::INFO << "@finalize()" << endreq;

  int samplingfailure = m_nbarTrk->samplingFailureN();
  int totalsampling = m_nbarTrk->totalSampling();


  log << MSG::INFO << "Total Sampling    :     " << totalsampling << endreq;
  if(samplingfailure != 0){
  		log << MSG::ERROR << "Sampling Failure  :     " << samplingfailure << endreq;
  		log << MSG::ERROR << "Failure  Ratio			   :     " << (double)samplingfailure/totalsampling << endreq;
  }

  if (m_nbarTrk) delete m_nbarTrk;

 // PartId2Name::release();

  StatusCode sc = Service::finalize();
  return sc;
}

void* AntiNeutronCorrectionSvc::setErrorMatrix(RecEmcShower *nbarTrk){
		m_nbarTrk->setErrorMatrix(nbarTrk);
	//	return m_nbarTrk;
}

RecEmcShower* AntiNeutronCorrectionSvc::getNbarShower(){
		return	m_nbarTrk->setNbarShower();
}

void* AntiNeutronCorrectionSvc::setAntiNeutronTrk(const HepLorentzVector initp4, const HepLorentzVector initposi, const HepLorentzVector IP,  int statistical_uncertainty ){
  		MsgStream log(messageService(), name());
		if(initp4.rho() == 0){
				log << MSG::ERROR << "In AntiNeutronTrk* setAntiNeutronTrk(const HepLorentzVector, const HepLorentzVector, const HepLorentzVector, int): wrong anti-Neutron truth four momentum" << endreq;
				abort();
		}
		if(IP.rho() == 0){
				log << MSG::ERROR << "In AntiNeutronTrk* setAntiNeutronTrk(const HepLorentzVector, const HepLorentzVector, const HepLorentzVector, int): wrong IP truth" << endreq;
				abort();
		}
		if(initposi.rho() == 0){
				log << MSG::ERROR << "In AntiNeutronTrk* setAntiNeutronTrk(const HepLorentzVector, const HepLorentzVector, const HepLorentzVector, int): wrong anti-Neutron initial position!" <<  endreq;
				abort();
		}
		m_nbarTrk->setAntiNeutronTrk(initp4, initposi, IP,  statistical_uncertainty);
	//	return m_nbarTrk;
}

AntiNeutronTrk* AntiNeutronCorrectionSvc::setEnergy(double nbar_energy){
		m_nbarTrk->setEnergy(nbar_energy);
		return m_nbarTrk;
}

AntiNeutronTrk* AntiNeutronCorrectionSvc::setJpsirunNo(const int runNo){
		m_nbarTrk->setJpsirunNo(runNo);
		return m_nbarTrk;
}

bool AntiNeutronCorrectionSvc::isAntiNeutronCorrectionValid(){
		return m_nbarTrk->validp4();

}

int AntiNeutronCorrectionSvc::EMC_module(){
		return m_nbarTrk->EMC_module();
}

