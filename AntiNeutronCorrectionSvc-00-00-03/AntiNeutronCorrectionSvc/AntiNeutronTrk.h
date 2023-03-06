#ifndef Anti_Neutron_Trk_H
#define Anti_Neutron_Trk_H
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>
#include "TGraph2D.h"
#include "TGraph.h"
#include <TString.h>
#include <TRandom.h>
#include <TH2D.h>
#include <vector>
#include "EvtRecEvent/EvtRecTrack.h"
#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "GaudiKernel/Service.h"
#include "EmcRecEventModel/RecEmcShower.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "AntiNeutronCorrectionSvc/AN_XYZ_ErrorMatrx.h"
using CLHEP::HepLorentzVector;

typedef HepGeom::Point3D<double> HepPoint3D;

typedef::vector<int> Vint;

using namespace std;

class AntiNeutronTrk{
		private:
				static AntiNeutronTrk *m_pointer;
				double init_px;
				double init_py;
				double init_pz;
				double init_e;

				TString theidx;
				TString p0validx;
				TString p1validx;
				TString p2validx;
				TString theidx1;
				TString pvalidx1;
				TString Evalidx;
				TString ntheidx;
				TString nphiidx;
				TString nangidx;
				TString nenergy_idx;
				TString ncosthe_idx;

				Hep3Vector nbarposi;
				Hep3Vector nbarposi3;  //  consider secondary vertex 

				TH2D *hnbar_ratio0912;
				TH2D *hnbar_ratio18;
				TH2D *hnbar_ratio;
				TH2D *hnbar_ratio0912E; // for the endcaps, re-divide into 300 bins in costheta
				TH2D *hnbar_ratio18E;
				TH2D *hnbar_ratioE;


				TH2D *hnbar_ratio0912_err;
				TH2D *hnbar_ratio18_err;
				TH2D *hnbar_ratio_err;
				TH2D *hnbar_ratio0912_errE;  // for the endcaps, re-divide into 300 bins in costheta
				TH2D *hnbar_ratio18_errE;
				TH2D *hnbar_ratio_errE;

				TH2D *hratio_multiTrk;
				TGraph2D *gnbar_ratio;

				int okreadfile;
				int m_EMC_module;
				int samplingFailure;
				int totalsampling;
				int okratio;
				int valid_count1;
				int valid_count2;
				int energy_flg;
				int efficiency_flg;
				int error_matrix_flg;
				struct nbar_xyz{
						double x;
						double y;
						double z;
						double theta;
						double phi;
						double dx;
						double dy;
						double dz;
						double dtheta;
						double dphi;
						double secmom;
						double energy;
						int hits;
				};
				nbar_xyz nbar_posi;	
				TFile *f1;
				TTree *t1;

				TFile *f2;
				TTree *t2;

				TFile *f3;
				TTree *t3;

				TFile *f_0912_1;
				TFile *f_0912_2;

				TFile *f_18_1;
				TFile *f_18_2;

				TRandom* ratiorandom;
				TRandom* rdm_nbardz;
				TRandom* rdm_nbardy;
				TRandom* rdm_nbardx;
				TRandom* rdm_nbarlatmom;
				TRandom* rdm_nbarz;
				TRandom* rdm_nbary;
				TRandom* rdm_nbarx;
				TRandom* rdm_nbar_err;


				TRandom* rdm_nbarhits;
				TRandom* rdm_nbarenergy;
				TRandom* rdm_nbarsecmom;
				TRandom* rdm_nbartheta;
				TRandom* rdm_nbarphi;
				TRandom* rdm_nbardtheta;
				TRandom* rdm_nbardphi;
				TRandom* rdm_nbarcosfrac;
				TRandom* rdm_nbarangle;
				HepSymMatrix errMatrix;
				AN_XYZ_ErrorMatrx *Dxyz;
				RecEmcShower *emcTrk;
				double Nbar_Energy;
				double Nbar_SecMom;
				int Nbar_Hits;

		public:

				static AntiNeutronTrk * instance();
				AntiNeutronTrk();
				virtual ~AntiNeutronTrk();

				int openFile(const char* path0912name, const char* angle0912path, const char* path18name, const char* angle18path, int rdmSeed);  // open efficiency files
				int openFile(const char* effpathname, const char* errpathname, int rdmSeed);  // open efficiency files
				void setJpsirunNo(const int runNo);
				void setEventCut(const double energy, const double secmom, const int hits);
				int setAntiNeutronTrkR(const HepLorentzVector initp4);
				int setAntiNeutronTrk(const HepLorentzVector initp4, const HepLorentzVector initposi, const HepLorentzVector IP, const int statistical_uncertainty = 0);
				void setErrorMatrix();
				void setErrorMatrix(RecEmcShower *nbarTrk);
				RecEmcShower* setNbarShower();
				void setEnergy(double nbar_energy);
				double energy();
				double secondMoment();
				int numHits();
				double theta();
				double phi();
				double dtheta();
				double dphi();
				double x();
				double y();
				double z();
				double dx();
				double dy();
				double dz();
				bool validp4();
				int samplingFailureN();
				int totalSampling();
				int EMC_module() {return m_EMC_module;};
};


#endif
