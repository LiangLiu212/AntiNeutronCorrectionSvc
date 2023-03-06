#ifndef Anti_Neutron_Trk
#define Anti_Neutron_Trk
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
#include "GaudiKernel/Service.h"
#include "EmcRecEventModel/RecEmcShower.h"
#include "CLHEP/Vector/LorentzVector.h"
using CLHEP::HepLorentzVector;

#include  "AntiNeutronCorrectionSvc/AntiNeutronTrk.h"
#endif
TString Itoa(int idx, int base);
double myfunc(double phi1, double theta1, double phi2, double theta2, double angle);

AntiNeutronTrk * AntiNeutronTrk::m_pointer = 0 ;

AntiNeutronTrk * AntiNeutronTrk::instance() {
		if(m_pointer) return m_pointer;
		m_pointer = new AntiNeutronTrk();
		return m_pointer;
}

AntiNeutronTrk::AntiNeutronTrk(){;}

AntiNeutronTrk::~AntiNeutronTrk() {
		delete ratiorandom;
		delete rdm_nbarlatmom;
		delete rdm_nbarhits;
		delete rdm_nbarenergy;
		delete rdm_nbarsecmom;
		delete rdm_nbartheta;
		delete rdm_nbarphi;
		delete rdm_nbardtheta;
		delete rdm_nbardphi;
		delete rdm_nbar_err;
		delete rdm_nbarcosfrac;
		delete rdm_nbarangle;
		delete hnbar_ratio_err;
		delete hnbar_ratio;

}



//+++++++++++++++++++++ initialize the AntineutronCorrection package +++++++++++++++++++++
// Open the root file including effciency surface, CDF, and error matrix
// The parameters of the function are the path of the root file 
// AntineutronCorrection::AntineutronCorrection(PATH of efficiency surface file, PATH of error matrix
// 													file, RANDOM seed)
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


int AntiNeutronTrk::openFile(const char* effpathname, const char* errpathname,  int rdmSeed)
{

		cout << "random seed for anti-neutron selection : " << rdmSeed << endl;

		//++++++++++++++ open root file +++++++++++++++++++++

		f1 = new TFile(effpathname, "read"); 		// open efficiency file
		if(NULL == f1){
				cout << "Cannot open the efficiency file!" << endl;
				abort();
		}
		f2 = new TFile(errpathname, "read");		// open error matrix file
		if(NULL == f2){
				cout << "can not open the error matrix file!" << endl;
				abort();
		}

		// read the efficiency surface histogram

		TString rationame = "hratio";
		hnbar_ratio = (TH2D*)f1->Get(rationame);
		if(NULL == hnbar_ratio) {
				cout << "can not read the 2D efficiency surface histogram!" << endl;
				abort();
		}
		TString rationameE = "hratioEndCap";
		hnbar_ratioE = (TH2D*)f1->Get(rationameE);
		if(NULL == hnbar_ratioE){
				cout << "can not read the 2D (EndCap) efficiency surface histogram!" << endl;
				abort();
		}

		rdm_nbar_err = new TRandom();
		rdm_nbar_err->SetSeed(rdmSeed + 11);
		hnbar_ratio_err = new TH2D("hnbar_ratio_err", "hnbar_ratio_err", 50, -1, 1, 50, -1, 1); 	// for anti-neutron correction 


		for(int i =1; i <= 50; i++){
				for(int j = 1; j <=50; j++){
						Double_t bin_err = hnbar_ratio->GetBinError(i, j);
						hnbar_ratio_err->SetBinContent(i, j, rdm_nbar_err->Gaus(0, bin_err));
				}
		}
		
	
		hnbar_ratio_errE = new TH2D("hnbar_ratio_errE", "hnbar_ratio_errE", 50, -1, 1, 300, -1, 1); 	// for anti-neutron correction 
		for(int i =1; i <= 50; i++){
				for(int j = 1; j <=300; j++){
						Double_t bin_errE = hnbar_ratioE->GetBinError(i, j);
						hnbar_ratio_errE->SetBinContent(i, j, rdm_nbar_err->Gaus(0, bin_errE));
				}
		}
		
		valid_count1=0;
		valid_count2=0;
		samplingFailure = 0;
		totalsampling = 0;
		okreadfile = 1;
		energy_flg = 0;
		efficiency_flg = 0;
		error_matrix_flg = 0;

		ratiorandom = new TRandom();
		ratiorandom->SetSeed(rdmSeed);
		rdm_nbardtheta = new TRandom();
		rdm_nbardtheta->SetSeed(rdmSeed+1);
		rdm_nbardphi = new TRandom();
		rdm_nbardphi->SetSeed(rdmSeed+2);
		rdm_nbartheta = new TRandom();
		rdm_nbartheta->SetSeed(rdmSeed+3);
		rdm_nbarphi = new TRandom();
		rdm_nbarphi->SetSeed(rdmSeed+4);
		rdm_nbarenergy = new TRandom();
		rdm_nbarenergy->SetSeed(rdmSeed+5);
		rdm_nbarhits = new TRandom();
		rdm_nbarhits->SetSeed(rdmSeed+6);
		rdm_nbarsecmom = new TRandom();
		rdm_nbarsecmom->SetSeed(rdmSeed+7);
		rdm_nbarlatmom = new TRandom();
		rdm_nbarlatmom->SetSeed(rdmSeed+8);
		rdm_nbarcosfrac = new TRandom();
		rdm_nbarcosfrac->SetSeed(rdmSeed+9);
		rdm_nbarangle = new TRandom();
		rdm_nbarangle->SetSeed(rdmSeed+10);

		emcTrk = new RecEmcShower();
		Dxyz = new AN_XYZ_ErrorMatrx();
		return okreadfile;
}


int AntiNeutronTrk::openFile(const char* path0912name, const char* angle0912path, const char* path18name, const char* angle18path,   int rdmSeed)
{

		cout << "random seed for anti-neutron selection : " << rdmSeed << endl;

		//++++++++++++++ open root file +++++++++++++++++++++

		f_0912_1 = new TFile(path0912name, "read"); 		// open efficiency file
		if(NULL == f_0912_1){
				cout << "Cannot open the efficiency file!" << endl;
				abort();
		}
		f_0912_2 = new TFile(angle0912path, "read");		// open error matrix file
		if(NULL == f_0912_2){
				cout << "can not open the error matrix file!" << endl;
				abort();
		}

		// read the efficiency surface histogram

		TString rationame = "hratio";
		hnbar_ratio0912 = (TH2D*)f_0912_1->Get(rationame);
		if(NULL == hnbar_ratio0912) {
				cout << "can not read the 2D efficiency surface histogram!" << endl;
				abort();
		}
		TString rationameE = "hratioEndCap";
		hnbar_ratio0912E = (TH2D*)f_0912_1->Get(rationameE);
		if(NULL == hnbar_ratio0912E){
				cout << "can not read the 2D (EndCap) efficiency surface histogram!" << endl;
				abort();
		}






		f_18_1 = new TFile(path18name, "read"); 		// open efficiency file
		if(NULL == f_18_1){
				cout << "Cannot open the efficiency file!" << endl;
				abort();
		}
		f_18_2 = new TFile(angle18path, "read");		// open error matrix file
		if(NULL == f_18_2){
				cout << "can not open the error matrix file!" << endl;
				abort();
		}

		// read the efficiency surface histogram

		hnbar_ratio18 = (TH2D*)f_18_1->Get(rationame);
		if(NULL == hnbar_ratio18) {
				cout << "can not read the 2D efficiency surface histogram!" << endl;
				abort();
		}


		hnbar_ratio18E = (TH2D*)f_18_1->Get(rationameE);
		if(NULL == hnbar_ratio18E) {
				cout << "can not read the 2D (EndCap) efficiency surface histogram!" << endl;
				abort();
		}



		rdm_nbar_err = new TRandom();
		rdm_nbar_err->SetSeed(rdmSeed + 11);

		hnbar_ratio18_err = new TH2D("hnbar_ratio18_err", "hnbar_ratio18_err", 50, -1, 1, 50, -1, 1); 	// for anti-neutron correction 
		for(int i =1; i <= 50; i++){
				for(int j = 1; j <=50; j++){
						Double_t bin_err = hnbar_ratio18->GetBinError(i, j);
						hnbar_ratio18_err->SetBinContent(i, j, rdm_nbar_err->Gaus(0, bin_err));
				}
		}
	
		hnbar_ratio0912_err = new TH2D("hnbar_ratio0912_err", "hnbar_ratio0912_err", 50, -1, 1, 50, -1, 1); 	// for anti-neutron correction 
		for(int i =1; i <= 50; i++){
				for(int j = 1; j <=50; j++){
						Double_t bin_err = hnbar_ratio0912->GetBinError(i, j);
						hnbar_ratio0912_err->SetBinContent(i, j, rdm_nbar_err->Gaus(0, bin_err));
				}
		}
		
	
		hnbar_ratio18_errE = new TH2D("hnbar_ratio18_errE", "hnbar_ratio18_errE", 50, -1, 1, 50, -1, 1); 	// for anti-neutron correction 
		for(int i =1; i <= 50; i++){
				for(int j = 1; j <=50; j++){
						Double_t bin_errE = hnbar_ratio18E->GetBinError(i, j);
						hnbar_ratio18_errE->SetBinContent(i, j, rdm_nbar_err->Gaus(0, bin_errE));
				}
		}
	
		hnbar_ratio0912_errE = new TH2D("hnbar_ratio0912_errE", "hnbar_ratio0912_errE", 50, -1, 1, 50, -1, 1); 	// for anti-neutron correction 
		for(int i =1; i <= 50; i++){
				for(int j = 1; j <=50; j++){
						Double_t bin_errE = hnbar_ratio0912E->GetBinError(i, j);
						hnbar_ratio0912_errE->SetBinContent(i, j, rdm_nbar_err->Gaus(0, bin_errE));
				}
		}
		




		valid_count1=0;
		valid_count2=0;
		samplingFailure = 0;
		totalsampling = 0;
		okreadfile = 1;
		energy_flg = 0;
		efficiency_flg = 0;
		error_matrix_flg = 0;

		ratiorandom = new TRandom();
		ratiorandom->SetSeed(rdmSeed);
		rdm_nbardtheta = new TRandom();
		rdm_nbardtheta->SetSeed(rdmSeed+1);
		rdm_nbardphi = new TRandom();
		rdm_nbardphi->SetSeed(rdmSeed+2);
		rdm_nbartheta = new TRandom();
		rdm_nbartheta->SetSeed(rdmSeed+3);
		rdm_nbarphi = new TRandom();
		rdm_nbarphi->SetSeed(rdmSeed+4);
		rdm_nbarenergy = new TRandom();
		rdm_nbarenergy->SetSeed(rdmSeed+5);
		rdm_nbarhits = new TRandom();
		rdm_nbarhits->SetSeed(rdmSeed+6);
		rdm_nbarsecmom = new TRandom();
		rdm_nbarsecmom->SetSeed(rdmSeed+7);
		rdm_nbarlatmom = new TRandom();
		rdm_nbarlatmom->SetSeed(rdmSeed+8);
		rdm_nbarcosfrac = new TRandom();
		rdm_nbarcosfrac->SetSeed(rdmSeed+9);
		rdm_nbarangle = new TRandom();
		rdm_nbarangle->SetSeed(rdmSeed+10);
		return okreadfile;
}


void AntiNeutronTrk::setJpsirunNo(const int runNo){
		energy_flg = 2;
	 	if(runNo > 0){
				efficiency_flg = -1;
		}

}
//+++++++++++++++++++++++++++++ set anti-neutron track ++++++++++++++++
// copy the anti-neutron shower parameters frome 
// RecEmcShower to AntiNeutronTrk
/*
void AntiNeutronTrk::setErrorMatrix(RecEmcShower *nbarTrk){
		nbar_posi.x = nbarTrk->x();
		nbar_posi.y = nbarTrk->y();
		nbar_posi.z = nbarTrk->z();
		nbar_posi.dx = nbarTrk->dx();
		nbar_posi.dy = nbarTrk->dy();
		nbar_posi.dz = nbarTrk->dz();
		nbar_posi.theta = nbarTrk->theta();
		nbar_posi.phi = nbarTrk->phi();
		nbar_posi.dtheta = nbarTrk->dtheta();
		nbar_posi.dphi = nbarTrk->dphi();
		nbar_posi.energy = nbarTrk->energy();
		nbar_posi.secmom = nbarTrk->secondMoment();
		nbar_posi.hits = nbarTrk->numHits();
		okratio = 2;
}
*/


RecEmcShower* AntiNeutronTrk::setNbarShower(){
		return emcTrk;
	/*	nbarTrk->setEnergy(nbar_posi.energy);
		nbarTrk->setSecondMoment(nbar_posi.secmom);
		nbarTrk->setNumHits(nbar_posi.hits);
		HepSymMatrix matrix(3);
		int err3[3];
		err3[0] = nbar_posi.dx;
		err3[1] = nbar_posi.dy;
		err3[2] = nbar_posi.dz;
		for(int i = 0; i < 3; i++){
				for(int j = 0; j < 3; j++){
						matrix[i][j]=err3[i]*err3[j];
				}
		}
		nbarTrk->setErrorMatrix(matrix);
		nbarTrk->setDtheta(nbar_posi.dtheta);
		nbarTrk->setDphi(nbar_posi.dphi);
		nbarTrk->setPosition(HepPoint3D(nbar_posi.x, nbar_posi.y, nbar_posi.z));
		nbarTrk->setTime(999); */
}

void AntiNeutronTrk::setErrorMatrix(RecEmcShower *nbarTrk){    // nbarTrk is the anti-neutron track from Emc shower, 
		if(energy_flg ==  -1){
				cout << "Cannot call AntiNeutronTrk::setErrorMatrix after AntiNeutronTrk::setEnergy" << endl;
				abort();
		}
	//	cout << " Assignment of the anti-neutron error matrix for data. " << endl;
//+++++++++++++++++++++++++++++ sampling a new error matrix for Emc shower. ++++++++++++++++++++++++++++++
		double costhe = cos(nbarTrk->theta());
		double nbar_Energy = nbarTrk->energy();
		int costhe_idx = floor((costhe+1)/0.04) + 1;    // index of cos theta, divide cos theta to 50 intervals.
		int energy_idx = floor((nbar_Energy)/0.1) + 1;	// index of energy, divide energy to 20 intervals.
		if(costhe_idx < 1) costhe_idx = 1;
		if(costhe_idx > 50) costhe_idx = 50;
		if(energy_idx  > 20) energy_idx = 20;
		nenergy_idx = "_e" + Itoa(energy_idx, 10);
		ncosthe_idx = "_t" + Itoa(costhe_idx, 10);
		TString emcnbar_the = "nbar_deltatheta" + nenergy_idx + ncosthe_idx;
		TString emcnbar_phi = "nbar_deltaphi" + nenergy_idx +  ncosthe_idx;
		TGraph *gnbar_delthe = (TGraph*)f2->Get(emcnbar_the);
		if(NULL == gnbar_delthe) {
				cout << "can not read the error matrix delta theta C.D.F file!" << endl;
		}
		double rtn_nbardelthe = gnbar_delthe->Eval(rdm_nbardtheta->Rndm());
		TGraph *gnbar_delphi = (TGraph*)f2->Get(emcnbar_phi);
		if(NULL == gnbar_delphi) {
				cout << "can not read the error matrix delta phi C.D.F file!" << endl;
		}
		double rtn_nbardelphi = gnbar_delphi->Eval(rdm_nbardphi->Rndm());
		
		Dxyz->DtheDphi_to_DxDyDz(nbarTrk->theta(), nbarTrk->phi(), rtn_nbardelthe, rtn_nbardelphi);
		nbarTrk->setErrorMatrix(Dxyz->getDxDyDz());
		nbarTrk->setDtheta(rtn_nbardelthe);
		nbarTrk->setDphi(rtn_nbardelphi);

//		nbar_posi.dtheta = rtn_nbardelthe;   			// assignment of delta theta
//		nbar_posi.dphi = rtn_nbardelphi;				// assignment of delta phi

		delete gnbar_delthe;
		delete gnbar_delphi;
		okratio = 2;
}
void AntiNeutronTrk::setErrorMatrix(){
		if(energy_flg ==  -1){
				cout << "Cannot call AntiNeutronTrk::setErrorMatrix after AntiNeutronTrk::setEnergy" << endl;
				abort();
		}
	//	cout << " Assignment of the anti-neutron error matrix for data. " << endl;

		emcTrk->setPosition(HepPoint3D(nbar_posi.x, nbar_posi.y, nbar_posi.z));
		emcTrk->setEnergy(nbar_posi.energy);
		emcTrk->setNumHits(nbar_posi.hits);
		emcTrk->setSecondMoment(nbar_posi.secmom);
//+++++++++++++++++++++++++++++ sampling a new error matrix for Antineutron Correction. ++++++++++++++++++++++++++++++
		double costhe = cos(nbar_posi.theta);
		double nbar_Energy = nbar_posi.energy;
		if(nbar_Energy < Nbar_Energy) nbar_Energy = Nbar_Energy;
		int costhe_idx = floor((costhe+1)/0.04) + 1;
		int energy_idx = floor((nbar_Energy)/0.1) + 1;
		if(costhe_idx < 1) costhe_idx = 1;
		if(costhe_idx > 50) costhe_idx = 50;
		if(energy_idx  > 20) energy_idx = 20;
		nenergy_idx = "_e" + Itoa(energy_idx, 10);
		ncosthe_idx = "_t" + Itoa(costhe_idx, 10);
		TString emcnbar_the = "nbar_deltatheta" + nenergy_idx + ncosthe_idx;
		TString emcnbar_phi = "nbar_deltaphi" + nenergy_idx +  ncosthe_idx;
		TGraph *gnbar_delthe = (TGraph*)f2->Get(emcnbar_the);
		if(NULL == gnbar_delthe) {
				cout << "can not read the error matrix delta theta C.D.F file!" << endl;
				abort();
		}
		double rtn_nbardelthe = gnbar_delthe->Eval(rdm_nbardtheta->Rndm());
		TGraph *gnbar_delphi = (TGraph*)f2->Get(emcnbar_phi);
		if(NULL == gnbar_delphi) {
				cout << "can not read the error matrix delta phi C.D.F file!" << endl;
				abort();
		}
		double rtn_nbardelphi = gnbar_delphi->Eval(rdm_nbardphi->Rndm());
		Dxyz->DtheDphi_to_DxDyDz(nbar_posi.theta, nbar_posi.phi, rtn_nbardelthe, rtn_nbardelphi);
		emcTrk->setErrorMatrix(Dxyz->getDxDyDz());
		emcTrk->setDtheta(rtn_nbardelthe);
		emcTrk->setDphi(rtn_nbardelphi);

		delete gnbar_delthe;
		delete gnbar_delphi;
		okratio = 2;
}

int AntiNeutronTrk::setAntiNeutronTrkR(const HepLorentzVector initp4){     //  using replace method to test
		cout << "Using data replace the Monte Carlo Simulation." << endl;
		if(okreadfile != 1){
				cout << "please set a correct root file in AntiNeutronTrk()" << endl;
				cout << "please set a correct root file in AntiNeutronTrk()!" << endl;
				abort();
		}
		init_px = initp4.px();
		init_py = initp4.py();
		init_pz = initp4.pz();
		init_e = initp4.e();
		Double_t  nbar_px, nbar_py, nbar_pz, nbar_e;

		Double_t nbar_energy;
		Int_t nbar_hits;
		Double_t nbar_secmom;
		Double_t nbar_theta;
		Double_t nbar_phi;
		Double_t nbar_dtheta;
		Double_t nbar_dphi;
		Double_t nbar_x;
		Double_t nbar_y;
		Double_t nbar_z;
		Double_t nbar_dx;
		Double_t nbar_dy;
		Double_t nbar_dz;
		Double_t nbar_latmom;

		t1->SetBranchAddress("vtxnbar_px", &nbar_px);
		t1->SetBranchAddress("vtxnbar_py", &nbar_py);
		t1->SetBranchAddress("vtxnbar_pz", &nbar_pz);
		t1->SetBranchAddress("vtxnbar_e", &nbar_e);
		t1->SetBranchAddress("nbar_energy", &nbar_energy);
		t1->SetBranchAddress("nbar_hits", &nbar_hits);
		t1->SetBranchAddress("nbar_secmom", &nbar_secmom);
		t1->SetBranchAddress("nbar_theta", &nbar_theta);
		t1->SetBranchAddress("nbar_phi", &nbar_phi);
		t1->SetBranchAddress("nbar_dtheta", &nbar_dtheta);
		t1->SetBranchAddress("nbar_dphi", &nbar_dphi);
		t1->SetBranchAddress("nbar_x", &nbar_x);
		t1->SetBranchAddress("nbar_y", &nbar_y);
		t1->SetBranchAddress("nbar_z", &nbar_z);
		t1->SetBranchAddress("nbar_dx", &nbar_dx);
		t1->SetBranchAddress("nbar_dy", &nbar_dy);
		t1->SetBranchAddress("nbar_dz", &nbar_dz);
		t1->SetBranchAddress("nbar_latmom", &nbar_latmom);

		HepLorentzVector nbar_p4;

		int match_count = 0;

		int startevent = ratiorandom->Integer(1000000);

		cout << startevent << endl;
	
		for(int i = startevent; i < t1->GetEntries(); i++){
				t1->GetEntry(i);
				double delta_px = init_px - nbar_px;
				double delta_py = init_py - nbar_py;
				double delta_pz = init_pz - nbar_pz;
				double delta_E = fabs(init_e - nbar_e);
				match_count++;
				if( sqrt(delta_px*delta_px + delta_py*delta_py + delta_pz*delta_pz + delta_E*delta_E) < 0.02 ){
						nbar_posi.x 	 = nbar_x;
						nbar_posi.y 	 = nbar_y;
						nbar_posi.z 	 = nbar_z;
						nbar_posi.dx 	 = nbar_dx;
						nbar_posi.dy 	 = nbar_dy;
						nbar_posi.dz 	 = nbar_dz;
						nbar_posi.theta  = nbar_theta;
						nbar_posi.phi 	 = nbar_phi;
						nbar_posi.dtheta = nbar_dtheta;
						nbar_posi.dphi 	 = nbar_dphi;
						nbar_posi.energy = nbar_energy;
						nbar_posi.secmom = nbar_secmom;
						nbar_posi.hits 	 = nbar_hits;
						cout <<" success : " <<  match_count << endl;
						okratio = 2;
						return okratio;

				}

		}

		okratio = 0;
		return okratio;

}



int AntiNeutronTrk::setAntiNeutronTrk(const HepLorentzVector initp4, const HepLorentzVector initposi, const HepLorentzVector IP, const int statistical_uncertainty){
		if(okreadfile != 1){
				cout << "please set a correct root file in AntiNeutronTrk()" << endl;
				cout << "please set a correct root file in AntiNeutronTrk()!" << endl;
				abort();
		}

		if(statistical_uncertainty != 0 && statistical_uncertainty != 1){
				cout << "please set a correct root file in AntiNeutronTrk()" << endl;
				cout << "The value of statistical_uncertainty must be 0 or 1!" << endl;
				abort();
		}

		if(efficiency_flg == -1){
				cout << "AntiNeutronTrk::setAntiNeutronTrk(const int runNo), ERROR! It is design for MC!" << endl;
				abort();
		}
		

		init_px = initp4.px();
		init_py = initp4.py();
		init_pz = initp4.pz();
		init_e = initp4.e();

		// consider the second vertex
		//
		//
		// initial position of anti-neutron
		double x1 = initposi.x();   		
		double y1 = initposi.y();
		double z1 = initposi.z();
		

		double tan_phi = tan(initp4.phi());
		double tan_theta = tan(initp4.theta());

		double x2;
		double y2;
		double z2;

		double tmp_x2;
		double tmp_y2;
		double tmp_z2;

		double d = 95.5; 		// the radius of the EMC, 95.5 cm


		double tmp_x21 = (-x1 - y1*tan_phi - sqrt(-y1*y1 + d*d + 2*x1*y1*tan_phi - x1*x1*tan_phi*tan_phi + d*d*tan_phi*tan_phi))/(1+tan_phi*tan_phi);
		double tmp_x22 = (-x1 - y1*tan_phi + sqrt(-y1*y1 + d*d + 2*x1*y1*tan_phi - x1*x1*tan_phi*tan_phi + d*d*tan_phi*tan_phi))/(1+tan_phi*tan_phi);
		if(init_px*tmp_x21 > 0){
				tmp_x2 = tmp_x21;
				tmp_y2 = tmp_x21*tan_phi;
		}
		else {
				tmp_x2 = tmp_x22;
				tmp_y2 = tmp_x22*tan_phi;
		}

		tmp_z2 = sqrt(tmp_x2*tmp_x2 + tmp_y2*tmp_y2)/tan_theta;
		if(tmp_z2 + z1 > 140.0 ){     			// anti-neutron hits on the Endcap.
				tmp_z2 = 140.0 - z1;
				z2 = tmp_z2;
				double l1 = fabs(tmp_z2) * fabs(tan_theta);
				x2 = tmp_x2*l1/(sqrt(tmp_x2*tmp_x2 + tmp_y2*tmp_y2));
				y2 = tmp_y2*l1/(sqrt(tmp_x2*tmp_x2 + tmp_y2*tmp_y2));
		}
		else if(tmp_z2 + z1 < -140.0){			// anti-neutron hits on the Endcap.
				tmp_z2 = -140.0 - z1;
				z2 = tmp_z2;
				double l1 = fabs(tmp_z2) * fabs(tan_theta);
				x2 = tmp_x2*l1/(sqrt(tmp_x2*tmp_x2 + tmp_y2*tmp_y2));
				y2 = tmp_y2*l1/(sqrt(tmp_x2*tmp_x2 + tmp_y2*tmp_y2));
		}
		else {									// anti-neutron hits on the barrel.
				z2 = tmp_z2;
				x2 = tmp_x2;
				y2 = tmp_y2;
		}

		double x3 = x1+x2 - IP.x();
		double y3 = y1+y2 - IP.y();
		double z3 = z1+z2 - IP.z();
		nbarposi3.setX(x3);
		nbarposi3.setY(y3);
		nbarposi3.setZ(z3);


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		double pmom_val = sqrt(init_px*init_px + init_py*init_py +init_pz*init_pz);   	// true momentum of the anti-neutron
		double costhe = nbarposi3.cosTheta();											// cos theta from the hit point on emc to the IP.
		double nthe = nbarposi3.theta();												// theta from the hit point on emc to the IP.
		double nphi = nbarposi3.phi();													// phi from the hit point on emc to the IP.
		int costhe_ratio_idx = floor((costhe+1)/0.04) + 1;								// index of the cos theta, divide the cos theta into 50 intervals
		int costhe_ratio_idx_E = floor((costhe+1)/(2.0/300)) + 1;								// index of the cos theta, divide the cos theta into 50 intervals
		int pmom_ratio_idx = floor(pmom_val/0.024) + 1;									// index of the momenutm, divide the momentum into 50 intervals

		Double_t ratio = 0;
		TGraph *eff_interpolation;
		if(statistical_uncertainty == 0){						
				ratio =  hnbar_ratio->GetBinContent(pmom_ratio_idx, costhe_ratio_idx);   						
				if(fabs(cos(nthe)) > 0.72){  		// a naive way to weight the efficiency when the cos theta out of the detection acceptance................
						ratio = hnbar_ratioE->GetBinContent(pmom_ratio_idx, costhe_ratio_idx_E);
				}
		}
		else if(statistical_uncertainty == 1){				// consider the statistical uncertainty of the Jpsi -> pnbarpi sample
				ratio =  hnbar_ratio->GetBinContent(pmom_ratio_idx, costhe_ratio_idx) + hnbar_ratio_err->GetBinContent(pmom_ratio_idx, costhe_ratio_idx);
				if(fabs(cos(nthe)) > 0.72){
						ratio = hnbar_ratioE->GetBinContent(pmom_ratio_idx, costhe_ratio_idx_E) + hnbar_ratio_errE->GetBinContent(pmom_ratio_idx, costhe_ratio_idx_E);
				}
		}
		double ratran = ratiorandom->Rndm();
		okratio = 0;
		if( ratran <= ratio) {
				totalsampling++;

				theidx1 = "_t" + Itoa(costhe_ratio_idx, 10);									
				pvalidx1 = "_p" + Itoa(pmom_ratio_idx, 10);


				TString nbar_ename = "nbar_energy" + pvalidx1 + theidx1;
				TGraph *gnbar_energy = (TGraph*)f1->Get(nbar_ename); 							// anti-neutron energy sampling
				if(NULL == gnbar_energy) {
						cout << "can not read the energy C.D.F TGraph!" << endl;
						abort();
				}
				double rtn_nbarenergy = gnbar_energy->Eval(rdm_nbarenergy->Rndm(), 0);

				TString nbar_hitsname = "nbar_hits" + pvalidx1 + theidx1;
				TGraph *gnbar_hits = (TGraph*)f1->Get(nbar_hitsname);
				if(NULL == gnbar_hits) {
						cout << "can not read the energy C.D.F TGraph!" << endl;
						abort();
				}
				double rtn_nbarhits = ceil(gnbar_hits->Eval(rdm_nbarhits->Rndm(), 0));

				TString nbar_secmomname = "nbar_secmom" + pvalidx1 + theidx1;
				TGraph *gnbar_secmom = (TGraph*)f1->Get(nbar_secmomname);
				if(NULL == gnbar_secmom) {
						cout << "can not read the second moment C.D.F TGraph!" << endl;
						abort();
				}
				double rtn_nbarsecmom = gnbar_secmom->Eval(rdm_nbarsecmom->Rndm(), 0);


				nbar_posi.energy = rtn_nbarenergy;
				nbar_posi.hits = rtn_nbarhits;
				nbar_posi.secmom = rtn_nbarsecmom;


				int ct_idx = floor((costhe+1)/0.1) + 1;
				int p_idx = floor(pmom_val/0.1) + 1;
				int E_idx2 = floor((rtn_nbarenergy - Nbar_Energy)/0.1) + 1;
				int E_idx;

				//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
				// for the anti-neutron position smear
				// we using energy, cos theta and momentum 3D function 
				//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
				
				if(E_idx2 < 1) 					E_idx = 1;
				if(E_idx2 >=1 && E_idx2 <=6 ) 	E_idx = E_idx2;
				if(E_idx2 == 7 || E_idx2 == 8) 	E_idx = 7;
				if(E_idx2 == 9 || E_idx2 == 10) E_idx = 8;
				if(E_idx2 > 10) 				E_idx = 9;
				if(ct_idx > 20) ct_idx = 20;
				if(p_idx > 12) p_idx = 12;



				theidx = "_t" + Itoa(ct_idx, 10);
				p1validx = "_p" + Itoa(p_idx, 10);
				Evalidx = "_E" + Itoa(E_idx, 10);

				TString nbar_angle0name = "nbar_angle0" + p1validx + theidx + Evalidx;				// angle C.D.F of barrel
				TGraph *gnbar_angle0 = (TGraph*)f1->Get(nbar_angle0name);
				if(NULL == gnbar_angle0) {
						cout << "cannot open file : angle0" << endl;
						abort();
				}

				TString nbar_theta0name = "nbar_theta0" + p1validx + theidx + Evalidx;				// theta C.D.F of barrel
				TGraph *gnbar_theta0 = (TGraph*)f1->Get(nbar_theta0name);
				if(NULL == gnbar_theta0) {
						cout << "cannot open file : theta0" << endl;
						abort();
				}

				TString nbar_angle1name = "nbar_angle1" + p1validx + theidx + Evalidx;				// angle C.D.F of Endcaps.
				TGraph *gnbar_angle1 = (TGraph*)f1->Get(nbar_angle1name);
				if(NULL == gnbar_angle1) {
						cout << "cannot open file : angle1" << endl;
						abort();
				}

				TString nbar_theta1name = "nbar_theta1" + p1validx + theidx + Evalidx;				// theta C.D.F of Endcaps.
				TGraph *gnbar_theta1 = (TGraph*)f1->Get(nbar_theta1name);
				if(NULL == gnbar_theta1) {
						cout << "cannot open file : theta1" << endl;
						abort();
				}



				TString nbar_thetaratio = "ratio_thetaP" + Itoa(pmom_ratio_idx, 10);
				TGraph *gnbar_thetaratio = (TGraph*)f1->Get(nbar_thetaratio);
				if(NULL == gnbar_thetaratio) {
						cout << "cannot open file : thetaP" << endl;
						abort();
				}

				double rtn_nbarangle;
				double theta_phi;
				double rdm_theta;
				double rat_theta;
				double rtn_nbardtheta;
				double rtn_nbartheta;
				double rtn_nbardphi;
				double rtn_nbarphi;

				double theta_ratio =  gnbar_thetaratio->Eval(cos(nthe));
				double rdm_theta_ratio = rdm_nbarcosfrac->Rndm();
				double rdm_angle = rdm_nbarangle->Rndm();
				valid_count1 = 0;
				valid_count2 = 0;

				if(theta_ratio > rdm_theta_ratio){
						m_EMC_module = 0;
						rtn_nbarangle = gnbar_angle0->Eval(rdm_angle);
						do{
								if(valid_count1 == 1000) break;
								valid_count1++;
								theta_phi = rdm_nbartheta->Rndm()*2*CLHEP::pi;
								rtn_nbardtheta = rtn_nbarangle*sin(theta_phi);							// delta theta
								rtn_nbardphi = rtn_nbarangle*cos(theta_phi)/sin(nthe);					// delta phi
								rtn_nbartheta = nthe + rtn_nbardtheta;									// theta
								rtn_nbarphi = nphi + rtn_nbardphi;										// phi

								rdm_theta = rdm_nbarphi->Rndm();
								rat_theta = gnbar_theta0->Eval(rtn_nbartheta);							// 
								if(valid_count1 > 500){
										rdm_theta = 0;
								}
						}while( fabs(cos(rtn_nbartheta)) > 0.8 || rtn_nbartheta < 0 || rtn_nbartheta > CLHEP::pi || rdm_theta > rat_theta );
						if(valid_count1 == 1000){   // considering the detection acceptance in a very naive way................
								cout << "barrel : " << nthe << " " << rtn_nbardtheta << endl;
								if(cos(rtn_nbartheta) > 0.8){
										rtn_nbartheta = nthe + 1.5*fabs(rtn_nbardtheta);
								}
								else if(cos(rtn_nbartheta) < -0.8){
										rtn_nbartheta = nthe - 1.5*fabs(rtn_nbardtheta);
								}

						}
				}
				else{
						m_EMC_module = 1;
						rtn_nbarangle = gnbar_angle1->Eval(rdm_angle);
						do{
								if(valid_count2 == 1000) break;
								valid_count2++;
								theta_phi = rdm_nbartheta->Rndm()*2*CLHEP::pi;
								rtn_nbardtheta = rtn_nbarangle*sin(theta_phi);
								rtn_nbardphi = rtn_nbarangle*cos(theta_phi)/sin(nthe);
								rtn_nbartheta = nthe + rtn_nbardtheta;
								rtn_nbarphi = nphi + rtn_nbardphi;
								rdm_theta = rdm_nbarphi->Rndm();
								rat_theta = gnbar_theta1->Eval(rtn_nbartheta);
								if(valid_count2 > 500){
										rdm_theta = 0;
								}
						}while(fabs(cos(rtn_nbartheta)) > 0.92 ||  fabs(cos(rtn_nbartheta)) < 0.86  || rtn_nbartheta < 0 || rtn_nbartheta > CLHEP::pi || rdm_theta > rat_theta );
						if(valid_count2 == 1000){    // considering the detection acceptance in a very naive way................
								cout << "EndCap : " << nthe << " " << rtn_nbardtheta << endl;
								if(cos(nthe) > 0){
										//			rtn_nbartheta = acos(rdm_nbarphi->Rndm()*0.06+0.86);
										if(cos(rtn_nbartheta) > 0.92){
												rtn_nbartheta = nthe + 1.5*fabs(rtn_nbardtheta);
										}
										else if(cos(rtn_nbartheta) < 0.86){
												rtn_nbartheta = nthe - 2.5*fabs(rtn_nbardtheta);
										}
								}
								else if(cos(nthe) < 0) {
										//	rtn_nbartheta = acos(-rdm_nbarphi->Rndm()*0.06 - 0.86);
										if(cos(rtn_nbartheta) < -0.92){
												rtn_nbartheta = nthe - 1.5*fabs(rtn_nbardtheta);
										}
										else if(cos(rtn_nbartheta) > -0.86){
												rtn_nbartheta = nthe + 2.5*fabs(rtn_nbardtheta);
										}
								}
								//		rtn_nbartheta = nthe + 3*rtn_nbardtheta;

						}
				}

				if(valid_count1 == 1000 || valid_count2 == 1000){
						samplingFailure++;
				}

				if(rtn_nbarphi > CLHEP::pi){
						rtn_nbarphi  =  rtn_nbarphi - 2*CLHEP::pi;
				}
				else if(rtn_nbarphi < -CLHEP::pi){
						rtn_nbarphi = 2*CLHEP::pi + rtn_nbarphi;
				}

				nbar_posi.theta = rtn_nbartheta;
				nbar_posi.phi = rtn_nbarphi;
				nbar_posi.dphi = rtn_nbardphi;
				nbar_posi.dtheta = rtn_nbardtheta;
				Dxyz->DtheDphi_to_DxDyDz(rtn_nbartheta, rtn_nbarphi, rtn_nbardtheta, rtn_nbardphi);
				nbar_posi.x = Dxyz->getPosition().x();
				nbar_posi.y = Dxyz->getPosition().y();
				nbar_posi.z = Dxyz->getPosition().z();
				nbar_posi.dx = sqrt(Dxyz->getDxDyDz()[0][0]);
				nbar_posi.dy = sqrt(Dxyz->getDxDyDz()[1][1]);
				nbar_posi.dz = sqrt(Dxyz->getDxDyDz()[2][2]);
				setErrorMatrix();
				okratio = 1;
				delete gnbar_angle0; gnbar_angle0= NULL;
				delete gnbar_theta0; gnbar_theta0=NULL;
				delete gnbar_energy; gnbar_energy=NULL;
				delete gnbar_angle1; gnbar_angle1=NULL;
				delete gnbar_theta1; gnbar_theta1=NULL;
				delete gnbar_thetaratio; gnbar_thetaratio=NULL;
				delete gnbar_hits; gnbar_hits= NULL;
				delete gnbar_secmom; gnbar_secmom=NULL;
		}
		else{
				okratio = 0;
		}
	//	delete eff_interpolation; eff_interpolation=NULL;
		energy_flg = 1;

		return okratio;
}



double  AntiNeutronTrk::energy(){
		if(okratio == 0 ) {
				return -1;
		}
		else if(okratio == 2 || okratio ==1){
				return nbar_posi.energy;
		}
		else{
				cout << "AntiNeutronTrk Wrong!" << endl;
				abort();
		}
}



int  AntiNeutronTrk::numHits(){
		if(okratio == 0 ) {
				return -1;
		}
		else if(okratio ==1 || okratio ==2){
				return nbar_posi.hits;
		}
		else{
				cout << "AntiNeutronTrk Wrong!" << endl;
				abort();
		}
}

double  AntiNeutronTrk::secondMoment(){
		if(okratio == 0 ) {
				return -1;
		}
		else if(okratio ==1 || okratio ==2){
				return nbar_posi.secmom;
		}
		else{
				cout << "AntiNeutronTrk Wrong!" << endl;
				abort();
		}
}

double  AntiNeutronTrk::theta(){
		if(okratio == 0 ) {
				return 999;
		}
		else if(okratio ==1 || okratio == 2){
				return nbar_posi.theta;
		}
		else{
				cout << "AntiNeutronTrk Wrong!" << endl;
				abort();
		}
}

double  AntiNeutronTrk::phi(){
		if(okratio == 0 ) {
				return 999;
		}
		else if(okratio ==1 || okratio == 2){
				return nbar_posi.phi;
		}
		else{
				cout << "AntiNeutronTrk Wrong!" << endl;
				abort();
		}
}

double  AntiNeutronTrk::dtheta(){
		if(okratio == 0 ) {
				return 999;
		}
		else if(okratio ==1 || okratio == 2){
				return nbar_posi.dtheta;
		}
		else{
				cout << "AntiNeutronTrk Wrong!" << endl;
				abort();
		}
}

double  AntiNeutronTrk::dphi(){
		if(okratio == 0 ) {
				return 999;
		}
		else if(okratio ==1 || okratio ==2){
				return nbar_posi.dphi;
		}
		else{
				cout << "AntiNeutronTrk Wrong!" << endl;
				abort();
		}
}

double  AntiNeutronTrk::x(){
		if(okratio == 0 ) {
				return 999;
		}
		else if(okratio ==1 || okratio ==2){
				return nbar_posi.x;
		}
		else{
				cout << "AntiNeutronTrk Wrong!" << endl;
				abort();
		}
}

double  AntiNeutronTrk::y(){
		if(okratio == 0 ) {
				return 999;
		}
		else if(okratio ==1 || okratio ==2){
				return nbar_posi.y;
		}
		else{
				cout << "AntiNeutronTrk Wrong!" << endl;
				abort();
		}
}

double  AntiNeutronTrk::z(){
		if(okratio == 0 ) {
				return 999;
		}
		else if(okratio ==1 || okratio ==2){
				return nbar_posi.z;
		}
		else{
				cout << "AntiNeutronTrk Wrong!" << endl;
				abort();
		}
}

double  AntiNeutronTrk::dx(){
		if(okratio == 0 ) {
				return 999;
		}
		else if(okratio ==1 || okratio ==2){
				return nbar_posi.dx;
		}
		else{
				cout << "AntiNeutronTrk Wrong!" << endl;
				abort();
		}
}
double  AntiNeutronTrk::dy(){
		if(okratio == 0 ) {
				return 999;
		}
		else if(okratio ==1 || okratio ==2){
				return nbar_posi.dy;
		}
		else{
				cout << "AntiNeutronTrk Wrong!" << endl;
				abort();
		}
}
double  AntiNeutronTrk::dz(){
		if(okratio == 0 ) {
				return 999;
		}
		else if(okratio ==1 || okratio ==2){
				return nbar_posi.dz;
		}
		else{
				cout << "AntiNeutronTrk Wrong!" << endl;
				abort();
		}
}

void AntiNeutronTrk::setEnergy(double nbar_energy){
		nbar_posi.energy = nbar_energy;
		energy_flg = -1;
}


bool AntiNeutronTrk::validp4(){
		if(okratio == 0 ) {
				return false;
		}
		else if(okratio ==1 || okratio ==2){
				return true;
		}
		else{
				return false;
		}
}

int AntiNeutronTrk::samplingFailureN(){
		return samplingFailure;
}
int AntiNeutronTrk::totalSampling(){
		return totalsampling;
}

double myfunc(double phi1, double theta1, double phi2, double theta2, double angle){
		double x2 = sin(theta2)*cos(phi2);
		double y2 = sin(theta2)*sin(phi2);
		double z2 = cos(theta2);
		double x1 = sin(theta1)*cos(phi1);
		double y1 = sin(theta1)*sin(phi1);
		double z1 = cos(theta1);
		TVector3 nbar_pos, nbar_mc;
		nbar_pos.SetXYZ(x2, y2, z2);
		nbar_mc.SetXYZ(x1, y1, z1);
		return nbar_pos.Angle(nbar_mc)  - angle;
}


TString Itoa(int idx, int base){
		TString str;
		do{
				int remainder = idx%base;
				char remainder_str = remainder + '0';
				str = remainder_str + str;
				idx = idx/base;
		}while(idx != 0);
		//	cout << str << endl;
		return str;

}

void AntiNeutronTrk::setEventCut(const double energy, const double secmom, const int hits){
		Nbar_Hits = hits;
		Nbar_SecMom = secmom;
		Nbar_Energy = energy;
}


