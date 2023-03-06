#ifndef IAnti_Neutron_Correction_Svc
#define IAnti_Neutron_Correction_Svc

#include <iostream>
#include "GaudiKernel/IService.h"

/* Declaration of the interface ID */
static const InterfaceID IID_IAntiNeutronCorrectionSvc("IAntiNeutronCorrectionSvc", 1, 0);


class IAntiNeutronCorrectionSvc : virtual public IService
{
public:
  virtual ~IAntiNeutronCorrectionSvc() {}

  static const InterfaceID& interfaceID() { return IID_IAntiNeutronCorrectionSvc; }
//  virtual AntiNeutronTrk *setErrorMatrix(RecEmcShower *nbarTrk ) = 0;
//  virtual AntiNeutronTrk *setErrorMatrix(AntiNeutronTrk *nbarTrk) = 0;
//  virtual AntiNeutronTrk *setAntiNeutronTrk(const HepLorentzVector initp4, const HepLorentzVector initposi, const HepLorentzVector IP,  int statistical_uncertainty = 0) = 0;

};

#endif
