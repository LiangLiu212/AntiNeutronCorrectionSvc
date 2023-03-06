#include "GaudiKernel/DeclareFactoryEntries.h"

#include "AntiNeutronCorrectionSvc/AntiNeutronCorrectionSvc.h"

DECLARE_SERVICE_FACTORY( AntiNeutronCorrectionSvc )

DECLARE_FACTORY_ENTRIES( AntiNeutronCorrectionSvc ) { 
  DECLARE_SERVICE( AntiNeutronCorrectionSvc );
}
