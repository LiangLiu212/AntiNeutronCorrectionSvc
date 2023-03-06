#ifndef AN_XYZ_ERRORMATRX_H
#define AN_XYZ_ERRORMATRX_H
#include "CLHEP/Geometry/Point3D.h"
#include "EmcRecEventModel/RecEmcShower.h"

typedef HepGeom::Point3D<double> HepPoint3D;

class AN_XYZ_ErrorMatrx{

		private:
				HepSymMatrix m_matrix;
				HepPoint3D m_position;
		public:
				AN_XYZ_ErrorMatrx();
				~AN_XYZ_ErrorMatrx(){;};
				void DtheDphi_to_DxDyDz(double theta, double phi, double Dtheta, double Dphi);
				void setDxDyDz(const HepSymMatrix matrix) {m_matrix = matrix;};
				void setPosition( const HepPoint3D position) {m_position = position;};
				HepSymMatrix getDxDyDz() {return m_matrix;};
				HepPoint3D getPosition() {return m_position;};

};
#endif
