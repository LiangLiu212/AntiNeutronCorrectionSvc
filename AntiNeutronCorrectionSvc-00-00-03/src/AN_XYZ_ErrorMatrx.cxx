#include "AntiNeutronCorrectionSvc/AN_XYZ_ErrorMatrx.h"
AN_XYZ_ErrorMatrx::AN_XYZ_ErrorMatrx(){
		setDxDyDz(HepSymMatrix(3,0));
		setPosition(HepPoint3D(0.0,0.0,0.0));

}

void AN_XYZ_ErrorMatrx::DtheDphi_to_DxDyDz(double theta, double phi, double Dtheta, double Dphi){
		double x, y, z;
		double dx, dy, dz;

		if(cos(theta) > 0.826){
				z =140.0;
				x = 140.0*tan(theta)*cos(phi);
				y = 140.0*tan(theta)*sin(phi);
				dx = 140.0*tan(theta + Dtheta) * Dphi *sin(phi);
				dy = 140.0*tan(theta + Dtheta) * Dphi *cos(phi);
				dz = 140.0*Dtheta/(cos(theta)*cos(theta));
		}
		else if(cos(theta) < -0.826){
				z = -140.0;
				x = -140.0*tan(theta)*cos(phi);
				y = -140.0*tan(theta)*sin(phi);
				dx = -140.0*tan(theta + Dtheta) * Dphi *sin(phi);
				dy = -140.0*tan(theta + Dtheta) * Dphi *cos(phi);
				dz = 140.0*Dtheta/(cos(theta)*cos(theta));
		}

		else{
				z = 95.5 / tan(theta);
				x = 95.5* cos(phi);
				y = 95.5 * sin(phi);
				dx = 95.5 * Dphi *sin(phi);
				dy = 95.5 * Dphi *cos(phi);
				dz = 95.5*Dtheta/(sin(theta)*sin(theta));
		}
		double d[3]; 
		d[0] = dx;
		d[1] = dy;
		d[2] = dz;

		for(int i = 0; i < 3; i++){
				for(int j = 0; j < 3; j++){
						m_matrix[i][j]=d[i]*d[j];
				}
		}
		m_position.set(x, y, z);
}
