#include "rotations.h"

namespace Rot {

  //**********************************************************************************************************************

    /* Rotation Matrix using euler angles z-y-z convention */

  //**********************************************************************************************************************

    itpp::mat EulerRot_t(double alpha, double beta, double gamma)
    {
      // to act upon cartesian 3-vectors
      // R(a,b,g) = Rz(a) Ry(b) Rz(g)
      
      itpp::mat Rzg(3,3); Rzg.zeros();
      Rzg(0,0) = cos(gamma); Rzg(0,1) = -1.0 * sin(gamma);
      Rzg(1,0) = -1.0 * Rzg(0,1) ; Rzg(1,1) = Rzg(0,0);
      Rzg(2,2) = 1.0;
      
      itpp::mat Ryb(3,3); Ryb.zeros();
      Ryb(0,0) = cos(beta); Ryb(0,2) = sin(beta);
      Ryb(1,1) = 1.0;
      Ryb(2,0) = -1.0 * Ryb(0,2); Ryb(2,2) = Ryb(0,0);
      
      itpp::mat Rza(3,3); Rza.zeros();
      Rza(0,0) = cos(alpha); Rza(0,1) = -1.0 * sin(alpha);
      Rza(1,0) = -1.0 * Rza(0,1); Rza(1,1) = Rza(0,0);
      Rza(2,2) = 1.0;
      
      return Rza * Ryb * Rzg;

    };

    XMLArray::Array<int> EulerRotVec_t(Hadron::CubicCanonicalRotation_t ref_angles_z, Hadron::CubicCanonicalRotation_t ref_angles, XMLArray::Array<int> in)
    {
        // to act upon cartesian 3-vectors
        // R(a,b,g) = Rz(a) Ry(b) Rz(g)

        itpp::vec vec(3), ans(3);
        itpp::mat R_z(3,3), R(3,3);
        XMLArray::Array<int> out(3);

        vec(0) = in[0]; vec(1) = in[1]; vec(2) = in[2];

        R_z = Rot::EulerRot_t(ref_angles_z.alpha, ref_angles_z.beta, ref_angles_z.gamma);
        R   = Rot::EulerRot_t(ref_angles.alpha, ref_angles.beta, ref_angles.gamma); 
    
        ans = R * itpp::inv(R_z) * vec ;
        out[0] = (int)round(ans(0)); out[1] = (int)round(ans(1)); out[2] = (int)round(ans(2));

        return out;
    };

  
};
