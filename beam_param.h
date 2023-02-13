#ifndef BEAM_PARAM_H_INCLUDED
#define BEAM_PARAM_H_INCLUDED

/*===================================================================
	Project:	Magister's diploma
	
	Author:		F. Soifer
	
	Date:		2023-02-12

	Description:	Class Beam_parameters (Beam preparings)

	Comments:	!!!!!!!!!!!!!!!!!!!!!!!

  ===================================================================*/


#include <cmath>


namespace ofvp {
namespace create_vortex{

// Class Beam_geometry

// Include:	Radiuses on both axes (circle or
//		ellips - circle is default version.
//		phi_0 is an initial angle from +x axe
//		turn up to the right (main diag of
//		an ellips)


	class Beam_geometry {
		double r_x_0;
		double r_y_0;
		double phi_0; // is not used
	public:
		Beam_geometry (double r_x_0, double r_y_0, double phi_0 = 0.0):
			r_x_0(r_x_0), r_y_0(r_y_0), phi_0(phi_0)
		{}
		double get_r_x_0() const {
			return r_x_0;
		}
		void set_r_x_0(double x) {
			r_x_0 = x;
		}
		double get_r_y_0() const {
			return r_y_0;
		}
		void set_r_y_0(double y) {
			r_y_0 = y;
		}
		double get_phi_0() const {
			return phi_0;
		}
		void set_phi_0(double ang) {
			phi_0 = ang;
		}
		double degr_of_ellipt() const {
			if (r_x_0 >= r_y_0) return r_x_0 / r_y_0;
			return r_y_0 / r_x_0;
		}
	};


// Class Beam_parameters

// Include:	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

// Notes:	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	class Beam_parameters: public Beam_geometry {
		double lambda; // wavelength
		double k; // wave number k
		double z_diffraction; // diffraction length of the beam
		unsigned char m; //topological charge
	public:
		Beam_parameters (double l, unsigned char M, double r_x_0, double r_y_0, double phi_0 = 0.0):
			lambda(l), m(M), Beam_geometry(r_x_0, r_y_0, phi_0)
		{
			k = 2.0 * M_PI / lambda;
			z_diffraction = k * (get_r_x_0()+get_r_y_0())/2 * (get_r_x_0()+get_r_y_0())/2; // is it correct with elleptic beams ??????
		}
		double get_lambda() const {
			return lambda;
		}
		void set_lambda(double x) {
			lambda = x;
			k = 2.0 * M_PI / lambda;
			z_diffraction = k * (get_r_x_0()+get_r_y_0())/2 * (get_r_x_0()+get_r_y_0())/2; // is it correct with elleptic beams ??????
		}
		double get_k() const {
			return k;
		}
		double get_z_diffraction() const {
			return z_diffraction;
		}
		void set_m(unsigned char M) {
			m = M;
		}
		unsigned char get_m() const {
			return m;
		}
	};


}
}


#endif // BEAM_PARAM_H_INCLUDED
