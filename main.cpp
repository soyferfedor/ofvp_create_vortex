#include "beam.hpp"
#include <memory>

using namespace ofvp::create_vortex;

/* ============================================================
	Project:	Magister's diploma
	
	Author:		F. Soifer
	
	Date:		2023-02-15

	Description:	MAIN file

	Comments:	!!!!!!!!!!!!!!!!!!!!!!!

AVAILABLE FUNCTIONS:
	FUNCTION												STATUS		UPD_DATE
	- initial_gauss()											ON		2023-02-15
	- initial_vortex()											ON		2023-02-15
	- step_diffraction(size_t num_of_calculating_steps_in_1_Z_diff, double z_finish)			ON		2023-02-15
	- phase_plate()												ON		2023-02-15
	- print_amplitude(std::string str, int num)								ON		2023-02-15
	- print_amplitude_center(double part_of_all_grid_in_1_D, std::string str, int num)			ON		2023-02-15
	- print_phase(std::string str, int num)									ON		2023-02-15
	- print_spectrum(std::string str, int num)								ON		2023-02-15
	- print_spectrum_shifted()										OFF		2023-02-15
	- initial_picture(std::string file_name, double N_x, double N_y)					OFF		2023-02-15

CLOSE FUTURE PLANS:
	- 2 or more peaks in one beam (2 or 3 mini-beams in one)
	- noise (phase and amplitude types)
	- interesting phase plates (2 singular points)
FUTURE PLANS:
	- elliptical beam (not just circular)
	- create an option to digitize an image of an arbitrary beam




 ============================================================ */

int main() {
	double wavelength = 8e-7;													 // in meters
	double initial_radius = 1.0;                                                 // r/r_0
	double max_amplitude_of_initial_Hauss = 1.0;								 // I/I_0
	unsigned char topological_charge = 1;
	size_t num_of_points_in_grid_in_1D = 2048;										 // N; total points = N*N
	double x_min = -8.0;                                  		                 // r/r_0
	double x_max = 8.0;                                      		             // r/r_0
	size_t num_of_calculating_steps_in_1_Z_diff = 1000;								 // NUM OF STEPS IN Z_diff
	double z_finish = 0.01;		

	std::unique_ptr<Beam> b1 (new Beam(num_of_points_in_grid_in_1D, x_min, x_max,
			    num_of_points_in_grid_in_1D, x_min, x_max,
			    wavelength, topological_charge, initial_radius, initial_radius));
	return 0;
}
