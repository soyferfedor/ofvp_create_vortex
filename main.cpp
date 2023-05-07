#include "beam.hpp"
#include "excep.hpp"
#include "log.hpp"
#include <vector>

using namespace ofvp::create_vortex;

/* ============================================================
	Project:	Magister's diploma
	
	Author:		F. Soifer
	
	Date:		2023-03-02

	Description:	MAIN file

	Comments:	!!!!!!!!!!!!!!!!!!!!!!!

AVAILABLE FUNCTIONS:
	FUNCTION												STATUS		UPD_DATE
	- initial_gauss()											ON		2023-02-15
	- initial_vortex()											ON		2023-03-02
	- initial_gaussvortex()											ON		2023-03-02
	- initial_ellips(double b_over_a)									ON		2023-03-01
	- initial_picture(std::string file_name, double N_x, double N_y)					OFF		2023-02-15
	- initial_2_max_points(double r_0_over_r_1)								ON		2023-02-19
	- initial_1_max_point(double r_0_over_r_1, double ang_0, double part_ampl_0)				ON		2023-04-24
	- initial_zero()											ON		2023-05-03
	- add_noise_phasescreen(std::string file_name, double sigma)						ON		2023-05-02
	- double find_r_beam()											ON		2023-05-03
	- step_diffraction(size_t num_of_calculating_steps_in_1_Z_diff, double z_finish)			ON		2023-02-15
	- phase_plate()												ON		2023-02-15
	- double_phase_plate(double num_r_new_center, char type_circle_1,
					char type_circle_2, char cut_1, char cut_2)				ON		2023-04-17
	- print_amplitude(std::string str, int num)								ON		2023-02-15
	- print_amplitude_center(double part_of_all_grid_in_1_D, std::string str, int num)			ON		2023-02-15
	- print_phase(std::string str, int num)									ON		2023-02-15
	- print_phase_center(double part_of_all_grid_in_1_D, std::string str, int num)				ON		2023-02-20
	- print_spectrum_x(std::string str, int num)								ON		2023-02-15		do not use 2d graphs (in work)
	- print_spectrum_center_x(double part_of_all_grid_in_1_D, std::string str, int num)			OFF		2023-02-27		do not use 2d graphs (in work)
	- print_spectrum(std::string str, int num)								ON		2023-02-27
	- print_spectrum_center(double part_of_all_grid_in_1_D, std::string str, int num)			OFF		2023-02-27
	- print_spectrum_shifted()										OFF		2023-02-15		now it is in print_spectrum()
	- double compare_angles_between_maxs (double ang_before)						ON		2023-04-24		returns double (!)

CLOSE FUTURE PLANS:
	- вихрь с флуктуацией на кольце-максимуме
	- interesting phase plates (2 singular points)
	- отраж от границ сетки (доработать)
	- добавить в инициализ. ф-ии возможность выбирать тип (фаз или ампл) и параметры шума (реализ шум парам в отд классе)
	- реализовать print_spectrum_center()
FUTURE PLANS:
	- elliptical beam (not just circular)
	- create an option to digitize an image of an arbitrary beam
	- реализ работу всех ф-ий для произвольной прямоуг сетки
	- созд функцию print которая будет корректно выводить лог и мэйк-инфо, и генер исключ
	- добавить генератор случайных флуктуаций с заданными парм (в т.ч. радиусом корреляции)




 ============================================================ */

int main() {
	double wavelength = 8e-7;													 // in meters
	double initial_radius = 1.0;                                                 // r/r_0
	double max_amplitude_of_initial_Hauss = 1.0;								 // I/I_0
	unsigned char topological_charge = 1;
	size_t num_of_points_in_grid_in_1D = 256;										 // N; total points = N*N
	double x_min = -8.0;                                  		                 // r/r_0
	double x_max = 8.0;                                      		             // r/r_0
	size_t num_of_calculating_steps_in_1_Z_diff = 10;								 // NUM OF STEPS IN Z_diff
	double z_finish = 1.0;

	Beam b1 (num_of_points_in_grid_in_1D, x_min, x_max,
			    num_of_points_in_grid_in_1D, x_min, x_max,
			    wavelength, topological_charge, initial_radius, initial_radius);
	//b1.initial_gauss().print_amplitude_center(0.6, "amplitude_init").step_diffraction(1, 1.0).print_amplitude_center(0.6, "amplitude_after_dif");
	//b1.initial_ellips(2.5).print_spectrum("spectrum").print_amplitude_center(0.6, "amplitude_init").step_diffraction(1, 3.0).print_amplitude_center(0.6, "amplitude_after_dif");
	//b1.initial_vortex().print_spectrum("spectrum").print_amplitude_center(0.8, "amplitude_init").step_diffraction(1, 1.0).print_amplitude_center(0.8, "amplitude_after_dif");
	//b1.initial_ellips(2.0).print_amplitude("amplitude_init").print_spectrum("spectrum").step_diffraction(1, 1.0);
	//b1.initial_2_max_points().print_phase("phase_init").print_amplitude("amplitude_init").print_spectrum("spectrum").step_diffraction(1, 0.1).print_phase("phase_0,1zd").print_amplitude("amplitude_after_0,1_z_dif").step_diffraction(1, 0.1).print_phase("phase_0,2zd").print_amplitude("amplitude_after_0,2_z_dif");
	
	
	//b1.initial_vortex().print_phase_center(0.6, "phase_ini");
	//b1.initial_gauss().phase_plate().print_phase_center(0.5);
	//b1.initial_gauss().print_amplitude_center(0.6, "amplitude_init_gauss").double_phase_plate(1.0, 0, 1).print_phase_center(0.2).print_amplitude_center(0.6, "amplitude_after_phase_plate").step_diffraction(1, 0.05).print_amplitude_center(0.6, "amplitude_after_dif_0_05");




	// double max vortex beams
	//b1.initial_2_max_points(1.0).print_spectrum("spectrum").print_amplitude_center(0.35, "amplitude_init").step_diffraction(1, 0.2).print_amplitude_center(0.35, "amplitude_after_dif_0,2_z_dif").step_diffraction(1, 0.3).print_amplitude_center(0.35, "amplitude_after_dif_0,5_z_dif").step_diffraction(1, 0.5).print_amplitude_center(0.35, "amplitude_after_dif_1,0_z_dif");
	//b1.initial_2_max_points(3.0).print_spectrum("spectrum").print_amplitude_center(0.25, "amplitude_init").step_diffraction(1, 0.2).print_amplitude_center(0.25, "amplitude_after_dif_0,2_z_dif").step_diffraction(1, 0.3).print_amplitude_center(0.25, "amplitude_after_dif_0,5_z_dif").step_diffraction(1, 0.5).print_amplitude_center(0.25, "amplitude_after_dif_1,0_z_dif");
	//b1.initial_2_max_points(3.0).print_spectrum("spectrum").print_amplitude_center(0.25, "amplitude_init").step_diffraction(1, 0.005).print_amplitude_center(0.25, "amplitude_after_dif_0,005_z_dif").step_diffraction(1, 0.005).print_amplitude_center(0.25, "amplitude_after_dif_0,01_z_dif").step_diffraction(1, 0.01).print_amplitude_center(0.25, "amplitude_after_dif_0,02_z_dif");


	// 1 max point example
	//b1.initial_1_max_point(1.1, 0.0, 0.8).print_amplitude_center(0.6).step_diffraction(1, 0.02).print_amplitude_center(0.6, "amplitude_after_dif_0_02").step_diffraction(1, 0.02).print_amplitude_center(0.6, "amplitude_after_dif_0_04").step_diffraction(1, 0.02).print_amplitude_center(0.6, "amplitude_after_dif_0_06").step_diffraction(1, 0.02).print_amplitude_center(0.6, "amplitude_after_dif_0_08");


	
/*	Beam *beam_arr = static_cast<Beam*>(operator new[] (100 * sizeof(Beam)));
	for (size_t i = 0; i < 100; ++i) {
		new (beam_arr + i) Beam(num_of_points_in_grid_in_1D, x_min, x_max,
					num_of_points_in_grid_in_1D, x_min, x_max,
					wavelength, topological_charge, initial_radius, initial_radius);
	}
*/

	//b1.initial_ellips(2.0).print_amplitude_center(0.7).phase_plate().step_diffraction(1, 0.5).print_amplitude_center(0.7, "amplitude_after_dif_0_5");

	//b1.initial_gauss().print_amplitude_center(0.6, "amplitude_init_gauss").double_phase_plate(1.0, 0, 1).print_phase_center(0.2).print_amplitude_center(0.6, "amplitude_after_phase_plate").step_diffraction(1, 0.05).print_amplitude_center(0.6, "amplitude_after_dif_0_05").print_phase_center(0.6);


	//b1.initial_vortex().print_amplitude();
	b1.initial_zero().add_noise_phasescreen("test/screen000.txt").print_amplitude();


	return 0;
}
