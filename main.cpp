#include "beam.hpp"
#include "excep.hpp"
#include "log.hpp"

using namespace ofvp::create_vortex;

/* ============================================================
	Project:	Magister's diploma
	
	Author:		F. Soifer
	
	Date:		2023-02-27

	Description:	MAIN file

	Comments:	!!!!!!!!!!!!!!!!!!!!!!!

AVAILABLE FUNCTIONS:
	FUNCTION												STATUS		UPD_DATE
	- initial_gauss()											ON		2023-02-15
	- initial_vortex()											ON		2023-02-15
	- initial_picture(std::string file_name, double N_x, double N_y)					OFF		2023-02-15
	- initial_2_max_points(double r_0_over_r_1)								ON		2023-02-19
	- step_diffraction(size_t num_of_calculating_steps_in_1_Z_diff, double z_finish)			ON		2023-02-15
	- phase_plate()												ON		2023-02-15
	- print_amplitude(std::string str, int num)								ON		2023-02-15
	- print_amplitude_center(double part_of_all_grid_in_1_D, std::string str, int num)			ON		2023-02-15
	- print_phase(std::string str, int num)									ON		2023-02-15
	- print_phase_center(double part_of_all_grid_in_1_D, std::string str, int num)				ON		2023-02-20
	- print_spectrum_x(std::string str, int num)								ON		2023-02-15		do not use 2d graphs (in work)
	- print_spectrum_center_x(double part_of_all_grid_in_1_D, std::string str, int num)			OFF		2023-02-27		do not use 2d graphs (in work)
	- print_spectrum(std::string str, int num)								ON		2023-02-27
	- print_spectrum_center(double part_of_all_grid_in_1_D, std::string str, int num)			OFF		2023-02-27
	- print_spectrum_shifted()										OFF		2023-02-15		now it is in print_spectrum()

CLOSE FUTURE PLANS:
	- interesting phase plates (2 singular points)
	- отраж от границ сетки (доработать)
	- добавить в инициализ. ф-ии возможность выбирать тип (фаз или ампл) и параметры шума (реализ шум парам в отд классе)
	- реализовать print_spectrum_center()
	- реализовать вывод графиков в 2д
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
	size_t num_of_points_in_grid_in_1D = 128;										 // N; total points = N*N
	double x_min = -8.0;                                  		                 // r/r_0
	double x_max = 8.0;                                      		             // r/r_0
	size_t num_of_calculating_steps_in_1_Z_diff = 10;								 // NUM OF STEPS IN Z_diff
	double z_finish = 1.0;		

	Beam b1 (num_of_points_in_grid_in_1D, x_min, x_max,
			    num_of_points_in_grid_in_1D, x_min, x_max,
			    wavelength, topological_charge, initial_radius, initial_radius);
	b1.initial_2_max_points().print_amplitude("amplitude_init").print_spectrum("spectrum").step_diffraction(1, 0.1).print_amplitude("amplitude_after_0,1_z_dif").step_diffraction(1, 0.1).print_amplitude("amplitude_after_0,2_z_dif");
	
	return 0;
}
