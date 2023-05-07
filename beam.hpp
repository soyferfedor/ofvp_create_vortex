#ifndef BEAM_H_INCLUDED
#define BEAM_H_INCLUDED

/*===================================================================
	Project:	Magister's diploma
	
	Author:		F. Soifer
	
	Date:		2023-02-27

	Description:	Class Beam

	Comments:	!!!!!!!!!!!!!!!!!!!!!!!

  ===================================================================*/


#include "grid_param.hpp"
#include "beam_param.hpp"
#include "work_info.hpp"
#include <fstream>

namespace ofvp {
namespace create_vortex{


	enum Initial_type {gauss = 1, vortex};

	template <typename T>
	std::string to_string(T val) {
		std::ostringstream oss;
		oss << val;
		return oss.str();
	}


/*	class Beam_type {
		Initial_type beam_t;
	public:
		Beam_type (unsigned char t) {
			try {
				if (t == 1)
					beam_t = gauss;
				else if (t == 2)
					beam_t = vortex;
				else
					throw "Code 2";
			}
			catch(...) {
				std::cout << "Incorrect beam type!" << std::endl;
				exit(2);
			}
		}
	};
*/

/*	enum Location_type {center = 1, diag, circle, square};
*/
	
/*	class Beam_numbers {
		unsigned char num;
		Location_type location_t;
	public:
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	};

	class Beam {
		Grid_parameters grid_p;
		Beam_parameters beam_p;
	public:
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	};
*/

	class Beam: public Grid_parameters, Beam_parameters {
		void multipl_by_exp(std::complex<double>& obj, double kx, double ky, double k, double dz, double lambda) {
			double fi = (kx*kx + ky*ky) * dz / (2.0 * k);
			obj *= std::exp(std::complex<double>(0.0, fi));
			//double real = obj.real(), imag = obj.imag();
			//obj = std::complex<double> (real * cos(fi) - imag * sin(fi), real * sin(fi) + imag * cos(fi));
		}
		Beam& propagation_in_1_step(double dz) {			// working just with N_x = N_y !!!!!!!!!!!!!!!!!!!!!! // rewrite more optimized
			std::complex<double>* in = set_ptr_in();
			std::complex<double>* out = set_ptr_out();
			int N = static_cast<int>(get_N_x());
			double dk = 2.0 * M_PI / (N * get_dx());
			fftw_execute_forward();
			for (int i = 0; i < N; ++i) {
				for (int j = 0; j < N; ++j) {
					if (i <= N/2-1) {
						if (j <= N/2-1) 
							multipl_by_exp(out[i*N + j],
											j * dk, i * dk,
											get_k(), dz, get_lambda());
						else 
							multipl_by_exp(out[i*N + j],
											(j-N) * dk, i * dk,
											get_k(), dz, get_lambda());
					} else {
						if (j <= N/2-1) 
							multipl_by_exp(out[i*N + j],
											j * dk, (i-N) * dk,
											get_k(), dz, get_lambda());
						else 
							multipl_by_exp(out[i*N + j],
											(j-N) * dk, (i-N) * dk,
											get_k(), dz, get_lambda());
					}
				}
			}
			fftw_execute_backward();
			for (size_t i = 0; i < get_num_points(); ++i)
				in[i] /= get_num_points();
			return *this;
		}
		std::complex<double> phase_plate_1_point(std::complex<double> E0, int i, int j) { // working just with N_x = N_y !!!!!!!!!!!!!!!!!!!!!!
			double xj, yi, fi, real, imag, Re, Im;
			xj = (j - (int)get_N_x()/2.0 + 0.5) * get_dx();
			yi = (i - (int)get_N_x()/2.0 + 0.5) * get_dx();
			fi = atan2(yi, xj);
			real = E0.real();
			imag = E0.imag();
			Re = (real*cos(get_m()*fi) - imag*sin(get_m()*fi))
						//* pow(sqrt(xj*xj + yi*yi)/beam.get_r0(), beam.get_m())
						* exp(- (xj*xj + yi*yi) /2.0/get_r_x_0()/get_r_x_0());
			Im = (imag*cos(get_m()*fi) + real*sin(get_m()*fi))
						//* pow(sqrt(xj*xj + yi*yi)/beam.get_r0(), beam.get_m())
						* exp(- (xj*xj + yi*yi) /2.0/get_r_x_0()/get_r_x_0());
			return std::complex<double> (Re, Im);
		}
		std::complex<double> double_phase_plate_1_point(std::complex<double> E0, int i, int j, double l, char type_circle_1, char type_circle_2) { // working just with N_x = N_y !!!!!!!!!!!!!!!!!!!!!!
			double xj, yi, fi, real, imag, Re, Im;
			int N_one_side = floor(double(get_N_x())*l/(get_x_max()-get_x_min())/2.0);
			int N_half = floor(double(get_N_x())/2.0);
			xj = (j - (int)get_N_x()/2.0 + 0.5) * get_dx();
			if (i < N_half) {
				yi = (i - N_half + 0.5 + N_one_side) * get_dx();
				fi = pow(-1, type_circle_1) * atan2(yi, xj);
			}
			else {
				yi = (i - N_half + 0.5 - N_one_side) * get_dx();
				fi = pow(-1, type_circle_2) * atan2(yi, xj);
			}
			real = E0.real();
			imag = E0.imag();
			Re = (real*cos(get_m()*fi) - imag*sin(get_m()*fi))
						//* pow(sqrt(xj*xj + yi*yi)/beam.get_r0(), beam.get_m())
						* exp(- (xj*xj + yi*yi) /2.0/get_r_x_0()/get_r_x_0());
			Im = (imag*cos(get_m()*fi) + real*sin(get_m()*fi))
						//* pow(sqrt(xj*xj + yi*yi)/beam.get_r0(), beam.get_m())
						* exp(- (xj*xj + yi*yi) /2.0/get_r_x_0()/get_r_x_0());
			return std::complex<double> (Re, Im);
		}
		size_t find_max() const{						// working just with N_x = N_y !!!!!!!!!!!!!!!!!!!!!!
			const std::complex<double>* in = get_ptr_in();
			size_t N = get_N_x();
			double A_max = 0.0;
			size_t n_max = 0;
			double A_now;
			for (size_t i = 0; i < N; ++i) {
				for (size_t j = 0; j < N; ++j) {
					A_now = sqrt(in[i*N + j].real()*in[i*N + j].real() + in[i*N + j].imag()*in[i*N + j].imag());
					if (A_now > A_max) {
						A_max = A_now;
						n_max = i*N + j;
					}
				}
			}
			return n_max;
		}
		size_t find_min() const{ // rewrite
			size_t N = get_N_x();
			return N/2*N + N/2;
		}
		double ampl_in_coord(double x, double y) const{						// working just with N_x = N_y !!!!!!!!!!!!!!!!!!!!!!
			const std::complex<double>* in = get_ptr_in();
			double dx = get_dx();
			size_t N = get_N_x();
			double dist_now;
			size_t i_close, j_close;
			for (int i = 0; i < N; ++i) {
				for (int j = 0; j < N; ++j) {
					dist_now = sqrt((coord_x(i)-x)*(coord_x(i)-x) + ((coord_y(j)-y)*(coord_y(j)-y)));
					if (dist_now < dx) {
						i_close = i;
						j_close = j;
						goto ret_res;
					}
				}
			}
			exit(9); // if goto didn't work
			ret_res:
			return sqrt(in[i_close*N + j_close].real()*in[i_close*N + j_close].real() + in[i_close*N + j_close].imag()*in[i_close*N + j_close].imag());
		}
		double find_r_beam() const {
			int N = (int)get_N_x();
			int idx_min = (int)find_min();
			int idx_max = (int)find_max();
			double x_min = coord_x(idx_min);
			double y_min = coord_y(idx_min);
			double x_max = coord_x(idx_max);
			double y_max = coord_y(idx_max);
			return std::sqrt((x_max-x_min)*(x_max-x_min) + (y_max-y_min)*(y_max-y_min));
		}
		void receive_A_of_alpha(double *arr, size_t num_points = 360, double angle_beg = 0.0) {
			size_t point_min = find_min();
			int j_min = (int)point_min % (int)get_N_x();
			int i_min = (int)point_min / (int)get_N_x();
			double x_min = coord_x(i_min);
			double y_min = coord_y(j_min);
			double radius = find_r_beam();
			double dalpha = 360.0 / num_points;
			double x, y;
			arr = new (std::nothrow) double[num_points];
			if (!arr) {
				std::cout << "Bad alloc memory in receive_A_of_alpha() function!";
				exit(10);
			}
			for (int i = 0; i < (int)num_points; i++) {
				x = x_min + radius * sin(angle_beg + i * dalpha);
				y = y_min + radius * cos(angle_beg + i * dalpha);
				arr[i] = ampl_in_coord(x, y);
			}
		}
	public:
		Beam (size_t N_x, double x_min, double x_max, size_t N_y, double y_min, double y_max, double l, unsigned char M, double r_x_0, double r_y_0, double phi_0 = 0.0):
			Grid_parameters(N_x, x_min, x_max, N_y, y_min, y_max), Beam_parameters(l, M, r_x_0, r_y_0, phi_0) {}
		Beam& initial_gauss() {						// working just with N_x = N_y !!!!!!!!!!!!!!!!!!!!!!
			std::complex<double>* in = this->set_ptr_in();
			size_t N = get_N_x();
			double min = get_x_min();
			double dx = get_dx();
			double x, y;
			for (size_t i = 0; i < N; ++i) {
				for (size_t j = 0; j < N; ++j) {
					x = min + i*dx;
					y = min + j*dx;
					in[i*N + j] = std::complex<double>(exp(-(x*x + y*y) / 2), 0.);
				}
			}
			return *this;
		}
		Beam& initial_vortex() { // working just with N_x = N_y !!!!!!!!!!!!!!!!!!!!!!
			double x, y, phi;
			std::complex<double>* in = set_ptr_in();
			size_t N = get_N_x();
			for (size_t i = 0; i < N; ++i) {
				for (size_t j = 0; j < N; ++j) {
					x = coord_x(i);
					y = coord_y(j);
					phi = atan2(y, x) * (1.0 /*+ bp_.noise_amplitude()*sin(bp_.noise_theta()*hypot(x,y))*/);
					in[i*N + j] = std::pow(x*x + y*y, get_m()/2.0) * std::exp(-(x*x + y*y) / 2 / std::pow(get_r_x_0(), 2)) * exp(std::complex<double>(0, get_m()*phi));
				}
			}
			return *this;
		}
		Beam& initial_gaussvortex() { // working just with N_x = N_y !!!!!!!!!!!!!!!!!!!!!!
			double x, y, phi;
			std::complex<double>* in = set_ptr_in();
			size_t N = get_N_x();
			for (size_t i = 0; i < N; ++i) {
				for (size_t j = 0; j < N; ++j) {
					x = coord_x(i);
					y = coord_y(j);
					phi = atan2(y, x) * (1.0 /*+ bp_.noise_amplitude()*sin(bp_.noise_theta()*hypot(x,y))*/);
					in[i*N + j] = std::exp(-(x*x + y*y) / 2 / std::pow(get_r_x_0(), 2)) * exp(std::complex<double>(0, get_m()*phi));
				}
			}
			return *this;
		}
		Beam& initial_ellips(double b_over_a = 1.5) {  // working just with N_x = N_y !!!!!!!!!!!!!!!!!!!!!!
			std::complex<double>* in = this->set_ptr_in();
			size_t N = get_N_x();
			double min = get_x_min();
			double dx = get_dx();
			double x, y;
			double a = 1.0;
			double b = a * b_over_a;
			for (size_t i = 0; i < N; ++i) {
				for (size_t j = 0; j < N; ++j) {
					x = min + i*dx;
					y = min + j*dx;
					in[i*N + j] = std::complex<double>(exp(-(x*x/a/a + y*y/b/b) / 2), 0.);
				}
			}
			return *this;
		}
		Beam& initial_2_max_points(double r_0_over_r_1 = 2.0) { // working just with N_x = N_y !!!!!!!!!!!!!!!!!!!!!!
			double x, y, phi;
			double r_0 = get_r_x_0();
			double r_1 = r_0 / r_0_over_r_1;
			std::complex<double>* in = set_ptr_in();
			size_t N = get_N_x();
			for (size_t i = 0; i < N; ++i) {
				for (size_t j = 0; j < N; ++j) {
					x = coord_x(i);
					y = coord_y(j);
					phi = atan2(y, x) * (1.0 /*+ bp_.noise_amplitude()*sin(bp_.noise_theta()*hypot(x,y))*/);
					in[i*N + j] = std::pow(x*x + y*y, get_m()/2.0) * (std::exp(-(x*x + (y-1)*(y-1))/2/r_1/r_1) + std::exp(-(x*x + (y+1)*(y+1))/2/r_1/r_1)) * exp(std::complex<double>(0, get_m()*phi));
				}
			}
			return *this;
		}
		Beam& initial_1_max_point(double r_0_over_r_1 = 1.5, double ang_0 = 0.0, double part_ampl_0 = 0.3) {
															// r_0_over_r_1 - how many times A_noise is more than A_0 
															// ang_0 - angle between OX in positive direction
															// part_ampl_0 - how many times A_noise is less than A_0
			double x, y, phi, x_g, y_g;
			double r_0 = get_r_x_0();
			double r_g = r_0 / r_0_over_r_1;
			std::complex<double>* in = set_ptr_in();
			size_t N = get_N_x();
			for (size_t i = 0; i < N; ++i) {
				for (size_t j = 0; j < N; ++j) {
					
					x = coord_x(i);
					y = coord_y(j);
					x_g = x - r_0 * cos(ang_0);
					y_g = y - r_0 * sin(ang_0);
					phi = atan2(y, x) * (1.0 /*+ bp_.noise_amplitude()*sin(bp_.noise_theta()*hypot(x,y))*/);
					in[i*N + j] = (std::pow(x*x + y*y, get_m()/2.0) * std::exp(-(x*x + y*y) / 2 / r_0 / r_0) +
						       part_ampl_0 * std::exp(-(x_g*x_g + y_g*y_g) / 2 / r_g / r_g))
						      * exp(std::complex<double>(0, get_m()*phi));
				}
			}
			return *this;
		}
		Beam& initial_picture(std::string file_name, double N_x, double N_y) {	//			WRITE!!!!!!
			return *this;
		}
		Beam& initial_zero() { // working just with N_x = N_y !!!!!!!!!!!!!!!!!!!!!!
			std::complex<double>* in = set_ptr_in();
			size_t N = get_N_x();
			for (size_t i = 0; i < N; ++i) {
				for (size_t j = 0; j < N; ++j) {
					in[i*N + j] = std::complex<double>(0.0,0.0);
				}
			}
			return *this;
		}
		Beam& add_noise_phasescreen(std::string file_name, double sigma = 1.0) {
			std::ifstream input(file_name);
			if (!input.is_open()) {
				std::cout << "Input file " << file_name << " with noise parameters can't be opened!" << std::endl;
				exit(7); 
			}
			int N = (int)get_N_x();
			int grid_size = (int)N*N*2;
			double *noise_arr = new double(grid_size);
			int counter = 0;
			double input_val;
			//while (getline(input, input_val) && counter < grid_size) {
			//	noise_arr[counter] = input_val;
			//	counter++;
			//}
			while (!input.eof() && counter < grid_size) {
				input >> input_val;
				noise_arr[counter] = input_val;
//std::cout << counter << '\t' << input_val << '\n';
				counter++;
			}
			if (counter != grid_size) {
				std::cout << "There are not so much values in " << file_name << "! Please, check the file!" << std::endl;
				exit(8);
			}
			std::complex<double>* in = set_ptr_in();
			for (int i = 0; i < N; ++i) {
				for (int j = 0; j < N; ++j) {
					in[i*N + j] = std::complex<double>(in[i*N + j].real() + noise_arr[(i*N + j) * 2]*sigma, in[i*N + j].imag() + noise_arr[(i*N + j) * 2 + 1]*sigma);
				}
			}
			input.close();
			return *this;
		}
		Beam& step_diffraction(size_t num_of_calculating_steps_in_1_Z_diff, double z_finish) {
			double dz = get_z_diffraction() / static_cast<double>(num_of_calculating_steps_in_1_Z_diff);
			for (double z_now = 0.0; z_now <= z_finish * get_z_diffraction(); z_now += dz)            // NUM OF Z_DIFFR (LENGH)
				propagation_in_1_step(dz);
			return *this;
		}
		Beam& phase_plate() { // working just with N_x = N_y !!!!!!!!!!!!!!!!!!!!!!
			std::complex<double>* in = set_ptr_in();
			size_t N = get_N_x();
			for (size_t i = 0; i < N; ++i) {												 // !!!!! put it in into the Grid_constructor
				for (size_t j = 0; j < N; ++j)
					in[i*N + j] = phase_plate_1_point(in[i*N + j], (int)i, (int)j);
			}
			return *this;
		}
		Beam& double_phase_plate(double num_r_new_center = 1.0, char type_circle_1 = 0, char type_circle_2 = 0, char cut_1 = 0, char cut_2 = 0) { // working just with N_x = N_y !!!!!!!!!! and 
																				//cut_1 and cut_2 are not working
			std::complex<double>* in = set_ptr_in();
			int N = (int)get_N_x();
			double r_0 = get_r_x_0();
			double l = num_r_new_center * r_0;
			int N_one_side = floor(double(N)*l/(get_x_max()-get_x_min())/2.0);
			int N_half = floor(double(N)/2.0);
			for (int i = N_half - 2*N_one_side; i < N_half + 2*N_one_side + 1; ++i) {						 // !!!!! put it in into the Grid_constructor
				for (int j = N_half - 2*N_one_side; j < N_half + 2*N_one_side + 1; ++j) {
					in[i*N + j] = double_phase_plate_1_point(in[i*N + j], i, j, l, type_circle_1, type_circle_2);
				}
			}
			return *this;
		}
		double compare_angles_between_maxs (double ang_before) { 		// working just with N_x = N_y !!!!!!!!!!!!!!!!!!!!!!
			size_t n_max =  find_max(); // write it in private
			int j = (int)n_max % (int)get_N_x();
			int i = (int)n_max / (int)get_N_x();
			double x_new = coord_x(i);
			double y_new = coord_y(j);
			double ang_new = atan2(y_new, x_new) + M_PI;
			return std::abs (ang_new - ang_before);
		}
		int find_num_of_max_on_circle() {
			int num = 0;
			double* A_of_alpha;
			int points = 360;
			double angle_beg = 0.0;
			receive_A_of_alpha(A_of_alpha, (size_t)points, angle_beg);
			// count num of maximums
			return num;
		}
		Beam& print_amplitude(std::string str = "output_ampl", int num = -1) { 		// working just with N_x = N_y !!!!!!!!!!!!!!!!!!!!!!
			str = "out/3D_" + str;
			if (num != -1)
				str += to_string(num);
			str += ".txt";
			std::ofstream fout(str);  // f.write
			if (!fout.is_open()) { // rewrite as exception
				std::cout << "Error in opening output file!";
				exit(3);
			}
			print_in_work_info(str);
			const std::complex<double>* in = set_ptr_in();
			size_t N = get_N_x();
			double res;
			for (size_t i = 0; i < N; ++i) {
				for (size_t j = 0; j < N; ++j) {
					res = sqrt(in[i*N + j].real()*in[i*N + j].real() + in[i*N + j].imag()*in[i*N + j].imag());
					fout << res << '\t' << coord_x(i) << '\t' << coord_y(j) << std::endl;
				}
			}
			fout.close();
			return *this;
		}
		Beam& print_amplitude_center(double part_of_all_grid_in_1_D = 0.125, std::string str = "output_ampl", int num = -1) { // working just with N_x = N_y !!!!!!!!!!!!!!!!!!!!!!
			str = "out/3D_" + str;
			if (part_of_all_grid_in_1_D > 1.0) { // rewrite as exception
				std::cout << "Part can't be more than 1.0!";
				exit(4);
			}
			if (num != -1)
				str += to_string(num);
			str += ".txt";
			std::ofstream fout(str);  // f.write
			if (!fout.is_open()) { // rewrite as exception
				std::cout << "Error in opening output file!";
				exit(3);
			}
			print_in_work_info(str);
			const std::complex<double>* in = set_ptr_in();
			size_t N = get_N_x();
			size_t N_close = floor(static_cast<double>(N) * static_cast<double>(part_of_all_grid_in_1_D) / 2.0);
			double res;
			for (size_t i = floor(N/2) - N_close; i < floor(N/2) + N_close; ++i) {
				for (size_t j = floor(N/2) - N_close; j < floor(N/2) + N_close; ++j) {
					res = sqrt(in[i*N + j].real()*in[i*N + j].real() + in[i*N + j].imag()*in[i*N + j].imag());
					fout << res << '\t' << coord_x(i) << '\t' << coord_y(j) << std::endl;
				}
			}
			fout.close();
			return *this;
		}
		Beam& print_phase(std::string str = "output_phas", int num = -1) { // working just with N_x = N_y !!!!!!!!!!!!!!!!!!!!!!
			str = "out/3D_" + str;
			if (num != -1)
				str += to_string(num);
			str += ".txt";
			std::ofstream fout(str);  // f.write
			if (!fout.is_open()) { // rewrite as exception
				std::cout << "Error in opening output file!";
				exit(3);
			}
			print_in_work_info(str);
			const std::complex<double>* in = set_ptr_in();
			size_t N = get_N_x();
			double res;
			for (size_t i = 0; i < N; ++i) {
				for (size_t j = 0; j < N; ++j) {
					res = arg(in[i*N + j]);
					fout << res << '\t' << coord_x(i) << '\t' << coord_y(j) << '\n';
				}
			}
			return *this;
		}
		Beam& print_phase_center(double part_of_all_grid_in_1_D = 0.125, std::string str = "output_phas", int num = -1) { // working just with N_x = N_y !!!!!!!!!!!!!!!!!!!!!!
			str = "out/3D_" + str;
			if (num != -1)
				str += to_string(num);
			str += ".txt";
			std::ofstream fout(str);  // f.write
			if (!fout.is_open()) { // rewrite as exception
				std::cout << "Error in opening output file!";
				exit(3);
			}
			print_in_work_info(str);
			const std::complex<double>* in = set_ptr_in();
			size_t N = get_N_x();
			size_t N_close = floor(static_cast<double>(N) * static_cast<double>(part_of_all_grid_in_1_D) / 2.0);
			double res;
			for (size_t i = floor(N/2) - N_close; i < floor(N/2) + N_close; ++i) {
				for (size_t j = floor(N/2) - N_close; j < floor(N/2) + N_close; ++j) {
					res = arg(in[i*N + j]);
					fout << res << '\t' << coord_x(i) << '\t' << coord_y(j) << '\n';
				}
			}
			return *this;
		}
		Beam& print_spectrum_x(std::string str = "output_spec_x", int num = -1) { // working just with N_x = N_y !!!!!!!!!!!!!!!!!!!!!!
			str = "out/2D_" + str;
			if (num != -1)
				str += to_string(num);
			str += ".txt";
			std::ofstream fout(str);  // f.write
			if (!fout.is_open()) { // rewrite as exception
				std::cout << "Error in opening output file!";
				exit(3);
			}
			print_in_work_info(str);
			fftw_execute_forward();
			std::complex<double>* in = set_ptr_in();
			const std::complex<double>* out = set_ptr_out();
			size_t N = get_N_x();
			double L = abs(get_x_max() - get_x_min());
			double dk_x = 2.0 * M_PI / L;
			double k_x_min = - dk_x * (N/2 - 1);
			double res;
			for (size_t i = 0; i < N; ++i) {
				res = out[N/2*N + i].real()*out[N/2*N + i].real() + out[N/2*N + i].imag()*out[N/2*N + i].imag();
				fout << res << '\t' << k_x_min + i*dk_x << '\n';
			}
			fftw_execute_backward();
			for (size_t i = 0; i < get_num_points(); ++i)
				in[i] /= get_num_points();
			fout.close();
			return *this;
		}
		Beam& print_spectrum(std::string str = "output_spec", int num = -1) {  // working just with N_x = N_y !!!!!!!!!!!!!!!!!!!!!!	
			str = "out/3D_" + str;
			if (num != -1)
				str += to_string(num);
			str += ".txt";
			std::ofstream fout(str);  // f.write
			if (!fout.is_open()) { // rewrite as exception
				std::cout << "Error in opening output file!";
				exit(3);
			}
			print_in_work_info(str);
			fftw_execute_forward();
			std::complex<double>* in = set_ptr_in();
			const std::complex<double>* out = set_ptr_out();
			size_t N = get_N_x();
			double res;
			for (size_t i = 0; i < N; ++i) {
				for (size_t j = 0; j < N; ++j) {
					res = abs(out[i*N + j]);
					fout << res << '\t' << spectr_coord_shift_x(i) << '\t' << spectr_coord_shift_y(j) << '\n';
				}
			}
			fftw_execute_backward();
			for (size_t i = 0; i < get_num_points(); ++i)
				in[i] /= get_num_points();
			fout.close();
			return *this;
		}		
		Beam& print_spectrum_shifted() { //				WRITE!!!!
			return *this;
		}
		~Beam() {}
	};




}
}


#endif // BEAM_H_INCLUDED
