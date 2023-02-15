#ifndef BEAM_H_INCLUDED
#define BEAM_H_INCLUDED

/*===================================================================
	Project:	Magister's diploma
	
	Author:		F. Soifer
	
	Date:		2023-02-15

	Description:	Class Grid_parameters (Beam preparings)

	Comments:	!!!!!!!!!!!!!!!!!!!!!!!

  ===================================================================*/


#include "grid_param.hpp"
#include "beam_param.hpp"
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
			double real = obj.real(), imag = obj.imag();
			obj = std::complex<double> (real * cos(fi) - imag * sin(fi), real * sin(fi) + imag * cos(fi));
		}
		Beam& propagation_in_1_step(double dz) {			// working just with N_x = N_y !!!!!!!!!!!!!!!!!!!!!! // rewrite more optimized
			std::complex<double>* in = set_ptr_in();
			std::complex<double>* out = set_ptr_out();
			size_t N = get_N_x();
			double dk = 2.0 * M_PI / (N * get_dx());
			fftw_execute_forward();
			for (size_t i = 0; i < N; ++i) {
				for (size_t j = 0; j < N; ++j) {
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
					in[i*N + j] = std::exp(-(x*x + y*y) / 2 / std::pow(get_r_x_0(), 2)) * exp(std::complex<double>(0, get_m()*phi));
				}
			}
			return *this;
		}
		Beam& step_diffraction(size_t num_of_calculating_steps_in_1_Z_diff, double z_finish) {
			double dz = get_z_diffraction() / static_cast<double>(num_of_calculating_steps_in_1_Z_diff);
			for (double z_now = 0.0; z_now <= z_finish * get_z_diffraction(); z_now += dz)            // NUM OF Z_DIFFR (LENGH)
				propagation_in_1_step(dz);
			return *this;
		}
		Beam& phase_plate() { //				WRITE!!!!
			return *this;
		}
		Beam& print_amplitude(std::string str = "output_ampl", int num = -1) { 		// working just with N_x = N_y !!!!!!!!!!!!!!!!!!!!!!
			if (num != -1)
				str += to_string(num);
			str += ".txt";
			std::ofstream fout(str);  // f.write
			if (!fout.is_open()) { // rewrite as exception
				std::cout << "Error in opening output file!";
				exit(3);
			}
			const std::complex<double>* in = set_ptr_in();
			size_t N = get_N_x();
			double res;
			for (size_t i = 0; i < N; ++i) {
				for (size_t j = 0; j < N; ++j) {
					res = sqrt(in[i*N + j].real()*in[i*N + j].real() + in[i*N + j].imag()*in[i*N + j].imag());
					fout << res << '\t' << get_x_min() + i*get_dx() << '\t' << get_x_min() + j*get_dx() << std::endl;
				}
			}
			fout.close();
			return *this;
		}
		Beam& print_amplitude_center(double part_of_all_grid_in_1_D = 0.125, std::string str = "output_ampl", int num = -1) { // working just with N_x = N_y !!!!!!!!!!!!!!!!!!!!!!
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
			const std::complex<double>* in = set_ptr_in();
			size_t N = get_N_x();
			size_t N_close = floor(static_cast<double>(N) * static_cast<double>(part_of_all_grid_in_1_D) / 2.0);
			double res;
			for (size_t i = floor(N/2) - N_close; i < floor(N/2) + N_close; ++i) {
				for (size_t j = floor(N/2) - N_close; j < floor(N/2) + N_close; ++j) {
					res = sqrt(in[i*N + j].real()*in[i*N + j].real() + in[i*N + j].imag()*in[i*N + j].imag());
					fout << res << '\t' << get_x_min() + i*get_dx() << '\t' << get_x_min() + j*get_dx() << std::endl;
				}
			}
			fout.close();
			return *this;
		}
		Beam& print_phase(std::string str = "output_phas", int num = -1) { // working just with N_x = N_y !!!!!!!!!!!!!!!!!!!!!!
			if (num != -1)
				str += to_string(num);
			str += ".txt";
			std::ofstream fout(str);  // f.write
			if (!fout.is_open()) { // rewrite as exception
				std::cout << "Error in opening output file!";
				exit(3);
			}
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
		Beam& print_spectrum(std::string str = "output_spec", int num = -1) { // working just with N_x = N_y !!!!!!!!!!!!!!!!!!!!!!
			if (num != -1)
				str += to_string(num);
			str += ".txt";
			std::ofstream fout(str);  // f.write
			if (!fout.is_open()) { // rewrite as exception
				std::cout << "Error in opening output file!";
				exit(3);
			}
			fftw_execute_forward();
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
