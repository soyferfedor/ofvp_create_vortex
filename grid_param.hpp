#ifndef GRID_PARAM_H_INCLUDED
#define GRID_PARAM_H_INCLUDED

/*===================================================================
	Project:	Magister's diploma
	
	Author:		F. Soifer
	
	Date:		2023-02-27

	Description:	Class Grid_parameters (Beam preparings)

	Comments:	!!!!!!!!!!!!!!!!!!!!!!!

  ===================================================================*/


#include "fftw3.h"
#include <cmath>
#include <complex>
#include <iostream>


namespace ofvp {
namespace create_vortex{

// Structure FFTW_parameters

// Include:	2 fftw plans from fftw3.h


	struct FFTW_parameters {
		fftw_plan plan_forward; 
		fftw_plan plan_backward;
	};


// Class Grid_arrays

// Include:	Num of points on both axes, min and
//		max values on every axe, distance
//		between points and 2 pointers
//		(before and after fftw execution)


	class Grid_arrays {
		size_t N_x;
		size_t N_y;
		double dx;
		double dy;
		double x_min;
		double x_max;
		double y_min;
		double y_max;
		std::complex<double> *in;  // unique_ptr ??
		std::complex<double> *out;
	public:
		Grid_arrays(size_t N_x, double x_min, double x_max, size_t N_y, double y_min, double y_max):
			N_x(N_x), x_min(x_min), x_max(x_max), N_y(N_y), y_min(y_min), y_max(y_max)
		{
			dx = (x_max-x_min)/(N_x-1);
			dy = (y_max-y_min)/(N_y-1);
			try {
				in = new std::complex<double> [N_x*N_y];  // in = new (nothrow) complex<double> [N*N]
				out = new std::complex<double> [N_x*N_y];
			}
			catch (std::bad_alloc &e) {
				std::cout << e.what() << std::endl;
				exit(1);
			}
		}
		Grid_arrays(const Grid_arrays&) =delete;
		Grid_arrays(Grid_arrays&&) =delete;  // move ctor
		Grid_arrays& operator=(const Grid_arrays&) =delete;
		Grid_arrays& operator=(Grid_arrays&&) =delete;
		const std::complex<double>* get_ptr_in() const {
			return in;
		}
		const std::complex<double>* get_ptr_out() const {
			return out;
		}
		std::complex<double>* set_ptr_in() {
			return in;
		}
		std::complex<double>* set_ptr_out() {
			return out;
		}
		size_t get_N_x() const {
			return N_x;
		}
		size_t get_N_y() const {
			return N_y;
		}
		size_t get_num_points() const {
			return N_x*N_y;
		}
		double get_x_min() const {
			return x_min;
		}
		double get_x_max() const {
			return x_max;
		}
		double get_y_min() const {
			return y_min;
		}
		double get_y_max() const {
			return y_max;
		}
		double get_dx() const {
			return dx;
		}
		double get_dy() const {
			return dy;
		}
		double coord_x(size_t idx) const {
			return x_min + dx * idx;
		}
		double coord_y(size_t idx) const {
			return y_min + dy * idx;
		}
		double spectr_coord_x(size_t idx) const {
			double dk = 2.0 * M_PI / (get_x_max() - get_x_min());
			double k_min = - dk * (get_N_x()/2 - 1);
			return k_min + idx * dk;
		}
		double spectr_coord_shift_x(size_t idx) const {
			double dk = 2.0 * M_PI / (get_x_max() - get_x_min());
			if (idx < get_N_x()/2)
				return idx * dk;
			else
				return - dk * (get_N_x() - idx);
		}
		double spectr_coord_y(size_t idx) const {
			double dk = 2.0 * M_PI / (get_y_max() - get_y_min());
			double k_min = - dk * (get_N_y()/2 - 1);
			return k_min + idx * dk;
		}
		double spectr_coord_shift_y(size_t idx) const {
			double dk = 2.0 * M_PI / (get_y_max() - get_y_min());
			if (idx < get_N_y()/2)
				return idx * dk;
			else
				return - dk * (get_N_y() - idx);
		}
		~Grid_arrays() {// dtor
			delete[] in;
			delete[] out;
		}
	};


// Class Grid_parameters

// Include:	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

// Notes:	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	class Grid_parameters: public  Grid_arrays{
		FFTW_parameters parameters;
	public:
		Grid_parameters (size_t N_x, double x_min, double x_max, size_t N_y, double y_min, double y_max):
			Grid_arrays(N_x, x_min, x_max, N_y, y_min, y_max)
		{
			parameters.plan_forward = fftw_plan_dft_2d(N_x, N_y,
							reinterpret_cast<fftw_complex*>(set_ptr_in()),
							reinterpret_cast<fftw_complex*>(set_ptr_out()),
							FFTW_FORWARD, FFTW_ESTIMATE);
			parameters.plan_backward = fftw_plan_dft_2d(N_x, N_y,
							reinterpret_cast<fftw_complex*>(set_ptr_out()),
							reinterpret_cast<fftw_complex*>(set_ptr_in()),
							FFTW_BACKWARD, FFTW_ESTIMATE);
		}
		Grid_parameters(const Grid_parameters&) =delete;
		Grid_parameters(Grid_parameters&&) =delete;  // move ctor
		Grid_parameters& operator=(const Grid_parameters&) =delete;
		Grid_parameters& operator=(Grid_parameters&&) =delete;

		void fftw_execute_forward() {
			fftw_execute(parameters.plan_forward);
		}
		void fftw_execute_backward() {
			fftw_execute(parameters.plan_backward);
		}
		~Grid_parameters() {
			fftw_destroy_plan(parameters.plan_forward);
			fftw_destroy_plan(parameters.plan_backward);
		}	
	};


}
}


#endif // GRID_PARAM_H_INCLUDED
