#ifndef NOISE_PARAM_H_INCLUDED
#define NOISE_PARAM_H_INCLUDED

/*===================================================================
	Project:	Magister's diploma
	
	Author:		F. Soifer
	
	Date:		2023-02-12

	Description:	Class Grid_parameters (Beam preparings)

	Comments:	!!!!!!!!!!!!!!!!!!!!!!!

  ===================================================================*/
  
  
namespace ofvp {
namespace create_vortex{

	enum Noise_type {
		JUST_AMPL = 1,
		JUST_PHAS,
		AMPL_PHAS
	};
	
	enum Noise_math_type {
		DELTA_CORR = 1,
		DELTA_CORR_BIG_GRID,
		WITH_RAD_COR
	};

	struct beam_noise_parameters {
		Noise_type noise_type;
		Noise_math_type noise_math_type;
		double	noise_ampl_perc;
		double	noise_ampl_abs;
		double	noise_phas_perc;
		double	noise_phas_abs;
	};




}
}


#endif // NOISE_PARAM_H_INCLUDED
