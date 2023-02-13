#ifndef BEAM_H_INCLUDED
#define BEAM_H_INCLUDED

/*===================================================================
	Project:	Magister's diploma
	
	Author:		F. Soifer
	
	Date:		2023-02-12

	Description:	Class Grid_parameters (Beam preparings)

	Comments:	!!!!!!!!!!!!!!!!!!!!!!!

  ===================================================================*/


#include "grid_param.h"
#include "beam_param.h"


namespace ofvp {
namespace create_vortex{


	enum Initial_type {gauss = 1, vortex};


	class Beam_type {
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


	enum Location_type {center = 1, diag, circle, square};
	
	
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

}
}


#endif // BEAM_H_INCLUDED
