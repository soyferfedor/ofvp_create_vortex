#ifndef WORK_INFO_H_INCLUDED
#define WORK_INFO_H_INCLUDED

/*===================================================================
	Project:	Magister's diploma
	
	Author:		F. Soifer
	
	Date:		2023-02-12

	Description:	Working info (for .py environment)

	Comments:	!!!!!!!!!!!!!!!!!!!!!!!

  ===================================================================*/


#include <iostream>
#include <fstream>


namespace ofvp {
namespace create_vortex{


	void print_in_work_info(std::string filename) {
		static bool flag_open = false;
		if (flag_open == false) {
			std::ofstream fout("out/work_info.txt");  // f.write
			if (!fout.is_open()) { // rewrite as exception
				std::cout << "C++: Error in opening work_info.txt file!";
				exit(51);
			}
			fout << filename << '\n';
			fout.close();
		}
		else {
			std::ofstream fout{ "out/work_info.txt", std::ios::app };
			if (!fout.is_open()) { // rewrite as exception
				std::cout << "C++: Error in opening work_info.txt file!";
				exit(51);
			}
			fout << filename << '\n';
			fout.close();
		}
		flag_open = true;
	}
  
}
}


#endif // WORK_INFO_H_INCLUDED
