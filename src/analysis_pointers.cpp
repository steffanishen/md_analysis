/*
    read_dcd : c++ class + main file example for reading a CHARMM dcd file
    Copyright (C) 2013  Florent Hedin
    
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

//***************** Partially contributed by Meng Shen since 06/10/2018 ***************************
//***************** Email: steffanishen@gmail.com ***************************************


#include <cstdlib>
#include <cstring>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>

#include <math.h>
#include "analysis_pointers.hpp"


#define PI 3.14159265

using namespace std;

ANALYSIS_POINTERS::ANALYSIS_POINTERS(PSF *system, vector<GROUP*> sels, INPUT* input, vector<PSF*> monomers)
{
    this->system = system;

    for (unsigned i = 0; i < sels.size(); i++) {
        this->sels.push_back(sels[i]);
    }
    this->input = input;

    for (unsigned i = 0; i < monomers.size(); i++) {
        this->monomers.push_back(monomers[i]);
    }
}

vector<ANALYSIS*> ANALYSIS_POINTERS::init() {
    for (auto &analysis_opt : input->analysis) {
        int vector1d;
 	int vector2d;
 	int voidf;
	string filename;
	if (analysis_opt[0] == "anglezs_rings") {
		int groupid;
		for (int argid = 1; argid < analysis_opt.size(); argid++) {
		    if (analysis_opt[argid] == "group") {
			groupid = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "vector1d") {
			vector1d = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "vector2d") {
			vector2d = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "voidf") {
			voidf = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "filename") {
			filename = analysis_opt[argid+1];
		    }  
		}
		analysis.push_back(new ANALYSIS_ANGLEZS_RINGS(system,sels[groupid],vector1d,vector2d,voidf,filename));

        } else if (analysis_opt[0] == "mindangles_seg_dislocated") {
		int groupid;
	 	int whichN;
		for (int argid = 1; argid < analysis_opt.size(); argid++) {
		    if (analysis_opt[argid] == "group") {
			groupid = stoi(analysis_opt[argid+1]);
		    } else if (analysis_opt[argid] == "whichN") {
			whichN = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "vector1d") {
			vector1d = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "vector2d") {
			vector2d = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "voidf") {
			voidf = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "filename") {
			filename = analysis_opt[argid+1];
		    }
		}
        	analysis.push_back(new ANALYSIS_MINDANGLES_SEG_DISLOCATED(system,sels[groupid],whichN,vector1d,vector2d,voidf,filename));
        } else if (analysis_opt[0] == "scale") {
		int groupid;
		float scaling;
		for (int argid = 1; argid < analysis_opt.size(); argid++) {
		    if (analysis_opt[argid] == "group") {
			groupid = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "vector1d") {
			vector1d = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "vector2d") {
			vector2d = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "voidf") {
			voidf = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "filename") {
			filename = analysis_opt[argid+1];
		    }  else if (analysis_opt[argid] == "scaling") {
			scaling = stoi(analysis_opt[argid+1]);
		    }  
		}
	        analysis.push_back(new ANALYSIS_SCALE(system,sels[groupid],vector1d,vector2d,voidf,filename,scaling));
        } else if (analysis_opt[0] == "mindangles_seg") {
		int groupid;
	 	int whichN;
		for (int argid = 1; argid < analysis_opt.size(); argid++) {
		    if (analysis_opt[argid] == "group") {
			groupid = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "vector1d") {
			vector1d = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "vector2d") {
			vector2d = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "voidf") {
			voidf = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "filename") {
			filename = analysis_opt[argid+1];
		    }
		}
        	analysis.push_back(new ANALYSIS_MINDANGLES_SEG(system,sels[groupid],vector1d,vector2d,voidf,filename));
        } else if (analysis_opt[0] == "avedangles_seg") {
		int groupid;
	 	int whichN;
		for (int argid = 1; argid < analysis_opt.size(); argid++) {
		    if (analysis_opt[argid] == "group") {
			groupid = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "vector1d") {
			vector1d = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "vector2d") {
			vector2d = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "voidf") {
			voidf = stoi(analysis_opt[argid+1]);
		    }  else if (analysis_opt[argid] == "filename") {
			filename = analysis_opt[argid+1];
		    }
		}
        	analysis.push_back(new ANALYSIS_AVEDANGLES_SEG(system,sels[groupid],vector1d,vector2d,voidf,filename));
        } else if (analysis_opt[0] == "random_walk") {
                int groupid;
	        string name1;
	        string name2;
	        string name3;
	        string name4;
	        string name_ref1;
	        string name_ref2;
		float fringe;
	        float bl;
                for (int argid = 1; argid < analysis_opt.size(); argid++) {
                    if (analysis_opt[argid] == "group") {
                        groupid = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "vector1d") {
                        vector1d = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "vector2d") {
                        vector2d = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "voidf") {
                        voidf = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "headname") {
                        name1 = analysis_opt[argid+1];
                    }  else if (analysis_opt[argid] == "headbondname") {
                        name2 = analysis_opt[argid+1];
                    }  else if (analysis_opt[argid] == "endname") {
                        name3 = analysis_opt[argid+1];
                    }  else if (analysis_opt[argid] == "endbondname") {
                        name4 = analysis_opt[argid+1];
                    }  else if (analysis_opt[argid] == "name_ref1") {
                        name_ref1 = analysis_opt[argid+1];
                    }  else if (analysis_opt[argid] == "name_ref2") {
                        name_ref2 = analysis_opt[argid+1];
                    }  else if (analysis_opt[argid] == "bl") {
                        bl = stof(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "filename") {
                        filename = analysis_opt[argid+1];
                    }  else if (analysis_opt[argid] == "fringe") {
                        fringe = stof(analysis_opt[argid+1]);
                    }
                }
        	analysis.push_back(new ANALYSIS_RANDOM_WALK(system,sels[groupid],vector1d,vector2d,voidf,bl,filename,monomers,fringe));
        } else if (analysis_opt[0] == "straight_chain") {
                int groupid;
	        string name1;
	        string name2;
	        string name3;
	        string name4;
	        string name_ref1;
	        string name_ref2;
		float fringe;
	        float bl;
                for (int argid = 1; argid < analysis_opt.size(); argid++) {
                    if (analysis_opt[argid] == "group") {
                        groupid = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "vector1d") {
                        vector1d = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "vector2d") {
                        vector2d = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "voidf") {
                        voidf = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "headname") {
                        name1 = analysis_opt[argid+1];
                    }  else if (analysis_opt[argid] == "headbondname") {
                        name2 = analysis_opt[argid+1];
                    }  else if (analysis_opt[argid] == "endname") {
                        name3 = analysis_opt[argid+1];
                    }  else if (analysis_opt[argid] == "endbondname") {
                        name4 = analysis_opt[argid+1];
                    }  else if (analysis_opt[argid] == "name_ref1") {
                        name_ref1 = analysis_opt[argid+1];
                    }  else if (analysis_opt[argid] == "name_ref2") {
                        name_ref2 = analysis_opt[argid+1];
                    }  else if (analysis_opt[argid] == "bl") {
                        bl = stof(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "filename") {
                        filename = analysis_opt[argid+1];
                    }  else if (analysis_opt[argid] == "fringe") {
                        fringe = stof(analysis_opt[argid+1]);
                    }
                }
        	analysis.push_back(new ANALYSIS_STRAIGHT_CHAIN(system,sels[groupid],vector1d,vector2d,voidf,bl,filename,monomers,fringe));
        } else if (analysis_opt[0] == "packing") {
                int groupid;
	        string name1;
	        string name2;
	        string name3;
	        string name4;
	        string name_ref1;
	        string name_ref2;
	        float lx;
	        float ly;
	        float lz;
                for (int argid = 1; argid < analysis_opt.size(); argid++) {
                    if (analysis_opt[argid] == "group") {
                        groupid = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "vector1d") {
                        vector1d = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "vector2d") {
                        vector2d = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "voidf") {
                        voidf = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "filename") {
                        filename = analysis_opt[argid+1];
                    }  else if (analysis_opt[argid] == "lx") {
                        lx = stof(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "ly") {
                        ly = stof(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "lz") {
                        lz = stof(analysis_opt[argid+1]);
                    }
                }
        	analysis.push_back(new ANALYSIS_PACKING(system,sels[groupid],vector1d,vector2d,voidf,filename,lx,ly,lz));
        } else if (analysis_opt[0] == "rotate") {
                int groupid;
	        string name1;
	        string name2;
	        string name3;
	        string name4;
	        string name_ref1;
	        string name_ref2;
	        float bl;
		int nbinsangle;
		int axisid;
                for (int argid = 1; argid < analysis_opt.size(); argid++) {
                    if (analysis_opt[argid] == "group") {
                        groupid = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "vector1d") {
                        vector1d = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "vector2d") {
                        vector2d = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "voidf") {
                        voidf = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "headname") {
                        name1 = analysis_opt[argid+1];
                    }  else if (analysis_opt[argid] == "headbondname") {
                        name2 = analysis_opt[argid+1];
                    }  else if (analysis_opt[argid] == "endname") {
                        name3 = analysis_opt[argid+1];
                    }  else if (analysis_opt[argid] == "endbondname") {
                        name4 = analysis_opt[argid+1];
                    }  else if (analysis_opt[argid] == "name_ref1") {
                        name_ref1 = analysis_opt[argid+1];
                    }  else if (analysis_opt[argid] == "name_ref2") {
                        name_ref2 = analysis_opt[argid+1];
                    }  else if (analysis_opt[argid] == "bl") {
                        bl = stof(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "filename") {
                        filename = analysis_opt[argid+1];
                    }  else if (analysis_opt[argid] == "nbinsangle") {
                        nbinsangle = stoi(analysis_opt[argid+1]);
                    }  else if (analysis_opt[argid] == "axisid") {
                        axisid = stoi(analysis_opt[argid+1]);
                    }
                }
        	analysis.push_back(new ANALYSIS_ROTATE(system,sels[groupid],vector1d,vector2d,voidf,bl,name1,name2,name3,name4,name_ref1,name_ref2,filename,nbinsangle,axisid));
        } else {
	   error1.error_exit("Wrong analysis option, please check the input script!") ;
        }
    }
    return analysis;
}


ANALYSIS_POINTERS::~ANALYSIS_POINTERS()
{
    for (unsigned i = 0; i < analysis.size(); i++) {
	analysis[i] = NULL;
    }
    analysis.clear() ;
    
    for (unsigned i = 0; i < sels.size(); i++) {
        sels[i] = NULL;
    }
    sels.clear();

    system = NULL;
    input = NULL;

}

