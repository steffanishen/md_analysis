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

//***************** Fully contributed by Meng Shen since 06/10/2018 ***************************
//***************** Email: steffanishen@gmail.com ***************************************
//****************** Note the x, y, z and pbc need to be carefully considered ******************************************

#include <cstdlib>
#include <cstring>

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "psf.hpp"

using namespace std;

PSF::PSF(string filename)
{
    
    psff.exceptions(std::ifstream::failbit);
    try
    {
        psff.open(filename,ios::in);
    }
    catch(std::ifstream::failure e)
    {
        cerr << "Exception opening/reading psf file '" << filename << "' : " << std::endl;
        cerr << "Please check the path of the psf file and if it exists." << endl;
    } 
    
    psf_first_read=true;
}

void PSF::alloc()
{
    atom_index = new int[NATOM];
    segid = new int[NATOM];
    resid = new int[NATOM];
    resname = new string[NATOM];
    atomname = new string[NATOM];
    atomtype = new string[NATOM];
    charge = new float[NATOM];
    mass = new float[NATOM];
    beta = new int[NATOM];
    nbonds_atom = new int[NATOM];
    ibond_atom.resize(NATOM);
    x = new float[NATOM];
    y = new float[NATOM];
    z = new float[NATOM];
    for (unsigned i = 0; i < 6; i++) {
       pbc[i] = 0.0; 
    }

    //atoms->resize(NATOM);
}

void PSF::alloc_bonds()
{
   ibond = new int*[NBONDS];
   for (int i = 0; i < NBONDS; i++) {
       ibond[i] = new int[2];
   }
}

void PSF::alloc_angles()
{
   iangle = new int*[NANGLES];
   for (int i = 0; i < NANGLES; i++) {
       iangle[i] = new int[3];
   }
}

void PSF::alloc_dihedrals()
{
   idihedral = new int*[NDIHEDRALS];
   for (int i = 0; i < NDIHEDRALS; i++) {
       idihedral[i] = new int[4];
   }
}

void PSF::alloc_impropers()
{
   improper = new int*[NIMPROPERS];
   for (int i = 0; i < NIMPROPERS; i++) {
       improper[i] = new int[4];
   }
}


void PSF::read_header()
{
    string lineStr; 
    for (int i = 0; i < 2; i++) {
        getline(psff,lineStr);
    }

    getline (psff,lineStr);
    istringstream iss(lineStr);
    iss >> NTITLE;
  
    cout << "PSF NTITLE: " << NTITLE << endl; 
    for (int i = 0; i < NTITLE; i++) {
        string temp;
        getline(psff,lineStr);
        iss.str(lineStr);
        iss >> temp;
    }

    getline (psff,lineStr);

    getline (psff,lineStr);
    iss.str(lineStr);
    iss >> NATOM;
 
    alloc();
}

void PSF::read_atoms()
{
    cout << "PSF NATOM: " << NATOM << endl;
    string lineStr; 
    int segid_run = 0;
    int resid_run = 0;
    vector<int> residues_item;
    vector<vector<int>> segments_item;  

    for (int i = 0; i < NATOM; i++) {
        getline (psff,lineStr);
        istringstream iss(lineStr);
        iss >> atom_index[i]; // comment: this is just a placeholder to release iss.
        atom_index[i] = i;
        iss >> segid[i];
        iss >> resid[i];
        iss >> resname[i];
        iss >> atomname[i];
        iss >> atomtype[i];
        iss >> charge[i];
        iss >> mass[i];
        iss >> beta[i];
        nbonds_atom[i] = 0;
        ATOM atom;
        atom.atom_index = atom_index[i];
        atom.segid = segid[i];
        atom.resid = resid[i];
        atom.resname = resname[i];
        atom.atomname = atomname[i];
        atom.atomtype = atomtype[i];
        atom.charge = charge[i];
        atom.mass = mass[i];
        atom.beta = beta[i];
        atoms.push_back(atom);
	if (i == 0) {
            residues_item.push_back(atom.atom_index);
	    resid_run = resid[i];
	    segid_run = segid[i];
        } else if (i == NATOM - 1) {
  	    if (resid[i] != resid_run || segid[i] != segid_run) {
	        residues.push_back(residues_item);
	        segments_item.push_back(residues_item);
	        resid_run = resid[i];
		residues_item.clear();
  	    }
	    if (segid[i] != segid_run) {
	        segments.push_back(segments_item);
		segid_run = segid[i];
		segments_item.clear();
	    }
	    residues_item.push_back(atom.atom_index);
	    residues.push_back(residues_item);
	    segments_item.push_back(residues_item);
	    segments.push_back(segments_item);
	    residues_item.clear();
	    segments_item.clear();
        } else {
  	    if (resid[i] != resid_run || segid[i] != segid_run) {
	        residues.push_back(residues_item);
	        segments_item.push_back(residues_item);
	        resid_run = resid[i];
		residues_item.clear();
	    }
	    if (segid[i] != segid_run) {
	        segments.push_back(segments_item);
		segid_run = segid[i];
		segments_item.clear();
	    }
            residues_item.push_back(atom.atom_index);
        }
    }

    getline (psff,lineStr);
    getline (psff,lineStr);
    istringstream iss1(lineStr);
    iss1 >> NBONDS; 
    alloc_bonds();
}

void PSF::read_bonds()
{
    cout << "PSF NBONDS: " << NBONDS << endl;
    string lineStr; 
    int i = 0;
    while (i < NBONDS) {
      getline(psff,lineStr);
      istringstream iss(lineStr);
      int itemp;
      vector<int> intline;
      while (iss >> itemp ) {
          intline.push_back(itemp); // Store atom indices for bonds of each line in intline
      }
      int npairs = 0;
      int ncol = intline.size();
      npairs = ncol/2;
      int icol = 0;
      for (int j = 0; j < npairs; j++) {
          for (int indb = 0; indb < 2; indb++) {
              ibond[i][indb] = intline[icol] - 1;
              icol++;
          }
          i++;
      }
      //cout << endl;
    }


    for (i = 0; i < NBONDS; i++) {
            int iatom = ibond[i][0];
            int jatom = ibond[i][1];
            ibond_atom[iatom].push_back(jatom);
            ibond_atom[jatom].push_back(iatom);
            nbonds_atom[iatom] = nbonds_atom[iatom] + 1; 
            nbonds_atom[jatom] = nbonds_atom[jatom] + 1; 
    }

    getline (psff,lineStr);
    getline (psff,lineStr);
    istringstream iss(lineStr);
    iss >> NANGLES; 
    alloc_angles();
}

void PSF::read_angles()
{
    cout << "PSF NANGLES: " << NANGLES << endl;
    string lineStr; 
    int i = 0;
    while (i < NANGLES) {
      getline(psff,lineStr);
      istringstream iss(lineStr);
      int itemp;
      vector<int> intline;
      while (iss >> itemp ) {
          intline.push_back(itemp);
      }
      int ntrimers = 0;
      int ncol = intline.size();
      ntrimers = ncol/3;
      int icol = 0;
      for (int j = 0; j < ntrimers; j++) {
          for (int indb = 0; indb < 3; indb++) {
              iangle[i][indb] = intline[icol] - 1;
              icol++;
          }
          i++;
      }
      //cout << endl;
    }
    

    getline (psff,lineStr);
    getline (psff,lineStr);
    istringstream iss(lineStr);
    iss >> NDIHEDRALS; 
    alloc_dihedrals();
}

void PSF::read_dihedrals()
{
    cout << "PSF NDIHEDRALS: " << NDIHEDRALS << endl;
    string lineStr; 
    int i = 0;
    while (i < NDIHEDRALS) {
      getline(psff,lineStr);
      istringstream iss(lineStr);
      int itemp;
      vector<int> intline;
      while (iss >> itemp ) {
          intline.push_back(itemp);
      }
      int nquadruples = 0;
      int ncol = intline.size();
      nquadruples = ncol/4;
      int icol = 0;
      for (int j = 0; j < nquadruples; j++) {
          for (int indb = 0; indb < 4; indb++) {
              idihedral[i][indb] = intline[icol] - 1;
              icol++;
          }
          i++;
      }
      //cout << endl;
    }
    
//    for (int i = 0; i < NDIHEDRALS; i++) {
//        cout << i << " " << idihedral[i][0] << " " << idihedral[i][1] << " " << idihedral[i][2] << " " << idihedral[i][3] << endl;
//    }

    getline (psff,lineStr);
    getline (psff,lineStr);
    istringstream iss(lineStr);
    iss >> NIMPROPERS; 
    alloc_impropers();
}

void PSF::read_impropers()
{
    cout << "PSF NIMPROPERS: " << NIMPROPERS<< endl;
    string lineStr; 
    int i = 0;
    while (i < NIMPROPERS) {
      getline(psff,lineStr);
      istringstream iss(lineStr);
      int itemp;
      vector<int> intline;
      while (iss >> itemp ) {
          intline.push_back(itemp);
      }
      int nquadruples = 0;
      int ncol = intline.size();
      nquadruples = ncol/4;
      int icol = 0;
      for (int j = 0; j < nquadruples; j++) {
          for (int indb = 0; indb < 4; indb++) {
              improper[i][indb] = intline[icol] - 1;
              icol++;
          }
          i++;
      }
      //cout << endl;
    }
    
    for (int i = 0; i < NIMPROPERS; i++) {
        //cout << i << " " << improper[i][0] << " " << improper[i][1] << " " << improper[i][2] << " " << improper[i][3] << endl;
    }

}

void PSF::read_psf() {
    read_header();
    read_atoms();
    read_bonds();
    read_angles();
    read_dihedrals();
    read_impropers();
/*
    for (int i = 0; i < NATOM; i++) {
        //cout << i << " " << nbonds_atom[i] << endl; 
        cout << i << " " << endl;
        for (int j = 0; j < nbonds_atom[i]; j++) {
            cout << ibond_atom[i][j] << " ";
        }
        cout << endl;
    }
*/
}


GROUP PSF:: select_atoms(vector<string> str) {
    GROUP atoms_select;
//    vector<int> flagarray(NATOM,0);
    string temp;
    vector<string> arg1;
    vector<string> arithmetic_analog;
    int segid_run = -1;
    int resid_run = -1;
    vector<int> residues_item;
    vector<vector<int>> segments_item;  

    for (auto &atom : atoms) {
        int flag;
        flag = atom.select_atoms(str);
        if (flag == 1) atoms_select.atoms.push_back(atom);
    }

    atoms_select.NATOM = atoms_select.atoms.size();

    for (size_t i = 0; i < atoms_select.NATOM; i++ ) {
        ATOM atom = atoms_select.atoms[i];
	if (i == 0) {
            residues_item.push_back(atom.atom_index);
	    resid_run = resid[atom.atom_index];
	    segid_run = segid[atom.atom_index];
        } else if (i == size_t(atoms_select.NATOM -1)) {
  	    if (resid[atom.atom_index] != resid_run || segid[atom.atom_index] != segid_run) {
	        atoms_select.residues.push_back(residues_item);
	        segments_item.push_back(residues_item);
	        resid_run = resid[atom.atom_index];
		residues_item.clear();
	    }
	    if (segid[atom.atom_index] != segid_run) {
	        atoms_select.segments.push_back(segments_item);
		segid_run = segid[atom.atom_index];
		segments_item.clear();
	    }
	    residues_item.push_back(atom.atom_index);
	    atoms_select.residues.push_back(residues_item);
	    segments_item.push_back(residues_item);
	    atoms_select.segments.push_back(segments_item);
	    residues_item.clear();
	    segments_item.clear();
        } else {
  	    if (resid[atom.atom_index] != resid_run || segid[atom.atom_index] != segid_run) {
	        atoms_select.residues.push_back(residues_item);
	        segments_item.push_back(residues_item);
	        resid_run = resid[atom.atom_index];
		residues_item.clear();
	    }
	    if (segid[atom.atom_index] != segid_run) {
	        atoms_select.segments.push_back(segments_item);
		segid_run = segid[atom.atom_index];
		segments_item.clear();
	    }
            residues_item.push_back(atom.atom_index);
        }
//        cout << "debug: we monitor atom " << atom.atom_index <<" the size of segments_item " <<  segments_item.size() << endl; 
    }

    for (auto &segment : atoms_select.segments) {
	vector<int> segment_temp;
	for (auto &residue : segment) {
	    for (auto &atomid : residue) {
	        segment_temp.push_back(atomid);
	    }
	}
	atoms_select.segments_ind.push_back(segment_temp);
    }


    return atoms_select;
}

vector<float> PSF::vector_subtract(vector<float> v1, vector<float> v2) {
    vector<float> v3(v1.size());
    for (int i = 0; i < v1.size(); i++) {
	v3[i] = v1[i] - v2[i];
    }
    return v3;
}

void PSF::determine_names() {

    Vch.clear();
    Vcc.clear();
    Vcc_pre.clear();
    Vcc_next.clear();

    int ind_a1_special;
    string name_a1_special;

// determine the monomer type: 1st, body, or last
    if (nbonds_atom[0] == 1 && ( resid[0] != resid[1] || resname[0] != resname[1] )) {
        if (nbonds_atom[NATOM-1] == 1 && ( resid[NATOM-1] != resid[NATOM-2] || resname[NATOM-1] != resname[NATOM-2])) {
            monomerflag = 2;
        } else {
	    monomerflag = 3;
	}
    } else {
	    monomerflag = 1;
    }
    if (monomerflag < 1 || monomerflag > 3) {
	cout << "Incomplete/wrong information for monomers!!!" << endl;
	exit(1);
    }

// determine the special atoms in the monomer
    if (monomerflag < 3) { // 1st or the body monomer
	ind_a2 = atom_index[ibond_atom[NATOM-1][0]];
	if (monomerflag == 2) {
            ind_a1 = atom_index[ibond_atom[0][0]];
	} else {
	    ind_a1 = atom_index[0];
        } 
	name1 = atomname[ind_a1];
	name2 = atomname[ind_a2];

        name_next = atomname[NATOM-1];

	if (ind_a1 != ind_a2) {
            name_pre = atomname[ind_a1];
	}

	int ibond_count = 0;
	for (int i = 0; i < nbonds_atom[ind_a2]; i++) {
	    int ind_temp = atom_index[ibond_atom[ind_a2][i]];
	    if (ind_temp != ind_a1 && ind_temp != atom_index[NATOM-1] && ind_temp != atom_index[0]) {
		if (ibond_count == 0) {
   		    ind_ref1 = ind_temp;
		    name_ref1 = atomname[ind_temp];
		}
		if (ibond_count == 1) {
		    ind_a1_special = ind_temp;
		    name_a1_special = atomname[ind_temp];
		    break;
		}
		ibond_count++;
	    }
	    if (ind_a1 == ind_a2) {
                name_pre = atomname[0];
	    }
	}

        vector<float> vec0 = {x[ind_a1],y[ind_a1],z[ind_a1]};
        vector<float> vec1 = {x[ind_a2],y[ind_a2],z[ind_a2]};
        vector<float> vec2 = {x[ind_ref1],y[ind_ref1],z[ind_ref1]};
        vector<float> vec3 = {x[NATOM-1],y[NATOM-1],z[NATOM-1]};


        Vch = vector_subtract(vec2,vec1); // the reference C2-H3 bond
        Vcc_next = vector_subtract(vec3,vec1); // the reference C2-C1 bond
	if (ind_a1 != ind_a2) {
            Vcc = vector_subtract(vec0,vec1); // the reference C1-C2 bond
        } else {
	    if (monomerflag == 2) {
                vector<float> vec4 = {x[0],y[0],z[0]};
                Vcc = vector_subtract(vec4,vec1); // the reference C1-C2 bond
	    } else {
	        cout << "Warning: First monomer contains only one backbone atom!" << endl;
		ind_a1 = ind_a1_special;
                name_pre = atomname[ind_a1];
                vector<float> vec4 = {x[ind_a1],y[ind_a1],z[ind_a1]};
                Vcc = vector_subtract(vec4,vec1); // the reference C1-C2 bond
	    //    exit(1);
	    }
        }

    } else {
        ind_a1 = atom_index[ibond_atom[0][0]];
	name1 = atomname[ind_a1];
    }

    if (monomerflag > 1) {
        vector<float> vec1 = {x[ind_a1],y[ind_a1],z[ind_a1]};
        vector<float> vec4 = {x[0],y[0],z[0]};
        Vcc_pre = vector_subtract(vec1,vec4); // the reference C1-C2 bond
    }

    monomer_resname = resname[ind_a1];
    cout << "NATOM: " << NATOM << " nbonds_atom: " << nbonds_atom[NATOM-1] << " ind_a2: " << ind_a2 << " monomer_resname: " << monomer_resname << " monomerflag: " << monomerflag <<endl; //debug 
}



PSF::~PSF()
{
    psff.close();
    
    delete[] atom_index;
    delete[] segid;
    delete[] resid;
    delete[] resname;
    delete[] atomname;
    delete[] atomtype;
    delete[] charge;
    delete[] mass;
    delete[] beta;
    delete[] nbonds_atom;
    for (int i = 0; i < NBONDS; i++) {
        delete[] ibond[i];
    }
    delete[] ibond;
    for (int i = 0; i < NANGLES; i++) {
        delete[] iangle[i];
    }
    delete[] iangle;
    for (int i = 0; i < NDIHEDRALS; i++) {
        delete[] idihedral[i];
    }
    delete[] idihedral;
    for (int i = 0; i < NIMPROPERS; i++) {
        delete[] improper[i];
    }
    delete[] improper;
    ibond_atom.clear();
    atoms.clear();
    residues.clear();
    segments.clear();

    if (sizeof(x) > 0) {
        delete[] x;
    }
    if (sizeof(y) > 0) {
        delete[] y;
    }
    if (sizeof(z) > 0) {
        delete[] z;
    }
}
