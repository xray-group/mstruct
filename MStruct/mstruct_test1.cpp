/* 
 * mstruct_test1.cpp
 * 
 * MStruct++ - Object-Oriented computer program/library for MicroStructure analysis
 * 					   from powder diffraction data.
 * 
 * Copyright (C) 2009  Zdenek Matej
 * 
 * This file is part of MStruct++.
 * 
 * MStruct++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MStruct++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MStruct++.  If not, see <http://www.gnu.org/licenses/>.
 * 
 */
 
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>

#include "mstruct_tests.h"
#include "MStruct.h"

#include "cctbx/sgtbx/space_group.h"
#include "cctbx/sgtbx/rot_mx_info.h"
#include "scitbx/vec3.h"

#include "ObjCryst/IO.h"
#include "Quirks/VFNStreamFormat.h"

#include "CrystVector/CrystVector.h"
#include "ObjCryst/DiffractionDataSingleCrystal.h"

#include <boost/algorithm/string/trim.hpp>

using namespace std;
using namespace ObjCryst;

namespace MStruct {

int mstruct_test1(int argc, char* argv[], std::istream &ccin)
{
	cout << "Running test no. 1" << endl;

	// create a TextureModelHKL object
	TextureModelHKL texture_model;
	texture_model.SetName("Texture Model 1");
	CrystVector_REAL hklz(3);
	cout << "texture z-axis:" << endl;
	ccin >> hklz(0) >> hklz(1) >> hklz(2); ccin.ignore(1024,'\n');// ignore_comments(ccin); 
	cout << "texture x-axis:" << endl;
	CrystVector_REAL hklx(3);
	ccin >> hklx(0) >> hklx(1) >> hklx(2); ccin.ignore(1024,'\n');// ignore_comments(ccin);
	cout << "Space group symbol, forceFriedelSymmetry:" << endl;
	string str;
	int friedelsymmetry;
	ccin >> str >> friedelsymmetry; ccin.ignore(1024,'\n');// ignore_comments(ccin);
	Crystal crystal = Crystal(5.,5.,5.,str.c_str());
	texture_model.SetCrystal(crystal,friedelsymmetry>0);
	CrystVector_REAL params(6);
	params(0) = 1.; // phi weight
	params(1) = 0.*DEG2RAD; // phi
	params(2) = 15.*DEG2RAD; // phi fwhm
	params(3) = 0.8; // phi2 weight
	params(4) = 30.*DEG2RAD; // phi2
	params(5) = 15.*DEG2RAD; // phi2 fwhm
	texture_model.AddComponent(TextureModelHKL::TEXTURE_MODEL_HKL_GAUSS,
														 TextureModelHKL::TEXTURE_MODEL_HKL_GAUSS,
														 TextureModelHKL::TEXTURE_MODEL_HKL_AXIAL,
														 1., hklz, hklx, params);

	// create a copy of the TextureModelHKL object
	TextureModelHKL texture_model_copy = texture_model;
	texture_model_copy.SetName("Texture Model 2");
	texture_model_copy.GetPar("Phi_fwhm_0").SetValue(20.*DEG2RAD);
	//hklz(0) = 1.; hklz(1) =  1.; hklz(2) = 1.;
	//hklx(0) = 1.; hklx(1) = -1.; hklx(2) = 0.;
	texture_model_copy.AddComponent(TextureModelHKL::TEXTURE_MODEL_HKL_SIMEK,
														 		  TextureModelHKL::TEXTURE_MODEL_HKL_AXIAL,
														 		  TextureModelHKL::TEXTURE_MODEL_HKL_AXIAL,
														 			1.,hklz,hklx);
	texture_model_copy.RemoveComponent(0);
	
	texture_model_copy.PrepareForOdfCalc();

	// create optimization object and set refinable object
	MStruct::LSQNumObj lsqOptObj("lsqOptObj");
	lsqOptObj.SetRefinedObj(texture_model);

	for(int i=0; i<gRefinableObjRegistry.GetNb(); i++) {
		const RefinableObj &obj = gRefinableObjRegistry.GetObj(i);
		cout<<i<<'\t'<<obj.GetClassName()<<'\t'<<obj.GetName()<<endl;
		for(int j=0; j<obj.GetNbPar(); j++)
			cout<<"\t\t\t"<<setw(14)<<obj.GetPar(j).GetName()<<":\t"<<obj.GetPar(j).GetHumanValue()<<'\n';
	}

	return 0;
}

// Cross product
CrystVector_REAL cross(const CrystVector_REAL &u, const CrystVector_REAL &v)
{
	CrystVector_REAL w(3);
	w(0) = u(1)*v(2) - u(2)*v(1);
	w(1) = u(2)*v(0) - u(0)*v(2);
	w(2) = u(0)*v(1) - u(1)*v(0);
	return w;
}

// Matrix multiplication
CrystMatrix_REAL mult(const CrystMatrix_REAL &A, const CrystMatrix_REAL &B)
{
	CrystMatrix_REAL C(3,3);
	for(int i=0; i<3; i++) {
		for(int j=0; j<3; j++) {
			C(i,j) = 0.;
				for(int k=0; k<3; k++) {
					C(i,j) += A(i,k)*B(k,j); // Cij = Aik Bkj
	}	}	}
	return C;
}

// This should be a simple program to simulatwe principal directions
// in HKL-pole figures. 
int mstruct_test2(int argc, char* argv[], std::istream &ccin)
{
	cout << "Running test no. 2" << endl;

	// Get the basic crystal orientation system 
	CrystVector_REAL hklz(3);
	cout << "texture z-axis" << endl;
	ccin >> hklz(0) >> hklz(1) >> hklz(2); ccin.ignore(1024,'\n');// ignore_comments(ccin); 
	cout << "texture x-axis" << endl;
	CrystVector_REAL hklx(3);
	ccin >> hklx(0) >> hklx(1) >> hklx(2); ccin.ignore(1024,'\n');// ignore_comments(ccin);
	// Crystal
	cout << "Space group symbol" << endl;
	string sg_symbol;
	ccin >> sg_symbol; ccin.ignore(1024,'\n');// ignore_comments(ccin);
	// Lattice parameters
	cout << "Lattice parameters: a, b, c, alpha, beta, gamma" << endl;
	REAL cell_a, cell_b, cell_c, cell_alpha, cell_beta, cell_gamma;
	ccin >> cell_a >> cell_b >> cell_c >> cell_alpha >> cell_beta >> cell_gamma;
	ccin.ignore(1024,'\n');// ignore_comments(ccin);
	cell_alpha *= DEG2RAD; cell_beta *= DEG2RAD; cell_gamma *= DEG2RAD;
	
	UnitCell uc = UnitCell(cell_a,cell_b,cell_c,cell_alpha,cell_beta,cell_gamma,sg_symbol.c_str());
	
	// Create transformation matrix rotating the given texture vector base to conventional
	// orthonormal coordinate axes (e1,e2,e3).
	CrystMatrix_REAL A(3,3);
	{
		CrystVector_REAL u = hklx;
		CrystVector_REAL w = hklz;
	
		uc.MillerToOrthonormalCoords(u(0),u(1),u(2));
		uc.MillerToOrthonormalCoords(w(0),w(1),w(2));
	
		u *= 1./sqrt(pow(u(0),2) + pow(u(1),2) + pow(u(2),2));
		w *= 1./sqrt(pow(w(0),2) + pow(w(1),2) + pow(w(2),2));
		
		// orthogonalise u
		CrystVector_REAL t = w;
		REAL s = u(0)*w(0) + u(1)*w(1) + u(2)*w(2);
		if (1.-fabs(s)<1.e-4) {
			cerr << "Error: Wrong base system." << endl;
			exit(-1);
		}
		t *= -s; u += t; u *= 1./sqrt(1.-s);
		
		CrystVector_REAL v = cross(w,u);
		
		t = cross(v,w);
		A(0,0) = t(0); A(0,1) = t(1); A(0,2) = t(2);
		t = cross(w,u);
		A(1,0) = t(0); A(1,1) = t(1); A(1,2) = t(2);
		t = cross(u,v);
		A(2,0) = t(0); A(2,1) = t(1); A(2,2) = t(2);

		A *= 1./(t(0)*w(0) + t(1)*w(1) + t(2)*w(2));
	}
	
	// Base system orietation
	cout << "Euler angles (phi1, phi, phi2: orientation of the base system)" << endl;
	REAL phi1, phi, phi2;
	ccin >> phi1 >> phi >> phi2; ccin.ignore(1024,'\n');// ignore_comments(ccin);
	phi1 *= DEG2RAD; phi *= DEG2RAD; phi2 *= DEG2RAD;
	// Create rotation matrix G
	CrystMatrix_REAL G(3,3);
	{
		// phi2
		G(0,0) =  cos(phi2); G(0,1) = sin(phi2); G(0,2) = 0.;
		G(1,0) = -sin(phi2); G(1,1) = cos(phi2); G(1,2) = 0.;
		G(2,0) =         0.; G(2,1) =        0.; G(2,2) = 1.;
		CrystMatrix_REAL t(3,3);
		// phi
		t(0,0) = 1.; t(0,1) =        0.; t(0,2) = 0.;
		t(1,0) = 0.; t(1,1) =  cos(phi); t(1,2) = sin(phi);
		t(2,0) = 0.; t(2,1) = -sin(phi); t(2,2) = cos(phi);
		//G *= t;
		G = mult(G,t);
		// phi1
		t(0,0) =  cos(phi1); t(0,1) = sin(phi1); t(0,2) = 0.;
		t(1,0) = -sin(phi1); t(1,1) = cos(phi1); t(1,2) = 0.;
		t(2,0) =         0.; t(2,1) =        0.; t(2,2) = 1.;
		//G *= t;
		G = mult(G,t);
	}
	
	{cout<<"A"<<endl;
	for(int i=0;i<3;i++){
		for(int j=0; j<3; j++)
			cout<<setw(10)<<A(i,j);
		cout<<"\n";
	}
	}

	{cout<<"G"<<endl;
	for(int i=0;i<3;i++){
		for(int j=0; j<3; j++)
			cout<<setw(10)<<G(i,j);
		cout<<"\n";
	}
	}

	// Calculate complete rotation matrix
	//G *= A;
	G = mult(G,A);

	{cout<<"G"<<endl;
	for(int i=0;i<3;i++){
		for(int j=0; j<3; j++)
			cout<<setw(10)<<G(i,j);
		cout<<"\n";
	}
	}

	// Get the simulated diffraction HKL vector indixes
	CrystVector_REAL hkl(3);
	cout << "Simulated diffraction (h k l)" << endl;
	ccin >> hkl(0) >> hkl(1) >> hkl(2); ccin.ignore(1024,'\n');// ignore_comments(ccin);
	// force Friedel Law
	cout << "Force Friedel Law - if 1 center of symmetry added (otherwise 0)" << endl;
	int friedelsymmetry;
	ccin >> friedelsymmetry; ccin.ignore(1024,'\n');// ignore_comments(ccin);
	
	// Generate list of equivalent reflections
	CrystMatrix_REAL HKL_list = uc.GetSpaceGroup().GetAllEquivRefl(hkl(0),hkl(1),hkl(2),
																																 false,friedelsymmetry>0);
	
	int nhkl = HKL_list.rows();
	
	CrystVector_REAL PSI_list(nhkl);
	CrystVector_REAL PHI_list(nhkl);

	// for all equvalent reflections
	for(int ihkl=0; ihkl<nhkl; ihkl++) {
		CrystVector_REAL HKL(3);
		HKL(0) = HKL_list(ihkl,0); HKL(1) = HKL_list(ihkl,1); HKL(2) = HKL_list(ihkl,2);
		uc.MillerToOrthonormalCoords(HKL(0),HKL(1),HKL(2));
		CrystVector_REAL Q(3);
		for(int i=0; i<3; i++) {
			Q(i) = 0.;
			for(int j=0; j<3; j++) Q(i) += G(i,j)*HKL(j); // j
		} // i
		Q *= 1./sqrt(pow(Q(0),2) + pow(Q(1),2) + pow(Q(2),2));
		PSI_list(ihkl) = acos(Q(2)); // [0,pi]
		REAL nQ = sqrt(pow(Q(0),2)+pow(Q(1),2));
		PHI_list(ihkl) = (nQ>1.e-4) ? atan2(Q(1)/nQ,Q(0)/nQ) : 0.; // [-pi,pi] 
	} // ihkl
	
	// print results
	cout << endl;
	cout<<"#"<<setw(3)<<"H"<<setw(4)<<"K"<<setw(4)<<"L"<<setw(10)<<"PSI"<<setw(10)<<"PHI"<<endl;
	cout<<fixed<<showpoint<<setprecision(2);
	for(int ihkl=0; ihkl<nhkl; ihkl++) {
		cout<<setw(4)<<int(HKL_list(ihkl,0));
		cout<<setw(4)<<int(HKL_list(ihkl,1));
		cout<<setw(4)<<int(HKL_list(ihkl,2));
		cout<<setw(10)<<PSI_list(ihkl)*RAD2DEG;
		cout<<setw(10)<<PHI_list(ihkl)*RAD2DEG;
		cout<<"\n";
	}

	return 0;
}

// This should be a simple program for texture calculations. 
int mstruct_test3(int argc, char* argv[], std::istream &ccin)
{
	cout << "Running test no. 3" << endl;
	
	Crystal crystal = Crystal(5.,5.,5.,"Fm-3m");
	
	MStruct::TextureModelFiber OdfModel;
	
	OdfModel.SetUnitCell(&crystal);
	OdfModel.SetFiberHKL(0,0,1,true);
	
	MStruct::TextureOdfNumCalculator OdfNumCalcultor;
	
	OdfNumCalcultor.SetUnitCell(&crystal);
	OdfNumCalcultor.SetOdfModel(&OdfModel);
	
	//std::ofstream fs("odf_3.txt");
	//OdfNumCalcultor.ExportOdfXPert(fs);
	//fs.close();
	
	// testing meaning of equivalent (hkl) directions for REAL (hkl) values
	CrystMatrix_REAL equivHKLs = crystal.GetSpaceGroup().GetAllEquivRefl(1.1,1.,1.,false,true);
	for(int ihkl=0; ihkl<equivHKLs.rows(); ihkl++)
		cout<<setw(10)<<equivHKLs(ihkl,0)<<setw(10)<<equivHKLs(ihkl,1)<<setw(10)<<equivHKLs(ihkl,2)<<"\n";
	
	cout << "ODFNorm: " << OdfNumCalcultor.GetOdfNorm() << endl;
	
	cout << "PF:" << endl;
	for(int i=0; i<=36; i++) {
		REAL psi = i*M_PI/72;
		REAL v1 = OdfNumCalcultor.CalcOdfProjection(psi,0.,0,0,1);
		REAL v2 = OdfNumCalcultor.CalcOdfProjection(psi,25.*DEG2RAD,0,0,1);
		REAL v3 = OdfNumCalcultor.CalcOdfProjection(psi,0.,1,1,1);
		REAL v4 = OdfNumCalcultor.CalcOdfProjection(psi,0.,1,0,1);
		cout << setw(10) << psi*RAD2DEG << setw(18) << v1 << setw(18) << v2;
		cout << setw(18) << v3 << setw(18) << v4 << "\n"; 
		//cout << setw(10) << psi*RAD2DEG << setw(18) << v3 << "\n"; 
	}
	
	// Pole Figure simulation
	{
		// define the Pole figure grid
		CrystVector_REAL psi(37);
		CrystVector_REAL phi(72);
		
		for(int ipsi=0; ipsi<=36; ipsi++) psi(ipsi)=ipsi*2.5*DEG2RAD;
		for(int iphi=0; iphi< 72; iphi++) phi(iphi)=iphi*5.*DEG2RAD;
		
		// calculate the pole figure
		CrystMatrix_REAL pf = OdfNumCalcultor.CalcPFProjection(psi,phi,1.,0.,0.);
		
		// export the calculated Pole Figure
		ofstream os("test1_pf.txt");
		OdfNumCalcultor.ExportPFProjectionXPert(os,psi,phi,pf);
		os.close();
		
		// check the normalisation
		REAL sum = 0.;
		int npsi = psi.numElements();
		for(int ipsi=0; ipsi<npsi; ipsi++) {
			REAL w = (ipsi>0 && ipsi<npsi-1) ? 1. : 0.5;
			REAL s = sin(psi(ipsi));
			for(int iphi=0; iphi<phi.numElements(); iphi++)
				sum += w*pf(ipsi,iphi)*s;
		}
		sum *= (psi(1)-psi(0))*(phi(1)-phi(0)) / (2.*M_PI);
		
		cout << "PF norm: " << sum << endl;
	}
			
	return 0;
}

// Testing of CalcPeakParams function
int mstruct_test4(int argc, char* argv[], std::istream &ccin)
{
	cout << "Running test no. 4" << endl;

	// Create a ReflectionProfilePseudoVoigt object
	ObjCryst::ReflectionProfilePseudoVoigt reflProfilePV;
	
	// Set profile parameters
	REAL intensity, center, fwhm, k;
	cout << "parameters of the pV-function: intensity, center, fwhm, k" << endl;
	ccin >> intensity >> center >> fwhm >> k; ccin.ignore(1024,'\n');// ignore_comments(ccin); 
	
	reflProfilePV.SetProfilePar (pow(fwhm*DEG2RAD,2), 0., 0., k, 0.);
	
	// Calculate profile
	REAL xmin, xmax, xstep;
	cout << "xmin, xmax, xstep" << endl;
	ccin >> xmin >> xmax >> xstep; ccin.ignore(1024,'\n');// ignore_comments(ccin); 
	
	int nbPoints = int( (xmax-xmin)/xstep ) + 1;
	
	cout << "    " << nbPoints << " points will be calcualted in the range: (" << xmin << ", " << xmin + (nbPoints-1)*xstep << ")" << endl;
	
	CrystVector_REAL x(nbPoints);
	for(int i=0; i<nbPoints; i++) x(i) = (xmin + i*xstep)*DEG2RAD;
	
	CrystVector_REAL y = reflProfilePV.GetProfile(x,center*DEG2RAD,0,0,0);
	
	y *= intensity;
	
	// Save calculated data into the file
	string filename;
	cout << "output file name" << endl;
	ccin >> filename; ccin.ignore(1024,'\n');// ignore_comments(ccin);
	
	ofstream s(filename.c_str());
	s<<"# intensity: "<< intensity<< ", center: "<<center<<", fwhm: "<<fwhm<<", k: "<<k<<"\n";
	for(int i=0; i<nbPoints; i++) {
		char buf[1024];
		sprintf(buf,"  %7.3f    %.3f\n",x(i)*RAD2DEG,y(i));
		s << buf;
	}
	s.close();
	
	x *= RAD2DEG;
	//x += -center;
	// Calculate peak parameters
	PeakParams p = CalcPeakParams (x,y);
	cout<<"Calculated peak parameters:\n";
	cout<<"xmax: "<<p.xmax<<"\n";
	cout<<"ymax: "<<p.ymax<<"\n";
	cout<<"intensity: "<<p.intensity*DEG2RAD<<"\n";
	cout<<"fwhm: "<<p.fwhm<<"\n";
	cout<<"asym: "<<p.asym<<"\n";
	
	return 0;
}

// Testing of Laue symmetry determination from a SpaceGroupId
int mstruct_test5(int argc, char* argv[], std::istream &is)
{
	cout << "Running test no. 5" << endl;

	// input-string-stream to which the input is translated to ignore comments etc. 
	istringstream ccin;
	
	// Create a UnitCell Object
	REAL a, b, c, alpha, beta, gamma;
	string spaceGroupId;
	
	cout << "Space group symbol" << endl;
	read_line (ccin, is); // read a line (ignoring all comments, etc.)
  ccin >> spaceGroupId;
  cout << "lattice parameters: a(A), b(A), c(A)" << endl;
  read_line (ccin, is); // read a line (ignoring all comments, etc.)
  ccin >> a >> b >> c;
  cout << "lattice parameters: alpha(deg), beta(deg), gamma(deg)" << endl;
  read_line (ccin, is); // read a line (ignoring all comments, etc.)
  ccin >> alpha >> beta >> gamma;
	alpha *= DEG2RAD; beta *= DEG2RAD; gamma *= DEG2RAD;
	
	const ObjCryst::UnitCell uc(a,b,c,alpha,beta,gamma,spaceGroupId);
	
	// Get cctbx::sgtbx::space_group
	const cctbx::sgtbx::space_group & sg = uc.GetSpaceGroup().GetCCTbxSpg();
	
	// Get Laue group code
	const cctbx::sgtbx::matrix_group::code & LaueCode = sg.laue_group_type();
	
	using namespace cctbx::sgtbx::matrix_group;
	
	// Determine the Laue Symmetry from the Laue Code
	cout << "Laue symmetry: ";
	
	if      ( LaueCode == code_1b )
		cout << "-1";
	else if ( LaueCode == code_2_m )
		cout << "2/m";
	else if ( LaueCode == code_mmm )
		cout << "mmm";
	else if ( LaueCode == code_4_m )
		cout << "4/m";
	else if ( LaueCode == code_4_mmm )
		cout << "4/mmm";
	else if ( LaueCode == code_3b )
		cout << "-3";
	else if ( LaueCode == code_3bm )
		cout << "-3m";
	else if ( LaueCode == code_6_m )
		cout << "6/m";
	else if ( LaueCode == code_6_mmm )
		cout << "6/mmm";
	else if ( LaueCode == code_m3b )
		cout << "m-3";
	else if ( LaueCode == code_m3bm )
		cout << "m-3m";
	else
		cout << "unknown";
		
	cout << endl;
	
	// If monoclinic, get the unique axis
	if ( LaueCode == code_2_m )
		cout << "\tResult of: GetUniqueAxis()" << endl;
		switch (uc.GetSpaceGroup().GetUniqueAxis()) {
			case 0 : cout << "Error: Unique axis 'a' not supported." << endl; break;
			case 1 : cout << "Unique axis 'b'." << endl; break;
			case 2 : cout << "Unique axis 'c'." << endl; break;
			default : cout << "Error: Can not identify unique axis." << endl; break;
		}
	
	// If monoclinic, get the unique axis
	if ( LaueCode == code_2_m ) {
		cout << "\tResult of: rot_mx analysis" << endl;
		int axis_id = -1;
		for(std::size_t i=0; i<sg.n_smx(); i++) {
			const cctbx::sgtbx::rot_mx_info & info = sg.smx(i).r().info();
			if ( abs(info.type()) == 2 ) {
				const scitbx::vec3< int > & v = info.ev();
				// cctbx returns strange values ('angles' in 'int' ?)
				const REAL vn = sqrt( pow(REAL(v[0]),2) + pow(REAL(v[1]),2) + pow(REAL(v[2]),2) );
				if      ( abs(abs(v[0])/vn-1.)<1e-4 )
					axis_id = 0;
				else if ( abs(abs(v[1])/vn-1.)<1e-4 )
					axis_id = 1;
				else if ( abs(abs(v[2])/vn-1.)<1e-4 )
					axis_id = 2;
				else
					cout << "Error: Can not identify unique axis." << endl;
				break;
			}
		}
		switch (axis_id) {
			case 0 : cout << "Unique axis 'a'." << endl; break;
			case 1 : cout << "Unique axis 'b'." << endl; break;
			case 2 : cout << "Unique axis 'c'." << endl; break;
			default : cout << "Error: Can not identify unique axis." << endl; break;
		}
	}
	
	cout << " End of program." << endl ;
	
	return 0;
}

// XECs Computation
int mstruct_test6(int argc, char* argv[], std::istream &is)
{
	cout << "Running test no. 6" << endl;

	cout << " XECs Computation" << endl;
	
	// input-string-stream to which the input is translated to ignore comments etc. 
	istringstream ccin;
	
	// Create a UnitCell Object
	REAL a, b, c, alpha, beta, gamma;
	string spaceGroupId;
	
	cout << "Space group symbol" << endl;
	read_line (ccin, is); // read a line (ignoring all comments, etc.)
  ccin >> spaceGroupId;
  cout << "lattice parameters: a(A), b(A), c(A)" << endl;
  read_line (ccin, is); // read a line (ignoring all comments, etc.)
  ccin >> a >> b >> c;
  cout << "lattice parameters: alpha(deg), beta(deg), gamma(deg)" << endl;
  read_line (ccin, is); // read a line (ignoring all comments, etc.)
  ccin >> alpha >> beta >> gamma;
	alpha *= DEG2RAD; beta *= DEG2RAD; gamma *= DEG2RAD;
	
	ObjCryst::UnitCell uc(a,b,c,alpha,beta,gamma,spaceGroupId);
	uc.SetName("Material_UnitCell");
	
	// Create a XECs Object representing the Reuss-Voigt model 
	XECsReussVoigt XECsCalc;
	XECsCalc.SetName("XECs_Calculator");
	
	// Initialise XECs Object with the material UnitCell
	XECsCalc.SetUnitCell(uc);
	
	// Get material Stiffness constants
		
		// Get Stiffness constants names
		const vector< string > & vCijNames = XECsCalc.GetStiffnessConstantsNames();
		// Number of requred Stiffness constants
		const int n = vCijNames.size();
		
		// Read  Stiffness constants values (as a map)
		map<string, REAL> mapCijValues;
		
		cout << "material ";
		for(int i=0; i<n; i++) cout << vCijNames[i] << " ";
		cout << "constants (in GPa) in the format: C11 value C12 value etc.)" << endl;
				
		int npar = 0; bool newline = true;
		while (npar<n) {
			if (newline) { read_line (ccin, is); newline = false; } // read a line (ignoring all comments, etc.)
			string name;
			ccin >> name;
			if ( !(name.length()==3 && (name[0]=='C' || name[0]=='c')) ) { newline = true; continue; }
			REAL value;
			ccin >> value;
    		for_each( name.begin() , name.end() , char2upper() ); // name to upper case
    		mapCijValues[name] = value;
			npar++;
		}
		
		cout << "Constants:" << endl;
		map< string, REAL>::const_iterator it = mapCijValues.begin();
		for(int i=0; i<n; i++) {
			cout << it->first << ": " << it->second << " (GPa)" << endl;
			it++;
		}
		
	// Set material Stiffness constants
	XECsCalc.SetStiffnessConstants(mapCijValues);
	
	// Get Cij matrix
	const CrystMatrix_REAL & matrixCij = XECsCalc.GetCijMatrix();
	
	cout << "Stiffness tensor in the Voigt notation (GPa):" << endl;
	for(int i=0; i<6; i++) {
		for(int j=0; j<6; j++)
			cout << setw(10) << matrixCij(i,j);
		cout << endl;
	}
	
	// Get model weigth
	REAL weight;
	cout << "model weight (0..Reuss,1..Voigt)" << endl;
	read_line (ccin, is); // read a line (ignoring all comments, etc.)
	ccin >> weight;
	
	// Set model weigth
	//XECsCalc.SetParams(weight);
	
	// Read number of reflections
	int nhkl;
	cout << "number of reflections" << endl;
	read_line (ccin, is); // read a line (ignoring all comments, etc.)
	ccin >> nhkl;
	
	// Get reflection indices
	CrystMatrix_long hkl(nhkl,3);
	
	for(int i=0; i<nhkl; i++) {
		// Get reflection indices
		cout << "h k l reflection indices" << endl;
		read_line (ccin, is); // read a line (ignoring all comments, etc.)
		ccin >> hkl(i,0) >> hkl(i,1) >> hkl(i,2);
	}
	
	// Calculation
	cout << " h k l   " << setw(12) << "S1 (1/TPa)" << setw(12) << "S2 (1/TPa)" << endl;
	for(int i=0; i<nhkl; i++) {
		// Get XECS
		REAL s1, s2;
		XECsCalc.GetXECs(s1,s2,hkl(i,0),hkl(i,1),hkl(i,2));
	
		cout << " " << hkl(i,0) << " " << hkl(i,1) << " " << hkl(i,2) << "   ";
		cout << setw(12) << s1 *1e3 << setw(12) << s2 *1e3 << endl;
	}

	cout << "The second testing calculation" << endl;

	// Calculation
	cout << " h k l   " << setw(12) << "S1 (1/TPa)" << setw(12) << "S2 (1/TPa)" << endl;
	for(int i=0; i<nhkl; i++) {
		// Get XECS
		REAL s1, s2;
		XECsCalc.GetXECs(s1,s2,hkl(i,0),hkl(i,1),hkl(i,2));
	
		cout << " " << hkl(i,0) << " " << hkl(i,1) << " " << hkl(i,2) << "   ";
		cout << setw(12) << s1 *1e3 << setw(12) << s2 *1e3 << endl;
	}
	
	cout << " End of program." << endl ;
	
	return 0;
}

// Refraction effect
int mstruct_test7(int argc, char* argv[], std::istream &is)
{
	cout << "Running test no. 7" << endl;

	cout << " Refraction Effect" << endl;
	
	// input-string-stream to which the input is translated to ignore comments etc. 
	istringstream ccin;
	
	// read crystal name (to load)
	string crystal_name;
	cout << "Crystal name (Crystal will be loaded from structures.xml file)" << endl;
	read_line (ccin, is); // read a line (ignoring all comments, etc.)
	ccin >> crystal_name;
	
	// load crystal from structures file
	ObjCryst::Crystal * crystal;
	XMLCrystFileLoadObject("structures.xml","Crystal",crystal_name,crystal);
	/* Unfortunatelly there is a clear shortage in ObjCryst:
	  T*obj in XMLCrystFileLoadObject function is not referenced as
	  a reference or double pointer, hence the value of created
	  T*obj si not returned to a calling scope */
	//if(crystal==NULL) throw ObjCrystException("Crystal Not Found!");
	{	
		RefinableObj &obj = gRefinableObjRegistry.GetObj(crystal_name,"Crystal");
		crystal = &(dynamic_cast<ObjCryst::Crystal&>(obj));
	}
 	 	 
	// create RefractionPositionCorr Object, which handles chi0 calculation
	MStruct::RefractionPositionCorr * refCorr = new MStruct::RefractionPositionCorr;
	
	// create all structures like in a normal powder pattern simulation
	
	// create PowderPatternDiffraction Object
	ObjCryst::PowderPatternDiffraction * scattDataPowder = new ObjCryst::PowderPatternDiffraction;
	// add Crystal
	scattDataPowder->SetCrystal(*crystal);
	
	// create PowderPattern Objcet
	MStruct::PowderPattern * pattern = new MStruct::PowderPattern;
	// set radiation
	pattern->SetWavelength("CuA1");
	// add Crystalline phase
	pattern->AddPowderPatternComponent(*scattDataPowder);
	
	// create ReflectionProfile Object
	MStruct::ReflectionProfile * profile =
			new MStruct::ReflectionProfile(*crystal,pattern->GetRadiation()); // (stupid constructor)
	
	// set profile
	profile->SetParentPowderPatternDiffraction(*scattDataPowder);
	scattDataPowder->SetProfile(profile); // this can be done in SetParentPowderPatternDiffraction ???
	
	// add profile-component
	profile->AddReflectionProfileComponent(*refCorr);
	refCorr->SetParentReflectionProfile(*profile);
	
	// set crystal
	refCorr->SetCrystal(*crystal,true);
	
	// test
	
	// create single crystal scattering data
	ObjCryst::DiffractionDataSingleCrystal scattData(*crystal);
	
	// set radiation
	scattData.SetWavelength("CuA1");
	
	// force switch on Imag-Scattering-Factor
	scattData.SetIsIgnoringImagScattFact(false);
	
	// Set (000) reflection
	CrystVector_REAL h(1), k(1), l(1); h = 0; k = 0; l = 0;
	scattData.SetHKL(h,k,l);
	
	// Get scattering factor
	const CrystVector_REAL & Fhkl_real = scattData.GetFhklCalcReal();
	const CrystVector_REAL & Fhkl_imag = scattData.GetFhklCalcImag();
	
	cout << " Fhkl_real: " << FormatVertVector<REAL>(Fhkl_real,8,4);
	cout << " Fhkl_imag: " << FormatVertVector<REAL>(Fhkl_imag,8,4) << flush;
		  	
	// calc chi0 (classical electron radius rel = 2.8179e-5 A)
	complex<REAL> chi0(Fhkl_real(0),Fhkl_imag(0));
	
	chi0 *= - 2.8179e-5 * pow(scattData.GetRadiation().GetWavelength()(0),2) / (M_PI * crystal->GetVolume());
	
	cout << "chi0 :" << scientific << setprecision(4) << chi0 << endl;
	
	// calculate refraction correction for an interval of incidence angles
	REAL min, max, step;
	cout << "incidence angle (deg): min, max, step" << endl;
	read_line (ccin, is); // read a line (ignoring all comments, etc.)
	ccin >> min >> max >> step;
	
	int n = int((max-min)/step+0.1) + 1;
	
	ostringstream os;
	os << scientific << setprecision(3) << fixed;
	for(int i=0; i<n; i++) {
		REAL alpha = min + step*i;
		pattern->SetIncidenceAngle(alpha*DEG2RAD);
		profile->SetIncidenceAngle(alpha*DEG2RAD);
		REAL dth2 = refCorr->GetPositionCorr(0.*DEG2RAD,0.,0.,0.)*RAD2DEG;
		os << setw(12) << alpha << setw(12) << dth2 << "\n";
	}
	cout << os.str() << flush;
	
	cout << "Chemical formula test: " << endl;
	
	string formula;
	cout << "chemical formula" << endl;
	read_line (ccin, is); // read a line (ignoring all comments, etc.)
	// remove leading and terminating blank spaces
	string::size_type ind;
	ind = ccin.str().find_first_not_of(" \t\n");
	if(ind!=string::npos) ccin.str( ccin.str().substr(ind) );
	// get chemical formula from the input (line) - 2spaces,\t,\n,/ are considered as separators
	ind = ccin.str().find_first_of("\t\n/");
	if(ind==string::npos) ind = ccin.str().find("  ");
	if(ind!=string::npos) formula = ccin.str().substr(0,ind); // use line up to the found separator
	else formula = ccin.str(); // use whole line
	
	RefractionPositionCorr::ChemicalFormula formulaObj;
	
	formulaObj.SetFormula(formula);
	
	cout << " End of program." << endl ;
	
	return 0;
}

int mstruct_test8(int argc, char* argv[], std::istream &iss)
{
  // ObjCryst::CubicSpline test

  cout << "Running test no. 8" << endl;

  cout << " CubicSpline Test" << endl;
	
  // input-string-stream to which the input is translated to ignore comments etc. 
  istringstream ccin;
  
  // input filename
  string filename;
  cout << "input file name" << endl;
  read_line (ccin, iss); // read a line (ignoring all comments, etc.)
  ccin >> filename;
  
  // x,y, vectors
  std::vector< REAL > tx, ty;

  ifstream fin(filename.c_str());
  string s;
  while( getline(fin,s) ) {
    boost::algorithm::trim(s); // trim the string
    if( s.size() == 0 ) continue; // empty line
    if( s.size() >= 1 && s[0] == '#' ) {
      cout << "  ignoring comment: " << s << endl;
    } else {
      istringstream ss(s);
      double x, y;
      while( ss >> x >> y ) {
	cout << "  got a numbers: " << x << ", " << y << endl;
	tx.push_back(x);
	ty.push_back(y);
      }
    }
  }
  fin.close();
  
  // set zero values in x=0 and x=1
  CrystVector_REAL x(tx.size()+2);
  CrystVector_REAL y(tx.size()+2);
  
  x(0) = 0.;           y(0) = 0.;
  x(tx.size()+1) = 1.; y(tx.size()+1) = 0.;

  // copy data
  for(int i=0; i<tx.size(); i++ )
    { x(i+1) = tx[i]; y(i+1) = ty[i]; }

  CubicSpline p(x,y,0.,0.);
  
  const int N = 1000;

  CrystVector_REAL yy  = p(0.,1./N,N+1);
  CrystVector_REAL yy1 = p.Derivative(0.,1./N,N+1);
  CrystVector_REAL yy2 = p.SecondDerivative(0.,1./N,N+1);

  cout << setprecision(4);

  for(int i=0; i<=N; i++) {
    cout << setw(14) << 1./N*i << setw(14) << yy(i);
    cout << setw(14) << yy1(i) << setw(14) << yy2(i) <<"\n"; 
  }

  return 0;

} // mstruct_test8

int mstruct_test9(int argc, char* argv[], std::istream &iss)
{
  // ObjCryst::CubicSpline test

  cout << "Running test no. 9" << endl;

  cout << "  SizeDistribBroadeningEffect::UniformiseDistributionMC1() Test" << endl;
	
  // input-string-stream to which the input is translated to ignore comments etc. 
  istringstream ccin;
  
  // input filename
  string filename;
  cout << "distribution-file name" << endl;
  read_line (ccin, iss); // read a line (ignoring all comments, etc.)
  ccin >> filename;

  REAL temp = 1.e-7;
  long Niter = 0;
  
  cout << "temperature, nb. iterations" << endl;
  read_line (ccin, iss); // read a line (ignoring all comments, etc.)
  ccin >> temp >> Niter;

  MStruct::SizeDistribBroadeningEffect effect;

  effect.ReadDistributionFromFile(filename.c_str());

  effect.UniformizeDistributionMC1( Niter, temp );

  effect.WriteDistributionToFile(filename.c_str());

  return 0;

} // mstruct_test9

int mstruct_test10(int argc, char* argv[], std::istream &iss)
{
  // testing class TurbostraticHexStructWB
  cout << "/*** This is TurbostraticHexStructWB TEST ***/\n";
  
  TurbostraticHexStructWB turboStructEffect;
  turboStructEffect.SetName("diffData_tCarbon");

  // input-string-stream to which the input is translated to ignore comments etc. 
  istringstream ccin;

  {
    std::ofstream s1("mstruct-TurbostraticHexStructWB-atoms.txt");
    turboStructEffect.mi0Calculator.CreateHexLayerAtoms( 2.461, 20.0 );
    turboStructEffect.mi0Calculator.PrintAtoms(s1);
    s1.close();
  }
  
  /*  {
    std::ofstream s1("mstruct-TurbostraticHexStructWB-atoms-30.txt");
    turboStructEffect.mi0Calculator.CreateHexLayerAtoms( 2.461, 30.0 );
    turboStructEffect.mi0Calculator.PrintAtoms(s1);
    s1.close();
    }*/
  
  {
    std::ofstream s1("mstruct-TurbostraticHexStructWB-hist.txt");
    turboStructEffect.mi0Calculator.AccumHist2d( 20.0, true, true );
    turboStructEffect.mi0Calculator.PrintHistogram(s1);
    s1.close();
  }
  
  {
    std::cout << "float(" << sizeof(float) << "): " << std::numeric_limits<float>::min() << ", " << std::numeric_limits<float>::max() << '\n';
    std::cout << "double(" << sizeof(double) << "): " << std::numeric_limits<double>::min() << ", " << std::numeric_limits<double>::max() << '\n';
    std::cout << "REAL(" << sizeof(REAL) <<"): " << std::numeric_limits<REAL>::min() << ", " << std::numeric_limits<REAL>::max() << '\n';

    unsigned int nQ = 1024;
    REAL Qmax = 4*M_PI*sin(140.*M_PI/360.)/1.54;
    REAL dQ = Qmax/nQ;
    CrystVector_REAL Q(nQ);
    
    for(int m=0; m<nQ; m++)
      Q(m) = m*dQ;

    std::ofstream s1("mstruct-TurbostraticHexStructWB-i0.txt");
    /*turboStructEffect.mi0Calculator.SetQ( 4*M_PI*sin(  0.*M_PI/360.)/1.54,
					  4*M_PI*sin(140.*M_PI/360.)/1.54,
					  20.0, 40);*/
    turboStructEffect.mi0Calculator.SetQ(Q);    
    //turboStructEffect.mi0Calculator.CalcI0();
    //turboStructEffect.mi0Calculator.PrintI0(s1);
    s1.close();
    /*for(int i=0; i<10; i++)
      turboStructEffect.mi0Calculator.CalcI0();*/
   
    turboStructEffect.mi00lCalculator.SetQ(Q);
    //std::cout << "i00lCalculator.CalcIq - start" << std::endl;
    //turboStructEffect.mi00lCalculator.CalcIq(2, 6.88, 20.0);
    CrystMatrix_REAL mIq = turboStructEffect.mi00lCalculator.CalcIq(5, 6.88, 20.0);
    std::ofstream sf("MStruct-turboStructEffect-iq.txt");
    turboStructEffect.mi00lCalculator.PrintIq(sf);
    sf.close();
    /*for(int i=0; i<mIq.rows(); i++) {
      std::ostringstream ss;
      ss << "MStruct-turboStructEffect-iq-" << (i+1) << ".txt";
      std::ofstream sf(ss.str().c_str());
      sf << "#    Q    iq\n";
      CrystVecto
      for(int j=0; j<mIq.cols(); j++) {
	cout << std::fixed << std::showpoint << std::setprecision(5) << turboStructEffect.mi00lCalculator.mQ(j);
	cout << "  ";
	cout << std::scientific << std::showpoint << std::setprecision(4) << mIq(i,j);
      }
      sf.close();
      }*/
    //std::cout << "i00lCalculator.CalcIq - end" << std::endl;
    
    s1.open("mstruct-TurbostraticHexStructWB-iq.txt");
    turboStructEffect.mi00lCalculator.PrintIq(s1);
    s1.close();
    
    //turboStructEffect.mi00lCalculator.SetQ(Q);
    //turboStructEffect.mi00lCalculator.CalcIq(5, 6.88, 20.0);
    
    s1.open("mstruct-TurbostraticHexStructWB-i00l.txt");
    turboStructEffect.mi00lCalculator.CalcI00l(5, 6.88, 20.0);
    //turboStructEffect.mi00lCalculator.CalcI00l(12, 6.88, 20.0);
    turboStructEffect.mi00lCalculator.PrintI00l(s1);
    s1.close();
  }

  {
    // create PowderPattern Object
    MStruct::PowderPattern * pattern = new MStruct::PowderPattern;
    // set radiation
    pattern->SetWavelength("CuA1");
    //pattern->GetRadiation().SetLinearPolarRate(0.0); // no monochromator - no polarization
    cout << "Lambda: " << std::fixed << std::setprecision(7) << pattern->GetRadiation().GetWavelength()(0) << "\n";
    REAL tthmono = 45.*DEG2RAD; // LiF (monochromator)
    REAL A = cos(tthmono)*cos(tthmono);
    REAL f = (1.-A)/(1.+A);
    pattern->GetRadiation().SetLinearPolarRate(f);
    // add Background/TotalScattering phase
    pattern->AddPowderPatternComponent(turboStructEffect);

    // force switch of Imag-Scattering-Factor
    //.SetIsIgnoringImagScattFact(false);

    //pattern->SetPowderPatternPar(0.0, 140.*M_PI/180./1024, 1025);
    pattern->SetPowderPatternPar(15.0*M_PI/180., 0.05*M_PI/180., int((140.-15.)/0.05)+1);
    {
      CrystVector_REAL t(pattern->GetNbPoint());
      t = 0.;
      pattern->SetPowderPatternObs(t);
    }
    pattern->SetWeightToUnit();
    pattern->SetScaleFactor(turboStructEffect, 1.0);

    pattern->Prepare();
    pattern->GetPowderPatternCalc();

    pattern->SavePowderPattern("powderPattern.dat");
  }

  cout << " End of program." << endl ;

  return 0;
} // mstruct_test10

} // namespace MStruct


