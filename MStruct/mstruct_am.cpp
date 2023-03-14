/* 
 * mstruct_am.cpp
 * 
 * MStruct++ - Object-Oriented computer program/library for MicroStructure analysis
 * 					   from powder diffraction data.
 * 
 * Copyright (C) 2009-2014  Zdenek Matej, Charles University in Prague
 * Copyright (C) 2014-2022  Zdenek Matej, MAX IV Laboratory, Lund University
 * Copyright (C) 2016-2019  Milan Dopita, Jan Endres, Charles University in Prague
 * Copyright (C) 2017-2018  Jiri Wollman, Charles University in Prague
 * 
 * This file is part of MStruct++.
 * 
 * MStruct++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * MStruct++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MStruct++.  If not, see <https://www.gnu.org/licenses/>.
 * 
 */
 
/*
 *  MStruct refinement program.
 *
 */


#include "config.h"
#include "MStruct.h"
#include "mstruct_tests.h"

#include <fstream>
#include <iomanip>
#include <cstring>

#include "ObjCryst/Atom.h"
#include "ObjCryst/IO.h"

#include <complex>
#include <limits>

#include <algorithm>
#include <functional>

#include <stdlib.h> // rand, srand
#include <time.h> // time

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/case_conv.hpp>

using namespace ObjCryst;

namespace MStruct {

#define BACKGROUND_INTERPOLATED        0
#define BACKGROUND_INVX                1
#define BACKGROUND_CHEBYSHEV           2
#define BACKGROUND_CHEBYSHEV_LOCAL     3
#define BACKGROUND_TURBOSTRATIC_CARBON 4

class PowderPatternBackground : public ObjCryst::PowderPatternBackground {
public:
  void SetInterpolationModel(const string &modelName);
};

void PowderPatternBackground::SetInterpolationModel(const string &modelName)
{
  mInterpolationModel.SetChoice(modelName);
} 
/*
class imstructstream: public istream {
	public:
		cout ();
};
*/
} // namespace MStruct

// this function reads all empty and commented (//) lines from the input
// stream until an uncommented and non-empty line is found 
void ignore_comments(istream & s);

// rename a refinable parameter of the refinble object (obj), parameter with name (old_name)
// must be a parameter of the object (obj) and is renamed (with new_name) only in the case
// the paremeter is used, function returns a reference to the renamed parameter,
// exceptions not handled 
RefinablePar& rename_par(RefinableObj &obj, const string &old_name, const string &new_name);

// rename crystal lattice parameters and Scattering powers Biso factor parameters
// to ensure better unique access to them. Parameters names are extended with
// 'Crystal name' suffix (eg. 'a' -> 'a_Anatase)
void rename_crystal_params(Crystal &crystal);

// extended get parameter function, if a parameter name 'extended_name'
// is a simple parameter name (not closer specified by ':' deliminator)
// function simply returns a reference to a parameter of the refinable object
// identifing paremeter by its name 'extended_name' in the refinable object
// 'ref_obj' (same as the RefinableObj method GetPar(const string&), in the case
// the 'extended_name' is a composed extended name containing a deliminator ':'
// (eg. Anatase:a), functions tries to find the specified parameter in the specified
// object in the global refinable object registry and returna reference to it.
// (eg. for 'Anatase:a' function scours the refinable object called 'Anatase'
// for a parameter called 'a') 
RefinablePar& get_par(const string &extended_name, RefinableObj &ref_obj);

// Global counter of handled SIGINT signals

namespace MStruct {
  extern unsigned int Global_SIGINThandled_counter;
} // namespace MStruct

int main (int argc, char *argv[])
{
	if(argc>=2 && (string(argv[1])==string("-v") || string(argv[1])==string("--version"))) //print version
   {
   		// print version and license information
      cout << "version: " << mstruct_version_str << "\n";
      cout << "mstruct Copyright\n";
      cout << "(C) 2009-2018 Charles University in Prague\n";
      cout << "(C) 2014-2022 Zdenek Matej, MAX IV Laboratory, Lund University\n";
      cout << "e-mail: <zdenek.matej@maxiv.lu.se>\n";
      cout << "web: <https://xray.cz/mstruct/>\n";
      cout << "License GNU GPL: <https://www.gnu.org/licenses/old-licenses/gpl-2.0.html>.\n";
      cout << "This program comes with ABSOLUTELY NO WARRANTY;\n";
      cout << "This is free software, and you are welcome to redistribute it.\n";
      cout << flush;
      return 0;
   }

   cout << " Beginning program ...." << endl ;

   int level =12;
   if(argc>=2) //debug level hase been supplied
   {
      level=atoi(argv[1]);
   }
   VFN_DEBUG_GLOBAL_LEVEL(level);
   
   if(argc>=3) //it was supllied if the calc should be saved
   {
     bsavecalc=atoi(argv[2])!=0;
   }
   
   
  // initialisation of the random number genarator
   if(0) {
     srand((unsigned) time(NULL));
   }

  // input stream (standard input or a supplied imput file)
  // (a short comedy for MSVC - not working simply without pointers)
	istream *pointer_to_ccin = &cin;
	bool using_file = false;
	ifstream supplied_input_file;
	if(argc>=4) //an input file was supplied
	{
		supplied_input_file.open(argv[3]);
		pointer_to_ccin = &supplied_input_file;
		using_file = true;
	}
	istream & imp_file = *pointer_to_ccin;
	// input-string-stream to which the input is translated to ignore comments etc. 
	istringstream ccin;

  // job type (job_type>10 - run test number job_type-1)
	int job_type = 0;
	cout << "job type (0-data refinement,1-grid refinement)" << endl;
	read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
	ccin >> job_type;
	
	if( job_type > 10 ) {
		// run a test
		int val = 0;
		switch (job_type) {
			case 11:
				val = MStruct::mstruct_test1(argc, argv, ccin);
				break;
			case 12:
				val = MStruct::mstruct_test2(argc, argv, ccin);
				break;
			case 13:
				val = MStruct::mstruct_test3(argc, argv, ccin);
				break;
			case 14:
				val = MStruct::mstruct_test4(argc, argv, ccin);
				break;
			case 15:
				val = MStruct::mstruct_test5(argc, argv, imp_file);
				break;
			case 16:
				val = MStruct::mstruct_test6(argc, argv, imp_file);
				break;
			case 17:
				val = MStruct::mstruct_test7(argc, argv, imp_file);
				break;
			case 18:
				val = MStruct::mstruct_test8(argc, argv, imp_file);
				break;
			case 19:
			        val = MStruct::mstruct_test9(argc, argv, imp_file);
				break;
			case 20:
			        val = MStruct::mstruct_test10(argc, argv, imp_file);
				break;
		        case 21:
			        val = MStruct::mstruct_test11(argc, argv, imp_file);
				break;
			case 22:
			        val = MStruct::mstruct_test12(argc, argv, imp_file);
				break;
			default:
				cout << "Job/test type not recognised!" << endl;
				break;
		}
		if(using_file == true) supplied_input_file.close();
		return val;
	}

  // data filename (2th, obs, sigma), data type (0-xy,1-xysigma)
   string data_filename;
   int data_type;
   cout << "data filename (2th,obs(,sigma)),data format type (0-xy,1-xysigma)" << endl;
   read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
   ccin >> data_filename >> data_type;
   
	// maximum sin(theta)/lambda
	 string maxsinthetalambda;
	 cout <<"maximum sin(theta)/lambda" << endl;
	 read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
   ccin >> maxsinthetalambda;
   
  // background fileaname (2th, bkg), background type (linear, cubic spline)
   string bkg_filename;
   bool settings_overwrite_bgr = false;
   int bkg_type;
   cout << "background filename (2th,bkg),background type (0-linear,1-cubic spline)" << endl;
   read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
   ccin >> bkg_filename >> bkg_type;
  // if "general" background file is given, the general (multicomponent background is assumend)
  // bkg_type is the number of background components in such a case
  // if a not-"general" background file is suppled, a (default) simple background speciafied in the input file
  // is assumed
   vector< int > v_bkg_type_ids; // bkg. comp. type-ids
   vector< string > v_bkg_strings; // bkg. comp. filenames or types
   vector< int > v_bkg_type_ints; // bkg. comp. integer atributes
   vector< CrystVector_REAL > v_bkg_params; // bkg. comp. REAL-parameters vectors
   vector< CrystVector_long > v_bkg_params_flags; // bkg. comp. parameters flags vectors
   {
   	 string str(bkg_filename);
   	 for_each( str.begin() , str.end() , char2lower() );
   	 if ( str!=string("general") ) { // (default) simple background specified in the input file
   	 	 
   	  // background points which should be refined
   		 int nb_bkg_points_ref = 0;
   		 CrystVector_long bkg_points_ref; 
   		 cout << "number of refined background points" << endl;
   		 read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
   		 ccin >> nb_bkg_points_ref;
   		 if(nb_bkg_points_ref > 0) {
     		 bkg_points_ref.resize(nb_bkg_points_ref);
     		 cout << "numbers of refined background points (first is 0)" << endl;
     		 read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
     		 for(int i=0; i<nb_bkg_points_ref; i++) {
     		 if ( ccin.eof() ) read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
     			 ccin >>  bkg_points_ref(i);
     		 }
   		 }
   		 
   		 v_bkg_type_ids.push_back( BACKGROUND_INTERPOLATED );
   		 v_bkg_strings.push_back( bkg_filename ); // original upper-lower case filename
   		 v_bkg_type_ints.push_back( bkg_type ); 
   		 v_bkg_params.push_back( CrystVector_REAL(0) ); // empty vector
   		 v_bkg_params_flags.push_back( bkg_points_ref );
 
   	 } else {
   	 	 // multicomponent (general) background
   	 	 const int ncomp = bkg_type;
   	 	 for (int icomp=0; icomp<ncomp; icomp++) {
   	 	 	// background component type-name
   	 	 	 str.clear();
   	 	 	 bkg_type = 0;
   	 	 	 cout << "background component type (chebyshev, invX, interpolated), X-func. type (0-X,1-sin(Th))" << endl;
   	 	 	 read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
   	 	 	 ccin >> str >> bkg_type;
   	 	 	// convert to lowercase
   	 	 	 for_each( str.begin() , str.end() , char2lower() );
   	 	 	 bool btype_not_found = true;
   	 	 	 
   	 	 	 if ( btype_not_found && str==string("interpolated") ) {
   	 	 	 	// background filename
   	 	 	 	 string bkg_filename;
   				 int bkg_type;
   				 cout << "background filename (2th,bkg),background type (0-linear,1-cubic spline)" << endl;
   				 read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
   				 ccin >> bkg_filename >> bkg_type;
   				// background points which should be refined
   		 		 int nb_bkg_points_ref = 0;
   		 		 CrystVector_long bkg_points_ref; 
   		 		 cout << "number of refined background points" << endl;
   		 		 read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
   		 		 ccin >> nb_bkg_points_ref;
   		 		 if(nb_bkg_points_ref > 0) {
     		 		 bkg_points_ref.resize(nb_bkg_points_ref);
     		 		 cout << "numbers of refined background points (first is 0)" << endl;
     		 		 read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
     		 		 for(int i=0; i<nb_bkg_points_ref; i++) {
     		 			 if ( ccin.eof() ) read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
     			 		 ccin >>  bkg_points_ref(i);
     		 		 }
   		 		 }
   		 		 
   		 		 v_bkg_type_ids.push_back( BACKGROUND_INTERPOLATED );
   		 		 v_bkg_strings.push_back( bkg_filename ); // original upper-lower case filename
   		 		 v_bkg_type_ints.push_back( bkg_type ); 
   		 		 v_bkg_params.push_back( CrystVector_REAL(0) ); // empty vector
   		 		 v_bkg_params_flags.push_back( bkg_points_ref ); 
   		 		 
   		 		 btype_not_found = false;
   	 	 	 } // interpoalted background 
   	 	 	 
   	 	 	 if ( btype_not_found && str==string("invx") ) {

   	 	 	 	// no parameters, no options 	 	 	
				   	 	 	 	
   	 	 	 	 v_bkg_type_ids.push_back( BACKGROUND_INVX );
   		 		 v_bkg_strings.push_back(""); // empty string
   		 		 v_bkg_type_ints.push_back( bkg_type ); 
   		 		 v_bkg_params.push_back( CrystVector_REAL(0) ); // empty vector
   		 		 v_bkg_params_flags.push_back( CrystVector_long(0) );  // empty vector
   	 	 	 	
   	 	 	 	 btype_not_found = false;
   	 	 	 } // invX background
   	 	 	 
   	 	 	 if ( btype_not_found && str==string("chebyshev") ) {
   	 	 	 	// polynomial degree (number of coefficients+1)
   	 	 	   	 int ndegree = 0;
   	 	 	   	 cout << "polynomial degree (number of coefficients+1)" << endl;
   	 	 	   	 read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
   		 	   	 ccin >> ndegree;
   		 	   	 CrystVector_REAL bcoefs;
   		 	   	 CrystVector_long bcoef_flags;
   		 		 if(ndegree+1 > 0) {
     		 		 bcoefs.resize(ndegree+1);
     		 		 bcoef_flags.resize(ndegree+1);
     		 		 cout << "values of coefficients (starting with zero-order), " << ndegree+1 << " values expected" << endl;
     		 		 read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
     		 		 for(int i=0; i<=ndegree; i++) {
     		 			 if ( ccin.eof() ) read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
     			 		 ccin >>  bcoefs(i);
     		 		 }
     		 		 cout << "coefficients refinement flags (0-refined,1-fixed), " << ndegree+1 << " integer values expected" << endl;
     		 		 read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
     		 		 for(int i=0; i<=ndegree; i++) {
     		 			 if ( ccin.eof() ) read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
     			 		 ccin >>  bcoef_flags(i);
     		 		 }
   		 		 }
   		 		 
   		 		 v_bkg_type_ids.push_back( BACKGROUND_CHEBYSHEV );
   		 		 v_bkg_strings.push_back(""); // empty string
   		 		 v_bkg_type_ints.push_back( bkg_type ); 
   		 		 v_bkg_params.push_back( bcoefs );
   		 		 v_bkg_params_flags.push_back( bcoef_flags );
   	 	 	   	
   	 	 	   	 btype_not_found = false;
   	 	 	 } // chebyshev background

			 if ( btype_not_found && str==string("chebyshev_local") ) {
   	 	 	 	// nothing to be read here   		 		 
   		 		 v_bkg_type_ids.push_back( BACKGROUND_CHEBYSHEV_LOCAL );
   		 		 v_bkg_strings.push_back(""); // empty string
   		 		 v_bkg_type_ints.push_back( bkg_type ); 
   		 		 v_bkg_params.push_back( CrystVector_REAL(0) );
   		 		 v_bkg_params_flags.push_back( CrystVector_long(0) );
   	 	 	   	
   	 	 	   	 btype_not_found = false;
   	 	 	 } // local chebyshev background
			 
			 if ( btype_not_found && str==string("noncrystallinephase") ) { // noncrystallinePhase
			   // effect type and name
			   string effect_type, effect_name;
			   cout << "effect type(tubostraticCarbonWB,...), effect name" << endl;
			   read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
			   ccin >> effect_type >> effect_name;
			   for_each( effect_type.begin() , effect_type.end() , char2lower() );
			   if ( effect_type.compare("tubostraticcarbonwb")==0 ) { // tubostraticCarbonWB
			     // graphitic layer parameters
			     cout << "La(nm), varLa(nm^2), a0(A), BISOa(A)" << endl;
			     REAL La, varLa, a0, BISOa;
			     read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
			     ccin >> La >> varLa >> a0 >> BISOa;
			     // inter-layer parameters
			     cout << "Lc(nm), varLc(nm^2), c0(A), BISOc(A)" << endl;
			     REAL Lc, varLc, c0, BISOc;
			     read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
			     ccin >> Lc >> varLc >> c0 >> BISOc;
			     // random fraction
			     cout << "fraction of non-organized atoms (kdiso=0...none)" << endl;
			     REAL kdiso;
			     read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
			     ccin >> kdiso;

			     // model parameters (La, varLa, a0, BISOa, Lc, varLc, c0, BISOc, kdiso)
			     CrystVector_REAL mparams(9);
			     mparams(0) = La; mparams(1) = varLa; mparams(2) = a0; mparams(3) = BISOa;
			     mparams(4) = Lc; mparams(5) = varLc; mparams(6) = c0; mparams(7) = BISOc;
			     mparams(8) = kdiso;
			     
			     // set model configuration
			     v_bkg_type_ids.push_back( BACKGROUND_TURBOSTRATIC_CARBON );
			     v_bkg_strings.push_back( effect_name ); // effect name
			     v_bkg_type_ints.push_back( 0 ); // model type 
			     v_bkg_params.push_back( mparams ); // basic model params
			     CrystVector_long flags( mparams.numElements() );
			     flags = 0;
			     v_bkg_params_flags.push_back( flags );
			   } else
			     cerr << "< main(...) (Warning) Unknown nonCrystalline effect type: "<<effect_type<<" >"<<endl;
			     
			   btype_not_found = false;
   	 	 	 } // turbostratic
   	 	 	 
   	 	 	// Type of the given background component not found
 			 if(btype_not_found) {
 			 	cout << "Warning: Type of the given background component not recognised!" << endl;
 			 	continue;
 			 }
 			 
   	 	 } // icomp
   	 } // multicomponent (general) background
   } // background block

  // incidence angle
   REAL omega;
   cout << "incidence angle (deg)-2Theta scan, -1.0: 2Theta/Theta scan, -2.0: 2Theta/Theta scan with variable slits" << endl;
   read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
   ccin >> omega;
   omega *= DEG2RAD;
	
  // instrumental profile params
   REAL instW, instU, instV;
   cout << "instrumental profile params (W,U,V)" << endl;
   read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
   ccin >> instW >> instU >> instV;
   instW *= DEG2RAD*DEG2RAD; instU *= DEG2RAD*DEG2RAD;
   instV *= DEG2RAD*DEG2RAD;

   REAL instEta0, instEta1;
   cout << "instrumental profile params (Eta0,Eta1)" << endl;
   read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
   ccin >> instEta0 >> instEta1;

   REAL instAsym0, instAsym1, instAsym2, instAsym2ThetaMax;
   cout << "instrumental profile params (Asym0,Asym1,Asym2,Asym2ThetaMax(deg))" << endl;
   read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
   ccin >> instAsym0 >> instAsym1 >> instAsym2 >> instAsym2ThetaMax;
   instAsym2ThetaMax *= DEG2RAD; 

  // Create Diffraction data object
   MStruct::PowderPattern data;
   data.SetName("pattern0");
   data.SetIncidenceAngle(omega);
   
  // wavelength type and monochromator factor
   string wavelength_type;
   REAL pol_rate = 0.;
   cout << "wavelength type (Cu,CuA1 or value), linear polarization rate (A=0.8,f=(1-A)/(1+A)=0.111 graphite mon.,f=0. unmonochromatized)" << endl;
   read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
   ccin >> wavelength_type >> pol_rate;
   // detect wavelength type (is the first character number or decimal?)
   boost::trim(wavelength_type);
   if ( wavelength_type.length()>0 && (isdigit(wavelength_type[0]) || wavelength_type[0]=='.') ) {
     // wavelength value specified, check if not energy (eV, keV)
     boost::algorithm::to_lower(wavelength_type);
     if (wavelength_type.length()>2 && wavelength_type.compare(wavelength_type.length()-2,string::npos,"ev")==0)
       if (wavelength_type.length()>3 && wavelength_type.compare(wavelength_type.length()-3,string::npos,"kev")==0)
	 data.SetEnergy( atof(wavelength_type.substr(0,wavelength_type.length()-3).c_str()) ); // [keV]
       else
	 data.SetEnergy( atof(wavelength_type.substr(0,wavelength_type.length()-2).c_str())/1.e3 ); // [eV]
     else
       data.SetWavelength( atof(wavelength_type.c_str()) ); // [A]
   } else
     data.SetWavelength(wavelength_type); // Cu, CuA1, ...
   // set polarization rate
   data.SetLinearPolarRate(pol_rate);
   
  // Import measured data
   if(data_type==0)
     data.ImportPowderPattern2ThetaObs(data_filename.c_str());
   else if(data_type==1)
     data.ImportPowderPattern2ThetaObsSigma(data_filename.c_str());
  // set maximum sin(theta)/lambda
   {
   	 REAL val;
   	 istringstream s; s.str(maxsinthetalambda);
   	 s >> val;
   	 if (s.fail()==true) {
   	 		val = 1.1*sin(data.GetPowderPatternXMax()/2)/data.GetWavelength();
   	 		cout << "Using default value for maximum sin(theta)/lambda: " << val << endl; 
   	 }
   	 data.SetMaxSinThetaOvLambda(val);
   }
  // Create background and import background points, set backgroud options
   vector< ObjCryst::PowderPatternComponent* > vBackgroundComponents;
   for(int icomp=0; icomp<(int)v_bkg_type_ids.size(); icomp++)
   	 switch (v_bkg_type_ids[icomp]) {
   	 	 case BACKGROUND_INTERPOLATED:
   	 	   {
   	 	   	 const string bkg_filename(v_bkg_strings[icomp]);
   	 	   	 const int bkg_type = v_bkg_type_ints[icomp];
   	 	   	 
   	 		 MStruct::PowderPatternBackground * bkgData= new MStruct::PowderPatternBackground;
   			 bkgData->SetName(string("bkgData.")+bkg_filename);
   			 data.AddPowderPatternComponent(*bkgData);
   			 bkgData->ImportUserBackground(bkg_filename.c_str());
			 
   			 if(bkg_type==0)
     		   bkgData->SetInterpolationModel(string("Linear"));
   			 else if(bkg_type==1)
     		   bkgData->SetInterpolationModel(string("Spline"));
     		 
     		 vBackgroundComponents.push_back(bkgData);   	 
   	 	   }
   	 	   break;
   	 	 case BACKGROUND_INVX:
   	 	   {
   	 	   	  MStruct::PowderPatternBackgroundInvX * bkgData= new MStruct::PowderPatternBackgroundInvX;
   			  bkgData->SetName("bkgData_InvX");
   			  data.AddPowderPatternComponent(*bkgData);
   			 
   			  bkgData->SetXFunctionType(v_bkg_type_ints[icomp]);
   			  bkgData->UseVariableSlitIntensityCorr(omega*RAD2DEG<=-1.999);
   			  vBackgroundComponents.push_back(bkgData);
   			}
   	 	   break;
   	 	 case BACKGROUND_CHEBYSHEV:
   	 	   {
   	 	   	 MStruct::PowderPatternBackgroundChebyshev * bkgData= new MStruct::PowderPatternBackgroundChebyshev;
   			 bkgData->SetName("bkgData_Chebyshev");
   			 data.AddPowderPatternComponent(*bkgData);
   			 
   			 bkgData->SetCoefficients(v_bkg_params[icomp]);
   			 
   			 bkgData->SetXFunctionType(v_bkg_type_ints[icomp]);
   			 bkgData->UseVariableSlitIntensityCorr(omega*RAD2DEG<=-1.999);
   			 vBackgroundComponents.push_back(bkgData);
   	 	   }
   	 	   break;
		 case BACKGROUND_CHEBYSHEV_LOCAL:
   	 	   {
   	 	   	 MStruct::LocalBackgroundChebyshev * bkgData= new MStruct::LocalBackgroundChebyshev;
   			 bkgData->SetName("localBkg_Chebyshev");
   			 data.AddPowderPatternComponent(*bkgData);   		      
   			 
   			 bkgData->SetXFunctionType(v_bkg_type_ints[icomp]);
   			 bkgData->UseVariableSlitIntensityCorr(omega*RAD2DEG<=-1.999);
   			 vBackgroundComponents.push_back(bkgData);
   	 	   }
   	 	   break; 
	         case BACKGROUND_TURBOSTRATIC_CARBON:
		   {
		     cout << "TODO:: carbonWB" << endl;
		     MStruct::TurbostraticHexStructWB * turboStructEffect = new MStruct::TurbostraticHexStructWB;
		     turboStructEffect->SetName(v_bkg_strings[icomp]);
		     data.AddPowderPatternComponent(*turboStructEffect);   		      
   			 
		     //turboStructEffect->SetXFunctionType(v_bkg_type_ints[icomp]);
		     turboStructEffect->UseVariableSlitIntensityCorr(omega*RAD2DEG<=-1.999);
		     vBackgroundComponents.push_back(turboStructEffect);

		     // set model parameters
		     REAL La = v_bkg_params[icomp](0), varLa = v_bkg_params[icomp](1);
		     REAL a0 = v_bkg_params[icomp](2), BISOa = v_bkg_params[icomp](3);
		     turboStructEffect->SetLayerParameters(La,varLa,a0,BISOa);
		     REAL Lc = v_bkg_params[icomp](4), varLc = v_bkg_params[icomp](5);
		     REAL c0 = v_bkg_params[icomp](6), BISOc = v_bkg_params[icomp](7);
		     turboStructEffect->SetInterLayerParameters(Lc,varLc,c0,BISOc);
		     REAL fracDisor = v_bkg_params[icomp](8);
		     turboStructEffect->SetFractionDisorder(fracDisor);

		     // set TotalScattering phase Polarization and Absorption corrections parameters
		     // thickness = 1.e6 nm (1cm), depth = 0, absfactor = 9.2 1/cm , omega = -1 deg
		     turboStructEffect->mAbsorptionCorr.SetAbsorptionCorrParams( 1.e6, 0., 9.2, -1.*DEG2RAD); // TODO:: set correct values
		   }
		   break;
   	 	 default:
   	 	   cerr << "< main(...)\n";
				 cerr << "Program logical error during creating background objects.\n >" << endl; 
				 throw ObjCrystException("main(...): Program error.");
   	 	   break;
   	 } // switch

  // Add excluded regions
   {
     int nexreg = 0;
     cout << "number of excluded regions" << endl;
     read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
     ccin >> nexreg;
     for(int n=0;n<nexreg;n++) {
       cout << "min2Theta, max2Theta (deg)" << endl;
       double min2Theta, max2Theta;
       read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
       ccin >> min2Theta >> max2Theta;
       min2Theta *= DEG2RAD; max2Theta *= DEG2RAD;
       if(min2Theta<data.GetPowderPatternXMin())
	  		 min2Theta = data.GetPowderPatternXMin();
       if (max2Theta>data.GetPowderPatternXMax())
				 max2Theta = data.GetPowderPatternXMax();
       data.AddExcludedRegion(min2Theta,max2Theta);
     }
   }

	 vector< int > vHKLChoice;
	 vector< int > vHKLChoicePrint;
   vector< MStruct::PowderPatternDiffraction* > vDiffData;	   
   vector< MStruct::HKLPseudoVoigtBroadeningEffectA* > vHKLEffect;
	 vector<string> vPhasesNames;
	 vector< MStruct::SizeDistribBroadeningEffect* > vSizeDistribEffect;
	 vector<string> vSizeDistribFileName;
	  
  // Number of Phases
   int nbphases = 0;
   cout << "number of phases" << endl;
   read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
   ccin >> nbphases;

   for(int iphase=0; iphase<nbphases; iphase++) {
    // Create crystal structure
     Crystal *crystal = NULL;
    // Load Crystal from a file
     string phase_name;
     string crystal_filename;
     string crystal_name;
     cout << "phase name (diffDataCrystal)" << endl;
     read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
     ccin >> phase_name;
     cout << "filename, name (filename-crystal xml file,name-crystal name)" << endl;
     read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
     ccin >> crystal_filename >> crystal_name;
     vPhasesNames.push_back(phase_name);
		 XMLCrystFileLoadObject(crystal_filename,"Crystal",crystal_name,crystal);
		/* Unfortunatelly there is a clear shortage in ObjCryst:
	     T*obj in XMLCrystFileLoadObject function is not referenced as
	     a reference or double pointer, hence the value of created
	     T*obj si not returned to a calling scope */
	   //if(crystal==NULL) throw ObjCrystException("Crystal Not Found!");
 	 	 {	
 	 	 	 RefinableObj &obj = gRefinableObjRegistry.GetObj(crystal_name,"Crystal");
 	 		 crystal = &(dynamic_cast<ObjCryst::Crystal&>(obj));
 	 	 }
		// rename crystal lattice parameters (to ensure better acces to them)
		 //rename_crystal_params(*crystal);
		 
    // Add Crystal as a Crystalline phase
     if(job_type==3)
       vDiffData.push_back(new MStruct::SizeDistribPowderPatternDiffraction);
     else
       vDiffData.push_back(new MStruct::PowderPatternDiffraction);
     vDiffData[iphase]->SetName(phase_name);
    // rename "Global Biso" -> "Global_Biso"
     rename_par(*vDiffData[iphase],"Global Biso","Global_Biso");
    // Set crystal
     vDiffData[iphase]->SetCrystal(*crystal);
    // Add crystalline phase
     data.AddPowderPatternComponent(*vDiffData[iphase]);
     
    // Set absorption and texture corr. params.
     REAL thickness, depth, afactor;
     cout << "absorp corr params: thickness(nm),depth(nm),abs.factor" << endl;
     read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
     ccin >> thickness >> depth >> afactor;

     vDiffData[iphase]->SetAbsorptionCorrParams(thickness,depth,afactor,omega);

	  // Texture
     vDiffData[iphase]->SetTextureCorrParams(omega);

     int nbtexturephases;
     cout << "nb of texture phases" << endl;
     read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
     ccin >> nbtexturephases;
     for(int i=0;i<nbtexturephases;i++) {
     	 REAL fraction;
     	 bool sflag = false;
     	 cout << "texture params: texture component fraction, force texture symmetry(0,1)" << endl; 
       read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
       ccin >> fraction >> sflag;

       int offset = 0;

			 CrystVector_REAL params(14);
			 params = 0.;
			 
       cout << "component base hz kz lz (space delimitaded)" << endl;
       read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
       ccin >> params(offset+6) >> params(offset+7) >> params(offset+8);
    
       cout << "component base hx kx lx (space delimitaded)" << endl;
       read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
       ccin >> params(offset+9) >> params(offset+10) >> params(offset+11);

       cout << "th component: weight (0.-1.), th0 (deg), dth0 (deg)"  << endl;
       read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
       ccin >> params(offset+0) >> params(offset+1) >> params(offset+2);
       params(offset+1) *= DEG2RAD; params(offset+2) *= DEG2RAD;
    
       cout << "psi component: weight (0.-1.), psi0 (deg), dpsi0 (deg)"  << endl;
       read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
       ccin >> params(offset+3) >> params(offset+4) >> params(offset+5);
       params(offset+4) *= DEG2RAD; params(offset+5) *= DEG2RAD;
      
       cout << "tilt: eth0 (deg), epsi (deg)" << endl;
       read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
       ccin >> params(offset+12) >> params(offset+13);
       params(offset+12) *= DEG2RAD; params(offset+13) *= DEG2RAD;
       
       vDiffData[iphase]->AddTextureCorrPhase(fraction,params,sflag);
     }

		 // not used now
     /*if(nbtexturephases>0) {
       REAL eth0, epsi0;
       cout << "tilt: eth0(deg),epsi0(deg)" << endl;
       read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
       ccin >> eth0 >> epsi0;
       eth0 *= DEG2RAD; epsi0 *= DEG2RAD;

       //vDiffData[iphase]->SetTextureTilt(eth0,epsi0);
     }*/
  
    // set HKL intensity corr
     cout << "hkl file(0-don't use,1-generate,2-free all,3-read)" << endl;
     vHKLChoice.push_back(0);
     vHKLChoicePrint.push_back(0);
     read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
     ccin >> vHKLChoice[iphase];
     vDiffData[iphase]->SetHKLIntensityCorrChoice(vHKLChoice[iphase]);
     if (vHKLChoice[iphase]>0) ccin >> vHKLChoicePrint[iphase];
     if(vHKLChoice[iphase]>0 && vHKLChoice[iphase]<3) {
       REAL level= (vHKLChoice[iphase]==2) ? 0.04 : -1.;
       vDiffData[iphase]->GenerateHKLIntensityCorrForAllReflections(level);
      // write IhklCorr params to file (later will be reloaded) 
       vDiffData[iphase]->WriteHKLIntensityCorrToFile(string(string("Ihkl_")+phase_name+string(".dat")).c_str());
     }

    // Profile (instrumental + physical profile)
     MStruct::ReflectionProfile * reflProfile
       = new MStruct::ReflectionProfile(*crystal,data.GetRadiation());
     reflProfile->SetParentPowderPatternDiffraction(*vDiffData[iphase]);
     reflProfile->SetIncidenceAngle(omega);
    // Number of additional broadening effects
     int nbbeffects = 0;
     cout << "number of additional broadening effects" << endl;
     read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
     ccin >> nbbeffects;
    // Instrumental profile
     MStruct::ReflectionProfileComponent * instReflProfileComp = NULL;
    // if U, V and W instrumental parameters are all very close to zero, than no instrumental profile is added
     if ( abs(instW) > numeric_limits<REAL>::epsilon() || abs(instV) > numeric_limits<REAL>::epsilon() ||
          abs(instU) > numeric_limits<REAL>::epsilon() )
   	 {
		 	 MStruct::ReflectionProfileComponent * reflProfileComp;
		 	 if(nbbeffects>0) {
			 	 MStruct::PseudoVoigtBroadeningEffectA * t = new MStruct::PseudoVoigtBroadeningEffectA;
			 	 t->SetProfilePar(instW,instU,instV,instEta0,instEta1);
			 	 t->SetAsymXMax(instAsym2ThetaMax);
	     	 reflProfileComp = t;
		 	 } else {
			 	 MStruct::PseudoVoigtBroadeningEffect * t = new MStruct::PseudoVoigtBroadeningEffect;
			 	 // PseudoVoigtBroadeningEffect ignores Asymmetry. Why ??? 
			 	 t->SetProfilePar(instW,instU,instV,instEta0,instEta1);
			 	 t->SetAsymXMax(instAsym2ThetaMax);
	     	 reflProfileComp = t;
		 	 }
		 	 reflProfileComp->SetName("instrProfile");
			 reflProfileComp->GetPar("Asym0").SetIsLimited(false);
			 reflProfileComp->GetPar("Asym1").SetIsLimited(false);
	     reflProfileComp->GetPar("Asym2").SetIsLimited(false);
			 reflProfileComp->GetPar("Asym0").SetValue(instAsym0);
	     reflProfileComp->GetPar("Asym1").SetValue(instAsym1);
	     reflProfileComp->GetPar("Asym2").SetValue(instAsym2);
	     reflProfileComp->SetParentReflectionProfile(*reflProfile);
	     reflProfile->AddReflectionProfileComponent(*reflProfileComp);
	    // save the pointer to the instrumental profile component
	     instReflProfileComp = reflProfileComp;
     }
    // Other broadening effects (including virtual and unused effects)
     vector<ObjCryst::ReflectionProfile*> vReflProfiles; // vector of ReflectionProfile Objects
     vector<MStruct::ReflectionProfileComponent*> vReflProfComponents; // vector of ReflectionProfileComponent Objects
     vector< vector<string> > vReflProfCompNames; // vector of ReflectionProfileComponent-s names for each ReflectionProfile Object
     int parentReflProfNb = -1; // position (number) of a parent ReflectionProfile Object in the vector "vReflProfiles"
     vHKLEffect.push_back(NULL);
     for(int ibeffect=0; ibeffect<nbbeffects; ibeffect++) {
     	// broadening component type
     	 string btype, bname;
     	 cout << "broadening component type (pVoigt(A),SizeLn,SizeDistrib,dislocSLvB,faultsVfcc,HKLpVoigtA),effect name,comp=" << ibeffect+1 << endl;
     	 read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
     	 ccin >> btype >> bname;
     	// broadening type has not been recognised yet
     	 bool btype_not_found = true;
     	// pseodoVoigt broadening
     	 if(btype_not_found && (btype==string("pVoigt") || btype==string("pVoigtA"))) {
    	 		MStruct::PseudoVoigtBroadeningEffectA * pVoigtEffect
       			= new MStruct::PseudoVoigtBroadeningEffectA;
       		pVoigtEffect->SetName(bname);
       		vReflProfComponents.push_back(pVoigtEffect);
       		
       	 // parameters
       	 	REAL W, U, V;
   				cout << "profile params (W,U,V)" << endl;
   				read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
   				ccin >> W >> U >> V;
   				W *= DEG2RAD*DEG2RAD; U *= DEG2RAD*DEG2RAD; V *= DEG2RAD*DEG2RAD; 

   				REAL Eta0, Eta1;
   				cout << "profile params (Eta0,Eta1)" << endl;
   				read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
   				ccin >> Eta0 >> Eta1;
       	 	
       	 	REAL Asym0, Asym1, Asym2, Asym2ThetaMax;
   				cout << "profile params (Asym0,Asym1,Asym2,Asym2ThetaMax(deg))" << endl;
   				read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
   				ccin >> Asym0 >> Asym1 >> Asym2 >> Asym2ThetaMax;
   				Asym2ThetaMax *= DEG2RAD; 
       			
     			pVoigtEffect->SetProfilePar(W,U,V,Eta0,Eta1);
     			pVoigtEffect->GetPar("Asym0").SetValue(Asym0);
     			pVoigtEffect->GetPar("Asym1").SetValue(Asym1);
     			pVoigtEffect->GetPar("Asym2").SetValue(Asym2);
     			pVoigtEffect->SetAsymXMax(Asym2ThetaMax);
     			//pVoigtEffect->SetParentReflectionProfile(*reflProfile);
     			//reflProfile->AddReflectionProfileComponent(*pVoigtEffect);
     			btype_not_found = false;
     	 }
     	// Size broadening
     	 if(btype_not_found && btype==string("SizeLn")) {
     	 	  MStruct::SizeBroadeningEffect * sizeEffect
       			 = new MStruct::SizeBroadeningEffect;
       		cout << "SizeBroadeningEffect (1):"<< sizeEffect->GetName() << ", bname:" << bname << endl;
       		sizeEffect->SetName(bname);
       		cout << "SizeBroadeningEffect (2):"<< sizeEffect->GetName() << endl;
       		vReflProfComponents.push_back(sizeEffect);
       		
       	 // parameters
       		REAL M, sigma;
       		cout << "M(nm),sigma" << endl;
       		read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
       		ccin >> M >> sigma;
     			sizeEffect->SetProfilePar(M,sigma);
     			//sizeEffect->SetParentReflectionProfile(*reflProfile);
     			//reflProfile->AddReflectionProfileComponent(*sizeEffect);
     			btype_not_found = false;
     	 }

	// Simple Interference affected Size broadening
	 if(btype_not_found && btype==string("SizeInterfSimple")) {
	   MStruct::InterferenceSimpleSizeBroadeningEffect * sizeEffect
	     = new MStruct::InterferenceSimpleSizeBroadeningEffect;
	   cout << "SizeBroadeningEffect (1): "<< sizeEffect->GetName() << ", bname: " << bname << endl;
	   sizeEffect->SetName(bname);
	   cout << "SizeBroadeningEffect (2): "<< sizeEffect->GetName() << endl;
	   vReflProfComponents.push_back(sizeEffect);
       		
       	 // parameters
	   REAL MG, sigmaG, MC, sigmaC, s0coherence;
	   cout << "MGrains(nm),sigmaGrains" << endl;
	   read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
	   ccin >> MG >> sigmaG;
	   cout << "MCrystallites(nm),sigmaCrystallites" << endl;
	   read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
	   ccin >> MC >> sigmaC;
	   cout << "s0=1/d (1/A) coherence limit" << endl;
	   read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
	   ccin >> s0coherence;
	   sizeEffect->SetProfilePar(MG,sigmaG,MC,sigmaC,s0coherence);
	   //sizeEffect->SetParentReflectionProfile(*reflProfile);
	   //reflProfile->AddReflectionProfileComponent(*sizeEffect);
	   btype_not_found = false;
     	 }
     	// Size broadening (refinable size distribution) 
     	 if(btype_not_found && btype==string("SizeDistrib")) {
	   MStruct::SizeDistribBroadeningEffect * sizeEffect
	                               = new MStruct::SizeDistribBroadeningEffect;
	   sizeEffect->SetName(bname);
	   vReflProfComponents.push_back(sizeEffect);
	   
	  // parameters
	   
	   // distribution type/source (file/generated)
	   string distribType;
	   REAL lsqConstraintScale;
	   
	   cout << "distribution type/source (file/create), LSQ constraint scale factor" << endl;
	   read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
	   ccin >> distribType >> lsqConstraintScale;
	   
	   // convert string to lowercase
	   for_each( distribType.begin() , distribType.end() , char2lower() );
	   
	   string filename;

	   if ( distribType==string("file") ) {
	    // distribution will be loaded from file
	     cout << "name of file with prescribed weighted distribution" << endl;
	     read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
	     ccin >> filename;
	   } else if ( distribType==string("create") ) {
	    // build the distribution and save the distribution in the file (to be reload later)
	     cout << "histogram spacing type (linear/log/sqrt), Dmin (nm), Dmax (nm), nb. of intervals" << endl;
	     string htype;
	     REAL Dmin, Dmax;
	     int nbIntervals;
	     read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
	     ccin >> htype >> Dmin >> Dmax >> nbIntervals;
	     // convert string to lowercase
	     for_each( htype.begin() , htype.end() , char2lower() );
	     sizeEffect->BuildDistribution(Dmin*10.,Dmax*10.,nbIntervals,htype);
	     // create distribution file name
	     filename = string("wD_") + crystal->GetName() + string(".dat");
	     // save data - will be reloaded later
	     //        (this means in addition that LSQRegType will also be altered later)
	     sizeEffect->WriteDistributionToFile( filename.c_str() );
	   }
	   
	   // set LSQ constraint scale
	   sizeEffect->SetLSQConstraintScale(lsqConstraintScale);

	   // get nb. of LSQ regularization methods used
	   int nbLSQRegOp = 0;
	   cout << "number of LSQ regularization methods used" << endl;
	   read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
	   ccin >> nbLSQRegOp;

	   // read all LSQ regularization operators
	   for(int iOp=0; iOp<nbLSQRegOp; iOp++) {

	     // get LSQ regularization type and weight
	     string LsqRegType;   
	     REAL alpha;

	     cout << "LSQ regularization method type (none/derivative/curvature), LSQ weight factor, [Niter]" << endl;
	     read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
	     ccin >> LsqRegType >> alpha;
	   
	     // convert string to lowercase
	     for_each( LsqRegType.begin() , LsqRegType.end() , char2lower() );

	     if ( LsqRegType==string("none") )
	       sizeEffect->AddLSQRegularizationMethod( MStruct::SizeDistribBroadeningEffect::LSQRegOpt_None );
	     else if ( LsqRegType==string("derivative") || LsqRegType==string("deriv") )
	       sizeEffect->AddLSQRegularizationMethod( MStruct::SizeDistribBroadeningEffect::LSQRegOpt_DistribDeriv,
						       alpha );
	     else if ( LsqRegType==string("volume-derivative") || LsqRegType==string("vol-deriv") )
	       sizeEffect->AddLSQRegularizationMethod( MStruct::SizeDistribBroadeningEffect::LSQRegOpt_VolumeDistribDeriv,
						       alpha );
	     else if ( LsqRegType==string("both-derivative") || LsqRegType==string("both-deriv") )
	       sizeEffect->AddLSQRegularizationMethod( MStruct::SizeDistribBroadeningEffect::LSQRegOpt_BothDistribDeriv,
						       alpha );
	     else if ( LsqRegType==string("curvature") || LsqRegType==string("curv") )
	       sizeEffect->AddLSQRegularizationMethod( MStruct::SizeDistribBroadeningEffect::LSQRegOpt_CurvIntegral,
						       alpha );
	     else if ( LsqRegType==string("uniformize") || LsqRegType==string("unif") || LsqRegType==string("uniformizeMC1") ) {
	       // this is strictly not a regularization method but rather a distrubition uniformization should be applied at
	       // the beginning of each optimization run
	       // ... load an additional parameter (nb of uniformization algorithm steps)
	       long Niter;
	       ccin >> Niter;
	       sizeEffect->SetUniformizeDistributionAtBeginOptimization ( true, Niter, alpha );
	     }
	     else {
	       cerr<<"Error: Unknow LSQ rgularization method: "<<LsqRegType;
	       cerr<<" for SizeDistribBroadeningEffect: "<<sizeEffect->GetName()<<endl;
	     }
	   } // iOp
	   
	   // do other common work to finish SizeDistribBroadeningEffect settings

	   //sizeEffect->ReadDistributionFromFile(filename.c_str());
	   vSizeDistribEffect.push_back(sizeEffect);
	   vSizeDistribFileName.push_back(filename);
	   // this effect is registered by the PowwderPattern as an additional LSQ func obj later
	   //sizeEffect->SetParentReflectionProfile(*reflProfile);
	   //reflProfile->AddReflectionProfileComponent(*sizeEffect);
	   btype_not_found = false;
     	 }
	// Size broadening (refinable size distribution - random version with same bits) 
     	 if(btype_not_found && btype==string("RandomSizeDistrib")) {
	   MStruct::RandomSizeDistribBroadeningEffect * sizeEffect
	     = new MStruct::RandomSizeDistribBroadeningEffect;
	   sizeEffect->SetName(bname);
	   vReflProfComponents.push_back(sizeEffect);
       		
	  // parameters
	   int Nbins;
	   REAL Dmax;
	   CrystVector_REAL distrib;

	   cout << "nb. of distribution bins, Dmax value (nm)" << endl;
	   read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
	   ccin >> Nbins >> Dmax;
	   
	   cout << Nbins << " - values of crystallite size distribution" << endl;
	   distrib.resize(Nbins);
	   read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
	   for(int i=0; i<Nbins; i++) ccin >> distrib(i);

	   sizeEffect->SetDistribution(Nbins,Dmax,distrib);
     			
	   //sizeEffect->SetParentReflectionProfile(*reflProfile);
	   //reflProfile->AddReflectionProfileComponent(*sizeEffect);
	   btype_not_found = false;
     	 }
	 
	 
	// Size broadening from ellipsoids with log-normal distribution
	 if(btype_not_found && ( btype==string("EllipSize") ) ) {
	   MStruct::EllipSizeBroadeningEffect * sizeEffect
	     = new MStruct::EllipSizeBroadeningEffect;
	   sizeEffect->SetName(bname);
	   vReflProfComponents.push_back(sizeEffect);

	  // rod axis
	   REAL hh, kk, ll;
	   cout << "(hkl) indices of the ellipsoid revolution axis (i.e. indices of the plane perpendiculat to revolution axis)" << endl;
	   read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
	   ccin >> hh >> kk >> ll;
	   sizeEffect->SetEllipAxis(hh,kk,ll);

	  // model parameters-set
	   string paramset;
	   cout << "choice of parameters set (AC, AE, CE)" << endl;
	   read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
	   ccin >> paramset;

	   boost::algorithm::to_lower(paramset);
	   if(paramset==string("ac")) {
	     sizeEffect->SetModelParSet(MStruct::EllipSizeBroadeningEffect::PARAM_SET_AC);
	     REAL diameterA, sigma, diameterC;
	     cout << "diameterA (nm), sigma, diameterC (nm)" << endl;
	     read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
	     ccin >> diameterA >> sigma >> diameterC;
	     sizeEffect->SetEllipDiameterA(diameterA, sigma);
	     sizeEffect->SetEllipDiameterC(diameterC, sigma);
	   }
	   else if(paramset==string("ae")) {
	     sizeEffect->SetModelParSet(MStruct::EllipSizeBroadeningEffect::PARAM_SET_AE);
	     REAL diameterA, sigma, epsilon;
	     cout << "diameterA (nm), sigma, epsilon" << endl;
	     read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
	     ccin >> diameterA >> sigma >> epsilon;
	     sizeEffect->SetEllipDiameterA(diameterA,sigma);
	     sizeEffect->SetEllipticity(epsilon);
	   }
	   else if(paramset==string("ce")) {
	     sizeEffect->SetModelParSet(MStruct::EllipSizeBroadeningEffect::PARAM_SET_CE);
	     REAL diameterC, sigma, epsilon;
	     cout << "diameterC (nm), sigma, epsilon" << endl;
	     read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
	     ccin >> diameterC >> sigma >> epsilon;
	     sizeEffect->SetEllipDiameterC(diameterC,sigma);
	     sizeEffect->SetEllipticity(epsilon);
	   }
	   else
	     throw ObjCrystException("Unknown model parameters-set option for CircRodsGamma broadening effect!");
	   
	   btype_not_found = false;
	 }
	 
	 
	// Size broadening from circular rods with Gamma distribution
	 if(btype_not_found && ( btype==string("CircRodsGamma") ) ) {
	   MStruct::CircRodsGammaBroadeningEffect * sizeEffect
	     = new MStruct::CircRodsGammaBroadeningEffect;
	   sizeEffect->SetName(bname);
	   vReflProfComponents.push_back(sizeEffect);

	  // rod axis
	   REAL hh, kk, ll;
	   cout << "(hkl) indices of the rod axis (i.e. indices of the basal plane)" << endl;
	   read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
	   ccin >> hh >> kk >> ll;
	   sizeEffect->SetRodAxis(hh,kk,ll);

	  // model parameters-set
	   string paramset;
	   cout << "choice of parameters set (DL, aD, AL)" << endl;
	   read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
	   ccin >> paramset;

	   boost::algorithm::to_lower(paramset);
	   if(paramset==string("dl") || paramset==string("d-l")) {
	     sizeEffect->SetModelParSet(MStruct::CircRodsGammaBroadeningEffect::PARAM_SET_DL);
	     REAL D, theta, L;
	     cout << "diameter (nm), theta, length (nm)" << endl;
	     read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
	     ccin >> D >> theta >> L;
	     sizeEffect->SetRodLength(L);
	     sizeEffect->SetRodDiameter(D,theta);
	   }
	   else if(paramset==string("ad") || paramset==string("a-d") || paramset==string("dlratio-d")) {
	     sizeEffect->SetModelParSet(MStruct::CircRodsGammaBroadeningEffect::PARAM_SET_aD);
	     REAL D, theta, a;
	     cout << "diameter (nm), theta, L/D ratio" << endl;
	     read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
	     ccin >> D >> theta >> a;
	     sizeEffect->SetRodDiameter(D,theta);
	     sizeEffect->SetRodShapePar(a);
	   }
	   else if(paramset==string("al") || paramset==string("a-l") || paramset==string("dlratio-l")) {
	     sizeEffect->SetModelParSet(MStruct::CircRodsGammaBroadeningEffect::PARAM_SET_aL);
	     REAL L, theta, a;
	     cout << "length (nm), theta, L/D ratio" << endl;
	     read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
	     ccin >> L >> theta >> a;
	     sizeEffect->SetRodLength(L,theta);
	     sizeEffect->SetRodShapePar(a);
	   }
	   else
	     throw ObjCrystException("Unknown model parameters-set option for CircRodsGamma broadening effect!");
	   
	   btype_not_found = false;
	 }
	// Size broadening from oval/elliptical rods with Gamma size-distribution
	 if(btype_not_found && (btype==string("EllipRodsGamma") || btype==string("OvalRodsGamma"))) {
	   MStruct::EllipRodsGammaBroadeningEffect * sizeEffect
	     = new MStruct::EllipRodsGammaBroadeningEffect;
	   sizeEffect->SetName(bname);
	   vReflProfComponents.push_back(sizeEffect);

	  // rod axis
	   REAL hh, kk, ll;
	   cout << "[hkl] fractional coordinates of the rod axis (direct space direction)" << endl;
	   read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
	   ccin >> hh >> kk >> ll;
	   sizeEffect->SetRodAxis(hh,kk,ll);
	   
	  // secondary axis
	   cout << "[hkl] fractional coordinates of the secondary axis" << endl;
	   read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
	   ccin >> hh >> kk >> ll;
	   sizeEffect->SetRodSecondaryAxis(hh,kk,ll);
	   
          // basal ellipse rotation
	   REAL psiD;
	   cout << "angle(deg) of rotation of basal ellipse (basal twist)" << endl;
	   read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
	   ccin >> psiD;
	   sizeEffect->SetRodBasalRotation(psiD*M_PI/180.);

	  // model parameters-set
	   string paramset;
	   cout << "choice of parameters set (DLf, aDf, aLf, DDL, aDD)" << endl;
	   read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
	   ccin >> paramset;

	   boost::algorithm::to_lower(paramset);
	   if(paramset==string("dlf") || paramset==string("d-l-f")) {
	     sizeEffect->SetModelParSet(MStruct::EllipRodsGammaBroadeningEffect::PARAM_SET_DLf);
	     REAL D, theta, L, f;
	     cout << "diameter (nm), theta, length (nm), flattening (0.-for circle)" << endl;
	     read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
	     ccin >> D >> theta >> L >> f;
	     sizeEffect->SetRodLength(L);
	     sizeEffect->SetRodMajorDiameter(D,theta);
	     sizeEffect->SetRodShapeFlattening(f);
	   }
	   else if(paramset==string("adf") || paramset==string("a-d-f") || paramset==string("dlratio-d-flatt")) {
	     sizeEffect->SetModelParSet(MStruct::EllipRodsGammaBroadeningEffect::PARAM_SET_aDf);
	     REAL D, theta, aL, f;
	     cout << "diameter (nm), theta, L/D ratio, flattening (0.-for circle)" << endl;
	     read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
	     ccin >> D >> theta >> aL >> f;
	     sizeEffect->SetRodMajorDiameter(D,theta);
	     sizeEffect->SetRodShapePar(aL);
	     sizeEffect->SetRodShapeFlattening(f);
	   }
	   else if(paramset==string("alf") || paramset==string("a-l-f") || paramset==string("dlratio-l-flatt")) {
	     sizeEffect->SetModelParSet(MStruct::EllipRodsGammaBroadeningEffect::PARAM_SET_aLf);
	     REAL L, theta, aL, f;
	     cout << "length (nm), theta, L/D ratio, flattening (0.-for circle)" << endl;
	     read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
	     ccin >> L >> theta >> aL >> f;
	     sizeEffect->SetRodLength(L,theta);
	     sizeEffect->SetRodShapePar(aL);
	     sizeEffect->SetRodShapeFlattening(f);
	   }
	   else if(paramset==string("ddl") || paramset==string("d-d-l") || paramset==string("da-db-l")) {
	     sizeEffect->SetModelParSet(MStruct::EllipRodsGammaBroadeningEffect::PARAM_SET_DDL);
	     REAL Da, theta, Db, L;
	     cout << "diameter-A (nm), theta, diameter-B (nm), length (nm)" << endl;
	     read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
	     ccin >> Da >> theta >> Db >> L;
	     sizeEffect->SetRodLength(L);
	     sizeEffect->SetRodDiameters(Da,Db,theta);
	   }
	   else if(paramset==string("add") || paramset==string("a-d-d") || paramset==string("dlratio-da-db")) {
	     sizeEffect->SetModelParSet(MStruct::EllipRodsGammaBroadeningEffect::PARAM_SET_aDD);
	     REAL Da, theta, Db, aL;
	     cout << "diameter-A (nm), theta, diameter-B (nm), L/D ratio" << endl;
	     read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
	     ccin >> Da >> theta >> Db >> aL;
	     sizeEffect->SetRodShapePar(aL);
	     sizeEffect->SetRodDiameters(Da,Db,theta);
	   }
	   else
	     throw ObjCrystException("Unknown model parameters-set option for EllipRodsGamma broadening effect!");
	   
	   btype_not_found = false;
	 }
	// Anisotropic microstrain broadening (simple model according to Popa, ref: TODO) 
	 if(btype_not_found && ( btype==string("StrainSimpleAniz") ) ) {
	   MStruct::StrainSimplePopaAniz * strainEffect
	     = new MStruct::StrainSimplePopaAniz;
	   strainEffect->SetName(bname);
	   vReflProfComponents.push_back(strainEffect);

	  // model anisotropy coefficients
	   int nbCoefs = 9;
	   CrystVector_REAL coefs(nbCoefs);
	   cout << "microstrain model coefficients (E1, E2, ...)" << endl;
	   read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
	   for(int k=0; k<nbCoefs; k++)
	     ccin >> coefs(k);
	  
          // model parameters
	   REAL ehhScale, eta0;
	   cout << "microstrain Scale, Eta0 (0...Gaussian, 1...Cauchy)" << endl;
	   read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
	   ccin >> ehhScale >> eta0;

	   strainEffect->SetPopaCoefs(coefs);
	   strainEffect->SetProfileParams(ehhScale, eta0);
	   
	   btype_not_found = false;
	 }
     	// Dislocation broadening dislocSLvN
     	 if(btype_not_found && ( btype==string("dislocSLvB") || btype==string("dislocSLvB+") ||
				 btype==string("dislocSLvBplus") || btype==string("dislocSLvBPlus"))) {
     	 	 	MStruct::DislocationBroadeningEffectSvB * strainEffect
       			 = new MStruct::DislocationBroadeningEffectSvB;
       	  strainEffect->SetName(bname);
       	  vReflProfComponents.push_back(strainEffect);
       	  
	 // use MWilk or Re and formula type for "-log(eta)"
	  int useMWilk = 0, formulaType = 0, argType = 0;
	  if(btype==string("dislocSLvB+") || btype==string("dislocSLvBplus") || btype==string("dislocSLvBPlus")) {
	    cout << "use MWilk instead of Re (0-No,1-Yes), formula (0-vanBerkum,1-fullWilkens,2-Kaganer), argument (0-x/Re,1-(b*g)*x/Re,2-(sqrt(C)*b*g)/Re" << endl;
	    read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
	    ccin >> useMWilk >> formulaType >> argType;
	  }

	  strainEffect->SetUseMWilk( useMWilk==1 );
	  strainEffect->SetChklChoiceABC( false );
	  strainEffect->SetFormula( formulaType , argType );

       	 // parameters
       	  REAL ReOrMWilk, Rou;
	  if( useMWilk==1 )
	    cout << "MWilk,rou(1/nm^2)" << endl;
	  else
	    cout << "Re(nm),rou(1/nm^2)" << endl;
	  read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
	  ccin >> ReOrMWilk >> Rou;

       	  REAL Cg0 = 1., Q1 = 0., Q2 = 0.;
       	  cout << "Cg0,q1,q2" << endl;
       	  read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
       	  ccin >> Cg0 >> Q1 >> Q2;
       	  strainEffect->SetProfilePar(ReOrMWilk,Rou,Cg0,Q1,Q2);
       	  //strainEffect->SetParentReflectionProfile(*reflProfile);
	  //reflProfile->AddReflectionProfileComponent(*strainEffect);
	  btype_not_found = false;
     	 }
     	// Deformation and Twin faults broadening for fcc material as desribed by Velterop2000 faultsVfcc
     	 if(btype_not_found && btype==string("faultsVfcc")) {
     	 	 	MStruct::FaultsBroadeningEffectVelteropFCC * faultsEffect
       			 = new MStruct::FaultsBroadeningEffectVelteropFCC;
       		faultsEffect->SetName(bname);
       		vReflProfComponents.push_back(faultsEffect);
       		
       	 // parameters
       	  REAL alpha, beta;
       	  cout << "alpha, beta" << endl;
       	  read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
       	  ccin >> alpha >> beta;
       	  
       	  faultsEffect->SetProfilePar(alpha,beta);
       	  //strainEffect->SetParentReflectionProfile(*reflProfile);
     			//reflProfile->AddReflectionProfileComponent(*strainEffect);
     	 		btype_not_found = false;
     	 }
     	// Deformation and Twin faults broadening for fcc material as desribed by Balogh&Ungar(JAP2006) faultsBfcc
     	 if(btype_not_found && btype==string("faultsBfcc")) {
     	 	 	MStruct::FaultsBroadeningEffectFCCBaloghUngar * faultsEffect
       			 = new MStruct::FaultsBroadeningEffectFCCBaloghUngar;
       		faultsEffect->SetName(bname);
       		vReflProfComponents.push_back(faultsEffect);
       		
       	 // parameters
       	  string type;
       	  REAL alpha;
       	  cout << "type(intrinsic/twins/extrinsic), alpha" << endl;
       	  read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
       	  ccin >> type >> alpha;
       	 
       	 // convert input type string to lowercase
       	  for_each( type.begin() , type.end() , char2lower() );
       	  
       	 // set profile parametrs
       	  if ( type==string("intrinsic") )
       	  		faultsEffect->SetProfilePar(MStruct::FaultsBroadeningEffectFCCBaloghUngar::INTRINSIC, alpha);
       	  else if ( type==string("twins") )
       	  		faultsEffect->SetProfilePar(MStruct::FaultsBroadeningEffectFCCBaloghUngar::TWINS, alpha);
       	  else if ( type==string("extrinsic") )
       	  		faultsEffect->SetProfilePar(MStruct::FaultsBroadeningEffectFCCBaloghUngar::EXTRINSIC, alpha);
       	  else {
       	  	cerr << "< Application: input argument mismatch\n";
						cerr << "\t"<<"Can not set the defect type ="<<type;
						cerr << " "<<"for object name: "<<bname;
						cerr << " and "<<"type: MStruct::FaultsBroadeningEffectFCCBaloghUngar.\n";
						cerr << "\t"<<"Unknown or unsupported defect type. >";
						throw ObjCrystException("Application: Wrong input argument.");
       	  }
     	 
       	  //strainEffect->SetParentReflectionProfile(*reflProfile);
     			//reflProfile->AddReflectionProfileComponent(*strainEffect);
     	 		btype_not_found = false;
     	 }

	// [11-23] displacement faults broadening in hexagonal WC faultsWC11m23
     	 if(btype_not_found && btype==string("faultsWC11m23")) {
     	 	 	MStruct::FaultsBroadeningEffectWC11m23 * faultsEffect
       			 = new MStruct::FaultsBroadeningEffectWC11m23;
       		faultsEffect->SetName(bname);
       		vReflProfComponents.push_back(faultsEffect);
	 // model options
	  string approx_type;
	  cout << "approximation type (simple,b_TPA_unaff)" << endl;
	  read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
       	  ccin >> approx_type;

	 // convert approximation type string to lowercase
       	  for_each( approx_type.begin() , approx_type.end() , char2lower() );
       	  
       	  if ( approx_type==string("simple") ) {
	    faultsEffect->SetModelApproxType(MStruct::FaultsBroadeningEffectWC11m23::APPROX_MODEL_SIMPLE);
	  } else if ( approx_type==string("b_tpa_unaff") || approx_type==string("b_ppt_unaff") ) {
	    faultsEffect->SetModelApproxType(MStruct::FaultsBroadeningEffectWC11m23::APPROX_MODEL_BROKEN_TPA_UNAFFECTED);
	  } else {
	    cerr << "< Application: unknown input option\n";
	    cerr << "\t"<<"Can not set model approximation type ="<<approx_type;
	    cerr << " "<<"for object name: "<<bname;
	    cerr << " and "<<"type: MStruct::FaultsBroadeningEffectWC11m23.\n";
	    cerr << "\t"<<"Unknown or unsupported model approximation type. >\n";
	    throw ObjCrystException("Application: Wrong input argument.");
	  }

       	 // parameters
       	  REAL alpha;
       	  cout << "alpha" << endl;
       	  read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
       	  ccin >> alpha;
       	  
       	  faultsEffect->SetProfilePar(alpha);
       	  //faultsEffect->SetParentReflectionProfile(*reflProfile);
     	  //reflProfile->AddReflectionProfileComponent(*faultsEffect);
     	  btype_not_found = false;
     	 }
     	// Phenomenological broadening
     	 if(btype_not_found && btype==string("HKLpVoigtA")) {
					if(vHKLEffect[iphase]==NULL)
     	 			vHKLEffect[iphase] = new MStruct::HKLPseudoVoigtBroadeningEffectA;
     	 		else {
     	 			cout << "Warning: Can not have multiple HKLpVoigtA effects for one phase!" << endl;
     	 			continue;
     	 		}
     	 		vHKLEffect[iphase]->SetName(bname);
     	 		vReflProfComponents.push_back(vHKLEffect[iphase]);
     	 		
     	 		//reflProfile->AddReflectionProfileComponent(*vHKLEffect[iphase]);
     	 		btype_not_found = false;
     	 }
     	// Common MStruct ReflectionProfile
     	 if(btype_not_found && btype==string("ReflProf")) {
					MStruct::ReflectionProfile * reflProfile 
												= new MStruct::ReflectionProfile(*crystal,data.GetRadiation());
					cout << "ReflProf (1):"<< reflProfile->GetName() << endl;
					reflProfile->SetName(bname);
					cout << "ReflProf (2):"<< reflProfile->GetName() << endl;
				 // add this ReflectionProfile to the list
				  vReflProfiles.push_back(reflProfile);
				 // number of additional broadening effects
				  int neffectnames = 0;
				  cout << "number of additional broadening effects (components only, instrumental effect already included)" << endl;
				  read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
				  ccin >> neffectnames;
				 // names of additional broadening effects
				 	vector< string > vCompNames;
				 	for(int ieffname=0; ieffname<neffectnames; ieffname++) {
				 		string compName;
				 		cout << "broadening component name" << endl;
				  	read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
				  	ccin >> compName;
				  	vCompNames.push_back(compName);
				 	}
				 // save components names
				 	vReflProfCompNames.push_back(vCompNames);
				 // flag if this ReflectionProfile is the the top-most parent ReflectionProfile for the given phase
				 	int parentFlag = -1;
				 	cout << "top parent effect (1-yes,0-no)" << endl;
				  read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
				  ccin >> parentFlag;
				  if (parentFlag>0) {
				  	if (parentReflProfNb<0)
				  		parentReflProfNb = vReflProfiles.size()-1;
				  	else {
				  		cerr << "< Application: input argument mismatch\n";
							cerr << "\t"<<"Can not set the ReflectionProfile Object="<<bname<<" as a top parent\n";
							cerr << "\t"<<"Reflection Profile="<<vReflProfiles[parentReflProfNb]->GetName()<<" already set";
							cerr << " as the top-most parent Reflection Profile >";
							throw ObjCrystException("Application: Wrong input argument.");
				  	}
				  }
				 // set some parameters specific for this ReflectionProfile Object
					reflProfile->SetParentPowderPatternDiffraction(*vDiffData[iphase]);
     			reflProfile->SetIncidenceAngle(omega);
     			cout << "ReflProf (3):"<< reflProfile->GetName() << endl;
     		 // 
     	 		btype_not_found = false;
     	 }
     	// Double component MStruct ReflectionProfile
     	 if(btype_not_found && btype==string("DoubleCompReflProf")) {
					MStruct::DoubleComponentReflectionProfile * reflProfile 
												= new MStruct::DoubleComponentReflectionProfile();
					reflProfile->SetName(bname);
				 // add this ReflectionProfile to the list
				  vReflProfiles.push_back(reflProfile);
				 // ReflectionProfile component names (have to be 2)
				 	vector< string > vCompNames;
				 	string compName;
				 	cout << "name of the first component" << endl;
				 	read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
				 	ccin >> compName;
				  vCompNames.push_back(compName);
				  cout << "name of the second component" << endl;
				 	read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
				 	ccin >> compName;
				 	vCompNames.push_back(compName);
				 // save components names
				 	vReflProfCompNames.push_back(vCompNames);
				 // parameters
				  REAL fraction = 0.;
				  cout << "fraction of the second component" << endl;
				 	read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
				 	ccin >> fraction;
				 	reflProfile->SetProfileParams(fraction);
				 // flag if this ReflectionProfile is the the top-most parent ReflectionProfile for the given phase
				 	int parentFlag = 0;
				 	cout << "top parent effect (1-yes,0-no)" << endl;
				  read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
				  ccin >> parentFlag;
				  if (parentFlag>0) {
				  	if (parentReflProfNb<0)
				  		parentReflProfNb = vReflProfiles.size()-1;
				  	else {
				  		cerr << "< Application: input argument mismatch.\n";
							cerr << "\t"<<"Can not set the ReflectionProfile Object="<<bname<<" as a top parent\n";
							cerr << "\t"<<"Reflection Profile="<<vReflProfiles[parentReflProfNb]->GetName()<<" already set";
							cerr << " as the top-most parent Reflection Profile >";
							throw ObjCrystException("Application: Wrong input argument.");
				  	}
				  }
     		 // 
     	 		btype_not_found = false;
     	 }
     	 // Residual Stress reflection position correction
	 if(btype_not_found && btype==string("StressSimple")) {
     	 	 	MStruct::ResidualStressPositionCorrection * stressCorr
       			 = new MStruct::ResidualStressPositionCorrection;
       		stressCorr->SetName(bname);
       		vReflProfComponents.push_back(stressCorr);
       		
       	 // parameters
       	 	string type;
       	  REAL stress;
       	  cout << "XECs model, stress (GPa)" << endl;
       	  read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
       	  ccin >> type >> stress;
       	  
       	 // convert XECs model type string to lowercase
       	  for_each( type.begin() , type.end() , char2lower() );
       	  
       	  if ( type==string("isotropic") ) {
       	   // isotropic XECs model
       	  	MStruct::XECsIsotropic * XECsModel = new MStruct::XECsIsotropic;
       	  	XECsModel->SetName(bname+string("XECs"));
       	   // read and set XECs model parameters (Young's modulus,Poisson's ratio)
       	    REAL E, ni;	
       	  	cout << "Young's modulus (GPa), Poisson's ratio" << endl;
       	  	read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
       	  	ccin >> E >> ni;
       	  	XECsModel->SetElasticParameters(E,ni);
       	   // set the XECs model
       	    stressCorr->SetXECsObj(*XECsModel);
       	  } else if ( type==string("reuss-voigt") ) {
       	   // anisotropic Reuss-Voigt XECs model
       	  	MStruct::XECsReussVoigt * XECsModel = new MStruct::XECsReussVoigt;
       	  	XECsModel->SetName(bname+string("XECs"));
       	   // set crystal for XECs Object
       	   	XECsModel->SetUnitCell(*crystal);
       	   // read single crystal stifness constants
       	     // get Stiffness constants names
							const vector< string > & vCijNames = XECsModel->GetStiffnessConstantsNames();
						 // number of requred Stiffness constants
							const int n = vCijNames.size();
							
						 // read  Stiffness constants values (as a map)
							map<string, REAL> mapCijValues;
							cout << "material ";
							for(int i=0; i<n; i++) cout << vCijNames[i] << " ";
							cout << "constants (in GPa) - in the format: C11 value C12 value etc." << endl;
							
							int npar = 0; bool newline = true;
							while (npar<n) {
								if (newline) { read_line (ccin, imp_file); newline = false; } // read a line (ignoring all comments, etc.)
								string name;
								ccin >> name;
								if ( !(name.length()==3 && (name[0]=='C' || name[0]=='c')) ) { newline = true; continue; }
								REAL value;
								ccin >> value;
    						for_each( name.begin() , name.end() , char2upper() ); // name to upper case
    						mapCijValues[name] = value;
								npar++;
							}
							
							cout << "Stiffness constants:" << endl;
							map< string, REAL>::const_iterator it = mapCijValues.begin();
							for(int i=0; i<n; i++) {
								cout << "\t" << it->first << ": " << it->second << " (GPa)" << endl;
								it++;
							}
							
						 // set material Stiffness constants
							XECsModel->SetStiffnessConstants(mapCijValues);
	
						 // get Cij matrix
							const CrystMatrix_REAL & matrixCij = XECsModel->GetCijMatrix();
	
							cout << "Stiffness tensor in the Voigt notation (GPa):" << endl;
							for(int i=0; i<6; i++) {
								cout << "\t";
								for(int j=0; j<6; j++)
									cout << setw(10) << matrixCij(i,j);
								cout << endl;
							}
							
       	   // read and set model parameters (Reuss-Voigt model weight)
       	    REAL weight;
						cout << "model weight (0..Reuss,1..Voigt)" << endl;
						read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
						ccin >> weight;
						XECsModel->SetParams(weight);
       	   // set the XECs model
       	    stressCorr->SetXECsObj(*XECsModel);
       	  } else {
       	  	cerr << "< Application: input argument mismatch\n";
						cerr << "\t"<<"Can not set the XECs model ="<<type;
						cerr << " "<<"for object name: "<<bname;
						cerr << " and "<<"type: MStruct::ResidualStressPositionCorrection.\n";
						cerr << "\t"<<"Unknown or unsupported XECs model. >";
						throw ObjCrystException("Application: Wrong input argument.");
       	  }
       	  
       	  stressCorr->SetParams(stress);
       	  //stressCorr->SetParentReflectionProfile(*reflProfile);
     			//reflProfile->AddReflectionProfileComponent(*stressCorr);
     	 		btype_not_found = false;
     	 }
     	// Refraction reflection position correction
     	 if(btype_not_found && btype==string("RefractionCorr")) {
     	 		MStruct::RefractionPositionCorr * refractionCorr = new MStruct::RefractionPositionCorr;
       		refractionCorr->SetName(bname);
       		vReflProfComponents.push_back(refractionCorr);
       		
       	 // parameters
       		string type;
       	  cout << "chi0 set directly-value,calculated from-crystal,calculated from chem.-formula" << endl;
       	  read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
       	  ccin >> type;
       	  
       	 // convert choice type string to lowercase
       	  for_each( type.begin() , type.end() , char2lower() );
       	  
       	  if ( type==string("value") ) {
       	  	cout << "chi0 real, imag" << endl;
       	  	read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
       	  	REAL rr, im;
       	  	ccin >> rr >> im;
       	  	refractionCorr->SetChi0(complex<REAL>(rr,im));
       	  }
       	  else if ( type==string("crystal") ) {
       	  	refractionCorr->SetCrystal(*crystal,true);
       	  } else if ( type==string("formula") ) {
       	  	
       	  	string formula;
       	  	REAL density;
       	  	
						cout << "chemical formula, density (g/cm3)" << endl;
						read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
						// get chemical formula from the input (line) - following strings are considered as separotors
						static const char* separators[23] = { " 0", " 1", " 2", " 3", " 4", " 5", " 6", " 7", " 8", " 9", " .", 
																						   		"\t0", "\t1", "\t2", "\t3", "\t4", "\t5", "\t6", "\t7", "\t8", "\t9", "\t.",
																								  "\n"};
						// remove leading blank spaces
						string::size_type ind;
						ind = ccin.str().find_first_not_of(" \t\n");
						if(ind!=string::npos) ccin.str( ccin.str().substr(ind) );
						// find formula-density separator
						for(int i=0; i<23; i++) {
							ind = ccin.str().find(separators[i]);
							if( ind != string::npos ) break;
						}
						formula = ccin.str().substr(0,ind); // use line up to the found separator for formula
						ccin.str( ccin.str().substr(ind) ); // cut formula part from the line
       	  	ccin >> density; // read density value
       	  	
       	  	refractionCorr->SetChemicalFormula(formula,density);
       	  } 
       	  
       	  cout << "relative density" << endl;
       	  read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
       	  REAL density;
       	  ccin >> density;
       	  refractionCorr->SetParams(density);
       	  
     	 		btype_not_found = false;
     	 }
	 // Geometrical broadening in parallel beam geometry with finite irradated area and slit acceptance
     	 if(btype_not_found && btype==string("instrGeomPBslit")) {
	   MStruct::ProfilePerspectiveA * perspectiveEffect = new MStruct::ProfilePerspectiveA;
	   perspectiveEffect->SetName(bname);
	   vReflProfComponents.push_back(perspectiveEffect);
       		
	   // parameters
	   REAL gamma, eta;
	   cout << "gamma(deg), eta (0/1...Gaussian/Loretzian)" << endl;
	   read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
	   ccin >> gamma >> eta;

	   // set parameters
	   perspectiveEffect->SetProfilePar(gamma, eta);

	   btype_not_found = false;
	 } // InstrGeomPBslit


	 
 			// Type of the given broadening not found
 			 if(btype_not_found) {
 			 		cout << "Warning: Type of the given broadening effect not recognised!" << endl;
 			 		continue;
 			 }
     } // ibeffect
    
    // Analyse the user given list of ReflectionProfiles and their Components to build a broadening model.
    
    // If no ReflectionProfile Object is specified by the user then the common MStruct::ReflectionProfile
    // is created automaticaly and all components in the user components list are included.
     if (vReflProfiles.empty() || parentReflProfNb<0) {
     	// ReflectionProfile Object was not set by user. Create and set it by default.
     	// create new ReflectionProfile Object
     	 MStruct::ReflectionProfile * reflProfile
     	 							= new MStruct::ReflectionProfile(*crystal,data.GetRadiation());
     	// set its specific parameters
     	 reflProfile->SetParentPowderPatternDiffraction(*vDiffData[iphase]);
     	 reflProfile->SetIncidenceAngle(omega);
     	// save it for the further processing
     	 vReflProfiles.push_back(reflProfile);
     	// set it as the parent reflection profile
     	 parentReflProfNb = vReflProfiles.size()-1;
     	// add to it all user supplied ReflectionProfileComponent objects (by default)
     	 vector< string > vCompNames;
     	 for(unsigned int ibeffect=0; ibeffect<vReflProfComponents.size(); ibeffect++)
     	 	 vCompNames.push_back(vReflProfComponents[ibeffect]->GetName());
     	 vReflProfCompNames.push_back(vCompNames);
     }
    
    // All ReflectionProfile and ReflectionProfileComponent Objects should be allocated now.
    // Just loop through the ReflectionProfile list and connect the Objects together.   
    
    // names of ReflectionProfile and ReflectionProfileComponent Objects will be helpful, becouse they will
    // be searched - prepare them
     vector<string> vReflProfNamesIndex (vReflProfiles.size());
     vector<string> vReflProfCompNamesIndex (vReflProfComponents.size());
     transform ( vReflProfiles.begin(), vReflProfiles.end(), vReflProfNamesIndex.begin(),
                                    	 	 			mem_fun (&ObjCryst::ReflectionProfile::GetName) );
     transform ( vReflProfComponents.begin(), vReflProfComponents.end(), vReflProfCompNamesIndex.begin(),
                                    	 	 			mem_fun (&MStruct::ReflectionProfileComponent::GetName) );
     {
     	 cout << endl;
     	 cout << "ReflectionProfile Objects aviable:\n";
     	 vector<string>::iterator iter;
     	 for(iter = vReflProfNamesIndex.begin(); iter != vReflProfNamesIndex.end(); iter++)
     	 	 cout << *iter << "\n";
     	 cout << "ReflectionProfileComponent Objects aviable:\n";
     	 for(iter = vReflProfCompNamesIndex.begin(); iter != vReflProfCompNamesIndex.end(); iter++)
     	 	 cout << *iter << "\n";
     	 cout << endl;
     }
     
     for(unsigned int ireflProf=0; ireflProf<vReflProfiles.size(); ireflProf++)
     {
     	 ObjCryst::ReflectionProfile * reflProfile = vReflProfiles[ireflProf];
    	// ReflectionProfile type has not been recognised yet
     	 bool type_not_found = true;
     	 
     	// common MStruct::ReflectionProfile
     	 if(type_not_found && reflProfile->GetClassName()==string("MStruct::ReflectionProfile")) {
     	 	// retype to the (right) inheritor class
     	 	 MStruct::ReflectionProfile * reflProfile =
     	 	 								dynamic_cast<MStruct::ReflectionProfile*>(vReflProfiles[ireflProf]);
     	 	// find and add all components for this reflection profile prescribed by the user
     	 	//   + instrumental profile (first)
     	 	 if( instReflProfileComp != NULL ) {
     	 	 	 MStruct::ReflectionProfileComponent * reflProfileComp = instReflProfileComp;
     	 	 	 reflProfileComp->SetParentReflectionProfile(*reflProfile);
     	 		 reflProfile->AddReflectionProfileComponent(*reflProfileComp);
     	 	 }
     	 	// component names for this ReflectionProfile Object
     	 	 vector<string> &vCompNames = vReflProfCompNames[ireflProf];			 
     	 	// loop over all ReflectionProfileComponents objects of this ReflectionProfile
     	 	 for(unsigned int ibeffect=0; ibeffect<vCompNames.size(); ibeffect++) {
     	 	 	// find a component with the given name
     	 	 	 MStruct::ReflectionProfileComponent * reflProfileComp = NULL;
     	 	 	 vector<string>::iterator iter = find_if ( vReflProfCompNamesIndex.begin(), vReflProfCompNamesIndex.end(),
     	 	 	 																						bind2nd( equal_to<string>(), vCompNames[ibeffect] ) );
     	 	 	 if( iter == vReflProfCompNamesIndex.end() ) {
     				 cerr << "< Application: input argument mismatch.\n";
						 cerr << "\t"<<"No ReflectionProfileComponent Object with name="<<vCompNames[ibeffect]<<" found\n";
						 cerr << "\t"<<"The Object is required by the ReflectionProfile Object name="<<reflProfile->GetName()<<"\n";
						 cerr << "\t"<<"Profile broadening model for Crystal Object name="<<crystal->GetName()<<" can not be created >\n";  
						 throw ObjCrystException("Application: Wrong input argument.");
   				 } else {
     				 reflProfileComp = *(vReflProfComponents.begin() + (iter-vReflProfCompNamesIndex.begin()));
   				 }
   			  // if a component found than add the component to this ReflectionProfile
   			   if( reflProfileComp != NULL ) {
     	 	 	 	 reflProfileComp->SetParentReflectionProfile(*reflProfile);
     	 		 	 reflProfile->AddReflectionProfileComponent(*reflProfileComp);
     	 	 	 }
     	 	 } // end of loop: for(int ibeffect=0, ...)
     	 	// 
     	 	 type_not_found = false;
     	 }
    
     	// Double component MStruct ReflectionProfile
     	 if(type_not_found && reflProfile->GetClassName()==string("MStruct::DoubleComponentReflectionProfile")) {
     	 	// retype to the (right) inheritor class
     	 	 MStruct::DoubleComponentReflectionProfile * reflProfile =
     	 	 								dynamic_cast<MStruct::DoubleComponentReflectionProfile*>(vReflProfiles[ireflProf]);
     	 	// find both components for this reflection profile prescribed by the user
     	 	 ObjCryst::ReflectionProfile * reflProfile1 = NULL;
     	 	 ObjCryst::ReflectionProfile * reflProfile2 = NULL;
     	 	 vector<string>::iterator iter = find_if ( vReflProfNamesIndex.begin(), vReflProfNamesIndex.end(),
     	 	  																				 bind2nd( equal_to<string>(), vReflProfCompNames[ireflProf][0] ) );
     	 	 if( iter == vReflProfNamesIndex.end() ) {
     			 cerr << "< Application: input argument mismatch\n";
					 cerr << "\t"<<"No ReflectionProfile Object with name="<<vReflProfCompNames[ireflProf][0]<<" found\n";
					 cerr << "\t"<<"The Object is required by the ReflectionProfile Object name="<<reflProfile->GetName()<<"\n";
					 cerr << "\t"<<"Profile broadening model for Crystal Object name="<<crystal->GetName()<<" can not be created >\n";  
					 throw ObjCrystException("Application: Wrong input argument.");
   			 } else {
     			 reflProfile1 = *(vReflProfiles.begin() + (iter-vReflProfNamesIndex.begin()));
   			 }
     	 	                          iter = find_if ( vReflProfNamesIndex.begin(), vReflProfNamesIndex.end(),
     	 	  																				 bind2nd( equal_to<string>(), vReflProfCompNames[ireflProf][1] ) );
     	 	 if( iter == vReflProfNamesIndex.end() ) {
     			 cerr << "< Application: input argument mismatch\n";
					 cerr << "\t"<<"No ReflectionProfile Object with name="<<vReflProfCompNames[ireflProf][1]<<" found\n";
					 cerr << "\t"<<"The Object is required by the ReflectionProfile Object name="<<reflProfile->GetName()<<"\n";
					 cerr << "\t"<<"Profile broadening model for Crystal Object name="<<crystal->GetName()<<" can not be created >\n";  
					 throw ObjCrystException("Application: Wrong input argument.");
   			 } else {
     			 reflProfile2 = *(vReflProfiles.begin() + (iter-vReflProfNamesIndex.begin()));
   			 }
   			// set ReflectionProfile sub-components
   			 reflProfile->SetComponents(*reflProfile1, *reflProfile2);
   			// 
     	 	 type_not_found = false;
     	 }
     	     	
     	// Type of the given ReflectionProfile Object not found
 			 if(type_not_found) {
 			 		cerr << "Warning: Type of the given broadening effect not recognised!" << endl;
 			 		continue;
 			 }
     }
  
    // set profile 
     vDiffData[iphase]->SetProfile(vReflProfiles[parentReflProfNb]);
    
   } // iphase
   
  // Prepare data
   data.Prepare();
   data.FitScaleFactorForRw();

   if (1) {
  // load IhklCorr params from files
  for(int iphase=0; iphase<nbphases; iphase++)
   if(vHKLChoice[iphase]>0)
     vDiffData[iphase]->ReadHKLIntensityCorrFromFile(string(string("Ihkl_")+vPhasesNames[iphase]+string(".dat")).c_str());

  // HKL peak profile params
   for(int iphase=0; iphase<nbphases; iphase++) {
     if(vHKLEffect[iphase]==NULL) continue;
     cout << "(nb of peaks with preset profiles)/(name of parameters file, flag, intensity) - phase: " << vPhasesNames[iphase] << endl;
     string str;
     read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
     ccin >> str; // read first word on the line (either number or string)
     char *endptr;
     long nbprofs = strtol(str.c_str(),&endptr,10);
     if((*endptr=='\0') && nbprofs>=0) { // first word is a nonegative number
       // read hkl diffraction params
       for(int iprof=0; iprof<nbprofs; iprof++) {
	 cout << "hkl,d2Theta(deg),fwhm(deg),eta,code(111)?" << endl;
	 REAL dx, fwhm, eta;
	 int h,k,l,code;
	 read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
	 ccin >> h >> k >> l >> dx >> fwhm >> eta >> code;
	 int b1=code/100, b2=(code-100*b1)/10, b3=(code-100*b1-10*b2);
	 vHKLEffect[iphase]->SetProfilePar(h,k,l,dx*DEG2RAD,fwhm*DEG2RAD,eta,1.0,b1==0,b2==0,b3==0,true);
       }
     } else { // first word is a string (or negative integer)
       vHKLEffect[iphase]->SetParametersFile( str.c_str() );
       int flag; // (0-don't use,1-generate,2-free strong,3-read)
       REAL relIntensity; // relative intensity level for generate option
       ccin >> flag >> relIntensity;
       if(flag==3) // option == read
	 vHKLEffect[iphase]->LoadParametersFile();
     }
   }
  
  // re-load all Size Distribututions for all effects (unfortunatelly distributions are earlier automatically fixed)
   for(unsigned int ieffect=0; ieffect<vSizeDistribEffect.size(); ieffect++)
   	 vSizeDistribEffect[ieffect]->ReadDistributionFromFile(vSizeDistribFileName[ieffect].c_str());

	// register all Size Distribututions as additional LSQ func. objects
	/*for(unsigned int ieffect=0; ieffect<vSizeDistribEffect.size(); ieffect++)
	  data.AddAdditionalLSQObj(*vSizeDistribEffect[ieffect]);*/ // hack
   	 
  // Set refinement status of background points and parameters
   for(int icomp=0; icomp<(int)v_bkg_type_ids.size(); icomp++)
   	 switch (v_bkg_type_ids[icomp]) {
   	 	case BACKGROUND_INTERPOLATED:
   	 	  {
   	 	   	const int nb_bkg_points_ref = v_bkg_params_flags[icomp].numElements();
   	 	   	const CrystVector_REAL &bkg_points_ref = v_bkg_params_flags[icomp];
   	 	   	
   	 	   	cout << "Interpolated background component: " << v_bkg_strings[icomp] << endl;
   	 	   	cout << "Number of background points set to be refined: " << nb_bkg_points_ref << endl;
   			if(nb_bkg_points_ref > 0) {
     		  for(int i=0; i<nb_bkg_points_ref; i++) {
       			string str; // param name
       		   // creating param name
       			ostringstream os(""); os << bkg_points_ref(i);
       			str = "Background_Point_" + os.str(); 
       			cout << "Setting param. " << str << " fixed=false" << endl;
       		  // set param
       		  //RefinablePar &param = lsqOptObj.GetFullRefinableObj().GetPar(str.c_str());
       			RefinablePar &param = get_par(str, *vBackgroundComponents[icomp]);
       			param.SetIsFixed(false);
       			param.Print();
     		  }
   			}
   	 	   }
   	 	   break;
		case BACKGROUND_CHEBYSHEV:
   	 	   {
   	 	   	 const int nb_coefs = v_bkg_params[icomp].numElements();
   	 	   	 const CrystVector_REAL &bkg_coefs = v_bkg_params[icomp];
   	 	   	 const CrystVector_REAL &bkg_coef_flags = v_bkg_params_flags[icomp];
   	 	   	 
   	 	   	 cout << "Chebyshev background component: " << vBackgroundComponents[icomp]->GetName() << endl;
   	 	   	 cout << "Chebyshev polynomials coefficients: " << endl;
   				 if(nb_coefs > 0) {
     			 	 for(int i=0; i<nb_coefs; i++) {
       				   string str; // param name
       				  // creating param name
       				   ostringstream os(""); os << i;
       				   str = "Background_Coef_" + os.str();
       				   bool bfixed = bkg_coef_flags(i)>0;
       				   cout << "Setting param. " << str << " fixed="<< bfixed << endl;
       				 // set param
       				  //RefinablePar &param = lsqOptObj.GetFullRefinableObj().GetPar(str.c_str());
       				   RefinablePar &param = get_par(str, *vBackgroundComponents[icomp]);
       				   param.SetIsFixed(bfixed);
       				   param.Print();
     				 }
   				 }
   	 	   }
   	 	   break;
	          case BACKGROUND_CHEBYSHEV_LOCAL:
   	 	   {
		     cout << "Local Chebyshev background component: " << vBackgroundComponents[icomp]->GetName() << endl;
		     cout << "This background component may not have any parameters.\n";
		     cout << "Use @local_bkg_chebyshev tag to define new local backgrounds." << endl;	 
   	 	   }
   	 	   break;
   	 	  case BACKGROUND_INVX:
   	 	   {
   	 	   	 cout << "InvX background component: " << vBackgroundComponents[icomp]->GetName() << endl;
   	 	   	 cout << "This background component has no parameters, but it is scalable indeed.\n";
   	 	   	 cout << "Use '" << string("Scale_"+vBackgroundComponents[icomp]->GetName()) << "' to adjust/fit data." << endl;
   	 	   }
   	 	   break;
	         case BACKGROUND_TURBOSTRATIC_CARBON:
		   {
		     cout << "tubostraticCarbonWB background component: " << vBackgroundComponents[icomp]->GetName() << endl;
		     cout << "This background component is scalable and has also other parameters. Use standard methods to set/refine them.\n";
		     cout << "Use '" << string("Scale_"+vBackgroundComponents[icomp]->GetName()) << "' to adjust/fit data." << endl;
		   }
		   break;
   	 	 default:
   	 	 	 cerr << "< main(...)\n";
				 cerr << "Program logical error during setting refinement status of background points.\n >" << endl; 
				 throw ObjCrystException("main(...): Program error.");
   	 	   break;
   	 }

     // Create the LSQ optimization object
   MStruct::LSQNumObj lsqOptObj("lsqOptObj");
   lsqOptObj.SetRefinedObj(data,0,true,true);
   lsqOptObj.AddRefinableObj(data);

   lsqOptObj.ObjCryst::LSQNumObj::PrepareRefParList(false);
   RefinableObj & lsqCompiledObj = lsqOptObj.GetCompiledRefinedObj();
   lsqCompiledObj.PrepareForRefinement();

   //lsqOptObj.FixAllPar();

	// for grid refinement (job=1) define a grid
   CrystVector_REAL gridx(0);
   CrystVector_REAL gridy(0);
   CrystMatrix_REAL gridRw(0,0);
	 string grid_xname("");
	 string grid_yname("");
	   
	 if(job_type==1) {
	 	 // refinement grid
	 	 REAL grid_xmin, grid_xmax, grid_xstep, grid_ymin, grid_ymax, grid_ystep;
	 	 int grid_dim, grid_xn, grid_yn;
	 	 
	   cout << "refinement grid dimensions (1 or 2)" << endl;
	   read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
	   ccin >> grid_dim;
	   
	   {
	     cout << "grid variable (x): name, min, max, step" << endl;
	     read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
	     ccin >> grid_xname >> grid_xmin >> grid_xmax >> grid_xstep;

	     grid_xn = int(1.01*(grid_xmax-grid_xmin)/grid_xstep)+1;
	     gridx.resize(grid_xn);
	     for(int ix=0; ix<grid_xn; ix++) {
	       gridx(ix) = grid_xmin + ix*grid_xstep;
	     }

	     gridRw.resize(1,grid_xn);
	   }

	   if(grid_dim==2) {
	     cout << "grid variable (y): name, min, max, step" << endl;
	     read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
	     ccin >> grid_yname >> grid_ymin >> grid_ymax >> grid_ystep;

	     grid_yn = int(1.01*(grid_ymax-grid_ymin)/grid_ystep)+1;
	     gridy.resize(grid_yn);
	     for(int iy=0; iy<grid_yn; iy++) {
	       gridy(iy) = grid_ymin + iy*grid_ystep;
	     }
	     
	     gridRw.resize(grid_yn,grid_xn);
	   }
	   
		 /*if(dim>1) {
		 	 cout << "grid variable: name, min, max, step" << endl;
	   	 read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
	   	 ccin >> grid_yname >> grid_ymin >> grid_ymax >> grid_ystep;
	   	 grid_yn = int((grid_ymax-grid_ymin)/grid_ystep)+1;
		 }*/
	 }
   string filename; // output filename
   int niter=0, nparams=0, niter1=0, niter2=0;
   REAL lsqTuneCoef=-1.;

  // Set PowderPattern params
   cout << "output filename" << endl;
   read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
   ccin >> filename;
   
   if(job_type!=3) {
     cout << "nb of interactions" << endl;
     read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
     ccin >> niter;
   } else {
     // SizeDistrib refinement
     cout << "nb of reduced,nb of full interactions,nb of cycles,tuning coef." << endl;
     read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
     ccin >> niter >> niter1 >> niter2 >> lsqTuneCoef;
   }

   cout << "nb of params to be set" << endl;
   read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
   ccin >> nparams;

   for(int i=0;i<nparams;i++) {
     // read param
     string str; // param name
     REAL scale, val, step; // scale, value, deriv step
     int limited; // limited
     REAL min, max;

     cout << "param name" << endl;
     read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
     ccin >> str;

     cout << "human scale" << endl;
     read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
     ccin >> scale;

     cout << "value, derivative step" << endl;
     read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
     ccin >> val >> step;

     cout << "limited (1,0), min, max" << endl;
     read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
     ccin >> limited; if (limited==1) ccin >> min >> max;

     // set param
     try {
     RefinablePar &param = (str.at(0)=='#') ?
       lsqCompiledObj.GetPar(atoi(str.substr(1).c_str())) :
       get_par(str, lsqOptObj.GetCompiledRefinedObj());
     param.Print();
     param.SetIsLimited(limited==1); // do first to allow user change default limits before setting value
     if (limited==1) {
       param.SetMin(min*scale);
       param.SetMax(max*scale);
     }
     param.SetValue(val*scale); // do first becouse of relative derivative step
     param.SetIsFixed(step<=0.); if (step>0.) param.SetDerivStep(step*scale);
     param.Print();
     } // try
		 catch (...) {
		 		cout << "Warning: Parameter: " << str << " not found and ignored!" << endl;
		 }
   }

  // Load additional Objects/Options/Constraints initiated by @-symbol
   if(!using_file)
     cout << "insert additional Objects/Options/Constraints initiated by @-symbol or type '@end'" << endl;

   bool need_refParList_update = false;
   read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)

   while( ccin.str().length()>0 ) {
     
     if( ccin.str()[0] != '@' ) {
       read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
       continue;
     }
      
     string keyword;
     ccin >> keyword;

     // @end
     if( keyword=="@end" || keyword=="@END" || keyword=="@End" )
       break;

     // @LSQConstraint
     if( strncmp(keyword.c_str(),"@LSQConstraint",14)==0 ) {

       std::string name = ( keyword.length()>16 && keyword.at(14)==':' ) ? keyword.substr(15) : "unnamed";

       int nb;
       ccin >> nb; // number of constrained params
       std::vector< const RefinablePar* > parList;
       CrystVector_REAL coef;

       // params list
       read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)

       parList.resize(nb);
       for(int i=0; i<nb; i++) {

	 std::string str;
	 ccin >> str;

	 // find param
	 try {
	   parList[i] = &( (str.at(0)=='#') ?
			   lsqCompiledObj.GetPar(atoi(str.substr(1).c_str())) :
			   get_par(str, lsqOptObj.GetCompiledRefinedObj()) );
	 } catch (...) {
	   parList[i] = NULL;
	   cout << "Warning: Parameter: " << str << " not found and ignored!" << endl;
	 }
       } // for(int i=0; i<nb; ...
       
       // constraint coeficients
       read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)

       coef.resize(nb);
       for(int i=0; i<nb; i++) {
	 ccin >> coef(i);
       }
       
       // add the constraint
       lsqCompiledObj.AddExternalLSQConstraint( parList, coef, name );
       need_refParList_update = true;

     } // if @LSQConstraint

     // @local_bkg_chebyshev
     if( strncmp(keyword.c_str(),"@local_bkg_chebyshev",19)==0 ) {
       // not implemented
     } // if @local_bkg_chebyshev

     // @PowderPattern:SetWavelength
     if( strncmp(keyword.c_str(),"@PowderPattern:SetWavelength",28)==0 ) {
       // change radiation for PowderPattern Object (string and ratio)
       // usage: @PowderPattern:SetWavelength Cu 0.13
       std::string str;
       REAL ratio;
       ccin >> str;
       ccin >> ratio;
       data.SetWavelength( str.c_str(), ratio );

       cout << keyword.c_str() << " ... " << data.GetName() << ":" << data.GetRadiation().GetClassName() << " is now ...\n";
       data.GetRadiation().XMLOutput( cout, 0);
       cout << "\n";
     } // if @PowderPattern:SetWavelength

    // @DislocationBroadeningEffectSvB:SetChklChoiceABC
     if( strncmp(keyword.c_str(),"@DislocationBroadeningEffectSvB:SetChklChoiceABC",48)==0 ) {
       // set parameters type (Cg0,q1,12 vs. A,B,C) for Chkl in DislocationBroadeningEffectSvB
       // usage: @DislocationBroadeningEffectSvB:SetChklChoiceABC strainProfWC true
       std::string effname; // effect name
       std::string choice; // choice
       ccin >> effname;
       ccin >> choice;
       // find dislocation broadening effect object with given name
       MStruct::DislocationBroadeningEffectSvB *dislEff = NULL;
       {
	 RefinableObj &obj = gRefinableObjRegistry.GetObj(effname.c_str(),"MStruct::DislocationBroadeningEffectSvB");
	 dislEff = &(dynamic_cast<MStruct::DislocationBroadeningEffectSvB&>(obj));
       }
       if( strncmp(choice.c_str(),"1",1)==0 || strncmp(choice.c_str(),"true",4)==0 ) {
	 cout << "Option: Setting Chkl choice to ABC for " << dislEff->GetName() << "\n";
	 dislEff->SetChklChoiceABC(true);
       } else {
	 cout << "Option: Setting Chkl choice to Cg0,q1,q2 for " << dislEff->GetName() << "\n";
	 dislEff->SetChklChoiceABC(false);
       }
     } // if @DislocationBroadeningEffectSvB:SetChklChoiceABC

    // @PowderPatternDiffraction:SetReflProfCalcParams
     if( strncmp(keyword.c_str(),"@PowderPatternDiffraction:SetReflProfCalcParams",47)==0 ) {
       // set parameters for reflection profile calculations (minRelativeIntensity and factor)
       // usage: @PowderPatternDiffraction:SetReflProfCalcParams diffObjectName 0.001  2.0
       std::string name; // phase name
       REAL minRelIntensity, factor;
       ccin >> name >> minRelIntensity >> factor;
       // find PowderPatternDiffraction object with given name
       MStruct::PowderPatternDiffraction *diffObj = NULL;
       try {
	 RefinableObj &obj = gRefinableObjRegistry.GetObj(name.c_str(),"MStruct::PowderPatternDiffraction");
	 diffObj = &(dynamic_cast<MStruct::PowderPatternDiffraction&>(obj));
       }
       catch (exception& e) {
	 cerr << "< Application: Can not find MStruct::PowderPatternDiffraction with name: " << name << "\n";
	 throw ObjCrystException("Application: Wrong option parameter.");
       }
       if (diffObj) {
	 cout << "Option: setting minRelIntensity = " << minRelIntensity << " and factor = " << factor;
	 cout << " for phase: " << diffObj->GetName() << "\n";
	 diffObj->SetReflProfCalcParams(minRelIntensity, factor);
       }
       
     } // if @PowderPattern:SetReflProfCalcParams

    // @overwrite_bgr_file
     if( strncmp(keyword.c_str(),"@overwrite_bgr_file",19)==0 ) {
       settings_overwrite_bgr = true;
     } // if @overwrite_bgr_file
     
     read_line (ccin, imp_file); // read a line (ignoring all comments, etc.)
   }
   
   // update Refined parameters list to include external constraints, etc.
   if(need_refParList_update)
     lsqOptObj.ObjCryst::LSQNumObj::PrepareRefParList(false);

  // Save the powder pattern in text format
   //data.FitScaleFactorForRw(); (don't refine scale factors)
   if( lsqCompiledObj.GetNbLSQConstraints() > 0 )
     lsqCompiledObj.PrintExternalConstraintsStatistics();
   REAL rw = data.GetRw();
   if (0){ // sV
     string::size_type loc = filename.find_last_of('.');
     string s = (loc==string::npos) ?
       filename+"-i" : filename.substr(0,loc)+"-i"+filename.substr(loc);
     data.SavePowderPattern(s.c_str());
   }
   cout << "Initial data (before fit-for analysis):" << endl;
   cout << "finished cycle #0/0. Rw="<<rw<<"->"<<rw<<endl;
   {
     const RefinableObj &obj = lsqCompiledObj;
     lsqCompiledObj.PrepareForRefinement();

     for(int i=0;i<obj.GetNbParNotFixed();i++) {
       const RefinablePar &par = obj.GetParNotFixed(i);
       cout << setw(29) << left << par.GetName() << right;
       REAL val = par.GetHumanValue();
       cout << setw(18) << val << setw(18) << val;
       cout << setw(18) << val << setw(18) << 0. << endl;;
     }
	 }
   
   bool useLevenbergMarquardt=true;
   bool silent=false;
   
   if(job_type==1) { // grid refinement
		// get param(s)
		 RefinablePar *pparamx = 0;
		 RefinablePar *pparamy = 0;
		 try {
		 	 string str = grid_xname;
			 pparamx = (str.at(0)=='#') ?
					&lsqOptObj.GetFullRefinableObj().GetPar(atoi(str.substr(1).c_str())) :
					&get_par(str, lsqOptObj.GetFullRefinableObj());
			 //lsqOptObj.GetFullRefinableObj().GetPar(str.c_str());
			 pparamx->SetIsFixed(true);
		 }
		 catch (...) {
			 cout << "Warning: Parameter: " << grid_xname << " not found and ignored!" << endl;
		 }
		 if(gridy.numElements()>0) {
		 try {
		 	 string str = grid_yname;
			 pparamy = (str.at(0)=='#') ?
					&lsqOptObj.GetFullRefinableObj().GetPar(atoi(str.substr(1).c_str())) :
					&get_par(str, lsqOptObj.GetFullRefinableObj());
			 //lsqOptObj.GetFullRefinableObj().GetPar(str.c_str());
			 pparamy->SetIsFixed(true);
		 }
		 catch (...) {
			 cout << "Warning: Parameter: " << grid_yname << " not found and ignored!" << endl;
		 }
		 }
		 if (pparamx!=0) {
		 	// the best configuration
		         REAL bestRw = numeric_limits<REAL>::max();
		 	 int best_set_id = lsqOptObj.GetFullRefinableObj().CreateParamSet("bestParamsSet");
		 	 int ibestvalx = -1;
			 int ibestvaly = -1;
		 	// save starting params. set
		 	 const int start_set_id = lsqOptObj.GetFullRefinableObj().CreateParamSet("startingParamsSet");

			 if(gridy.numElements()==0) {
			   for(int ix=0; ix<gridx.numElements(); ix++) {
			    // restore starting parmas. set
			     lsqOptObj.GetFullRefinableObj().RestoreParamSet(start_set_id);
			    // set param
			     pparamx->SetValue(gridx(ix)/pparamx->GetHumanScale());
			    // run the refinement
			     lsqOptObj.Refine(niter,useLevenbergMarquardt,silent);
			    // store the Rw value
			     gridRw(0,ix) = lsqOptObj.RwFactor();
			    // save the best configuration
			     if (gridRw(ix)<=bestRw) { // store the best configuration
			       bestRw = gridRw(0,ix);
			       lsqOptObj.GetFullRefinableObj().SaveParamSet(best_set_id);
			       ibestvalx = ix;
			     }
			   }
			  // print grid refinement results
			   cout << " ------------------------------------------------ " << endl;
			   cout << "GRID REFINEMENT results:" << endl;
			   cout << setw(25) << grid_xname << setw(15) << "Rw" << endl;
			   for(int ix=0; ix<gridx.numElements(); ix++)
			     cout << setw(25) << gridx(ix) << setw(15) << gridRw(0,ix) << endl;
			   if (ibestvalx>=0) {
			     cout << "Best parameters set: " << grid_xname << " = " << gridx(ibestvalx) << endl;
			    // restore the best configuration
			     lsqOptObj.GetFullRefinableObj().RestoreParamSet(best_set_id);
			     cout << "Best parameters set restored." << endl;
			   }
			   else
			     cerr << "Error: grid - refinement - best parameters set not found!" << endl;
			 } // if(gridy.numElements()==0)

			 if(gridy.numElements()>0) {			 
			   for(int iy=0; iy<gridy.numElements(); iy++)
			     for(int ix=0; ix<gridx.numElements(); ix++) {
			      // restore starting parmas. set
			       lsqOptObj.GetFullRefinableObj().RestoreParamSet(start_set_id);
			      // set params
			       pparamx->SetValue(gridx(ix)/pparamx->GetHumanScale());
			       pparamy->SetValue(gridy(iy)/pparamy->GetHumanScale());
			      // run the refinement
			       lsqOptObj.Refine(niter,useLevenbergMarquardt,silent);
			      // store the Rw value
			       gridRw(iy,ix) = lsqOptObj.RwFactor();
			      // save the best configuration
			       if (gridRw(iy,ix)<=bestRw) { // store the best configuration
				 bestRw = gridRw(iy,ix);
				 lsqOptObj.GetFullRefinableObj().SaveParamSet(best_set_id);
				 ibestvalx = ix;
				 ibestvaly = iy;
			       }
			     }
			  // print grid refinement results
			   cout << " ------------------------------------------------ " << endl;
			   cout << "GRID REFINEMENT results:" << endl;
			   cout << setw(20) << grid_xname << setw(20) << grid_yname << setw(20) << "Rw" << endl;
			   cout << setw(20) << "";
			   for(int iy=0; iy<gridy.numElements(); iy++)
			       cout << setw(20) << setw(20) << gridy(iy);
			   cout << endl;
			   for(int ix=0; ix<gridx.numElements(); ix++) {
			     cout << setw(20) << gridx(ix);
			     for(int iy=0; iy<gridy.numElements(); iy++)
			       cout << setw(20) << gridRw(iy,ix);
			     cout << endl;
			   }
			   if (ibestvalx>=0 && ibestvaly>=0) {
			     cout << "Best parameters set: ( " << grid_xname << " , " << grid_yname << " ) = ( ";
			     cout << gridx(ibestvalx) << " , " << gridy(ibestvaly) << " )" << endl;
			    // restore the best configuration
			     lsqOptObj.GetFullRefinableObj().RestoreParamSet(best_set_id);
			     cout << "Best parameters set restored." << endl;
			   }
			   else
			     cerr << "Error: grid - refinement - best parameters set not found!" << endl;
		 } // if(gridy.numElements()>0)

   } // if (pparamx!=0)
   }
	 else if(job_type==3) {
	   // Combined full and reduced SizeDistrib refinement

	   // list of all parameters of the optimised RefinableObj and their fixed/refined status
	   vector< pair< long, bool > > refObjParams;
	   for(long ipar=0; ipar<lsqCompiledObj.GetNbPar(); ipar++)
	     refObjParams.push_back( make_pair( ipar, lsqCompiledObj.GetPar(ipar).IsFixed() ) );

	   // saved values of regularization operators weights set before they are tuned
	   vector< REAL > userRegOpWeight;

	   for(int irep=0; irep<niter2; irep++) {
	     
	     // Tune LSQ regularization weights

	     if(lsqTuneCoef>0.) {

	       // calculate data Chi2
	       REAL chi2data = 0.;
	       {
		 int nfun = 0;
		 CrystVector_REAL t = data.GetLSQCalc(nfun);
		 t -= data.GetLSQObs(nfun);
		 t *= t;
		 t *= data.GetLSQWeight(nfun);
		 chi2data = t.sum();
	       }
 
	       vector< std::pair< RefinableObj *, unsigned int > > refinedObjMap = lsqOptObj.GetRefinedObjMap();
	       vector< std::pair< RefinableObj *, unsigned int > >::const_iterator itRefObj;
	       // the preparation cycle
	       itRefObj = refinedObjMap.begin();
	       REAL chi2regSumSq = 0.;
	       for( ; itRefObj!=refinedObjMap.end(); ++itRefObj ) {
		 for(int iOp=0; iOp<itRefObj->first->GetNbLSQRegularizationOperator(itRefObj->second); iOp++) {
		   long ind = gLSQRegularizationOperatorRegistry.Find(&(itRefObj->first->GetLSQRegularizationOperator(iOp,itRefObj->second)));
		   LSQRegularizationOperator & regOp = gLSQRegularizationOperatorRegistry.GetObj(ind);
		   if(irep==0)
		     userRegOpWeight.push_back(regOp.GetRegularizationOperatorWeight()); // save user's weight
		   chi2regSumSq += userRegOpWeight[iOp]*pow(regOp.GetValue(),2);
		 } // iOp
	       } // itRefObj
	       // the init cycle
	       itRefObj = refinedObjMap.begin();
	       for( ; itRefObj!=refinedObjMap.end(); ++itRefObj ) {
		 for(int iOp=0; iOp<itRefObj->first->GetNbLSQRegularizationOperator(itRefObj->second); iOp++) {
		   long ind = gLSQRegularizationOperatorRegistry.Find(&(itRefObj->first->GetLSQRegularizationOperator(iOp,itRefObj->second)));
		   LSQRegularizationOperator & regOp = gLSQRegularizationOperatorRegistry.GetObj(ind);
		   REAL chi2 = regOp.GetValue();
		   if(chi2regSumSq>1.e-7) {
		     REAL lambda = lsqTuneCoef*chi2data/chi2regSumSq *userRegOpWeight[iOp]*chi2;
		     regOp.SetRegularizationOperatorWeight(lambda);
		   }
		 } // iOp
	       } // itRefObj	       
	     }

	     // Reduced refinement

	     // fix parameters with some exceptions
	     vector< pair< long, bool > >::const_iterator it = refObjParams.begin();
	     for( ; it!=refObjParams.end(); ++it) {
	       ObjCryst::RefinablePar & par = lsqCompiledObj.GetPar(it->first);
	       if( par.GetType()!=ObjCryst::gpRefParTypeScattDataScale &&
		   par.GetType()!=ObjCryst::gpRefParTypeObjCryst &&
		   !par.GetType()->IsDescendantFromOrSameAs(MStruct::gpRefParTypeScattDataProfileSizeDistrib) )
		 par.SetIsFixed(true);
	     }
	     // run
	     if(niter!=0 && MStruct::Global_SIGINThandled_counter==0)
	       lsqOptObj.Refine(-niter,useLevenbergMarquardt,silent,true,0.001);

	     // unfix parameters after reduced refinement
	     it = refObjParams.begin();
	     for( ; it!=refObjParams.end(); ++it) {
	       ObjCryst::RefinablePar & par = lsqCompiledObj.GetPar(it->first);
	       par.SetIsFixed(it->second);
	     }
	     
	     // Full refinement
	     if(niter1!=0 && MStruct::Global_SIGINThandled_counter==0)
	       lsqOptObj.Refine(-niter1,useLevenbergMarquardt,silent,true,0.001);
	   }
	 }
	 else {
	   // Run the refinement (normal or Random-search?)
	   ObjCryst::RefinableObj * pRefObjRand = NULL;
	   {
	     // try to find an object of MStruct::RandomSizeDistribBroadeningEffect class
	     //ObjCryst::ObjRegistry< ObjCryst::RefinableObj > &objRegistry = lsqOptObj.GetFullRefinableObj().GetSubObjRegistry();
	     ObjCryst::ObjRegistry< ObjCryst::RefinableObj > &objRegistry = ObjCryst::gRefinableObjRegistry;
	     for(int i=0; i<objRegistry.GetNb(); i++) {
	       ObjCryst::RefinableObj &obj = objRegistry.GetObj(i);
	       //cout << "ClassName:" << obj.GetClassName() << endl;
	       if( obj.GetClassName()==string("MStruct::RandomSizeDistribBroadeningEffect") 
		|| obj.GetClassName()==string("MStruct::SizeDistribBroadeningEffect-X") ) {
		 pRefObjRand = &obj;
		 break; // the first one is used
	       }
	     }
	   }
	   
	   if(pRefObjRand==NULL) // normal LSQ refinement
	     lsqOptObj.Refine(-niter,useLevenbergMarquardt,silent,true,0.001);
	   else { // random averaging with LSQ refinement
	     // create new random optimization object
	     MStruct::RandomOptimizationObj randOptObj("RandomOptObj");
	     // set the random optimization object
	     randOptObj.SetLSQNumObj(&lsqOptObj,-niter,useLevenbergMarquardt,silent,true,0.001);
	     randOptObj.SetRefinableObjRandomised(pRefObjRand);
	     // run optimization/refinement
	     //long nbSteps = 1500;
	     long nbSteps = 400;
	     randOptObj.Optimize(nbSteps,false);
	     // print results into the file
	     randOptObj.WriteResultsToFile("randOptSet1.dat");
	   }
	 }

  // save IhklCorr params to file
   for(int iphase=0; iphase<nbphases; iphase++) { 
     if(vHKLChoice[iphase]>0) {
       vDiffData[iphase]->WriteHKLIntensityCorrToFile(string(string("Ihkl_")+vPhasesNames[iphase]+string(".dat")).c_str());
       if(vHKLChoicePrint[iphase]>0) {
       	 cout << "HKL Intensity Corrections for phase: " << vPhasesNames[iphase] << '\n';
       	 cout << "*********************************************************************\n";
       	 vDiffData[iphase]->PrintHKLIntensityCorr(cout);
       }
     }
   	
     // save size distributions
     // loop over all stored SizeDistribBroadeningEffect objects
     for(unsigned int ieffect=0; ieffect<vSizeDistribEffect.size(); ieffect++)
       vSizeDistribEffect[ieffect]->WriteDistributionToFile(vSizeDistribFileName[ieffect].c_str());

     // save HKLPseudoVoigtBroadeningEffect parameters
     if(vHKLEffect[iphase]!=NULL)
       vHKLEffect[iphase]->SaveParametersFile();
   }

  // Save the powder pattern in text format
   data.SavePowderPattern(filename.c_str());

   if( lsqCompiledObj.GetNbLSQConstraints() > 0 )
     lsqCompiledObj.PrintExternalConstraintsStatistics();
   
   } // if(1)

  // print values of Scale factors
   if(1) {
     // vector of data scale factors
     const CrystVector_REAL &scaleFactor = data.GetScaleFactor();
     cout << "\n";
     cout << setw(50) << "Scale factors" << setw(30) << "(Mcell*Vcell 1e-4)" << "\n";
     cout << "********************************************************************************\n";
     cout << scientific << setprecision(3);
     for(int icomp=0; icomp<data.GetNbPowderPatternComponent(); icomp++) {
       // check if this is a PowderPatternDiffraction component (for unscalable comps. scale=1)
       if( data.GetPowderPatternComponent(icomp).IsScalable()==false ) continue;
       // print name and human value
       RefinablePar &par = data.GetPar(scaleFactor.data()+icomp);
       cout << setw(30) << par.GetName();
       cout << setw(20) << par.GetHumanValue();
       // check if this is a diffraction componenet
       const string &cname = data.GetPowderPatternComponent(icomp).GetClassName();
       if( cname=="PowderPatternDiffraction" || cname=="MStruct::PowderPatternDiffraction" || 
	   cname=="MStruct::SizeDistribPowderPatternDiffraction" ) {
	 const ObjCryst::PowderPatternDiffraction &diff =
	   dynamic_cast<const ObjCryst::PowderPatternDiffraction&>(data.GetPowderPatternComponent(icomp));
	 // get reference to crystal
	 const ObjCryst::Crystal &crystal = diff.GetCrystal();
	 // print Mcell * Vcell (10-4)
	 ostringstream ss;
	 ss << "(" << fixed << setprecision(3)
	    << MStruct::CalcUnitCellMass(crystal)*crystal.GetVolume()*1e-4 << ")";
	 cout << setw(30) << ss.str();
	   } else {
		  // empty space printed
		   cout << setw(30) << "       ";
	   }
       cout << setw(20) << par.GetHumanSigma();
       cout << "\n";
     }

	  /*
			cout << "\n";
		cout << "Scale factors\n";
		cout << "*********************************************************************\n";
		cout << scientific << setprecision(3);
		const RefinableObj &obj = data;
    for(int iparam=0; iparam<obj.GetNbPar(); iparam++) {
    	const RefinablePar &par = obj.GetPar(iparam);
    	if (par.GetName().substr(0,5)=="Scale") {
    		cout << "\t" << par.GetName();
    		cout << "\t\t" << par.GetHumanValue() << "\n";
    	}
    }
    */
   }

   if(1 && vDiffData.size()>0) { // sV
  // Save calculated intensities
   ofstream f("phase1_par.txt");
   vDiffData[0]->PrintHKLInfo(f);
   f.close();
   /*ofstream f1("phase2_par.txt");
   vDiffData[1]->PrintHKLInfo(f1);
   f1.close();*/
   //
   /*ofstream f2("phase1_par_2.txt");
     vDiffData[0]->PrintHKLInfo2(f2,0.001);
   //vDiffData[1]->PrintHKLInfo(f2);
   f2.close();*/
  // Save the powder pattern in text format
   //data.SavePowderPattern("tio2.dat");
  // Save everything in xml so that we can reload it later
   XMLCrystFileSaveGlobal("xray.xml");
   } // sV

   if(using_file == true) supplied_input_file.close();

  // Save interpolated background
   if( settings_overwrite_bgr )
     for(int icomp=0; icomp<(int)v_bkg_type_ids.size(); icomp++)
       if (v_bkg_type_ids[icomp]==BACKGROUND_INTERPOLATED) {
	 string bkg_filename = v_bkg_strings[icomp];
	 cout << "Saving interpolated background: " << bkg_filename << "\n";
	 // get diffraction background component
	 const MStruct::PowderPatternBackground * bkgData = dynamic_cast<const MStruct::PowderPatternBackground*>(vBackgroundComponents[icomp]);
	 std::ofstream f(bkg_filename.c_str());
	 f << std::fixed;
	 const std::pair< const CrystVector_REAL *, const CrystVector_REAL * > 	points = bkgData->GetInterpPoints();
	 for(int ipoint=0; ipoint<points.first->numElements(); ipoint++) {
	   f << std::setw(10) << std::setprecision(3) << *(points.first->data()+ipoint) * RAD2DEG << "  ";
	   f << std::setw(12) << std::setprecision(3) << *(points.second->data()+ipoint) << "\n";
	 }
	 f.close();
       }
   
   cout << " End of program." << endl ;
   /*for(int i=0; i<gRefinableObjRegistry.GetNb(); i++) {
		const RefinableObj &obj = gRefinableObjRegistry.GetObj(i);
		cout<<i<<'\t'<<obj.GetClassName()<<'\t'<<obj.GetName()<<endl;
   }*/

   return 0;
}

// The function reads all empty and commented (//) lines from the input
// stream 'is' until an uncommented and non-empty line is found. This last read line is returned
// as content of the 'iss' input-string-stream which can be used for a standard input.
// If the 'load_objects' parameter is true, than all found 'objects' are loaded and poiters to
// them are returned in the output vector.
std::vector< void* > read_line ( istringstream &iss, istream &is, const bool load_objects)
{
	std::vector< void* > objects;
	iss.clear();
	
	string s;
	
	while( getline(is,s) ) {
		
		// trim the leading spaces
		size_t startpos = s.find_first_not_of(" \t"); // Find the first character position after excluding leading blank spaces
		size_t endpos = s.find_last_not_of(" \t\r\n"); // Find the first character position from reverse
		if ( ( string::npos == startpos ) || ( string::npos == endpos) )
			s = "";
		else
			s = s.substr( startpos, endpos-startpos+1 );
		
		// ignoring empty lines
		if ( s.empty() ) continue;
		
		// ignoring comments (//) 
		if ( s.size() >= 2 && s[0] == '/' && s[1] == '/' ) continue; 
		
		// if there is something on the line break the cycle
		if ( !s.empty() ) break;
	}
	
	iss.str(s);
	//cout << "iss: " << iss.str().c_str() << endl;
	
	return objects;
}

// this function reads all empty and commented (//) lines from the input
// stream until an uncommented and non-empty line is found 
void ignore_comments(istream & s)
{
		static int fn = 0;
		
		fn++;
		cout << "IgnoreComments START:" << fn << endl;
				
		//if (s.peek()==EOF) return;

		bool cycle = (s.peek()!=EOF);
		
		cout << "cycle: " << cycle << endl;
		
		while (cycle) {
			
			char t;
			int c = s.peek();
			int c1;
			
			switch (c) {
				case EOF:
							s.get(t);
							cycle = false;
							break;
				case ' ':
				case '\t':
				case '\n':
				case '\r':
							s.get(t);
							break;
				case '/':
							s.get(t);
							// We need to know another char
							c1 = s.peek();
							if (c1=='/') { // this is a comment line
								s.get(t);
								s.ignore(1024,'\n');
							} else if (c!=EOF) { // this is a sigle '/' or EOF
								// put back the read '/'
								s.putback(c);
								cycle = false;
							}
							else { // this is a '/' followed by EOF
								s.get(t);
								cerr << "Warning: Unmatched '/' at the end of the file." << endl;
								cycle = false;
							}
							break;
				default:
							cycle = false;
							break;
			}			
		}
		cout << "IgnoreComments END:" << fn << endl;
}

// rename a refinable parameter of the refinble object (obj), parameter with name (old_name)
// must be a parameter of the object (obj) and is renamed (with new_name) only in the case
// the paremeter is used, function returns a reference to the renamed parameter,
// exceptions not handled 
RefinablePar& rename_par(RefinableObj &obj, const string &old_name, const string &new_name)
{
		RefinablePar &par = obj.GetPar(old_name);
		if (par.IsUsed()==true) par.SetName(new_name);
		return par;
}

// Rename crystal lattice parameters and Scattering powers Biso factor parameters
// to ensure better unique access to them. Parameters names are extended with
// 'Crystal name' suffix (eg. 'a' -> 'a_Anatase)
void rename_crystal_params(Crystal &crystal)
{ 
		string crystal_name = crystal.GetName();
		
		// lattice parameters
		{
			rename_par(crystal,"a",string("a_")+crystal_name);
			rename_par(crystal,"b",string("b_")+crystal_name);
			rename_par(crystal,"c",string("c_")+crystal_name);
			rename_par(crystal,"alpha",string("alpha_")+crystal_name);
			rename_par(crystal,"beta",string("beta_")+crystal_name);
			rename_par(crystal,"gamma",string("gamma_")+crystal_name);
		}
		
		// scat. powers Biso factors
		{
			ObjRegistry< ScatteringPower > & registry = crystal.GetScatteringPowerRegistry();
			for(int iscat=0; iscat<registry.GetNb(); iscat++) 
				rename_par(registry.GetObj(iscat),"Biso",string("Biso_")+crystal_name);
		}
}

// extended get parameter function, if a parameter name 'extended_name'
// is a simple parameter name (not closer specified by ':' deliminator)
// function simply returns a reference to a parameter of the refinable object
// identifing paremeter by its name 'extended_name' in the refinable object
// 'ref_obj' (same as the RefinableObj method GetPar(const string&), in the case
// the 'extended_name' is a composed extended name containing a deliminator ':'
// (eg. Anatase:a), functions tries to find the specified parameter in the specified
// object in the global refinable object registry and returna reference to it.
// (eg. for 'Anatase:a' function scours the refinable object called 'Anatase'
// for a parameter called 'a') 
RefinablePar& get_par(const string &extended_name, RefinableObj &ref_obj)
{
		// separete object and parameter name
		string::size_type loc = extended_name.find( ':', 0 );
   	if( loc == string::npos ) {
    	// this is not a composed extended name - simple case
    	ref_obj.GetPar(extended_name.c_str()).SetExtName(extended_name);
    cout << "ext_name = " << extended_name << " : " << ref_obj.GetPar(extended_name.c_str()).GetExtName() << endl;
    return ref_obj.GetPar(extended_name.c_str());
   	} else {
   		// a composed extended name 
   	
   		// split the name
   		string par_name;
   		vector<string> obj_names;
   		
   		obj_names.push_back( extended_name.substr(0,loc) );
   		
   		string::size_type loc1 = extended_name.find( ':', loc+1 );
   		
   		while( loc1 != string::npos ) {
   			
   			obj_names.push_back( extended_name.substr(loc+1,loc1-loc-1) );
   			loc = loc1;
   			loc1 = extended_name.find( ':', loc+1 );
   			
   		}
   		
   		par_name.assign( extended_name, loc+1, extended_name.length() - loc -1);
   		
   		// get the first (topmost) object 'obj_name' from a global refinable object registry by its name
   		RefinableObj * obj = & gRefinableObjRegistry.GetObj(obj_names[0]);
   		// and continue getting objects recursively from its SubObjRegistry
   		for(unsigned int i=1; i<obj_names.size(); i++) obj = & obj->GetSubObjRegistry().GetObj(obj_names[i]);
   		
   		// get the parameter 'name' from the found object by its name
   		RefinablePar &par = obj->GetPar(par_name);
   		
    par.SetExtName(extended_name);
    cout << "ext_name = " << extended_name << " : " << par.GetExtName() << "  ->  " << par_name << endl;
			// the found parameter 'par' will be returned as a result
			return par;
   	}
}


