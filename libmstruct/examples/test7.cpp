// Test program for checking the installation of ObjCryst shared library.
// Compile and run this code using:
//
//      c++ testlib.cpp -lObjCryst
//      ./a.out


#include <iostream>
#include <iomanip>
#include <ObjCryst/RefinableObj/RefinableObj.h>
#include <ObjCryst/ObjCryst/IO.h>
#include <ObjCryst/ObjCryst/Crystal.h>
#include <ObjCryst/ObjCryst/PowderPattern.h>
#include <MStruct/MStruct.h>
#include <ObjCryst/Quirks/VFNStreamFormat.h>

using namespace std;

// Refraction effect
int mstruct_test7(int argc, char* argv[], std::istream &is)
{
  cout << "Running test no. 7.0" << endl;

  // crystal name (to load)
  string crystal_name = "Strontium Dodecairon(III) Oxide";

  // load crystal from structures file
  ObjCryst::Crystal * crystal;
  XMLCrystFileLoadObject("HexaFerit01.xml","Crystal", crystal_name, crystal);
  /* Unfortunatelly there is a clear shortage in ObjCryst:
     T*obj in XMLCrystFileLoadObject function is not referenced as
     a reference or double pointer, hence the value of created
     T*obj si not returned to a calling scope */
  //if(crystal==NULL) throw ObjCrystException("Crystal Not Found!");
  {
    ObjCryst::RefinableObj &obj = ObjCryst::gRefinableObjRegistry.GetObj(crystal_name, "Crystal");
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
  cout << "  Compare with http://sergey.gmca.aps.anl.gov/cgi/www_form.exe?template=x0h_form.htm,\n";
  cout << "  setting CuKa1, SrFe12O19 and density below.\n";
  // calculate refraction correction for an interval of incidence angles
  REAL min=0.00, max=1.0, step=0.05;
  cout << "incidence angle (deg): min, max, step" << endl;

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

  string formula = "TiO2";
  MStruct::RefractionPositionCorr::ChemicalFormula formulaObj;

  formulaObj.SetFormula(formula);
  
  // --- cleanup ---  
  ///  if(pattern != 0) delete pattern;
  ///if(scattDataPowder != 0) delete scattDataPowder;
  ///if(crystal != 0) delete crystal;
        
  cout << " End of test." << endl ;
	
  return 0;
}

int main(int argc, char* argv[])
{
  mstruct_test7(argc, argv, cin);
  std::cout << "Using of ObjCryst/MStruct shared library works!\n";
  return 0;
}

