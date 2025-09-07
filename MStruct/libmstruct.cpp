#include "config.h"
#ifdef __PYMSTRUCT_TEST__
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
//Includes for Python API
#include <ObjCryst/ObjCryst/IO.h>
#include "MStruct.h"
#include "IO.h"

namespace bp = boost::python;

string version()
{
  return mstruct_version_str + string("-python") + python_version_str;
}

template<typename s,typename T>
void _SetWavelength(s& self, const T wavelength)
{
  self.SetWavelength(wavelength);
}

template<typename s, typename T>
void _SetParentProfile(s& self, T& refProfile)
{
  refProfile.AddReflectionProfileComponent(self);
  self.SetParentReflectionProfile(refProfile);
}

template<typename T>
ObjCryst::RefinablePar& _GetPar(T& self, const string & name)
{
    return self.GetPar(name);
}

template<typename T>
ObjCryst::Radiation _Get_Radiation(T& self)
{
  return self.GetRadiation();
}

CrystVector_REAL* NumpyArrayToCrystVector_REAL(const boost::python::numpy::ndarray& vector)
{
  const long array_length = vector.get_nd();
  CrystVector_REAL* crystvector = new CrystVector_REAL(array_length);
  for (int i=0; i<array_length; i++)
  {
    crystvector[i] = boost::python::extract<REAL>(vector[i]);
  }

  return crystvector;
}

boost::python::numpy::ndarray CrystVector_REAL_to_NumpyArray(const CrystVector_REAL vector)
{
  const unsigned int length = vector.numElements();
  boost::python::tuple shape = boost::python::make_tuple(length);
  boost::python::numpy::dtype dt = boost::python::numpy::dtype::get_builtin<REAL>();
  boost::python::numpy::ndarray array = boost::python::numpy::zeros(shape, dt);
  for (int i=0; i<length; i++)
  {
    array[i] = vector(i);
  }
  return array;
}

void _SetPowderPatternObs(MStruct::PowderPattern& self, const boost::python::numpy::ndarray& vector)
{
  CrystVector_REAL* crystalvector = NumpyArrayToCrystVector_REAL(vector);
  self.SetPowderPatternObs(*crystalvector);
}

void _AddPowderPatternComponent(MStruct::PowderPattern& self, MStruct::PowderPatternDiffraction& powdif)
{
  self.AddPowderPatternComponent(powdif);
}

void _SetProfile(MStruct::PowderPatternDiffraction& self, MStruct::ReflectionProfile& profile)
{
  profile.SetParentPowderPatternDiffraction(self);
  self.SetProfile(&profile);
}

// Component should be general
void _AddComponent_Refraction(MStruct::ReflectionProfile& self, MStruct::RefractionPositionCorr& refCorr)
{
  self.AddReflectionProfileComponent(refCorr);
  refCorr.SetParentReflectionProfile(self);
}

void _SavePowderPattern(MStruct::PowderPattern& self, const string & filename)
{
    self.SavePowderPattern(filename);
}

boost::python::numpy::ndarray _GetPowderPatternX(MStruct::PowderPattern& self)
{
  return CrystVector_REAL_to_NumpyArray(self.GetPowderPatternX());
}

boost::python::numpy::ndarray _GetPowderPatternObs(MStruct::PowderPattern& self)
{
  return CrystVector_REAL_to_NumpyArray(self.GetPowderPatternObs());
}

ObjCryst::Crystal * _Create_Crystal()
{
    return new ObjCryst::Crystal;
}

ObjCryst::Crystal * _XMLLoadCrystal(const std::string& file_name, const std::string& crystal_name)
{
  ObjCryst::Crystal * crystal;
  ObjCryst::XMLCrystFileLoadObject(file_name,"Crystal", crystal_name, crystal);
  ObjCryst::RefinableObj &obj = ObjCryst::gRefinableObjRegistry.GetObj(crystal_name, "Crystal");
  crystal = &(dynamic_cast<ObjCryst::Crystal&>(obj));
  return crystal;
}

void _PowderPattern_SetObsToZero(MStruct::PowderPattern& self)
{
    CrystVector_REAL obs(self.GetNbPoint());
    obs = 0.;
    self.SetPowderPatternObs(obs);
}

boost::python::numpy::ndarray _GetPowderPatternCalc(MStruct::PowderPattern& self){
    CrystVector_REAL crystvector = self.GetPowderPatternCalc();
    return CrystVector_REAL_to_NumpyArray(crystvector);
}

void _ReflectionProfile_SetParentPowderPatternDiffraction(MStruct::ReflectionProfile& self,
                                                          MStruct::PowderPatternDiffraction& scattData)
{
    self.SetParentPowderPatternDiffraction(scattData);
}




/*
MStruct::ReflectionProfile * _Create_ReflectionProfile(ObjCryst::Crystal* crystal, PowderPattern* pattern)
{
  MStruct::ReflectionProfile * profile =
    new MStruct::ReflectionProfile(*crystal, pattern->GetRadiation());
    return profile;
}
*/

MStruct::PowderPattern& _gRefinableObjRegistry_GetPowderPattern(const std::string& name)
{
  // Get PowderPattern object
  ObjCryst::RefinableObj &obj = ObjCryst::gRefinableObjRegistry.GetObj(name, "MStruct::PowderPattern");
  MStruct::PowderPattern &data = dynamic_cast<MStruct::PowderPattern&>(obj);
  return data;
}

ObjCryst::Crystal& _gRefinableObjRegistry_GetCrystal(const std::string& name)
{
  // Get Crystal object
  ObjCryst::RefinableObj &obj = ObjCryst::gRefinableObjRegistry.GetObj(name, "Crystal");
  ObjCryst::Crystal &crystal = dynamic_cast<ObjCryst::Crystal&>(obj);
  return crystal;
}

// credit: @vincefn, https://github.com/diffpy/pyobjcryst/blob/main/src/extensions/lsq_ext.cpp

/*bool _LSQ_SafeRefine(MStruct::LSQNumObj & lsq, REAL maxChi2factor, int nbCycle, bool useLevenbergMarquardt,
                 const bool silent, const bool callBeginEndOptimization,const float minChi2var)
{
    //CaptureStdOut gag;
    //if(!silent) gag.release();

    std::list<ObjCryst::RefinablePar*> vnewpar;
    std::list<const ObjCryst::RefParType*> vnewpartype;
    return lsq.SafeRefine(vnewpar, vnewpartype, nbCycle, useLevenbergMarquardt, silent,
                          callBeginEndOptimization, minChi2var);
}*/

ObjCryst::RefinablePar& _RefinableObj_GetParString(ObjCryst::RefinableObj& obj, const string& s)
{
    obj.SetDeleteRefParInDestructor(0);
    return obj.GetPar(s);
}

// extended get parameter function: the 'extended_name' must be a composed
// extended name. i.e. containi a deliminator ':' as e.g. Anatase:a.
// Functions tries to find the specified parameter in the specified
// object in the global refinable object registry and return a reference to it.
// (e.g. for 'Anatase:a' function scours the refinable object called 'Anatase'
// for a parameter called 'a') 
ObjCryst::RefinablePar& _GetParExtString(const string &extended_name)
{
		// separete object and parameter name
		string::size_type loc = extended_name.find( ':', 0 );
   	if( loc == string::npos ) {
    	// this is not a composed extended name - simple case
       throw std::invalid_argument("this is not a composed extended name");
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
   		ObjCryst::RefinableObj * obj = & ObjCryst::gRefinableObjRegistry.GetObj(obj_names[0]);
   		// and continue getting objects recursively from its SubObjRegistry
   		for(unsigned int i=1; i<obj_names.size(); i++) obj = & obj->GetSubObjRegistry().GetObj(obj_names[i]);
   		
   		// get the parameter 'name' from the found object by its name
   		ObjCryst::RefinablePar &par = obj->GetPar(par_name);
   		
		return par;
   	}
}

void test_numpy(const boost::python::numpy::ndarray& test_array)
{
  std::cout << "In test_numpy" << endl;;
  int var = boost::python::extract<int>(test_array[0]);
  std::cout << var << endl;
}


BOOST_PYTHON_MODULE(libMStruct)
{
  using namespace boost::python;
  Py_Initialize();
  numpy::initialize();
  def("version", version);
  def("test_numpy", test_numpy);
  def("CreateCrystalFromXML", &_XMLLoadCrystal, return_value_policy<manage_new_object>());
  def("XMLCrystFileLoadAllObject", (void (*)(const std::string&)) &MStruct::XMLCrystFileLoadAllObject, (bp::arg("name")));
  def("XMLCrystFileSaveGlobal", (void (*)(const std::string&)) &ObjCryst::XMLCrystFileSaveGlobal, (bp::arg("name")));
  def("GetPowderPattern", _gRefinableObjRegistry_GetPowderPattern, return_value_policy<reference_existing_object>());
  def("GetCrystal", _gRefinableObjRegistry_GetCrystal, return_value_policy<reference_existing_object>());
  def("GetPar", _GetParExtString, return_value_policy<reference_existing_object>());
  def("CalcUnitCellMass", (REAL (*)(const ObjCryst::Crystal&)) &MStruct::CalcUnitCellMass, (bp::arg("crystal")));
  //def("Create_ReflectionProfile", &_Create_ReflectionProfile);

  class_<ObjCryst::Restraint, boost::noncopyable>("Restraint");

  class_<ObjCryst::RefinableObjClock>("RefinableObjClock")
    .def("Reset", &ObjCryst::RefinableObjClock::Reset, "Reset a Clock to 0, to force an update")
    ;
  
  // credit: @vincefn, https://github.com/diffpy/pyobjcryst/blob/main/src/extensions/lsq_ext.cpp
  // Class for holding RefinablePar instances created in c++. This should not
  // be exported to the pyobjcryst namespace
  class_<ObjCryst::RefinablePar, bases<ObjCryst::Restraint> > ("_RefinablePar", no_init)
      .def("SetHumanValue", &ObjCryst::RefinablePar::SetHumanValue)
      .def("Print", &ObjCryst::RefinablePar::Print)
      // Python-only attributes
      //.def("__str__", &__str__<ObjCryst::RefinablePar>)
      .add_property("value", &ObjCryst::RefinablePar::GetValue, &ObjCryst::RefinablePar::SetValue)
      ;

  class_<ObjCryst::RefinableObj, boost::noncopyable>("RefinableObj")
      .def(init<bool>())
      .def("GetPar", &_RefinableObj_GetParString,
            return_internal_reference<>())
      .def("GetClockMaster", &ObjCryst::RefinableObj::GetClockMaster,
            return_value_policy<copy_const_reference>())
      .def("PrepareForRefinement", &ObjCryst::RefinableObj::PrepareForRefinement)
      .def("Print", &ObjCryst::RefinableObj::Print)
       // Virtual
      .def("GetClassName", &ObjCryst::RefinableObj::GetClassName,
            &ObjCryst::RefinableObj::GetClassName,
            return_value_policy<copy_const_reference>())
      .def("GetName", &ObjCryst::RefinableObj::GetName,
            &ObjCryst::RefinableObj::GetName,
            return_value_policy<copy_const_reference>());

  class_<ObjCryst::ScatteringData, bases<ObjCryst::RefinableObj>, boost::noncopyable>(
            "ScatteringData", no_init)
      .def("GetCrystal", (ObjCryst::Crystal& (ObjCryst::ScatteringData::*)()) &ObjCryst::ScatteringData::GetCrystal,
        return_internal_reference<>())
      .def("HasCrystal", &ObjCryst::ScatteringData::HasCrystal);

  class_<ObjCryst::PowderPatternComponent, bases<ObjCryst::RefinableObj>, boost::noncopyable>
        ("_ObjCryst::PowderPatternComponent", no_init)
      .def("GetParentPowderPattern",
        (ObjCryst::PowderPattern& (ObjCryst::PowderPatternComponent::*)())
        &ObjCryst::PowderPatternComponent::GetParentPowderPattern,
          return_internal_reference<>());

  class_<ObjCryst::PowderPatternDiffraction, bases<ObjCryst::PowderPatternComponent, ObjCryst::ScatteringData> >(
                "_Obj_Cryst_PowderPatternDiffraction", no_init);

  class_<MStruct::PowderPattern, bases<ObjCryst::RefinableObj> >("PowderPattern")
      .def(init<>())
      .def("SetPowderPatternObs", &_SetPowderPatternObs)
      .def("SetWeightToUnit", &MStruct::PowderPattern::SetWeightToUnit)
      .def("SetSigmaToSqrtIobs", &MStruct::PowderPattern::SetSigmaToSqrtIobs)
      .def("SetWeightToInvSigmaSq", &MStruct::PowderPattern::SetWeightToInvSigmaSq)
      .def("Prepare", &MStruct::PowderPattern::Prepare)
      .def("SetObsToZero", &_PowderPattern_SetObsToZero)
      .def("SavePowderPattern", &_SavePowderPattern)
      .def("SetPar", &MStruct::PowderPattern::SetPowderPatternPar)
      .def("SetIncidenceAngle", &MStruct::PowderPattern::SetIncidenceAngle)
      .def("SetWavelength", &_SetWavelength<MStruct::PowderPattern, string>)
      .def("SetWavelength", &_SetWavelength<MStruct::PowderPattern, REAL>)
      .def("GetRadiation", &_Get_Radiation<MStruct::PowderPattern>)
      .def("GetCalc", &_GetPowderPatternCalc)
      .def("GetPowderPatternX", &_GetPowderPatternX)
      .def("GetPowderPatternObs", &_GetPowderPatternObs)
      .def("GetPar", &_GetPar<MStruct::PowderPattern>, return_value_policy<reference_existing_object>())
      .def("Print", &MStruct::PowderPattern::Print) 
      .def("AddComponent", &MStruct::PowderPattern::AddPowderPatternComponent)
      .def("AddComponent", &_AddPowderPatternComponent)
      .def("GetNbPowderPatternComponent", &ObjCryst::PowderPattern::GetNbPowderPatternComponent)
      .def("GetPowderPatternComponent", (ObjCryst::PowderPatternComponent& (ObjCryst::PowderPattern::*) (const int))
        &ObjCryst::PowderPattern::GetPowderPatternComponent, return_internal_reference<>())
      .def("FitScaleFactorForRw", &MStruct::PowderPattern::FitScaleFactorForRw)
      .def("GetScaleFactor", (REAL (MStruct::PowderPattern::*) (const int) const) &ObjCryst::PowderPattern::GetScaleFactor)
      .def("GetScaleFactor", (REAL (MStruct::PowderPattern::*) (const ObjCryst::PowderPatternComponent&) const) &ObjCryst::PowderPattern::GetScaleFactor)
      .def("SetScaleFactor", (void (MStruct::PowderPattern::*) (const int, REAL)) &ObjCryst::PowderPattern::SetScaleFactor);

  class_<MStruct::PseudoVoigtBroadeningEffect>("PseudoVoigtBroadeningEffect")
      .def(init<>())
      .def("SetParentProfile", &_SetParentProfile<MStruct::PseudoVoigtBroadeningEffect, MStruct::ReflectionProfile>);

  class_<ObjCryst::Crystal, bases<ObjCryst::RefinableObj> >("Crystal")
      .def("GetVolume", (REAL (ObjCryst::Crystal::*) () const) &ObjCryst::Crystal::GetVolume);

  //class_<ReflectionProfile, boost::noncopyable>("ReflectionProfile", no_init);
  class_<MStruct::ReflectionProfile>("ReflectionProfile", init<ObjCryst::Crystal&, ObjCryst::Radiation&>())
      .def("AddComponent", &_AddComponent_Refraction)
      .def("GetPar", &_GetPar<MStruct::ReflectionProfile>, return_value_policy<reference_existing_object>())
      .def("SetParentPowderPatternDiffraction", &_ReflectionProfile_SetParentPowderPatternDiffraction)
      .def("SetIncidenceAngle", &MStruct::ReflectionProfile::SetIncidenceAngle);

  class_<MStruct::RefractionPositionCorr>("RefractionPositionCorr")  
      .def(init<>())
      .def("GetPositionCorr", &MStruct::RefractionPositionCorr::GetPositionCorr)
      .def("SetParams", &MStruct::RefractionPositionCorr::SetParams)
      .def("GetPar", &_GetPar<MStruct::RefractionPositionCorr>, return_value_policy<reference_existing_object>())
      .def("SetParentProfile", &_SetParentProfile<MStruct::RefractionPositionCorr, MStruct::ReflectionProfile>)
      .def("SetCrystal", &MStruct::RefractionPositionCorr::SetCrystal);

  class_<ObjCryst::Radiation>("Radiation", init<const std::string, REAL>())
      //.def("GetWavelength", &ObjCryst::Radiation::GetWavelength)
      .def("Print", &ObjCryst::Radiation::Print);
  
  class_<MStruct::PowderPatternDiffraction, bases<ObjCryst::PowderPatternDiffraction> >("PowderPatternDiffraction")
      .def(init<>())
      .def("SetProfile", &_SetProfile)
      .def("SetIsIgnoringImagScattFact", &MStruct::PowderPatternDiffraction::SetIsIgnoringImagScattFact)
      .def("SetCrystal", &MStruct::PowderPatternDiffraction::SetCrystal);

  class_<ObjCryst::DiffractionDataSingleCrystal>("DiffractionDataSingleCrystal", init<ObjCryst::Crystal*>())
      .def("SetWavelength", &_SetWavelength<ObjCryst::DiffractionDataSingleCrystal, string>)
      .def("SetWavelength", &_SetWavelength<ObjCryst::DiffractionDataSingleCrystal, REAL>)
      .def("SetHKL", &ObjCryst::DiffractionDataSingleCrystal::SetHKL)
      .def("GetRadiation", &_Get_Radiation<ObjCryst::DiffractionDataSingleCrystal>)
      .def("SetIsIgnoringImagScattFact", &ObjCryst::DiffractionDataSingleCrystal::SetIsIgnoringImagScattFact)
    ;

  // credit: @vincefn, https://github.com/diffpy/pyobjcryst/blob/main/src/extensions/lsq_ext.cpp
  class_<MStruct::LSQNumObj>("LSQ")
      /// LSQNumObj::PrepareRefParList() must be called first!
      .def("SetParIsFixed",
              (void (MStruct::LSQNumObj::*)(const std::string&, const bool))
              &ObjCryst::LSQNumObj::SetParIsFixed,
              (bp::arg("parName"), bp::arg("fix")))
      .def("SetParIsFixed",
              (void (MStruct::LSQNumObj::*)(const ObjCryst::RefParType *, const bool))
              &ObjCryst::LSQNumObj::SetParIsFixed,
              (bp::arg("type"), bp::arg("fix")))
      .def("SetParIsFixed",
              (void (MStruct::LSQNumObj::*)(ObjCryst::RefinablePar &, const bool))
              &ObjCryst::LSQNumObj::SetParIsFixed,
              (bp::arg("par"), bp::arg("fix")))
      //void SetParIsFixed(RefinableObj &obj, const bool fix);
      .def("UnFixAllPar", &ObjCryst::LSQNumObj::UnFixAllPar)
      //void SetParIsUsed(const std::string& parName, const bool use);
      //void SetParIsUsed(const RefParType *type, const bool use);
      .def("Refine", &MStruct::LSQNumObj::Refine,
              (bp::arg("nbCycle")=1,
               bp::arg("useLevenbergMarquardt")=false,
               bp::arg("silent")=false,
               bp::arg("callBeginEndOptimization")=true,
               bp::arg("minChi2var")=0.01))
      /*.def("SafeRefine", &_LSQ_SafeRefine,
              (bp::arg("maxChi2factor")=1.01,bp::arg("nbCycle")=1,
               bp::arg("useLevenbergMarquardt")=false,
               bp::arg("silent")=false,
               bp::arg("callBeginEndOptimization")=true,
               bp::arg("minChi2var")=0.01))*/
      .def("Rfactor", &MStruct::LSQNumObj::Rfactor)
      .def("RwFactor", &MStruct::LSQNumObj::RwFactor)
      .def("ChiSquare", &MStruct::LSQNumObj::ChiSquare)
      .def("SetRefinedObj", &MStruct::LSQNumObj::SetRefinedObj,
              (bp::arg("obj"), bp::arg("LSQFuncIndex")=0,
               bp::arg("init")=true, bp::arg("recursive")=false),
               with_custodian_and_ward<1,2>())
      .def("AddRefinableObj", &ObjCryst::OptimizationObj::AddRefinableObj,
              (bp::arg("obj")))             
      .def("GetCompiledRefinedObj",
              (ObjCryst::RefinableObj& (MStruct::LSQNumObj::*)())
              &MStruct::LSQNumObj::GetCompiledRefinedObj,
              with_custodian_and_ward_postcall<0,1,return_internal_reference<> >())
      .def("PrintRefResults", &MStruct::LSQNumObj::PrintRefResults)
      .def("PrepareRefParList", &ObjCryst::LSQNumObj::PrepareRefParList,
              (bp::arg("copy_param")=false))
      .def("GetLSQCalc", &MStruct::LSQNumObj::GetLSQCalc,
              return_value_policy<copy_const_reference>())
      .def("GetLSQObs", &MStruct::LSQNumObj::GetLSQObs,
              return_value_policy<copy_const_reference>())
      .def("GetLSQWeight", &MStruct::LSQNumObj::GetLSQWeight,
              return_value_policy<copy_const_reference>())
      .def("GetLSQDeriv", &MStruct::LSQNumObj::GetLSQDeriv,
              bp::arg("par"),
              return_value_policy<copy_const_reference>())
      .def("BeginOptimization", &ObjCryst::LSQNumObj::BeginOptimization,
              (bp::arg("allowApproximations")=false,
               bp::arg("enableRestraints")=false))
      .def("EndOptimization", &ObjCryst::LSQNumObj::EndOptimization)
      ;
}

#endif /* #ifdef __PYMSTRUCT_TEST__ */
