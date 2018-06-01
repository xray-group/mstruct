#include "MStruct.h"
#ifdef __PYMSTRUCT_TEST__
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
//Includes for Python API
#include <ObjCryst/ObjCryst/IO.h>

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
  const uint length = vector.numElements();
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

void test_numpy(const boost::python::numpy::ndarray& test_array)
{
  std::cout << "In test_numpy" << endl;;
  int var = boost::python::extract<int>(test_array[0]);
  std::cout << var << endl;
}


BOOST_PYTHON_MODULE(libMStruct)
{
  using namespace boost::python;
  numpy::initialize();
  def("test_numpy", test_numpy);
  def("CreateCrystalFromXML", &_XMLLoadCrystal, return_value_policy<manage_new_object>());
  //def("Create_ReflectionProfile", &_Create_ReflectionProfile);

  class_<MStruct::PowderPattern>("PowderPattern")
      .def(init<>())
      .def("SetPowderPatternObs", &_SetPowderPatternObs)
      .def("SetWeightToUnit", &MStruct::PowderPattern::SetWeightToUnit)
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
      .def("GetPar", &_GetPar<MStruct::PowderPattern>, return_value_policy<reference_existing_object>())
      .def("Print", &MStruct::PowderPattern::Print) 
      .def("AddComponent", &MStruct::PowderPattern::AddPowderPatternComponent)
      .def("AddComponent", &_AddPowderPatternComponent);

  class_<ObjCryst::RefinablePar>("RefinablePar")
      .def("SetHumanValue", &ObjCryst::RefinablePar::SetHumanValue)
      .def("Print", &ObjCryst::RefinablePar::Print); 

  class_<MStruct::PseudoVoigtBroadeningEffect>("PseudoVoigtBroadeningEffect")
      .def(init<>())
      .def("SetParentProfile", &_SetParentProfile<MStruct::PseudoVoigtBroadeningEffect, MStruct::ReflectionProfile>);

  class_<ObjCryst::Crystal>("Crystal");

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

  class_<MStruct::PowderPatternDiffraction>("PowderPatternDiffraction")
      .def(init<>())
      .def("SetProfile", &_SetProfile)
      .def("SetIsIgnoringImagScattFact", &MStruct::PowderPatternDiffraction::SetIsIgnoringImagScattFact)
      .def("SetCrystal", &MStruct::PowderPatternDiffraction::SetCrystal);

  class_<ObjCryst::DiffractionDataSingleCrystal>("DiffractionDataSingleCrystal", init<ObjCryst::Crystal*>())
      .def("SetWavelength", &_SetWavelength<ObjCryst::DiffractionDataSingleCrystal, string>)
      .def("SetWavelength", &_SetWavelength<ObjCryst::DiffractionDataSingleCrystal, REAL>)
      .def("SetHKL", &ObjCryst::DiffractionDataSingleCrystal::SetHKL)
      .def("GetRadiation", &_Get_Radiation<ObjCryst::DiffractionDataSingleCrystal>)
      .def("SetIsIgnoringImagScattFact", &ObjCryst::DiffractionDataSingleCrystal::SetIsIgnoringImagScattFact);
}

#endif /* #ifdef __PYMSTRUCT_TEST__ */
