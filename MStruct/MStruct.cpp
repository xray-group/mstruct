/* 
 * MStruct.cpp
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
 
#ifdef MSVC
#pragma warning(disable : 4996)
#endif // MSVC

#include "MStruct.h"

#include <iostream>
#include <fstream>
#include <iomanip>

#include "Quirks/VFNStreamFormat.h"

#include "cctbx/sgtbx/space_group.h"
#include "cctbx/sgtbx/lattice_symmetry.h"
#include "cctbx/sgtbx/rot_mx_info.h"
#include "scitbx/vec3.h"
#include "cctbx/eltbx/tiny_pse.h"

#include "newmat/newmatap.h" //for matrix equation solution
#include "newmat/newmatio.h" //in CalcParams subroutine

#define _USE_MATH_DEFINES

#include <math.h> // <cmath>
#include <float.h>
#include <fftw3.h>
#include <limits>
#include <ctime>

#include <stdlib.h> // rand, srand
//#include <time.h> // time

#include <assert.h>

#include <stdexcept>

#include <signal.h> // To catch CTRL+C signal during LSQNumObj refinement

bool bsavecalc = false;
REAL xcenterlimits[2] = {25.2*DEG2RAD, 25.5*DEG2RAD}; // 14 16, 144 145 102 107, 135 138, 72 77, 15 22, 23 27, 5 27, 87 93

#define absorption_corr_factor 1.e4

//extern "C" double expei_(double *x); from ei.f
//extern "C" double daw_(double *x); from dw.f

// --- REAL RANDOM VALUES GENERATOR ---

#define real_rand() ( (REAL) ((double) rand() / (RAND_MAX + 1.)) )

double func_ei(const double x);
double func_daw(const double x);

// already defined in the file ReflectionProfile.cpp
#if defined(_MSC_VER) || defined(__BORLANDC__)
#undef min // Predefined macros.... (wx?)
#undef max

double erfc(const double x);// in C99, but not in VC++....
   // defined later in this file
#endif

//#define __DEBUG__ZDENEK__

using namespace ObjCryst;

namespace RotationTB {

void GetEulerAngles(REAL& phi1, REAL& phi, REAL& phi2,
										const CrystMatrix_REAL& A)
{
  phi = acos(A(2,2)); // a33 = cos(phi)
  if (fabs(phi)>FLT_EPSILON && fabs(phi-M_PI)>FLT_EPSILON) {
    REAL sinphi = sin(phi);
    // a31 =  sin(phi) sin(phi1)
    // a32 = -sin(phi) cos(phi1)
    phi1 = atan2(A(2,0)/sinphi,-A(2,1)/sinphi);
    // a13 = sin(phi) sin(phi2)
    // a23 = sin(phi) cos(phi2)
    phi2 = atan2(A(0,2)/sinphi,A(1,2)/sinphi); }
  else {
    // phi2 == 0, a11 = cos(phi1), a12 = sin(phi1)
    phi1 = atan2(A(0,1),A(0,0));
    phi2 = 0.0; }
}

CrystMatrix_REAL GetEulerMatrix(const REAL phi1, const REAL phi, const REAL phi2)
{
  CrystMatrix_REAL A(3,3);
  A(0,0) = cos(phi1)*cos(phi2) - cos(phi)*sin(phi1)*sin(phi2);
  A(0,1) = cos(phi2)*sin(phi1) + cos(phi1)*cos(phi)*sin(phi2);
  A(0,2) = sin(phi2)*sin(phi);
  A(1,0) = -(cos(phi2)*cos(phi)*sin(phi1)) - cos(phi1)*sin(phi2);
  A(1,1) = cos(phi1)*cos(phi2)*cos(phi) - sin(phi1)*sin(phi2);
  A(1,2) = cos(phi2)*sin(phi);
  A(2,0) = sin(phi1)*sin(phi);
  A(2,1) = -(cos(phi1)*sin(phi));
  A(2,2) = cos(phi);
  return A;
}

CrystMatrix_REAL MatrixMult(const CrystMatrix_REAL& A, const CrystMatrix_REAL& B)
{
  if(A.cols()!=B.rows())
    throw ObjCrystException("RotationTB::MatrixMult: Matrix dimension mismatch!");
  
  CrystMatrix_REAL C(A.rows(),B.cols());
  
  C = 0.;
  
  for(int i=0; i<C.rows(); i++)
    for(int j=0; j<C.cols(); j++)
      for(int k=0; k<A.cols(); k++) C(i,j) += A(i,k)*B(k,j);
		
  return C;
}

CrystMatrix_REAL MatrixTranspose(const CrystMatrix_REAL& A)
{
  CrystMatrix_REAL B(A.cols(),A.rows());
  
  for(int i=0; i<A.rows(); i++)
    for(int j=0; j<A.cols(); j++) B(j,i) = A(i,j);
  
  return B;
}

CrystVector_REAL VectorTransf(const CrystMatrix_REAL& A, const CrystVector_REAL& b)
{
  if(A.cols()!=b.numElements())
    throw ObjCrystException("RotationTB::VectorTransf: Matrix dimension mismatch!");
  
  CrystVector_REAL a(A.rows());
  
  for(int i=0; i<a.numElements(); i++) {
    a(i) = 0.;
    for(int k=0; k<A.cols(); k++) a(i) += A(i,k)*b(k);
  }
  
  return a; 
}

} // namespace RotationTB 

namespace MStruct {

const ObjCryst::RefParType *gpRefParTypeScattDataCorrHKLIntensity=0;
const ObjCryst::RefParType *gpRefParTypeMaterialData=0;
const ObjCryst::RefParType *gpRefParTypeScattDataProfileSizeDistrib=0;

long NiftyStaticGlobalObjectsInitializer_MStruct::mCount=0;

#define MSTRUCT_NONE	        0
#define MSTRUCT_SC		1
#define MSTRUCT_BCC		2
#define MSTRUCT_FCC		4
#define MSTRUCT_HCP		6

// ReflStore
int ReflStore::find(const ReflData& data)const {
  int ind = -1;
  REAL dxmin = FLT_MAX;
  std::vector<ReflData>::const_iterator it=mvData.begin();
  for(unsigned int i=0;i<mvData.size();i++,it++)
    if((*it).H==data.H && (*it).K==data.K && (*it).L==data.L) {
      REAL dx = fabs((*it).x-data.x);
      if(dx<=dxmin) {dxmin=dx;ind=i;}
    }
  return ind;
}

void ReflStore::print(std::ostream& s)const
{
  s<<"   **** ReflStore ****"<<endl;
  std::vector<ReflData>::const_iterator it=mvData.begin();
  int i=0;
  for(std::vector<ReflData>::const_iterator it=mvData.begin();it!=mvData.end();it++,i++) {
    s<<setw(12)<<i<<":";
    s<<setw(4)<<(*it).H<<setw(4)<<(*it).K<<setw(4)<<(*it).L;
    s<<setw(12)<<(*it).x<<setw(12)<<(*it).data<<endl;
  }
}

////////////////////////////////////////////////////////////////////////
//
//    LSQNumObj
//
////////////////////////////////////////////////////////////////////////

ObjCryst::ObjRegistry< ObjCryst::LSQNumObj > GlobalLSQNumObj_SIGINT_Registry;
/// Global counter of handled SIGINT signals
unsigned int Global_SIGINThandled_counter = 0;

// CTRL+C signal handler for LSQNumObjects finds 
void LSQNumObj_SIGINT_Handler (int param)
{
	cout << "\t...CTRL+C signal received..."<<endl;
	
	Global_SIGINThandled_counter++;
	
	// Stop refinement of all registred LSQNumObj objects
	for(long i=0; i<GlobalLSQNumObj_SIGINT_Registry.GetNb(); i++) {
		ObjCryst::LSQNumObj &obj = GlobalLSQNumObj_SIGINT_Registry.GetObj(i);
		obj.StopAfterCycle();
	}
}

void LSQNumObj::Refine (int nbCycle,bool useLevenbergMarquardt,
                        const bool silent, const bool callBeginEndOptimization,
                        const float minChi2var)
{
	// register CTRL+C signal handler
	void (*prev_fn)(int);	
	prev_fn = signal (SIGINT, LSQNumObj_SIGINT_Handler);
	if (prev_fn==SIG_IGN) signal (SIGINT, SIG_IGN);
	// Register this object to enable SIGINT handler functionality  
	GlobalLSQNumObj_SIGINT_Registry.Register(*this);
	
	// run refinement
	ObjCryst::LSQNumObj::Refine (nbCycle, useLevenbergMarquardt, silent,
				     callBeginEndOptimization, minChi2var);
	
	// de-register CTRL+C signal handler
	signal (SIGINT, prev_fn);
	// DeRegister this object to disable SIGINT handler functionality
	GlobalLSQNumObj_SIGINT_Registry.DeRegister(*this);
}

////////////////////////////////////////////////////////////////////////
//
//    RandomOptimizationObj
//
////////////////////////////////////////////////////////////////////////
RandomOptimizationObj::RandomOptimizationObj(std::string objName)
  :ObjCryst::OptimizationObj(objName),mpLSQNumObj(NULL),mLSQnbCycle(1),
   mLSQuseLevenbergMarquardt(false),mLSQsilent(false),mLSQcallBeginEndOptimization(true),
   mLSQminChi2var(0.01),mpRefinableObjRandomised(NULL)
{}

void RandomOptimizationObj::SetLSQNumObj(ObjCryst::LSQNumObj *pLSQNumObj, int nbCycle,
					 bool useLevenbergMarquardt, const bool silent,
					 const bool callBeginEndOptimization, const float minChi2var)
{
  mpLSQNumObj = pLSQNumObj;
  mLSQnbCycle= nbCycle;
  mLSQuseLevenbergMarquardt = useLevenbergMarquardt;
  mLSQsilent = silent;
  mLSQcallBeginEndOptimization = callBeginEndOptimization;
  mLSQminChi2var = minChi2var;
}

void RandomOptimizationObj::SetRefinableObjRandomised(ObjCryst::RefinableObj *pRefinableObj)
{
  mpRefinableObjRandomised = pRefinableObj;
}

void RandomOptimizationObj::Optimize(long &nbSteps, const bool silent, const REAL finalcost, const REAL maxTime)
{
  if(mpLSQNumObj==NULL) {
    cerr << "MStruct::RandomOptimizationObj::Optimize(): No LSQNumObject set." << endl;
    return;
  }

  if(mpRefinableObjRandomised==NULL) {
    cerr << "MStruct::RandomOptimizationObj::Optimize(): No RefinableObj to randomise set." << endl;
    return;
  }

  // clear old results
  mpRefinableObjRandomised->EraseAllParamSet();
  mSetInfo.clear();
  
  // save the fixed-refined state for all params of the randomised object
  vector<bool> vFixed(false,mpRefinableObjRandomised->GetNbPar());
  for(int ipar=0; ipar<mpRefinableObjRandomised->GetNbPar(); ipar++)
    vFixed[ipar] = mpRefinableObjRandomised->GetPar(ipar).IsFixed();
  
  // results will be consecutively saved in a file
  WriteCurrentParamSetToFile("randOptSet1-current.dat", false);

  // run optimization - random search for the best sets of configurations
  for(int istep=0; istep<nbSteps; istep++) {
    
    if(!silent) cout <<"--- Running trial step nb. "<<(istep+1)<<"/"<<nbSteps<<" ---\n";

    // randomise configuration
    mpRefinableObjRandomised->GlobalOptRandomMove(1.);
    
    // fix all parameters in the randomised RefinableObj before refinement
    //    mpRefinableObjRandomised->FixAllPar();

    // run LSQ refinement
    mpLSQNumObj->Refine(mLSQnbCycle,mLSQuseLevenbergMarquardt,mLSQsilent,
			mLSQcallBeginEndOptimization,mLSQminChi2var);
    
    // during the refinement some parameters can be fixed by chance, restore their fixed-refined state
    for(int ipar=0; ipar<mpRefinableObjRandomised->GetNbPar(); ipar++)
      mpRefinableObjRandomised->GetPar(ipar).SetIsFixed(vFixed[ipar]);

    // unfix all parameters in the randomised RefinableObj before a par. set is created
    //    mpRefinableObjRandomised->UnFixAllPar();

    // save configuration
    unsigned long id = mpRefinableObjRandomised->CreateParamSet();

    // get ChiSquare value
    REAL chiSq = mpLSQNumObj->ChiSquare();
    
    // save Param.Set ID and chiSq value
    mSetInfo.push_back( std::make_pair( id, chiSq ) );

    // write current parameter set into the file
    WriteCurrentParamSetToFile("randOptSet1-current.dat", true);

  } // for istep

}

void RandomOptimizationObj::WriteResultsToFile(const string &filename) const
{
  std::ofstream os(filename.c_str());
  
  // nb of parameters in the set
  const int nbpar = mpRefinableObjRandomised->GetNbPar();

  // print header
  os<<"#"<<setw(8)<<"id"<<setw(14)<<"chiSq";
  for(int i=0; i<nbpar; i++)
    os<<setw(14)<<mpRefinableObjRandomised->GetPar(i).GetName();
  os<<"\n";
      
  // print values
  for(int iset=0; iset<mSetInfo.size(); iset++) {
    
    unsigned long id = mSetInfo[iset].first;
    REAL chiSq = mSetInfo[iset].second;

    // unfix all parameters in the randomised RefinableObj before a par. set is restored
    //mpRefinableObjRandomised->UnFixAllPar();

    // restore parameter set
    mpRefinableObjRandomised->RestoreParamSet(id);
    
    os<<setw(9)<<id<<setw(14)<<scientific<<showpoint<<setprecision(3)<<chiSq;

    for(int i=0; i<nbpar; i++)
      os<<setw(14)<<mpRefinableObjRandomised->GetPar(i).GetHumanValue();
	
    os<<"\n";

  }

  // fix all parameters in the randomised RefinableObj after we finished
  //mpRefinableObjRandomised->FixAllPar();

  os.close();
}

void RandomOptimizationObj::WriteCurrentParamSetToFile(const string &filename, const bool append)
{
  std::ofstream os;

  if (append)
    os.open(filename.c_str(),ios_base::app);
  else
    os.open(filename.c_str());
  
  // nb of parameters in the set
  const int nbpar = mpRefinableObjRandomised->GetNbPar();

  // print header
  if (!append) {
    os<<"#"<<setw(8)<<"id"<<setw(14)<<"chiSq";
    for(int i=0; i<nbpar; i++)
      os<<setw(14)<<mpRefinableObjRandomised->GetPar(i).GetName();
    os<<"\n";
  }
      
  // print values
  
  if( mSetInfo.size()>0 ) {
    
    unsigned long id = mSetInfo[mSetInfo.size()-1].first;
    REAL chiSq = mpLSQNumObj->ChiSquare();
    
    os<<setw(9)<<id<<setw(14)<<scientific<<showpoint<<setprecision(4)<<chiSq;

    for(int i=0; i<nbpar; i++)
      os<<setw(14)<<mpRefinableObjRandomised->GetPar(i).GetHumanValue();
	
    os<<"\n";

  }

  os.close();
}

////////////////////////////////////////////////////////////////////////
//
//    PowderPatternBackgroundBase
//
////////////////////////////////////////////////////////////////////////

PowderPatternBackgroundBase::PowderPatternBackgroundBase()
:mMaxSinThetaOvLambda(10.),mUseVariableSlitIntensityCorr(false),mPowderPatternSinTheta(0)
{}

PowderPatternBackgroundBase::PowderPatternBackgroundBase(const PowderPatternBackgroundBase &old)
:ObjCryst::PowderPatternComponent(old),
mMaxSinThetaOvLambda(old.mMaxSinThetaOvLambda),
mUseVariableSlitIntensityCorr(old.mUseVariableSlitIntensityCorr),
mPowderPatternSinTheta(0)
{}

const string& PowderPatternBackgroundBase::GetClassName()const
{
	static const string className = "MStruct::PowderPatternBackgroundBase";
	return className;
}

void PowderPatternBackgroundBase::SetParentPowderPattern(ObjCryst::PowderPattern &s)
{
  if(mpParentPowderPattern!=0) 
    mClockMaster.RemoveChild(mpParentPowderPattern->GetIntegratedProfileLimitsClock());
  mpParentPowderPattern = &s;
  mClockMaster.AddChild(mpParentPowderPattern->GetIntegratedProfileLimitsClock());
  mClockMaster.AddChild(mpParentPowderPattern->GetClockPowderPatternPar());
  //mClockMaster.AddChild(mpParentPowderPattern->GetClockPowderPatternXCorr());
  //mClockMaster.AddChild(mpParentPowderPattern->GetClockPowderPatternRadiation());
}

const CrystVector_REAL&	PowderPatternBackgroundBase::GetPowderPatternCalc()const
{
	this->CalcPowderPattern();
	return mPowderPatternCalc;
}

pair< const CrystVector_REAL *, const ObjCryst::RefinableObjClock * > PowderPatternBackgroundBase::GetPowderPatternIntegratedCalc()const
{
	VFN_DEBUG_MESSAGE("MStruct::PowderPatternBackgroundBase::GetPowderPatternIntegratedCalc()",3)
	this->CalcPowderPatternIntegrated();
	return make_pair(&mPowderPatternIntegratedCalc,&mClockPowderPatternIntegratedCalc);
}

const CrystVector_REAL& PowderPatternBackgroundBase::GetPowderPatternCalcVariance()const
{
	this->CalcPowderPattern();
	return mPowderPatternCalcVariance;
}

pair< const CrystVector_REAL *, const RefinableObjClock * > PowderPatternBackgroundBase::GetPowderPatternIntegratedCalcVariance()const
{
	VFN_DEBUG_MESSAGE("MStruct::PowderPatternBackgroundBase::GetPowderPatternIntegratedCalcVariance()",3)
	this->CalcPowderPatternIntegrated();
	return make_pair(&mPowderPatternIntegratedCalcVariance,
                   &mClockPowderPatternIntegratedVarianceCalc);
}

bool PowderPatternBackgroundBase::HasPowderPatternCalcVariance()const
{
   #ifdef USE_BACKGROUND_MAXLIKE_ERROR
   return true;
   #else
   return false;
   #endif
}

void PowderPatternBackgroundBase::UseVariableSlitIntensityCorr(const bool b)
{
	mUseVariableSlitIntensityCorr = b;
	mClockMaster.Click();
}

const CrystVector_REAL& PowderPatternBackgroundBase::GetPowderPatternSinTheta()const
{
	VFN_DEBUG_MESSAGE("MStruct::PowderPatternBackgroundBase::GetPowderPatternSinTheta():Begin",3)
	TAU_PROFILE("MStruct::PowderPatternBackgroundBase::GetPowderPatternSinTheta()","void ()",TAU_DEFAULT);
	
	// TODO:: TOF data
	
	if(mpParentPowderPattern==NULL) {
		cerr << "< MStruct::PowderPatternBackgroundBase::GetPowderPatternSinTheta()\n";
		cerr << "Logical error: No parent PowderPattern object.\n >" << endl; 
		throw ObjCrystException("MStruct::PowderPatternBackgroundBase::GetPowderPatternSinTheta(): Program error.");
	}
	
	// check for changes
	if( mClockPowderPatternSinTheta<mpParentPowderPattern->GetClockPowderPatternPar() ||
		  mClockPowderPatternSinTheta<mpParentPowderPattern->GetClockPowderPatternRadiation() ) {
		
		VFN_DEBUG_MESSAGE("MStruct::PowderPatternBackgroundBase::GetPowderPatternSinTheta():Recalculating",3)
		// recalculate sin(theta) values
		const unsigned long nb = mpParentPowderPattern->GetNbPoint();
		
		mPowderPatternSinTheta.resize(nb);
		
		const REAL *p1 = mpParentPowderPattern->GetPowderPatternX().data();
		REAL *p2 = mPowderPatternSinTheta.data();
		
		for(unsigned long i=0; i<nb; i++) { *p2 = sin(0.5*(*p1)); p1++; p2++; }
		
		mClockPowderPatternSinTheta.Click();
	}
	
	VFN_DEBUG_MESSAGE("MStruct::PowderPatternBackgroundBase::GetPowderPatternSinTheta():End",3)
	return mPowderPatternSinTheta;
}

const ObjCryst::RefinableObjClock& PowderPatternBackgroundBase::GetClockPowderPatternSinTheta()const 
{
	return mClockPowderPatternSinTheta;
}	
	
void PowderPatternBackgroundBase::CalcPowderPatternIntegrated()const
{
	// NOTE:: This is only almost a copy-paste from the original source
	if(mClockPowderPatternIntegratedCalc>mClockMaster) return; // TODO:: Integrated ??? See original code in ObjCryst::PowderPatternBackground::CalcPowderPatternIntegrated()

	this->CalcPowderPattern();// :TODO: Optimize
  if( (mClockPowderPatternIntegratedCalc>mClockPowderPatternCalc)
     	 &&(mClockPowderPatternIntegratedCalc>mpParentPowderPattern->GetIntegratedProfileLimitsClock()))
         return;

	VFN_DEBUG_ENTRY("MStruct::PowderPatternBackgroundBase::CalcPowderPatternIntegrated()",3)
	TAU_PROFILE("MStruct::PowderPatternBackgroundBase::CalcPowderPatternIntegrated()","void ()",TAU_DEFAULT);
	const CrystVector_long *pMin=&(mpParentPowderPattern->GetIntegratedProfileMin());
	const CrystVector_long *pMax=&(mpParentPowderPattern->GetIntegratedProfileMax());

	const long numInterval=pMin->numElements();
	mPowderPatternIntegratedCalc.resize(numInterval);
	REAL * RESTRICT p2=mPowderPatternIntegratedCalc.data();
	for(int j=0;j<numInterval;j++)
	{
		const long max=(*pMax)(j);
		const REAL * RESTRICT p1=mPowderPatternCalc.data()+(*pMin)(j);
		*p2=0;           
		for(int k=(*pMin)(j);k<=max;k++) *p2 += *p1++;
		p2++;
  }
	#ifdef USE_BACKGROUND_MAXLIKE_ERROR
	mPowderPatternIntegratedCalcVariance.resize(numInterval);
	p2=mPowderPatternIntegratedCalcVariance.data();
	for(int j=0;j<numInterval;j++)
	{
		const long max=(*pMax)(j);
		const REAL *p1=mPowderPatternCalcVariance.data()+(*pMin)(j);
		*p2=0;           
		for(int k=(*pMin)(j);k<=max;k++) *p2 += *p1++;
		p2++;
	}
	mClockPowderPatternIntegratedVarianceCalc.Click();
	#endif
	mClockPowderPatternIntegratedCalc.Click();
	VFN_DEBUG_EXIT("MStruct::PowderPatternBackgroundBase::CalcPowderPatternIntegrated():End",3)
}

const CrystVector_long& PowderPatternBackgroundBase::GetBraggLimits()const
{
  // no integration interval for the background
  mIntegratedReflLimits.resize(0);
  return mIntegratedReflLimits;
}

void PowderPatternBackgroundBase::SetMaxSinThetaOvLambda(const REAL max)
{
	mMaxSinThetaOvLambda = max;
	mClockMaster.Click();
}

void PowderPatternBackgroundBase::Prepare()
{}

////////////////////////////////////////////////////////////////////////
//
//    PowderPatternBackgroundInvX
//
////////////////////////////////////////////////////////////////////////

PowderPatternBackgroundInvX::PowderPatternBackgroundInvX()
:mXFunctionType(FUNCTION_OF_X)
{
	mIsScalable = true;
}

PowderPatternBackgroundInvX::PowderPatternBackgroundInvX(const PowderPatternBackgroundInvX &old)
:PowderPatternBackgroundBase(old)
{
	mIsScalable = true;
}

const string& PowderPatternBackgroundInvX::GetClassName()const
{
	static const string className = "MStruct::PowderPatternBackgroundInvX";
	return className;
}

void PowderPatternBackgroundInvX::SetXFunctionType(const int type)
{
	if( (type!=FUNCTION_OF_X) && (type!=FUNCTION_OF_SIN_TH) ) {
		cerr << "< MStruct::PowderPatternBackgroundInvX::SetXFunctionType(int)\n";
		cerr << "Bad input argument 'type': " << type << ".\n >" << endl; 
		throw ObjCrystException("MStruct::PowderPatternBackgroundInvX::SetXFunctionType(int): Bad input argument.");
	}
	
	mXFunctionType = type;
	mClockMaster.Click();
}

void PowderPatternBackgroundInvX::CalcPowderPattern()const
{
	if (mClockPowderPatternCalc>mClockMaster) return;
   
  // this powder pattern component has no parameters (it is only scalable).
  // All necessary clocks are included in the master clocks
  
	TAU_PROFILE("MStruct::PowderPatternBackgroundInvX::CalcPowderPattern()","void ()",TAU_DEFAULT);
	VFN_DEBUG_MESSAGE("MStruct::PowderPatternBackgroundInvX::CalcPowderPattern()",3);
   
	const unsigned long nb = mpParentPowderPattern->GetNbPoint();
	mPowderPatternCalc.resize(nb);
	
	REAL *p2 = mPowderPatternCalc.data();
	
	try {
		if(mXFunctionType==FUNCTION_OF_X) {
			const REAL *p1 = mpParentPowderPattern->GetPowderPatternX().data();
			if(!mUseVariableSlitIntensityCorr)
				for(unsigned long i=0; i<nb; i++) { *p2 = 1./(*p1); p1++; p2++; }
			else {
				const REAL *p3 = GetPowderPatternSinTheta().data();
				for(unsigned long i=0; i<nb; i++) { *p2 = (abs(*p1)<1e-6) ? 0.5*(1.-pow(*p1,2)/6.) : (*p3)/(*p1); p1++; p2++; p3++; }
			}
		} else if(mXFunctionType==FUNCTION_OF_SIN_TH) {
			if(!mUseVariableSlitIntensityCorr) {
				const REAL *p3 = GetPowderPatternSinTheta().data();
				for(unsigned long i=0; i<nb; i++) { *p2 = 1./(*p3); p2++; p3++; }
			}	else
				for(unsigned long i=0; i<nb; i++) { *p2 = 1.; p2++; }
		} else {
			cerr << "< MStruct::PowderPatternBackgroundInvX::CalcPowderPattern()\n";
			cerr << "Unexpected/unsupported x-function type: "<< mXFunctionType <<".\n >" << endl; 
			throw ObjCrystException("MStruct::PowderPatternBackgroundInvX::CalcPowderPattern(): Program error.");
		}
	}
	catch (std::exception &e) {
		cerr << "< MStruct::PowderPatternBackgroundInvX::CalcPowderPattern()\n";
		cerr << "Unexpected exception: " << e.what() << "\n";
		cerr << "Unexpected exception thrown during calcualtion of the powder pattern background.\n >" << endl; 
		throw ObjCrystException("MStruct::PowderPatternBackgroundInvX::CalcPowderPattern(): Program error.");
	}
	
	VFN_DEBUG_MESSAGE("MStruct::PowderPatternBackgroundInvX::CalcPowderPattern()",3);
	#ifdef USE_BACKGROUND_MAXLIKE_ERROR
	{
		mPowderPatternCalcVariance.resize(nb);
		const REAL step=mModelVariance*mModelVariance/(REAL)nbPoint;
		REAL var=0;
		REAL *p=mPowderPatternCalcVariance.data();
		for(long i=0;i<nb;i++) {*p++ = var;var +=step;}
	}
	mClockPowderPatternVarianceCalc.Click();
	#endif
	mClockPowderPatternCalc.Click();
	VFN_DEBUG_MESSAGE("MStruct::PowderPatternBackgroundInvX::CalcPowderPattern():End",3);	
}

////////////////////////////////////////////////////////////////////////
//
//    PowderPatternBackgroundChebyshev
//
////////////////////////////////////////////////////////////////////////

PowderPatternBackgroundChebyshev::PowderPatternBackgroundChebyshev()
:mChebyshevCoef(0), mChebyshevPolynomials(0,0)
{
	mIsScalable = false,	
	Init();
}

PowderPatternBackgroundChebyshev::PowderPatternBackgroundChebyshev(const PowderPatternBackgroundChebyshev &old)
:PowderPatternBackgroundBase(old),mChebyshevCoef(old.mChebyshevCoef),mChebyshevPolynomials(0,0)
{
	mIsScalable = false,
	Init();
}

const string& PowderPatternBackgroundChebyshev::GetClassName()const
{
	static const string className = "MStruct::PowderPatternBackgroundChebyshev";
	return className;
}

void PowderPatternBackgroundChebyshev::SetCoefficients(const CrystVector_REAL &coef)
{
	VFN_DEBUG_MESSAGE("MStruct::PowderPatternBackgroundChebyshev::SetCoefficients(...):Begin",11);
	
	bool rebuild = ( coef.numElements() != mChebyshevCoef.numElements() );
	
	if (rebuild) {
		// Number of params. has changed. The old params. have to be removed.
		try {
			for(int i=0; i<mChebyshevCoef.numElements(); i++) {
				ObjCryst::RefinablePar &par = this->GetPar(&mChebyshevCoef(i));
				this->RemovePar(&par);
			}
		}
		catch (std::exception &e) {
			cerr << "< MStruct::PowderPatternBackgroundChebyshev::SetCoefficients()\n";
			cerr << "Unexpected exception: " << e.what() << "\n";
			cerr << "Unexpected exception thrown during removing old parameters from the object.\n >" << endl; 
			throw ObjCrystException("MStruct::PowderPatternBackgroundChebyshev::SetCoefficients(): Program error.");
		}
	}
	
	// set new coeficients
	mChebyshevCoef = coef;
	
	mClockMaster.Click(); // Values of the used Chebyshev polynomials has to be recomputed.
	
	if (rebuild) {
		// Number of params. has changed. Object has to be reinitialsied.
		Init();
		// Order of the Chebyshev polynomials changed - force Chebyshev polynomials values recalc
		mClockChebyshevPolynomialsCalc.Reset();
	}
	
	VFN_DEBUG_MESSAGE("MStruct::PowderPatternBackgroundChebyshev::SetCoefficients(...):End",11);
}

const CrystVector_REAL& PowderPatternBackgroundChebyshev::GetCoefficients()const
{
	return mChebyshevCoef;
}

void PowderPatternBackgroundChebyshev::SetXFunctionType(const int type)
{
	if( (type!=FUNCTION_OF_X) && (type!=FUNCTION_OF_SIN_TH) ) {
		cerr << "< MStruct::PowderPatternBackgroundChebyshev::SetXFunctionType(int)\n";
		cerr << "Bad input argument 'type': " << type << ".\n >" << endl; 
		throw ObjCrystException("MStruct::PowderPatternBackgroundChebyshev::SetXFunctionType(int): Bad input argument.");
	}
	
	mXFunctionType = type;
	mClockMaster.Click();
	// Order of the Chebyshev polynomials changed - force Chebyshev polynomials values recalc
	mClockChebyshevPolynomialsCalc.Reset();
}

const CrystVector_REAL& PowderPatternBackgroundChebyshev::GetLSQDeriv(const unsigned int nfunc, RefinablePar& par)
{
	// check LSQ func. nb. 
	if(nfunc!=0) {
		cerr << "< MStruct::PowderPatternBackgroundChebyshev::GetLSQDeriv(...)\n";
		cerr << "LSQ func. nb.: " << nfunc << ", Parameter: " << par.GetName() << "\n";
		par.Print();
		cerr << "Unexpected/unsupported LSQ func. nb.\n >" << endl; 
		throw ObjCrystException("MStruct::PowderPatternBackgroundChebyshev::GetLSQDeriv(...): Bad argument.");
  }
  
  // find parameter
  vector< ObjCryst::RefinablePar* >::iterator result;
  result = find( mvpRefPar.begin(), mvpRefPar.end(), &par );
 
  if( result == mvpRefPar.end() ) {
  	cerr << "< MStruct::PowderPatternBackgroundChebyshev::GetLSQDeriv(...)\n";
		cerr << "LSQ func. nb.: " << nfunc << ", Parameter: " << par.GetName() << "\n";
		par.Print();
		cerr << "Specified RefinableParameter was not found in the list of parameters.\n >" << endl; 
		throw ObjCrystException("MStruct::PowderPatternBackgroundChebyshev::GetLSQDeriv(...): Bad argument.");
  }
  
  // identify parameter (Background_Coef_nb), (nb)?
  string str = string(par.GetName());
  string::size_type loc;
  loc = str.find("Background_Coef_", 0); // Background_Coef_*
  if( loc == string::npos ) {
  	cerr << "< MStruct::PowderPatternBackgroundChebyshev::GetLSQDeriv(...)\n";
		cerr << "LSQ func. nb.: " << nfunc << ", Parameter: " << par.GetName() << "\n";
		par.Print();
		cerr << "Specified RefinableParameter was not identified.\n >" << endl; 
		throw ObjCrystException("MStruct::PowderPatternBackgroundChebyshev::GetLSQDeriv(...): Bad argument.");
  }
  
  // using rest of the parameter name
  istringstream ss(str.substr(loc+16));
  // Chebyshev coefficient number
  int coef_nb = -1;
  ss >> coef_nb;
  
  // check the Chebyshev coefficient number
  if( (coef_nb<0) || (coef_nb>=mChebyshevCoef.numElements()) ) {
  	cerr << "< MStruct::PowderPatternBackgroundChebyshev::GetLSQDeriv(...)\n";
		cerr << "LSQ func. nb.: " << nfunc << ", Parameter: " << par.GetName() << "\n";
		par.Print();
		cerr << "Identified Coefficient number ("<< coef_nb <<") cound not be accepted.\n";
		cerr << "Only numbers 0 ... " << mChebyshevCoef.numElements()-1 << " are accepted.\n >" << endl; 
		throw ObjCrystException("MStruct::PowderPatternBackgroundChebyshev::GetLSQDeriv(...): Bad argument.");
  }
  
	VFN_DEBUG_MESSAGE("PowderPatternBackgroundChebyshev::GetLSQDeriv(...): "<<
		      "Background coefficient nb: "<< coef_nb,11)
		      
  // calculate derivative (return the Chebyshev polynomial)
  const unsigned long nb = mpParentPowderPattern->GetNbPoint();
  mLSQDeriv.resize(nb);
  
  // zero order coef.
  if(coef_nb==0) mLSQDeriv = 1.; // not stored in mChebyshevPolynomials
  // higher order coef.
  else {
  	CalcChebyshevPolynomials(); // update (if necessary)
  	const REAL *p1 = mChebyshevPolynomials.data() + (coef_nb-1)*nb; // row with data
  	REAL *p2 = mLSQDeriv.data();
  	for(unsigned long i=0; i<nb; i++) { *p2 = (*p1); p1++; p2++; }
  }
  
  // variable slit intesity corr.
	if( mUseVariableSlitIntensityCorr )
		mLSQDeriv *= GetPowderPatternSinTheta();
		
	return mLSQDeriv;
}

void PowderPatternBackgroundChebyshev::CalcChebyshevPolynomials()const
{
	TAU_PROFILE("MStruct::PowderPatternBackgroundChebyshev::CalcChebyshevPolynomials()","void ()",TAU_DEFAULT);
	VFN_DEBUG_MESSAGE("MStruct::PowderPatternBackgroundChebyshev::CalcChebyshevPolynomials():Begin",3);
	
	// Recalaculate Chebyshev polynomials values if necessary.
	//   ChebyshevPolynomials have to be recalculated also in the case that an order of the highest
	//   used polynomial or the 'X'-function type is changed. It is assumend that
	//   recalculation is forced by mClockChebyshevPolynomialsCalc-reset in such cases.
  if ( (mClockChebyshevPolynomialsCalc < mpParentPowderPattern->GetClockPowderPatternPar()) ||
  	   (mClockChebyshevPolynomialsCalc < mpParentPowderPattern->GetClockPowderPatternRadiation()) ) {
  	VFN_DEBUG_MESSAGE("MStruct::PowderPatternBackgroundChebyshev::CalcChebyshevPolynomials():Recalculating Chebyshev polynomials matrix",3);
  	
  	const unsigned long nb = mpParentPowderPattern->GetNbPoint();
		const unsigned long M = mChebyshevCoef.numElements();
  	
  	// values of the zero-order Chybyshev polynomial are not stored
  	mChebyshevPolynomials.resize(M-1,nb);
  		
  	if (M>1) {
  		// the first order (x)
  		{
  			if(mXFunctionType==FUNCTION_OF_X) {
  				const REAL min = 0.; //mpParentPowderPattern->GetPowderPatternX().GetXMin()/2;
  				const REAL max = 2*M_PI; //mpParentPowderPattern->GetPowderPatternX().GetXMax()/2;
  				// TODO:: make usable also for TOF data
  				const REAL s = 0.5*(max-min);
  				const REAL x0 = 0.5*(max+min);
  				const REAL *p1 = mpParentPowderPattern->GetPowderPatternX().data();
 					REAL *p2 = mChebyshevPolynomials.data();
  				for(unsigned long i=0; i<nb; i++) { *p2 = ((*p1)-x0)/s; p1++; p2++; }
  			}	else if(mXFunctionType==FUNCTION_OF_SIN_TH) {
  				const REAL *p1 = GetPowderPatternSinTheta().data();
 					REAL *p2 = mChebyshevPolynomials.data();
  				for(unsigned long i=0; i<nb; i++) { *p2 = (*p1); p1++; p2++; }
  			} else {
  				cerr << "< MStruct::PowderPatternBackgroundChebyshev::CalcChebyshevPolynomials()\n";
					cerr << "Unexpected/unsupported x-function type: "<< mXFunctionType <<".\n >" << endl; 
					throw ObjCrystException("MStruct::PowderPatternBackgroundChebyshev::CalcChebyshevPolynomials(): Program error.");
  			}
  		}
  		// the second order (2 * x^2 - 1)
  		if (M>2) {
  			const REAL *p1 = mChebyshevPolynomials.data();
 				REAL *p2 = mChebyshevPolynomials.data() + nb;
  			for(unsigned long i=0; i<nb; i++) { *p2 = 2.*pow((*p1),2) - 1.; p1++; p2++; }
  		}
  		// higher orders (recurrence formula: Tn+1(x) = 2 * x * Tn(x) - Tn-1(x))
  		for (unsigned long m=3; m<M; m++) {
  			const REAL *p0 = mChebyshevPolynomials.data();
  			const REAL *p1 = mChebyshevPolynomials.data()+(m-3)*nb;
  			const REAL *p2 = mChebyshevPolynomials.data()+(m-2)*nb;
  			REAL *p3 = mChebyshevPolynomials.data()+(m-1)*nb;
  			for(unsigned long i=0; i<nb; i++) { *p3 = 2*(*p0)*(*p2) - (*p1); p1++; p2++; p3++; }
  		}
  	}
  		
  	mClockChebyshevPolynomialsCalc.Click();
  }
	
	VFN_DEBUG_MESSAGE("MStruct::PowderPatternBackgroundChebyshev::CalcChebyshevPolynomials():End",3);
}

void PowderPatternBackgroundChebyshev::CalcPowderPattern()const
{
	if (mClockPowderPatternCalc>mClockMaster) return;
  
	TAU_PROFILE("MStruct::PowderPatternBackgroundChebyshev::CalcPowderPattern()","void ()",TAU_DEFAULT);
	VFN_DEBUG_MESSAGE("MStruct::PowderPatternBackgroundChebyshev::CalcPowderPattern():Begin",3);

	try {
	  
  	const unsigned long nb = mpParentPowderPattern->GetNbPoint();
		mPowderPatternCalc.resize(nb);

		const unsigned long M = mChebyshevCoef.numElements();
		
		// Recalaculate Chebyshev polynomials values if necessary.
  	CalcChebyshevPolynomials();
	
		// calcualtion
	
		mPowderPatternCalc = 0.;
		
		// zero-order
		if (M>0) {
			REAL *p2 = mPowderPatternCalc.data();
			if (abs(mChebyshevCoef(0))>1e-6)
				for(unsigned long i=0; i<nb; i++) {	*p2 += mChebyshevCoef(0); p2++; }
		}
		
		// higher orders
		for (unsigned long m=1; m<M; m++) {
			const REAL *p1 = mChebyshevPolynomials.data()+(m-1)*nb;
			REAL *p2 = mPowderPatternCalc.data();
			if (abs(mChebyshevCoef(m))>1e-6)
				for(unsigned long i=0; i<nb; i++) {	*p2 += (*p1)*mChebyshevCoef(m); p1++; p2++; }
		}
		
		// variable slid intesity corr.
		if ( mUseVariableSlitIntensityCorr && (M>0) )
			mPowderPatternCalc *= GetPowderPatternSinTheta();
		
	}
	catch (std::exception &e) {
		cerr << "< MStruct::PowderPatternBackgroundChebyshev::CalcPowderPattern()\n";
		cerr << "Unexpected exception: " << e.what() << "\n";
		cerr << "Unexpected exception thrown during calcualtion of the powder pattern background.\n >" << endl; 
		throw;
	}
	
	VFN_DEBUG_MESSAGE("MStruct::PowderPatternBackgroundChebyshev::CalcPowderPattern()",3);
	#ifdef USE_BACKGROUND_MAXLIKE_ERROR
	{
		mPowderPatternCalcVariance.resize(nb);
		const REAL step=mModelVariance*mModelVariance/(REAL)nbPoint;
		REAL var=0;
		REAL *p=mPowderPatternCalcVariance.data();
		for(long i=0;i<nb;i++) {*p++ = var;var +=step;}
	}
	mClockPowderPatternVarianceCalc.Click();
	#endif
	mClockPowderPatternCalc.Click();
	VFN_DEBUG_MESSAGE("MStruct::PowderPatternBackgroundChebyshev::CalcPowderPattern():End",3);	
}

void PowderPatternBackgroundChebyshev::Init()
{
	// add parameters (coefficients of the Chebyshev polynomials)
	for(int i=0; i<mChebyshevCoef.numElements(); i++) {
		// create parameter name
		ostringstream ss;
		ss << "Background_Coef_" << i;
		// create a new parameter
		RefinablePar tmp(ss.str(),&mChebyshevCoef(i),-100000.,100000.,
										 gpRefParTypeObjCryst,REFPAR_DERIV_STEP_ABSOLUTE,
                     false,true,true,false,1.);
		tmp.AssignClock(mClockMaster);
		tmp.SetDerivStep(1.0);
		//tmp.SetGlobalOptimStep(.05);
		this->AddPar(tmp);
  }
}

// PowderPattern
PowderPattern::PowderPattern()
:mOmega(-1.)
{
  // Add scales to relevant params
  mClockMaster.AddChild(mClockScaleFactor);
  // Set Integration Option - not Integrated
  mOptProfileIntegration.SetChoice(1);
}

const CrystVector_REAL& PowderPattern::GetLSQCalc(const unsigned int n) const
{
  const long nbExclude=mExcludedRegionMinX.numElements();

  if(nbExclude==0 && mAdditionalLSQObjRegistry.GetNb()==0) {
  	return ObjCryst::PowderPattern::GetLSQCalc(n);
  } else if(nbExclude==0) {
 		mLSQCalcNotExcluded = ObjCryst::PowderPattern::GetLSQCalc(n);
  } else {
  	
	  const REAL *p1=ObjCryst::PowderPattern::GetLSQCalc(n).data();
	
	  if(mClockPowderPatternCalc>mClockLSQCalcNotExcluded) {
	
	    unsigned long maxPoints=mNbPointUsed;
	    mLSQCalcNotExcluded.resize(maxPoints);
	    
	    // copy data from not excluded regions
	    REAL *p2= mLSQCalcNotExcluded.data();
	    unsigned long min,max;
	    unsigned long i=0,nbPointsUsedNotExcluded=0;
	    for(int j=0;j<nbExclude;j++) {
	      min=(unsigned long)floor(this->X2Pixel(mExcludedRegionMinX(j)));
	      max=(unsigned long)ceil (this->X2Pixel(mExcludedRegionMaxX(j)));
	      if(min>maxPoints) break;
	      if(max>maxPoints)max=maxPoints;
	      for(;i<min;i++)//! min is the *beginning* of the excluded region !
	      {
					*p2 = *p1; p1++; p2++; nbPointsUsedNotExcluded++;
	      }
	      p1 += max-i;
	      i  += max-i;
	    }
	    for(;i<maxPoints;i++) {
	      *p2 = *p1; p1++; p2++; nbPointsUsedNotExcluded++;
	    }
	    mLSQCalcNotExcluded.resizeAndPreserve(nbPointsUsedNotExcluded);
	    
	    mClockLSQCalcNotExcluded.Click();
	  }
	 
 	} // else (nbExclude==0)

	// contributions of additional LSQ func. objects
	
	std::vector< const CrystVector_REAL* > vLSQCalc(mAdditionalLSQObjRegistry.GetNb());
	long nn = 0;
	for(int iobj=0; iobj<mAdditionalLSQObjRegistry.GetNb(); iobj++) {
		vLSQCalc[iobj] = &mAdditionalLSQObjRegistry.GetObj(iobj).GetLSQCalc(n);
		nn += vLSQCalc[iobj]->numElements();
	}
	
	long nn0 = mLSQCalcNotExcluded.numElements();
	
	mLSQCalcNotExcluded.resizeAndPreserve(nn0+nn);
	
	// copy data
	REAL *p2 = mLSQCalcNotExcluded.data() + nn0;
	for(int iobj=0; iobj<mAdditionalLSQObjRegistry.GetNb(); iobj++) {
		const REAL *p1 = vLSQCalc[iobj]->data();
		for(int i=0; i<vLSQCalc[iobj]->numElements(); i++) *p2++ = *p1++;
	}
	
  return mLSQCalcNotExcluded;
}

const CrystVector_REAL& PowderPattern::GetLSQObs(const unsigned int n) const
{
  const long nbExclude=mExcludedRegionMinX.numElements();

  if(nbExclude==0 && mAdditionalLSQObjRegistry.GetNb()==0) {
    return ObjCryst::PowderPattern::GetLSQObs(n);
  } if(nbExclude==0) {
  	mLSQObsNotExcluded = ObjCryst::PowderPattern::GetLSQObs(n);
  } else {
  
	  const REAL *p1=ObjCryst::PowderPattern::GetLSQObs(n).data();
	
	  if(mClockPowderPatternPar>mClockLSQObsNotExcluded) {
	
	    unsigned long maxPoints=mNbPointUsed;
	    mLSQObsNotExcluded.resize(maxPoints);
	    
	    // copy data from not excluded regions
	    REAL *p2=mLSQObsNotExcluded.data();
	    unsigned long min,max;
	    unsigned long i=0,nbPointsUsedNotExcluded=0;
	    for(int j=0;j<nbExclude;j++) {
	      min=(unsigned long)floor(this->X2Pixel(mExcludedRegionMinX(j)));
	      max=(unsigned long)ceil (this->X2Pixel(mExcludedRegionMaxX(j)));
	      if(min>maxPoints) break;
	      if(max>maxPoints)max=maxPoints;
	      for(;i<min;i++)//! min is the *beginning* of the excluded region !
	      {
					*p2 = *p1; p1++; p2++; nbPointsUsedNotExcluded++;
	      }
	      p1 += max-i;
	      i  += max-i;
	    }
	    for(;i<maxPoints;i++) {
	      *p2 = *p1; p1++; p2++; nbPointsUsedNotExcluded++;
	    }
	    mLSQObsNotExcluded.resizeAndPreserve(nbPointsUsedNotExcluded);
	    
	    mClockLSQObsNotExcluded.Click();
	  }
  
  } // else(nbExclude==0)
  
  // contributions of additional LSQ func. objects
	
	std::vector< const CrystVector_REAL* > vLSQObs(mAdditionalLSQObjRegistry.GetNb());
	long nn = 0;
	for(int iobj=0; iobj<mAdditionalLSQObjRegistry.GetNb(); iobj++) {
		vLSQObs[iobj] = &mAdditionalLSQObjRegistry.GetObj(iobj).GetLSQObs(n);
		nn += vLSQObs[iobj]->numElements();
	}
	
	long nn0 = mLSQObsNotExcluded.numElements();
	
	mLSQObsNotExcluded.resizeAndPreserve(nn0+nn);
	
	// copy data
	REAL *p2 = mLSQObsNotExcluded.data() + nn0;
	for(int iobj=0; iobj<mAdditionalLSQObjRegistry.GetNb(); iobj++) {
		const REAL *p1 = vLSQObs[iobj]->data();
		for(int i=0; i<vLSQObs[iobj]->numElements(); i++) *p2++ = *p1++;
	}
	
  return mLSQObsNotExcluded;
}

const CrystVector_REAL& PowderPattern::GetLSQWeight(const unsigned int n) const
{
  const long nbExclude=mExcludedRegionMinX.numElements();

  if(nbExclude==0 && mAdditionalLSQObjRegistry.GetNb()==0) {
    return ObjCryst::PowderPattern::GetLSQWeight(n);
  } else if(nbExclude==0) {
  	mLSQWeightNotExcluded = ObjCryst::PowderPattern::GetLSQWeight(n);
  } else {
	  const REAL *p1=ObjCryst::PowderPattern::GetLSQWeight(n).data();
	
	  if(mClockPowderPatternPar>mClockLSQWeightNotExcluded) {
	
	    unsigned long maxPoints=mNbPointUsed;
	    mLSQWeightNotExcluded.resize(maxPoints);
	    
	    // copy data from not excluded regions
	    REAL *p2=mLSQWeightNotExcluded.data();
	    unsigned long min,max;
	    unsigned long i=0,nbPointsUsedNotExcluded=0;
	    for(int j=0;j<nbExclude;j++) {
	      min=(unsigned long)floor(this->X2Pixel(mExcludedRegionMinX(j)));
	      max=(unsigned long)ceil (this->X2Pixel(mExcludedRegionMaxX(j)));
	      if(min>maxPoints) break;
	      if(max>maxPoints)max=maxPoints;
	      for(;i<min;i++)//! min is the *beginning* of the excluded region !
	      {
					*p2 = *p1; p1++; p2++; nbPointsUsedNotExcluded++;
	      }
	      p1 += max-i;
	      i  += max-i;
	    }
	    for(;i<maxPoints;i++) {
	      *p2 = *p1; p1++; p2++; nbPointsUsedNotExcluded++;
	    }
	    mLSQWeightNotExcluded.resizeAndPreserve(nbPointsUsedNotExcluded);
	    
	    mClockLSQWeightNotExcluded.Click();
	  }
  
  } // else(nbExclude==0)
  
  // contributions of additional LSQ func. objects
	
	std::vector< const CrystVector_REAL* > vLSQWeight(mAdditionalLSQObjRegistry.GetNb());
	long nn = 0;
	for(int iobj=0; iobj<mAdditionalLSQObjRegistry.GetNb(); iobj++) {
		vLSQWeight[iobj] = &mAdditionalLSQObjRegistry.GetObj(iobj).GetLSQWeight(n);
		nn += vLSQWeight[iobj]->numElements();
	}
	
	long nn0 = mLSQWeightNotExcluded.numElements();
	
	mLSQWeightNotExcluded.resizeAndPreserve(nn0+nn);
	
	// copy data
	REAL *p2 = mLSQWeightNotExcluded.data() + nn0;
	for(int iobj=0; iobj<mAdditionalLSQObjRegistry.GetNb(); iobj++) {
		const REAL *p1 = vLSQWeight[iobj]->data();
		for(int i=0; i<vLSQWeight[iobj]->numElements(); i++) *p2++ = *p1++;
	}
	
  return mLSQWeightNotExcluded;
}

const CrystVector_REAL& PowderPattern::GetLSQDeriv(const unsigned int n, RefinablePar &par) {
  
  VFN_DEBUG_MESSAGE("MStruct::PowderPattern::GetLSQDerivative: "<<
		    "calc derivatives for param: "<<par.GetName()<<endl,11)

  bool excluderegions = true;

  ofstream S;
  string filename("LSQDeriv_");
  if (bsavecalc) {
    filename += par.GetName() + string(".txt");
    S.open(filename.c_str());
  }

  bool bnotcalculated = true;
	
  string str = string(par.GetName());
  string::size_type loc;
	
  // "*_Ihkl_*"
  
  loc = str.find("_Ihkl_", 0); // CompName_Ihkl_h_k_l
  if(loc != string::npos) {
    VFN_DEBUG_MESSAGE("MStruct::PowderPattern::GetLSQDerivative: "<<
		      "Derivative for Ihkl detected: "<< str,11)
    mLSQDeriv = GetPowderPatternComponent(str.substr(0,loc)).GetLSQDeriv(n,par);
    // get scale factor (I can also find it in the registry)
    REAL val = GetPar("Scale_"+str.substr(0,loc)).GetValue();
    mLSQDeriv *= val;
    
    bnotcalculated = false;
  }

  // "*pD_*" (SizeDistrib)

  loc = str.find("pD_", 0); // pD_ProfileName_nb
  // If it is SizeDistrib we need to
  //   1) verify this is really derivative with respect to SizeDistrib (histogram) coef.
  //   2) identify the PowderPatternDiffraction component
  //   3) verify that the PowderPatternDiffraction component is of a SizeDistrib type
  if( loc != string::npos && // TODO:: create SizeDistribRefPartype
      par.GetType()->IsDescendantFromOrSameAs(gpRefParTypeScattDataProfileSizeDistrib) ) {
    VFN_DEBUG_MESSAGE("MStruct::PowderPattern::GetLSQDerivative: "<<
		      "Possible Derivative for SizeDistrib detected: "<< str,11)
    
    // Try to find the corresponding ReflectionProfileComponent from the "ProfileName"
    string::size_type loc1 = str.find_last_of("_");
    string profileName = str.substr(loc+3,loc1-(loc+3));
    MStruct::SizeDistribPowderPatternDiffraction *pDiffObj = NULL;
    MStruct::SizeDistribBroadeningEffect *pSizeDistribEffect = NULL;
    {
      unsigned int iComp=0;
      for( ; iComp<this->GetNbPowderPatternComponent(); iComp++) {
	long ind = -1;
	if(this->GetPowderPatternComponent(iComp).GetClassName()==string("MStruct::SizeDistribPowderPatternDiffraction"))
	  ind = dynamic_cast<const MStruct::SizeDistribPowderPatternDiffraction &>(this->GetPowderPatternComponent(iComp)).
	    GetProfile().GetSubObjRegistry().Find(profileName,"MStruct::SizeDistribBroadeningEffect",true);
	if(ind>=0) {
	  pDiffObj =
	    & dynamic_cast<MStruct::SizeDistribPowderPatternDiffraction &>(this->GetPowderPatternComponent(iComp));
	  pSizeDistribEffect =
	    & dynamic_cast<MStruct::SizeDistribBroadeningEffect &>(pDiffObj->GetProfile().GetSubObjRegistry().GetObj(ind));
	  // Confirm finally that the SizeDistribBroadeningEffect really owns the parameter
	  if(pSizeDistribEffect->GetParIndex(par.GetPointer(),true)==-1) {
	    // Finally we have not found what we wanted
	    pDiffObj = NULL;
	    pSizeDistribEffect = NULL;
	  } else
	    break;
	}
      }
    }

    if(pDiffObj!=NULL && pSizeDistribEffect!=NULL) {
      VFN_DEBUG_MESSAGE("MStruct::PowderPattern::GetLSQDerivative: "<<
			"Calculation is forwarded to the Object name: "<< pDiffObj->GetName(),11)
      mLSQDeriv = pDiffObj->GetLSQDeriv(n,par);
      // apply scale factor
      mLSQDeriv *= this->GetScaleFactor( *pDiffObj );

      bnotcalculated = false; 
    }

  }
    
  
  // "Background_Coef_*"
  if(bnotcalculated) {
    loc = str.find("Background_Coef_", 0);
    if (loc != string::npos) {
      VFN_DEBUG_MESSAGE("MStruct::PowderPattern::GetLSQDerivative: "<<
			"Derivative for Chebyshev polynomial coef. detected: "<< str,11)
      mLSQDeriv = GetPowderPatternComponent("bkgData_Chebyshev").GetLSQDeriv(n,par);
      bnotcalculated = false;
    }
  }
  
  // "Scale_"
  if(bnotcalculated) {
    
    loc = str.find("Scale_", 0);
 
    if (loc != string::npos) {
      VFN_DEBUG_MESSAGE("MStruct::PowderPattern::GetLSQDerivative: "<<
			"Derivative for Scale detected: "<< str.data()+loc+6,11)
      mLSQDeriv = GetPowderPatternComponent(str.data()+loc+6).GetPowderPatternCalc();
      bnotcalculated = false;
    }
  }
  	
  // standard calculation
  if(bnotcalculated) {	
    ObjCryst::PowderPattern::GetLSQDeriv(n,par);
    bnotcalculated = false;
    excluderegions = false;
  }

  if (bsavecalc) {
    for(int i=0;i<mLSQDeriv.numElements();i++) S << mLSQDeriv(i) << endl;
    S.close();
    VFN_DEBUG_MESSAGE("MStruct::PowderPattern::GetLSQDerivative: "<<
		      "Derivative saved into the file: "<<filename,11);
  }
  
  // excluded regions
  const long nbExclude=mExcludedRegionMinX.numElements();
  
  const REAL *p1=mLSQDeriv.data();

  if(excluderegions==true && nbExclude>0) {
    unsigned long maxPoints=mNbPointUsed;
    mLSQDerivNotExcluded.resize(maxPoints);
    
    // copy data from not excluded regions
    REAL *p2=mLSQDerivNotExcluded.data();
    unsigned long min,max;
    unsigned long i=0,nbPointsUsedNotExcluded=0;
    for(int j=0;j<nbExclude;j++) {
      min=(unsigned long)floor(this->X2Pixel(mExcludedRegionMinX(j)));
      max=(unsigned long)ceil (this->X2Pixel(mExcludedRegionMaxX(j)));
      if(min>maxPoints) break;
      if(max>maxPoints)max=maxPoints;
      for(;i<min;i++)//! min is the *beginning* of the excluded region !
      {
				*p2 = *p1; p1++; p2++; nbPointsUsedNotExcluded++;
      }
      p1 += max-i;
      i  += max-i;
    }
    for(;i<maxPoints;i++) {
      *p2 = *p1; p1++; p2++; nbPointsUsedNotExcluded++;
    }
    mLSQDerivNotExcluded.resizeAndPreserve(nbPointsUsedNotExcluded);
  } else {
  	mLSQDerivNotExcluded = mLSQDeriv;
  }
  
  // contributions of additional LSQ func. objects
	
  std::vector< const CrystVector_REAL* > vLSQDeriv(mAdditionalLSQObjRegistry.GetNb());
  long nn = 0;
  for(int iobj=0; iobj<mAdditionalLSQObjRegistry.GetNb(); iobj++) {
    vLSQDeriv[iobj] = &mAdditionalLSQObjRegistry.GetObj(iobj).GetLSQDeriv(n,par);
    nn += vLSQDeriv[iobj]->numElements();
  }
	
  long nn0 = mLSQDerivNotExcluded.numElements();
  
  mLSQDerivNotExcluded.resizeAndPreserve(nn0+nn);
  
  // copy data
  REAL *p2 = mLSQDerivNotExcluded.data() + nn0;
  for(int iobj=0; iobj<mAdditionalLSQObjRegistry.GetNb(); iobj++) {
    const REAL *p1 = vLSQDeriv[iobj]->data();
    for(int i=0; i<vLSQDeriv[iobj]->numElements(); i++) *p2++ = *p1++;
  }
  
  return mLSQDerivNotExcluded;  
}

void PowderPattern::SetLinearPolarRate(const REAL f)
{
  mRadiation.SetLinearPolarRate(f);
}

void PowderPattern::SetIncidenceAngle(const REAL omega)
{
	mOmega = omega;
	mClockMaster.Click();
}

REAL PowderPattern::GetIncidenceAngle(const REAL omega) const
{
	return mOmega;
}

void PowderPattern::BeginOptimization(const bool allowApproximations,
                                      const bool enableRestraints)
{
  // Scales of scalable components could be modified during call
  // of ObjCryst::PowderPattern BeginOptimization method. This is in our case
  // sometimes unsuitable. We would like to preserve scales of components here
  // if the derivative step is zero or positive.
  CrystVector_REAL scaleFactors = mScaleFactor;  
  ObjCryst::PowderPattern::BeginOptimization(allowApproximations, enableRestraints);
  for(int icomp=0; icomp<mScalableComponentIndex.numElements(); icomp++) {
    RefinablePar & param = GetPar(&mScaleFactor(mScalableComponentIndex(icomp)));
    if(param.IsUsed() && abs(param.GetDerivStep())>-1.e-7) 
      param.SetValue(scaleFactors(mScalableComponentIndex(icomp)));
  }
}

void PowderPattern::AddAdditionalLSQObj(ObjCryst::RefinableObj& obj)
{
	mAdditionalLSQObjRegistry.Register(obj);
}

void PowderPattern::RemoveAdditionalLSQObj(ObjCryst::RefinableObj& obj)
{
	mAdditionalLSQObjRegistry.DeRegister(obj);
}

// Texture Correction

/* vector cross product (c = a x b) */
void cross(REAL* a, REAL* b, REAL* c)
{
  c[0] = a[1]*b[2]-a[2]*b[1];
  c[1] = a[2]*b[0]-a[0]*b[2];
  c[2] = a[0]*b[1]-a[1]*b[0];
}

/* normalized Gauss function */
REAL fgauss(REAL* a, REAL x)
{
  return 1./a[1]/sqrt(2.*M_PI)*exp(-(x-a[0])*(x-a[0])/2./a[1]/a[1]);
}

/* calculate rotation matrix for given Euler angles,
   the so-called x-convention is used (see: Euler angles,
   www.mathworld.com */
void rotationMatrix(REAL* A, REAL phi, REAL th, REAL psi)
{
  A[0] = cos(phi)*cos(psi) - cos(th)*sin(phi)*sin(psi);
  A[1] = cos(psi)*sin(phi) + cos(phi)*cos(th)*sin(psi);
  A[2] = sin(psi)*sin(th);
  A[3] = -(cos(psi)*cos(th)*sin(phi)) - cos(phi)*sin(psi);
  A[4] = cos(phi)*cos(psi)*cos(th) - sin(phi)*sin(psi);
  A[5] = cos(psi)*sin(th);
  A[6] = sin(phi)*sin(th);
  A[7] = -(cos(phi)*sin(th));
  A[8] = cos(th);
}

/* calculate Euler angle from a given rotation matrix,
   the so-called x-convention is used (see: Euler angles,
   www.mathworld.com */
void getAngles(REAL* A, REAL& phi, REAL& th, REAL& psi)
{
  th = acos(A[8]);
  if (fabs(th)>FLT_EPSILON && fabs(th-M_PI)>FLT_EPSILON) {
    REAL sinth = sin(th);
    phi = atan2(A[6]/sinth,-A[7]/sinth);
    psi = atan2(A[2]/sinth,A[5]/sinth); }
  else {
    phi = atan2(A[1],A[0]);
    psi = 0.0; }
}

/* matrix multiplication A = B*C */
void matrixMult(REAL* B, REAL* C, REAL* A)
{
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++) {
      REAL t = 0.;
      for(int k=0; k<3; k++) t += B[i*3+k]*C[k*3+j];
      A[3*i+j] = t;
    }
} 

// TextureCalculator
TextureCalculator::TextureCalculator()
  :mpCrystal(0),mMultiplicity(0)
{
  mParams.resize(14);
  mvri.resize(0);
  mAAMatrix.resize(9);
  mAMatrix.resize(9);
  mTMatrix.resize(9);
  //mMAMatrix.resize(9);
  mAVector.resize(0);

  // set dafault ODF function params
  mParams = 0.; mParams(8) = 1.; mParams(9) = 1.;
  //SetTextureParams(params);

  /* define integration limits and an integration grid */
  phimin = 0.*DEG2RAD; phimax = 360.*DEG2RAD; phistep = 10.*DEG2RAD;
  nphi = int((phimax-phimin)/phistep)+1;

  thmin = 0.*DEG2RAD; thmax = 180.*DEG2RAD; thstep = 2.5*DEG2RAD;
  nth = int((thmax-thmin)/thstep)+1;

  psimin = 0.*DEG2RAD; psimax = 360.*DEG2RAD; psistep = 10.*DEG2RAD;
  npsi = int((psimax-psimin)/psistep)+1;

  // calc norm. factor of the ODF function
  //mOdfNFactor = CalcOdfNFactor();
  
  // set dafault ODF function params
  //CrystVector_REAL params(14);
  //params = 0.; params(8) = 1.; params(9) = 1.;
  //SetTextureParams(params);
}

void TextureCalculator::SetTextureParams(const CrystVector_REAL& params,
                                         bool bForceTextureSymmetry) {
  // if the texture is just rotated - only two last params are modified,
  // it is not necessary to recalculate the normalization factor
  
  bool recalcOdfNFactor = (fabs(params(0)-mParams(0))>1.e-7)
                       || (fabs(params(1)-mParams(1))>1.e-7)
                       || (fabs(params(2)-mParams(2))>1.e-7)
                       || (fabs(params(3)-mParams(3))>1.e-7)
                       || (fabs(params(4)-mParams(4))>1.e-7)
                       || (fabs(params(5)-mParams(5))>1.e-7)
                       || (fabs(params(6)-mParams(6))>1.e-7)
                       || (fabs(params(7)-mParams(7))>1.e-7)
                       || (fabs(params(8)-mParams(8))>1.e-7)
                       || (fabs(params(9)-mParams(9))>1.e-7)
                       || (fabs(params(10)-mParams(10))>1.e-7)
                       || (fabs(params(11)-mParams(11))>1.e-7);

  // set parameters
  mParams = params;
  
  // set flac for texture symmetry conditions
  mbForceTextureSymmetry = bForceTextureSymmetry;

  // Force Texture Symmetry !
  //    Generate list of all orientations of crystallities which are equivalent
  //    with supplied main and secondary HKL texture system
  if(mbForceTextureSymmetry==true)
  {
    // Generate list of all equivalent reflections to the main {HKL} texture axis
    mnmtaHKL =
      mpCrystal->GetSpaceGroup().GetAllEquivRefl(params(6),params(7),params(8),false,true);
  
    int tamultiplicity = mnmtaHKL.rows();

    REAL ac = mpCrystal->GetLatticePar(0)/
              mpCrystal->GetLatticePar(2); // a/c

    // Normalise the list of the main texture axis
    const REAL sm = sqrt(params(6)*params(6)+params(7)*params(7)+ac*ac*params(8)*params(8));
    mnmtaHKL /= sm;
    
    // Create list of all equivalent reflections to the secondary {HKL} texture axis
    mnstaHKL = CrystMatrix_REAL(mnmtaHKL.rows(),mnmtaHKL.cols());

    // Calculate the normalisation factor for the list of the secondary texture axis
    const REAL ss = sqrt(params(9)*params(9)+params(10)*params(10)+ac*ac*params(11)*params(11));
 
    mnstaHKL(0,0) = params(9)/ss; mnstaHKL(0,1) = params(10)/ss; mnstaHKL(0,2) = params(11)/ss;
    
    // Rotation matrix of the original main texture axis
    CrystMatrix_REAL AA1(3,3);
    REAL st, cosphi, sinphi, cospsi, sinpsi;

    {
      st = sqrt(mnmtaHKL(0,0)*mnmtaHKL(0,0)+mnmtaHKL(0,1)*mnmtaHKL(0,1));
      cosphi =  (st>1.e-4) ? mnmtaHKL(0,1)/st : 1.;
      sinphi =  (st>1.e-4) ? mnmtaHKL(0,0)/st : 0.;
      cospsi =  ac*mnmtaHKL(0,2);
      sinpsi =  st;

      AA1(0,0) = cosphi; AA1(0,1) = -sinphi; AA1(0,2) = 0.;
      AA1(1,0) = cospsi*sinphi; AA1(1,1) = cosphi*cospsi; AA1(1,2) = -sinpsi;
      AA1(2,0) = sinphi*sinpsi; AA1(2,1) = cosphi*sinpsi; AA1(2,2) = cospsi; 
    }

    for(int ihkl=1; ihkl<tamultiplicity; ihkl++) {
      
      // Create inversion of an appropriate rotation matrix of qn equivalent reflection
      CrystMatrix_REAL AA2i(3,3);

      {
	st = sqrt(mnmtaHKL(ihkl,0)*mnmtaHKL(ihkl,0)+mnmtaHKL(ihkl,1)*mnmtaHKL(ihkl,1));
	cosphi =  (st>1.e-4) ? mnmtaHKL(ihkl,1)/st : 1.;
	sinphi =  (st>1.e-4) ? mnmtaHKL(ihkl,0)/st : 0.;
	cospsi =  ac*mnmtaHKL(ihkl,2);
	sinpsi =  st;
	
	AA2i(0,0) = cosphi; AA2i(0,1) = cospsi*sinphi; AA2i(0,2) = sinphi*sinpsi;
	AA2i(1,0) = -sinphi; AA2i(1,1) = cosphi*cospsi; AA2i(1,2) = cosphi*sinpsi;
	AA2i(2,0) = 0.; AA2i(2,1) = -sinpsi; AA2i(2,2) = cospsi;
      }

      // Create a rotation matrix from the original main texture direction
      // to the current equivalent main texture direction: AA = AA2i*AA1
      CrystMatrix_REAL AA(3,3);
      AA = 0.;
      for(int i=0; i<AA.rows(); i++)
	for(int j=0; j<AA.cols(); j++)
	  for(int k=0; k<AA2i.cols(); k++)
	    AA(i,j) += AA2i(i,k)*AA1(k,j);
      
      // Create an appropriate secondary texture axis direction - transform the original main
      // secondary axis by the matrix AA
      cout << fixed;
      mnstaHKL(ihkl,0) = AA(0,0)*mnstaHKL(0,0)+AA(0,1)*mnstaHKL(0,1)+AA(0,2)*mnstaHKL(0,2);
      mnstaHKL(ihkl,1) = AA(1,0)*mnstaHKL(0,0)+AA(1,1)*mnstaHKL(0,1)+AA(1,2)*mnstaHKL(0,2);
      mnstaHKL(ihkl,2) = AA(2,0)*mnstaHKL(0,0)+AA(2,1)*mnstaHKL(0,1)+AA(2,2)*mnstaHKL(0,2);

    } // ihkl
    
    // The list of secondary texture axis is already normalised - we are dome

    // print
    
    #ifdef __DEBUG__ZDENEK__
    cout<<"List of equivalent texture systems for: (";
    cout<<params(6)<<" "<<params(7)<<" "<<params(8)<<")x(";
    cout<<params(9)<<" "<<params(10)<<" "<<params(11)<<")\n";
    
    for(int i=0; i< mnmtaHKL.rows(); i++) {
      cout<<"("<<mnmtaHKL(i,0)<<" "<<mnmtaHKL(i,1)<<" "<<mnmtaHKL(i,2)<<")x(";
      cout<<mnstaHKL(i,0)<<" "<<mnstaHKL(i,1)<<" "<<mnstaHKL(i,2)<<")\n";
    }
    cout<<endl;
    #endif // __DEBUG__ZDENEK__
  }
  else { // not mbForthTextureSymmetry==true

    // Lists of the main and secondary texture axis contains only one record
    mnmtaHKL.resize(1,3);
    mnstaHKL.resize(1,3);

    REAL ac = mpCrystal->GetLatticePar(0)/
              mpCrystal->GetLatticePar(2); // a/c

    // Normalization facor
    const REAL sm = sqrt(params(6)*params(6)+params(7)*params(7)+ac*ac*params(8)*params(8));
 
    mnmtaHKL(0,0) = params(6)/sm; mnmtaHKL(0,1) = params(7)/sm; mnmtaHKL(0,2) = params(8)/sm;

    // Normalization facor
    const REAL ss = sqrt(params(9)*params(9)+params(10)*params(10)+ac*ac*params(11)*params(11));
 
    mnstaHKL(0,0) = params(9)/ss; mnstaHKL(0,1) = params(10)/ss; mnstaHKL(0,2) = params(11)/ss;
    
  }

  // calc norm. factor of the ODF function
  if(recalcOdfNFactor==true) mOdfNFactor = CalcOdfNFactor();
}

void TextureCalculator::SetCrystal(const Crystal& crystal) {
  mpCrystal = &crystal;
}

/* ODF - orientation distribution function */
REAL TextureCalculator::funcodf(REAL phi,REAL th,REAL psi,REAL* aa) const
{
  // parameters
  if (mParams.numElements()<14) throw ObjCrystException("MStruct::TextureCalculator::funcodf(...)\
             Not enough texture parameters set.");
  const CrystVector_REAL &a = mParams;

  #ifdef __LOISA_123__
  cout<<"a:\n";
  cout<<a(0)<<" "<<a(1)*RAD2DEG<<" "<<a(2)*RAD2DEG<<"\n";
  cout<<a(3)<<" "<<a(4)*RAD2DEG<<" "<<a(5)*RAD2DEG<<"\n";
  cout<<a(6)<<" "<<a(7)<<" "<<a(8)<<"\n";
  cout<<a(9)<<" "<<a(10)<<" "<<a(11)<<"\n";
  cout<<a(12)*RAD2DEG<<" "<<a(13)*RAD2DEG<<"\n";
  #endif
	
  // if rotation matrix(phi,th,psi) not suplied calc it
  if (aa==0) {
    aa = mAAMatrix.data();
    rotationMatrix(aa,phi,th,psi);
  }
  else { // if supplied we need also angles (they are right?) so calc them
    getAngles(aa,phi,th,psi);
  }
  
  // tilt
  REAL *ta = mTMatrix.data();
  rotationMatrix(ta,0.,a(12),a(13));
  
  REAL tresult = 0.;

  const int nsymcomponents = (mbForceTextureSymmetry==true) ? mnmtaHKL.rows() : 1;

  for(int ita=0; ita<nsymcomponents; ita++) {
  
    // normalized HKLx and HKLz vectors
    CrystVector_REAL HKLz(3), HKLx(3);
  
    //HKLz(0) = a(6); HKLz(1) = a(7); HKLz(2) = a(8);
    //HKLz /= sqrt(a(6)*a(6)+ a(7)*a(7)+ a(8)*a(8));
    HKLz(0) = mnmtaHKL(ita,0); HKLz(1) = mnmtaHKL(ita,1); HKLz(2) = mnmtaHKL(ita,2);
    
    //HKLx(0) = a(9); HKLx(1) = a(10); HKLx(2) = a(11);
    //HKLx /= sqrt(a(9)*a(9)+ a(10)*a(10)+ a(11)*a(11));
    HKLx(0) = mnstaHKL(ita,0); HKLx(1) = mnstaHKL(ita,1); HKLx(2) = mnstaHKL(ita,2);

    REAL t = (HKLz(0)*aa[0]+HKLz(1)*aa[1]+HKLz(2)*aa[2])*ta[2] +
      (HKLz(0)*aa[3]+HKLz(1)*aa[4]+HKLz(2)*aa[5])*ta[5] +
      (HKLz(0)*aa[6]+HKLz(1)*aa[7]+HKLz(2)*aa[8])*ta[8];
    //REAL t = HKLz(0)*aa[6]+HKLz(1)*aa[7]+HKLz(2)*aa[8];
    //REAL t = 1/sqrt(3.)*(aa[6]+aa[7]+aa[8]);
    REAL th0  = (fabs(t)<1.) ? acos(t) : ((t>0.) ? 0. : M_PI); 
    //t = HKLx(0)*aa[0]+HKLx(1)*aa[1]+HKLx(2)*aa[2];
    t = (HKLx(0)*aa[0]+HKLx(1)*aa[1]+HKLx(2)*aa[2])*ta[0] +
      (HKLx(0)*aa[3]+HKLx(1)*aa[4]+HKLx(2)*aa[5])*ta[3] +
      (HKLx(0)*aa[6]+HKLx(1)*aa[7]+HKLx(2)*aa[8])*ta[6];
    REAL psi0 = (fabs(t)<1.) ? acos(t) : ((t>0.) ? 0. : M_PI);
    // we need to know also the sign of the psi0 angle
    REAL tt1[3], tt2[3], tt3[3];
    tt1[0] = HKLx(0)*aa[0]+HKLx(1)*aa[1]+HKLx(2)*aa[2]; tt2[0] = ta[0];
    tt1[1] = HKLx(0)*aa[3]+HKLx(1)*aa[4]+HKLx(2)*aa[5]; tt2[1] = ta[3];
    tt1[2] = HKLx(0)*aa[6]+HKLx(1)*aa[7]+HKLx(2)*aa[8]; tt2[2] = ta[6];
    cross(tt1,tt2,tt3);
    t = tt3[0]*ta[2]+tt3[1]*ta[5]+tt3[2]*ta[8];
    psi0 = (t>=0.) ? psi0 : 2.*M_PI-psi0;
    
    REAL result = 1.;
    CrystVector_REAL &a1 = mAVector;

	  if (fabs(a(0))>FLT_EPSILON) {
      REAL dx = (1.-cos(th0))/(1.-cos(a(1)/2.));
      //result *= exp(-dx/0.1);
      result *= exp(-M_LN2*pow(dx,REAL(a(2)/M_PI*180)));
	  }
	  if(0) {
    // the 1st component
    if (fabs(a(0))>FLT_EPSILON) {
      REAL dx = fabs(th0-a(1));
      dx = (fabs(th0+M_PI-a(1))>dx) ? dx : fabs(th0+M_PI-a(1));
      dx = (fabs(th0-M_PI-a(1))>dx) ? dx : fabs(th0-M_PI-a(1));
      a1 = a; a1(1) = 0.;
      result *= 1.-a(0) + a(0)*fgauss(a1.data()+1,dx);
    }
    if (fabs(a(3))>FLT_EPSILON) {
      REAL dx = fabs(psi0-a(4));
      dx = (fabs(psi0+2.*M_PI-a(4))>dx) ? dx : fabs(psi0+2.*M_PI-a(4));
      dx = (fabs(psi0-2.*M_PI-a(4))>dx) ? dx : fabs(psi0-2.*M_PI-a(4));
      a1 = a; a1(4) = 0.;
      result *= 1.-a(3) + a(3)*fgauss(a1.data()+4,dx);
    }
    }
    
    tresult += result;
  } // ihkl

  return tresult/nsymcomponents;
}

/* ODF - orientation distribution function */
REAL TextureCalculator::funcodf0(REAL phi,REAL th,REAL psi,REAL* aa) const
{
  // parameters
  if (mParams.numElements()<14) throw ObjCrystException("MStruct::TextureCalculator::funcodf(...)\
             Not enough texture parameters set.");
  const CrystVector_REAL &a = mParams;

  #ifdef __LOISA_123__
  cout<<"a:\n";
  cout<<a(0)<<" "<<a(1)*RAD2DEG<<" "<<a(2)*RAD2DEG<<"\n";
  cout<<a(3)<<" "<<a(4)*RAD2DEG<<" "<<a(5)*RAD2DEG<<"\n";
  cout<<a(6)<<" "<<a(7)<<" "<<a(8)<<"\n";
  cout<<a(9)<<" "<<a(10)<<" "<<a(11)<<"\n";
  cout<<a(12)*RAD2DEG<<" "<<a(13)*RAD2DEG<<"\n";
  #endif

  // if rotation matrix(phi,th,psi) not suplied calc it
  if (aa==0) {
    aa = mAAMatrix.data();
    rotationMatrix(aa,phi,th,psi);
  }
  else { // if supplied we need also angles (they are right?) so calc them
    getAngles(aa,phi,th,psi);
  }
  
  // tilt
  REAL *ta = mTMatrix.data();
  rotationMatrix(ta,0.,a(12),a(13));
  
  REAL tresult = 0.;

  const int nsymcomponents = (mbForceTextureSymmetry==true) ? mnmtaHKL.rows() : 1;

  for(int ita=0; ita<nsymcomponents; ita++) {
  
    // normalized HKLx and HKLz vectors
    CrystVector_REAL HKLz(3), HKLx(3);
  
    //HKLz(0) = a(6); HKLz(1) = a(7); HKLz(2) = a(8);
    //HKLz /= sqrt(a(6)*a(6)+ a(7)*a(7)+ a(8)*a(8));
    HKLz(0) = mnmtaHKL(ita,0); HKLz(1) = mnmtaHKL(ita,1); HKLz(2) = mnmtaHKL(ita,2);
    
    //HKLx(0) = a(9); HKLx(1) = a(10); HKLx(2) = a(11);
    //HKLx /= sqrt(a(9)*a(9)+ a(10)*a(10)+ a(11)*a(11));
    HKLx(0) = mnstaHKL(ita,0); HKLx(1) = mnstaHKL(ita,1); HKLx(2) = mnstaHKL(ita,2);

    REAL t = (HKLz(0)*aa[0]+HKLz(1)*aa[1]+HKLz(2)*aa[2])*ta[2] +
      (HKLz(0)*aa[3]+HKLz(1)*aa[4]+HKLz(2)*aa[5])*ta[5] +
      (HKLz(0)*aa[6]+HKLz(1)*aa[7]+HKLz(2)*aa[8])*ta[8];
    //REAL t = HKLz(0)*aa[6]+HKLz(1)*aa[7]+HKLz(2)*aa[8];
    //REAL t = 1/sqrt(3.)*(aa[6]+aa[7]+aa[8]);
    REAL th0  = (fabs(t)<1.) ? acos(t) : ((t>0.) ? 0. : M_PI); 
    //t = HKLx(0)*aa[0]+HKLx(1)*aa[1]+HKLx(2)*aa[2];
    t = (HKLx(0)*aa[0]+HKLx(1)*aa[1]+HKLx(2)*aa[2])*ta[0] +
      (HKLx(0)*aa[3]+HKLx(1)*aa[4]+HKLx(2)*aa[5])*ta[3] +
      (HKLx(0)*aa[6]+HKLx(1)*aa[7]+HKLx(2)*aa[8])*ta[6];
    REAL psi0 = (fabs(t)<1.) ? acos(t) : ((t>0.) ? 0. : M_PI);
    // we need to know also the sign of the psi0 angle
    REAL tt1[3], tt2[3], tt3[3];
    tt1[0] = HKLx(0)*aa[0]+HKLx(1)*aa[1]+HKLx(2)*aa[2]; tt2[0] = ta[0];
    tt1[1] = HKLx(0)*aa[3]+HKLx(1)*aa[4]+HKLx(2)*aa[5]; tt2[1] = ta[3];
    tt1[2] = HKLx(0)*aa[6]+HKLx(1)*aa[7]+HKLx(2)*aa[8]; tt2[2] = ta[6];
    cross(tt1,tt2,tt3);
    t = tt3[0]*ta[2]+tt3[1]*ta[5]+tt3[2]*ta[8];
    psi0 = (t>=0.) ? psi0 : 2.*M_PI-psi0;
    
    REAL result = 1.;
    CrystVector_REAL &a1 = mAVector;

    // the 1st component
    if (fabs(a(0))>FLT_EPSILON) {
      REAL dx = fabs(th0-a(1));
      dx = (fabs(th0+M_PI-a(1))>dx) ? dx : fabs(th0+M_PI-a(1));
      dx = (fabs(th0-M_PI-a(1))>dx) ? dx : fabs(th0-M_PI-a(1));
      a1 = a; a1(1) = 0.;
      result *= 1.-a(0) + a(0)*fgauss(a1.data()+1,dx);
    }
    if (fabs(a(3))>FLT_EPSILON) {
      REAL dx = fabs(psi0-a(4));
      dx = (fabs(psi0+2.*M_PI-a(4))>dx) ? dx : fabs(psi0+2.*M_PI-a(4));
      dx = (fabs(psi0-2.*M_PI-a(4))>dx) ? dx : fabs(psi0-2.*M_PI-a(4));
      a1 = a; a1(4) = 0.;
      result *= 1.-a(3) + a(3)*fgauss(a1.data()+4,dx);
    }
    
    tresult += result;
  } // ihkl

  return tresult/nsymcomponents;
}

/* ODF - orientation distribution function */
REAL TextureCalculator::funcodf1(REAL phi,REAL th,REAL psi,REAL* aa) const
{
  // parameters
  if (mParams.numElements()<8) throw ObjCrystException("MStruct::TextureCalculator::funcodf(...)\
             Not enough texture parameters set.");
  const CrystVector_REAL &a = mParams;

  // if rotation matrix(phi,th,psi) not suplied calc it
  if (aa==0) {
    aa = mAAMatrix.data();
    rotationMatrix(aa,phi,th,psi);
  }
  else { // if supplied we need also angles (they are right?) so calc them
    getAngles(aa,phi,th,psi);
  }
  
  // tilt
  REAL *ta = mTMatrix.data();
  rotationMatrix(ta,0.,a(6),a(7));

  REAL t = ta[2]*aa[2]+ta[5]*aa[5]+ta[8]*aa[8];
  REAL th0  = (fabs(t)<1.) ? acos(t) : ((t>0.) ? 0. : M_PI); 
  t = ta[0]*aa[0]+ta[3]*aa[3]+ta[6]*aa[6];
  REAL psi0 = (fabs(t)<1.) ? acos(t) : ((t>0.) ? 0. : M_PI); 

  REAL result = 1.;
  CrystVector_REAL &a1 = mAVector;

  // the 1st component
  if (fabs(a(0))>FLT_EPSILON) {
    REAL dx = fabs(th0-a(1));
    //dx = (fabs(th0+M_PI-a(1))>dx) ? dx : fabs(th0+M_PI-a(1));
    //dx = (fabs(th0-M_PI-a(1))>dx) ? dx : fabs(th0-M_PI-a(1));
    a1 = a; a1(1) = 0.;
    result *= 1.-a(0) + a(0)*fgauss(a1.data()+1,dx);
  }
  if (fabs(a(3))>FLT_EPSILON) {
    REAL dx = fabs(psi0-a(4));
    a1 = a; a1(4) = 0.;
    result *= 1.-a(3) + a(3)*fgauss(a1.data()+4,dx);
  }
  
  return result;
}

/* ODF - orientation distribution function */
REAL TextureCalculator::funcodf2(REAL phi,REAL th,REAL psi,REAL* aa) const
{
  // parameters
  if (mParams.numElements()<8) throw ObjCrystException("MStruct::TextureCalculator::funcodf(...)\
             Not enough texture parameters set.");
  const CrystVector_REAL &a = mParams;

  // if rotation matrix(phi,th,psi) not suplied calc it
  if (aa==0) {
    aa = mAAMatrix.data();
    rotationMatrix(aa,phi,th,psi);
  }
  else { // if supplied we need also angles (they are right?) so calc them
    getAngles(aa,phi,th,psi);
  }
  
  // tilt
  REAL *ta = mTMatrix.data();
  rotationMatrix(ta,a(13),a(12),0.);

  // Generate list of all equivalent reflections
  CrystMatrix_REAL HKLz =
    mpCrystal->GetSpaceGroup().GetAllEquivRefl(a(6),a(7),a(8),true);
  
  const int multiplicity = HKLz.rows();

  // their equivalent HKLx
  CrystMatrix_REAL HKLx(HKLz.rows(),HKLz.cols());

  // normalized texture base vectors
  CrystVector_REAL qz(3),qx(3);
  qz(0) = a(6); qz(1) = a(7); qz(2) = a(8);
  qz /= sqrt(qz(0)*qz(0)+qz(1)*qz(1)+qz(2)*qz(2));
  qx(0) = a(9); qx(1) = a(10); qx(2) = a(11);
  qx /= sqrt(qx(0)*qx(0)+qx(1)*qx(1)+qx(2)*qx(2));

  for(int ihkl=0; ihkl<multiplicity; ihkl++) {
    // normalized HKLz - hz
    CrystVector_REAL hz(3);
    hz(0) = HKLz(ihkl,0); hz(1) = HKLz(ihkl,1); hz(2) = HKLz(ihkl,2);
    hz /= sqrt(hz(0)*hz(0)+hz(1)*hz(1)+hz(2)*hz(2));
    
    const REAL nqz = sqrt(qz(0)*qz(0)+qz(1)*qz(1));
    const REAL nhz = sqrt(hz(0)*hz(0)+hz(1)*hz(1));

    if ((nqz>1.e-4) && (nhz>1.e-4)) {
      REAL qphi = acos((qz(0)*hz(0)+qz(1)*hz(1))/nqz/nhz) * 
	                   ((qz(0)*hz(1)-qz(1)*hz(0)>0) ? 1.: -1.);
      
    }
  }

  // tilted texture base vectors
  CrystVector_REAL bz(3);
  CrystVector_REAL bx(3);

  bz(0) = ta[0]*a(6) + ta[1]*a(7) + ta[2]*a(8);
  bz(1) = ta[3]*a(6) + ta[4]*a(7) + ta[5]*a(8);
  bz(2) = ta[6]*a(6) + ta[7]*a(7) + ta[8]*a(8);

  bx(0) = ta[0]*a(9) + ta[1]*a(10) + ta[2]*a(11);
  bx(1) = ta[3]*a(9) + ta[4]*a(10) + ta[5]*a(11);
  bx(2) = ta[6]*a(9) + ta[7]*a(10) + ta[8]*a(11);

  // rotated(aa) texture base vectors
  CrystVector_REAL cz(3);
  CrystVector_REAL cx(3);
  
  cz(0) = aa[0]*a(6) + aa[1]*a(7) + aa[2]*a(8);
  cz(1) = aa[3]*a(6) + aa[4]*a(7) + aa[5]*a(8);
  cz(2) = aa[6]*a(6) + aa[7]*a(7) + aa[8]*a(8);

  cx(0) = aa[0]*a(9) + aa[1]*a(10) + aa[2]*a(11);
  cx(1) = aa[3]*a(9) + aa[4]*a(10) + aa[5]*a(11);
  cx(2) = aa[6]*a(9) + aa[7]*a(10) + aa[8]*a(11);

  //REAL t = a(6)*aa[6]+a(7)*aa[7]+a(8)*aa[8];
  //t /= sqrt(a(6)*a(6)+ a(7)*a(7)+ a(8)*a(8));
  REAL t = bz(0)*cz(0)+bz(1)*cz(1)+bz(2)*cz(2);
  t /= a(6)*a(6)+ a(7)*a(7)+ a(8)*a(8);
  REAL th0  = (fabs(t)<1.) ? acos(t) : ((t>0.) ? 0. : M_PI); 
  //t =  a(9)*aa[0]+a(10)*aa[1]+a(11)*aa[2];
  //t /= sqrt(a(9)*a(9)+ a(10)*a(10)+ a(11)*a(11));
  t =  bx(0)*cx(0)+bx(1)*cx(1)+bx(2)*cx(2);
  t /= a(9)*a(9)+ a(10)*a(10)+ a(11)*a(11);
  REAL psi0 = (fabs(t)<1.) ? acos(t) : ((t>0.) ? 0. : M_PI); 

  REAL result = 1.;
  CrystVector_REAL &a1 = mAVector;

  // the 1st component
  if (fabs(a(0))>FLT_EPSILON) {
    REAL dx = fabs(th0-a(1));
    //dx = (fabs(th0+M_PI-a(1))>dx) ? dx : fabs(th0+M_PI-a(1));
    //dx = (fabs(th0-M_PI-a(1))>dx) ? dx : fabs(th0-M_PI-a(1));
    a1 = a; a1(1) = 0.;
    result *= 1.-a(0) + a(0)*fgauss(a1.data()+1,dx);
  }
  if (fabs(a(3))>FLT_EPSILON) {
    REAL dx = fabs(psi0-a(4));
    a1 = a; a1(4) = 0.;
    result *= 1.-a(3) + a(3)*fgauss(a1.data()+4,dx);
  }
  
  return result;
}

REAL TextureCalculator::CalcOdfNFactor() const
{
  /* calc an ODF normalization factor */
  REAL odfnfactor = 0.;

  /* integration of the ODF function over whole Euler angles space */
  {
    // integration over th, psi and phi (psi and phi are 2PI periodic)
    REAL sum3 = 0.;
    for(int ith=0; ith<nth-1; ith++) {
      REAL th = thmin + ith*thstep;
      REAL domega = sin(th)/8./M_PI/M_PI;
      REAL sum2 = 0.;
      for(int ipsi=0; ipsi<npsi-1; ipsi++) {
	REAL psi = psimin + ipsi*psistep;
	REAL sum1 = 0.;
	for(int iphi=0; iphi<nphi-1; iphi++) {
	  REAL phi = phimin + iphi*phistep;
	  sum1 += funcodf(phi,th,psi)*domega;	}
	sum2 += sum1*phistep; }
      sum3 += sum2*psistep; }
    odfnfactor = sum3*thstep;
    cout << "ODF normalization factor: " << odfnfactor << endl;	
  }

  return odfnfactor;
}

void TextureCalculator::PrepareForCalc(REAL h,REAL k,REAL l) const {

  if (mpCrystal==0) throw ObjCrystException("MStruct::TextureCalculator::PrepareForCalc(...)\
             Can not find Crystal object. This object was not properly initialized.");

  int &multiplicity = mMultiplicity;
  CrystVector_REAL &vri = mvri;

  REAL hh[3] = {h, k, l};

  // Generate list of all equivalent reflections
  CrystMatrix_REAL HKL =
    mpCrystal->GetSpaceGroup().GetAllEquivRefl(h,k,l,false,true);
  
  multiplicity = HKL.rows();

  #ifdef __DEBUG__ZDENEK__
 // print the HKL list
  cout << "GetAllEquivRefl(h,k,l,true,false)" << endl;
  cout << "Symmetry ecvivalent reflections which will be used:" << endl;
  for(int i=0;i<multiplicity;i++) {
    cout<<setw(12)<<HKL(i,0)<<setw(12)<<HKL(i,1)<<setw(12)<<HKL(i,2)<<endl;
  }
  #endif /*  __DEBUG__ZDENEK__ */

  REAL ac = mpCrystal->GetLatticePar(0)/
            mpCrystal->GetLatticePar(2); // a/c

  /* create matrix transforming given diffraction vector
     into the direction (0,0,1) */
  vri.resize(9*multiplicity);
  for(int irefl=0;irefl<multiplicity;irefl++) /* for each reflection */
  {
    /* generate an orthogonal system - diffraction vector
       (h,k,l) should be parallel with z-axes */

    // normalization of the diffraction vector
    REAL h0 = sqrt(h*h+k*k+l*l*ac*ac);
    hh[0]=HKL(irefl,0)/h0; hh[1]=HKL(irefl,1)/h0; hh[2]=ac*HKL(irefl,2)/h0;

    REAL a1[3] = {1., 0., 0.};
    REAL a2[3] = {0., 0., 0.};

    if (fabs(1.-fabs(a1[0]*hh[0]+a1[1]*hh[1]+a1[2]*hh[2]))<FLT_EPSILON) {
      a1[0] = 0.; a1[1] = 1.; a1[2] = 0.;
    }

    if (fabs(1.-fabs(a1[0]*hh[0]+a1[1]*hh[1]+a1[2]*hh[2]))<FLT_EPSILON) {
      ostringstream s;
      s << "MStruct::TextureCalculator::PrepareForCalc(...): ";
      s << "Given (hkl)=("<<h<<","<<k<<","<<l<<") can not be used."<<endl;
      throw ObjCrystException(s.str());
    }

    cross(hh,a1,a2);
    {
      REAL a0 = sqrt(a2[0]*a2[0]+a2[1]*a2[1]+a2[2]*a2[2]);
      a2[0] /= a0;  a2[1] /= a0;  a2[2] /= a0;
    } 
    cross(a2,hh,a1);
    {
      REAL a0 = sqrt(a1[0]*a1[0]+a1[1]*a1[1]+a1[2]*a1[2]);
      a1[0] /= a0;  a1[1] /= a0;  a1[2] /= a0;
    }

    /* create a transformation matrix */
    
    REAL r[9] = {a1[0], a2[0], hh[0],
        	 a1[1], a2[1], hh[1],
                 a1[2], a2[2], hh[2]};

    /* estimate angles */
    REAL phi, th, psi;
    getAngles(r,phi,th,psi);
    
    /* create an inversion of the transformation matrix */
    REAL *ri = &vri(9*irefl);
    rotationMatrix(ri,-psi,-th,-phi);
  }

}

/* Calc Euler orientation angles of crystals that are rotated by angle vphi0(i) around
   the given diffraction vector and all its equvalent vectors. Diffraction vector, which
   is considered, is inclened in sample space by angles th0 and psi0.
   Out matrix storage:
   phi((hkl)(0)) th((hkl)(0)) psi((hkl)(0)) ... phi((hkl)(mult)) th((hkl)(mult)) psi((hkl)(mult)) for vphi0(0)
   phi((hkl)(0)) th((hkl)(0)) psi((hkl)(0)) ... phi((hkl)(mult)) th((hkl)(mult)) psi((hkl)(mult)) for vphi0(1)
   ...
   size: vphi0.numElements() x 3*multiplicity */   
CrystMatrix_REAL TextureCalculator::GetCrystalOrientations(REAL th0,REAL psi0,CrystVector_REAL vphi0,
					  REAL h,REAL k,REAL l) const
{
  PrepareForCalc(h,k,l);

  int &multiplicity = mMultiplicity;
  CrystVector_REAL &vri = mvri;
  
  /* Prepare Output Matrix for Euler orientation angles for all equivalent
     diffractions (hkl) and all given angles vphi0(i) of rotation around
     appropriate diffraction vector. */
  // Storage:
  // phi((hkl)(0)) th((hkl)(0)) psi((hkl)(0)) ... phi((hkl)(mult)) th((hkl)(mult)) psi((hkl)(mult)) for vphi0(0)
  // phi((hkl)(0)) th((hkl)(0)) psi((hkl)(0)) ... phi((hkl)(mult)) th((hkl)(mult)) psi((hkl)(mult)) for vphi0(1)
  // ...
  CrystMatrix_REAL out(vphi0.numElements(),3*multiplicity);

  /* Calc Euler orientation angles of crystals that are rotated by angle vphi0(i) around
     the given diffraction vector and all its equvalent vectors. Diffraction vector, which
     is considered, is inclened in sample space by angles th0 and psi0. */
  for(int iphi=0; iphi<vphi0.numElements(); iphi++) {
    REAL phi0 = vphi0(iphi);
    /* calculete Euler angles in the sample coordinates */
    REAL* a = mAMatrix.data();
    REAL* aa = mAAMatrix.data();
    rotationMatrix(a,phi0,th0,psi0); // TODO::can be faster
    /* loop over all symmetry equivalent reflections */
    for(int irefl=0;irefl<multiplicity;irefl++) {
      REAL *ri = &vri(9*irefl);
      matrixMult(a,ri,aa);
      // we need angles not only rotation matrix
      getAngles(aa,out(iphi,3*irefl+0),out(iphi,3*irefl+1),out(iphi,3*irefl+2));
    }
  }

  return out;
}

REAL TextureCalculator::CalcCorr(REAL th0,REAL psi0,REAL h,REAL k,REAL l) const
{
  PrepareForCalc(h,k,l);

  int &multiplicity = mMultiplicity;
  CrystVector_REAL &vri = mvri; 
  REAL &odfnfactor = mOdfNFactor;
  
  /* calculate texture correction */

  #ifdef __DEBUG__ZDENEK__ 
  {
    cout<<"  TextureCalculator::CalcCorr: ";
    cout<<setprecision(2);
    cout<<"th0: "<<setw(8)<<th0*RAD2DEG;
    cout<<",psi0: "<<setw(8)<<psi0*RAD2DEG<<endl;
  }
  #endif /* __DEBUG__ZDENEK__ */

  /* integrate over phi0 */
  REAL sum = 0.;
  for(int iphi=0; iphi<nphi-1; iphi++) {
    REAL phi0 = phimin + iphi*phistep;
    #ifdef __DEBUG__ZDENEK__
    cout<<"phi0: "<<phi0<<endl;
    #endif /* __DEBUG__ZDENEK__ */
    /* calculete Euler angles in the sample coordinates */
    REAL* a = mAMatrix.data();
    REAL* aa = mAAMatrix.data();
    rotationMatrix(a,phi0,th0,psi0); // TODO::can be faster
    #ifdef __DEBUG__ZDENEK__
    cout<<"  a: "<<a<<endl;
    cout<<setprecision(6);
    for(int i=0;i<3;i++) {
      cout<<"    ";
      for(int j=0;j<3;j++) cout<<setw(15)<<a[3*i+j];
      cout<<endl;
    }
    #endif /* __DEBUG__ZDENEK__ */
    /* sum over all symmetry equivalent reflections */
    for(int irefl=0;irefl<multiplicity;irefl++) {
      REAL *ri = &vri(9*irefl);
      matrixMult(a,ri,aa);
      #ifdef __DEBUG__ZDENEK__ 
      {
	cout<<"  ri("<<irefl<<"): "<<ri<<endl;
	cout<<setprecision(6);
	for(int i=0;i<3;i++) {
	  cout<<"    ";
	  for(int j=0;j<3;j++) cout<<setw(15)<<ri[3*i+j];
	  cout<<endl;
	}
	cout<<"  aa: "<<aa<<endl;
	cout<<setprecision(6);
	for(int i=0;i<3;i++) {
	  cout<<"    ";
	  for(int j=0;j<3;j++) cout<<setw(15)<<aa[3*i+j];
	  cout<<endl;
	}
	REAL phi, psi, th;
	getAngles(aa,phi,th,psi);
	cout<<"->phi:"<<phi*RAD2DEG;
	cout<<",th:"<<th*RAD2DEG;
	cout<<",psi:"<<psi*RAD2DEG<<endl;
      }
      #endif /* __DEBUG__ZDENEK__ */
      sum += funcodf(0.,0.,0.,aa); // matrix used, angles omitted
    }
  }
  /* normalize integration sum */
  sum *= phistep;
  sum /= 2.*M_PI*odfnfactor*multiplicity;

  return sum;
}

CrystMatrix_REAL TextureCalculator::CalcCorr(const CrystVector_REAL &vth0,
					     const CrystVector_REAL &vpsi0,
					     REAL h,REAL k,REAL l) const
{
  CrystMatrix_REAL result(vth0.numElements(),vpsi0.numElements());
  
  PrepareForCalc(h,k,l);

  int &multiplicity = mMultiplicity;
  CrystVector_REAL &vri = mvri; 
  REAL &odfnfactor = mOdfNFactor;
  
  /* calculate texture correction */

  for(int i=0;i<vth0.numElements();i++) {
    REAL th0 = vth0(i);
    for(int j=0;j<vpsi0.numElements();j++) {
      REAL psi0 = vpsi0(j);

      /* integrate over phi0 */
      REAL sum = 0.;
      for(int iphi=0; iphi<nphi-1; iphi++) {
	REAL phi0 = phimin + iphi*phistep;
	/* calculete Euler angles in the sample coordinates */
	REAL* a = mAMatrix.data();
	REAL* aa = mAAMatrix.data();
	rotationMatrix(a,phi0,th0,psi0); // TODO::can be faster
	/* sum over all symmetry equivalent reflections */
	for(int irefl=0;irefl<multiplicity;irefl++) {
	  REAL *ri = &vri(9*irefl);
	  matrixMult(a,ri,aa);
	  sum += funcodf(0.,0.,0.,aa); // matrix used, angles omitted
	}
      }
      /* normalize integration sum */
      sum *= phistep;
      sum /= 2.*M_PI*odfnfactor*multiplicity;

      result(i,j) = sum;
    }
  }

  return result;
}

void TextureCalculator::ExportOdf(const string filename) const
{
  ofstream f(filename.c_str());
  
  REAL phi, th, psi, odfvalue;

  for(int iphi=0; iphi<nphi; iphi++) {
    phi = phimin + iphi*phistep;
    for(int ith=0; ith<nth; ith++) {
      th = thmin + ith*thstep;
      for(int ipsi=0; ipsi<npsi; ipsi++) {
	psi = psimin + ipsi*psistep;
	odfvalue = funcodf(phi,th,psi,NULL);
	f << fixed << showpoint << setprecision(2) << setw(10) << phi*RAD2DEG << setw(10) << th*RAD2DEG;
	f << setw(10) << psi*RAD2DEG<< scientific <<  setprecision(3) << setw(14) << odfvalue << "\n";
      }
    }
  }

  f.close();
}

// ScatteringCorr
ScatteringCorr::ScatteringCorr(const ScatteringData & data):
ObjCryst::ScatteringCorr(data)
{}

const CrystVector_REAL&  ScatteringCorr::GetCorr(bool needRecalc)const
{
  if(needRecalc==true) this->CalcCorr();
  return mCorr;
}

// AbsorptionCorr
AbsorptionCorr::AbsorptionCorr(const ScatteringData & data):
ScatteringCorr(data),mOmega(-1.),mAbsFactor(0.),mThickness(-1.),mDepth(0.)
{}

AbsorptionCorr::~AbsorptionCorr()
{}

const string & AbsorptionCorr::GetName() const
{
   //So far, we do not need a personalized name...
   const static string mName="AbsorptionCorr";
   return mName;
}

const string & AbsorptionCorr::GetClassName() const
{
   const static string className="AbsorptionCorr";
   return className;
}

void AbsorptionCorr::SetAbsorptionCorrParams(REAL thickness, REAL depth,
					     REAL absfactor, REAL omega)
{
  mOmega = omega;
  mThickness = thickness*10.; // nm -> A
  mDepth = depth*10.; // nm -> A
  mAbsFactor = absfactor*1.e-8; // cm -> A
  mClockCorrCalc.Reset();
}

void AbsorptionCorr::CalcCorr() const
{
   const CrystVector_REAL *theta=&(mpData->GetTheta());
   if(mpData->GetClockTheta()<mClockCorrCalc) return;
   VFN_DEBUG_MESSAGE("AbsorptionCorr::CalcCorr()",10)
   TAU_PROFILE("AbsorptionCorr::CalcCorr()","void ()",TAU_DEFAULT);
   mCorr.resize(mpData->GetNbRefl());
   #ifdef __DEBUG__
   stringstream s;
   s<< "   AbsorptionCorr: " << endl;
   s<<setw(10)<<"omega"<<setw(15)<<"pen.depth(tp)";
   s<<setw(15)<<"exp(-depth/tp)"<<" "<<"(1.-exp(-thickness/tp))";
   s<<" corr"<<endl;
   #endif
   for(long i=0;i<mpData->GetNbRefl();i++) {
     REAL omega = (mOmega>0.) ? mOmega : (*theta)(i);
     REAL tp = 1./(1./sin(omega) + 1./sin(2.*((*theta)(i))-omega))/mAbsFactor;
     REAL corr = (mThickness>=0.) ? tp*(1.-exp(-mThickness/tp))*exp(-mDepth/tp)
     															: tp*                         exp(-mDepth/tp);
     #ifdef __DEBUG__
     s<<setw(10)<<omega*RAD2DEG<<setw(15)<<tp/10.;
     s<<setw(15)<<exp(-mDepth/tp)<<setw(15)<<( (mThickness>=0.) ? (1.-exp(-mThickness/tp)) : 0. );
     if(mOmega*RAD2DEG<=-1.999) // Bragg-Brentano with variable slits
     	 s<<setw(15)<<corr/absorption_corr_factor<<endl;
     else
     	 s<<setw(15)<<corr/sin(omega)/absorption_corr_factor<<endl;
     #endif
     //mCorr(i) = corr/sin(omega)/absorption_corr_factor;
     if (mOmega*RAD2DEG<=-1.999) // Bragg-Brentano with variable slits
       mCorr(i) = corr/absorption_corr_factor;
     else
     	 mCorr(i) = corr/sin(omega)/absorption_corr_factor;
   }
   #ifdef __DEBUG__
   VFN_DEBUG_MESSAGE("AbsorptionCorr::CalcCorr(): "<<endl<<s.str(),11)
   #endif
   mClockCorrCalc.Click();
}

// TextureCorr
TextureCorr::TextureCorr(const ScatteringData & data):
  ScatteringCorr(data),mOmega(-1.),mDivergenceOmega(0.),mDivergence2Theta(0.),
  mDivergencePsi_i(0.04),mDivergencePsi_f(5.*DEG2RAD),mTiltTh(0.),mTiltPsi(0.)
{
  this->InitParameters();
}

TextureCorr::TextureCorr(const TextureCorr & old):
ScatteringCorr(*old.mpData)
{
  this->InitParameters();
}

TextureCorr::~TextureCorr()
{}

const string & TextureCorr::GetName() const
{
   //So far, we do not need a personalized name...
   const static string mName="TextureCorr";
   return mName;
}

const string & TextureCorr::GetClassName() const
{
   const static string className="TextureCorr";
   return className;
}

void TextureCorr::AddPhase(const REAL fraction,const REAL thweight,
			   const REAL th0,const REAL thwidth,
			   const REAL psiweight,const REAL psi0,
			   const REAL psiwidth)
{
  TexturePhase *phase = new TexturePhase;
  phase->mTextureObj.SetCrystal(mpData->GetCrystal());
  // Set parameters
  phase->fraction=fraction;
  phase->params(0)=thweight;phase->params(1)=th0;phase->params(2)=thwidth;
  phase->params(3)=psiweight;phase->params(4)=psi0;phase->params(5)=psiwidth;
  const int nbPhase=this->GetNbPhase();
  // Add parameters
  char buf [5];
  sprintf(buf,"%d",nbPhase);
  {
    RefinablePar tmp("Fraction_"+(string)buf,&(phase->fraction),0.,1.,
		     gpRefParTypeScattDataCorrInt,REFPAR_DERIV_STEP_ABSOLUTE,
		     true,true,true,false,1.);
    tmp.AssignClock(phase->mClockParams);
    tmp.SetDerivStep(1e-2);
    //tmp.SetGlobalOptimStep(.05);
    this->AddPar(tmp);
  }
  {
    RefinablePar tmp("ThWeight_"+(string)buf,&(phase->params(0)),0.,1.,
		     gpRefParTypeScattDataCorrInt,REFPAR_DERIV_STEP_ABSOLUTE,
		     true,true,true,false,1.);
    tmp.AssignClock(phase->mClockParams);
    tmp.SetDerivStep(1e-2);
    //tmp.SetGlobalOptimStep(.05);
    this->AddPar(tmp);
  }
  {
    RefinablePar tmp("Th0_"+(string)buf,&(phase->params(1)),0.,180.*DEG2RAD,
		     gpRefParTypeScattDataCorrInt,REFPAR_DERIV_STEP_ABSOLUTE,
		     true,true,true,true,RAD2DEG,180.*DEG2RAD);
    tmp.AssignClock(phase->mClockParams);
    tmp.SetDerivStep(2.5*DEG2RAD);
    //tmp.SetGlobalOptimStep(.05);
    this->AddPar(tmp);
  }
  {
    RefinablePar tmp("ThWidth_"+(string)buf,&(phase->params(2)),2.5*DEG2RAD,360.*DEG2RAD,
		     gpRefParTypeScattDataCorrInt,REFPAR_DERIV_STEP_ABSOLUTE,
		     true,true,true,false,RAD2DEG);
    tmp.AssignClock(phase->mClockParams);
    tmp.SetDerivStep(2.5*DEG2RAD);
    //tmp.SetGlobalOptimStep(.05);
    this->AddPar(tmp);
  }
  {
    RefinablePar tmp("PsiWeight_"+(string)buf,&(phase->params(3)),0.,1.,
		     gpRefParTypeScattDataCorrInt,REFPAR_DERIV_STEP_ABSOLUTE,
		     true,true,true,false,1.);
    tmp.AssignClock(phase->mClockParams);
    tmp.SetDerivStep(1e-2);
    //tmp.SetGlobalOptimStep(.05);
    this->AddPar(tmp);
  }
  {
    RefinablePar tmp("Psi0_"+(string)buf,&(phase->params(4)),-180.*DEG2RAD,180.*DEG2RAD,
		     gpRefParTypeScattDataCorrInt,REFPAR_DERIV_STEP_ABSOLUTE,
		     true,true,true,true,RAD2DEG,360.*DEG2RAD);
    tmp.AssignClock(phase->mClockParams);
    tmp.SetDerivStep(2.5*DEG2RAD);
    //tmp.SetGlobalOptimStep(.05);
    this->AddPar(tmp);
  }
  {
    RefinablePar tmp("PsiWidth_"+(string)buf,&(phase->params(5)),2.5*DEG2RAD,360.*DEG2RAD,
		     gpRefParTypeScattDataCorrInt,REFPAR_DERIV_STEP_ABSOLUTE,
		     true,true,true,false,RAD2DEG);
    tmp.AssignClock(phase->mClockParams);
    tmp.SetDerivStep(2.5*DEG2RAD);
    //tmp.SetGlobalOptimStep(.05);
    this->AddPar(tmp);
  }

  this->GetPar(&mTiltTh).SetIsUsed(true);
  this->GetPar(&mTiltPsi).SetIsUsed(true);

  mPhaseRegistry.Register(*phase);
  mClockMaster.AddChild(phase->mClockParams);
  mClockMaster.Click();
}

void TextureCorr::AddPhase(const REAL fraction,
							const CrystVector_REAL params,bool bForceTextureSymmetry)
{
  TexturePhase *phase = new TexturePhase;
  phase->mTextureObj.SetCrystal(mpData->GetCrystal());
  // Set parameters
  phase->fraction=fraction;
  phase->params = params;
  phase->bForceTextureSymmetry = bForceTextureSymmetry;
  
  const int nbPhase=this->GetNbPhase();
  // Add parameters
  char buf [5];
  sprintf(buf,"%d",nbPhase);
  {
    RefinablePar tmp("Fraction_"+(string)buf,&(phase->fraction),0.,1.,
		     gpRefParTypeScattDataCorrInt,REFPAR_DERIV_STEP_ABSOLUTE,
		     true,true,true,false,1.);
    tmp.AssignClock(phase->mClockParams);
    tmp.SetDerivStep(1e-2);
    //tmp.SetGlobalOptimStep(.05);
    this->AddPar(tmp);
  }
  {
    RefinablePar tmp("ThWeight_"+(string)buf,&(phase->params(0)),0.,1.,
		     gpRefParTypeScattDataCorrInt,REFPAR_DERIV_STEP_ABSOLUTE,
		     true,true,true,false,1.);
    tmp.AssignClock(phase->mClockParams);
    tmp.SetDerivStep(1e-2);
    //tmp.SetGlobalOptimStep(.05);
    this->AddPar(tmp);
  }
  {
    RefinablePar tmp("Th0_"+(string)buf,&(phase->params(1)),0.,180.*DEG2RAD,
		     gpRefParTypeScattDataCorrInt,REFPAR_DERIV_STEP_ABSOLUTE,
		     true,true,true,true,RAD2DEG,180.*DEG2RAD);
    tmp.AssignClock(phase->mClockParams);
    tmp.SetDerivStep(2.5*DEG2RAD);
    //tmp.SetGlobalOptimStep(.05);
    this->AddPar(tmp);
  }
  {
    RefinablePar tmp("ThWidth_"+(string)buf,&(phase->params(2)),2.5*DEG2RAD,360.*DEG2RAD,
		     gpRefParTypeScattDataCorrInt,REFPAR_DERIV_STEP_ABSOLUTE,
		     true,true,true,false,RAD2DEG);
    tmp.AssignClock(phase->mClockParams);
    tmp.SetDerivStep(2.5*DEG2RAD);
    //tmp.SetGlobalOptimStep(.05);
    this->AddPar(tmp);
  }
  {
    RefinablePar tmp("PsiWeight_"+(string)buf,&(phase->params(3)),0.,1.,
		     gpRefParTypeScattDataCorrInt,REFPAR_DERIV_STEP_ABSOLUTE,
		     true,true,true,false,1.);
    tmp.AssignClock(phase->mClockParams);
    tmp.SetDerivStep(1e-2);
    //tmp.SetGlobalOptimStep(.05);
    this->AddPar(tmp);
  }
  {
    RefinablePar tmp("Psi0_"+(string)buf,&(phase->params(4)),-180.*DEG2RAD,180.*DEG2RAD,
		     gpRefParTypeScattDataCorrInt,REFPAR_DERIV_STEP_ABSOLUTE,
		     true,true,true,true,RAD2DEG,360.*DEG2RAD);
    tmp.AssignClock(phase->mClockParams);
    tmp.SetDerivStep(2.5*DEG2RAD);
    //tmp.SetGlobalOptimStep(.05);
    this->AddPar(tmp);
  }
  {
    RefinablePar tmp("PsiWidth_"+(string)buf,&(phase->params(5)),2.5*DEG2RAD,360.*DEG2RAD,
		     gpRefParTypeScattDataCorrInt,REFPAR_DERIV_STEP_ABSOLUTE,
		     true,true,true,false,RAD2DEG);
    tmp.AssignClock(phase->mClockParams);
    tmp.SetDerivStep(2.5*DEG2RAD);
    //tmp.SetGlobalOptimStep(.05);
    this->AddPar(tmp);
  }
  
  // hklx, hklz and tilt are not refinaable

  //this->GetPar(&mTiltTh).SetIsUsed(true);
  //this->GetPar(&mTiltPsi).SetIsUsed(true);

  mPhaseRegistry.Register(*phase);
  mClockMaster.AddChild(phase->mClockParams);
  mClockMaster.Click();
}

int TextureCorr::GetNbPhase() const {return mPhaseRegistry.GetNb();}

void TextureCorr::SetTextureCorrParams(REAL omega)
{
  mOmega = omega;
  mClockMaster.Click();
  //mClockCorrCalc.Reset();
}

void TextureCorr::SetCrystal(ObjCryst::Crystal &crystal)
{
  for(int i=0;i<GetNbPhase();i++) {
    TexturePhase &phase=mPhaseRegistry.GetObj(i);
    phase.mTextureObj.SetCrystal(crystal);
  }
  mClockMaster.Click();
  //mClockCorrCalc.Reset();
}

void TextureCorr::CalcCorr() const
{
  // TODO:: can be much faster (I dont't have to recalc all texture obj)
   const CrystVector_REAL *theta=&(mpData->GetTheta());
   if(mpData->GetClockTheta()<mClockCorrCalc &&
      mClockMaster<mClockCorrCalc) return;
   VFN_DEBUG_MESSAGE("TextureCorr::CalcCorr()",11)
   TAU_PROFILE("TextureCorr::CalcCorr()","void ()",TAU_DEFAULT);
   const REAL lambda = mpData->GetRadiation().GetWavelength()(0);
   const CrystVector_REAL *sinthovlam=&(mpData->GetSinThetaOverLambda());
   // update params in TextureCalculators objects
   REAL fractionsum = 0.;
   for(int i=0;i<GetNbPhase();i++) {
     TexturePhase &phase=mPhaseRegistry.GetObj(i);
     if(fabs(phase.fraction)<FLT_MIN) continue;
     if(mClockCorrCalc<phase.mClockParams || mClockCorrCalc<mClockTiltParams) {
       //CrystVector_REAL params(8);
       //for(int j=0;j<6;j++) params(j)=phase.params(j);
       //params(6)=mTiltTh;params(7)=mTiltPsi;
       //phase.mTextureObj.SetTextureParams(params);
       phase.mTextureObj.SetTextureParams(phase.params,phase.bForceTextureSymmetry);
     }
     fractionsum += phase.fraction; 
   }
   // calc
   mCorr.resize(mpData->GetNbRefl());
   const CrystVector_REAL &H = mpData->GetH();
   const CrystVector_REAL &K = mpData->GetK();
   const CrystVector_REAL &L = mpData->GetL();
   #ifdef __DEBUG__
   stringstream s;
   s<<"   TextureCorr: "<<endl;
   s<<setw(10)<<"2theta"<<setw(6)<<"h"<<setw(4)<<"k"<<setw(4)<<"l";
   s<<setw(12)<<"psi"<<setw(12)<<"corr"<<endl;
   #endif
   for(long i=0;i<mpData->GetNbRefl();i++) {
     REAL omega = (mOmega>0.) ? mOmega : (*theta)(i);
     REAL psi = (*theta)(i)-omega;
     REAL corr = 1.;
     if(fabs(fractionsum)>FLT_MIN) {
       REAL sum = 0.;
       for(int j=0;j<GetNbPhase();j++) {
	 const TexturePhase &phase=mPhaseRegistry.GetObj(j);
	 if(fabs(phase.fraction)>FLT_MIN) {
	   // Integrate over axial divergence
	   const REAL sc = 1./2./((*sinthovlam)(i)*lambda);
	   const REAL phi_min  = -0.5*(mDivergencePsi_f+mDivergencePsi_i)*sc;
	   const REAL phi_step = (mDivergencePsi_f+mDivergencePsi_i)*sc/10;
	   for(int iax=1; iax<10; iax++) {
	     const REAL phi = phi_min+iax*phi_step;
	     const REAL interval_lenght = max(min(0.5*mDivergencePsi_f*sc-phi,0.5*mDivergencePsi_i*sc),
					      -0.5*mDivergencePsi_i*sc)
	                                 -max(-0.5*mDivergencePsi_f*sc-phi,-0.5*mDivergencePsi_i*sc);
	     // for this divergence idea
	     //sum += phase.fraction
	     //          * phase.mTextureObj.CalcCorr(psi,phi,H(i),K(i),L(i))
	     //          * interval_lenght * phi_step * RAD2DEG;
	   }
	   	 // else simply
	   	 sum += phase.fraction*phase.mTextureObj.CalcCorr(psi,0.,H(i),K(i),L(i));
	 }
       }
       corr *= sum/fractionsum;
     }
     #ifdef __DEBUG__
     s<<setw(10)<<2*(*theta)(i)*RAD2DEG;
     s<<setw(6)<<H(i)<<setw(4)<<K(i)<<setw(4)<<L(i);
     s<<setw(12)<<psi*RAD2DEG<<setw(12)<<corr<<endl;
     #endif
     mCorr(i) = corr;
   }
   #ifdef __DEBUG__
   VFN_DEBUG_MESSAGE("TextureCorr::CalcCorr(): "<<endl<<s.str(),11)
   #endif
   mClockCorrCalc.Click();
}

void TextureCorr::InitParameters()
{
  /*{
    RefinablePar tmp("tiltTh0",&mTiltTh,-180.*DEG2RAD,180.*DEG2RAD,
		     gpRefParTypeScattDataCorrInt,REFPAR_DERIV_STEP_ABSOLUTE,
		     true,true,false,true,RAD2DEG,360.*DEG2RAD);
    tmp.AssignClock(mClockTiltParams);
    tmp.SetDerivStep(2.5*DEG2RAD);
    this->AddPar(tmp);
  }
  {
    RefinablePar tmp("tiltPsi0",&mTiltPsi,0.*DEG2RAD,180.*DEG2RAD,
		     gpRefParTypeScattDataCorrInt,REFPAR_DERIV_STEP_ABSOLUTE,
		     true,true,false,true,RAD2DEG,180.*DEG2RAD);
    tmp.AssignClock(mClockTiltParams);
    tmp.SetDerivStep(2.5*DEG2RAD);
    this->AddPar(tmp);
  }*/
  mClockMaster.AddChild(mClockTiltParams);
}

// HKLIntensityCorr
HKLIntensityCorr::HKLIntensityCorr(const ScatteringData & data):
ScatteringCorr(data)
{}

HKLIntensityCorr::HKLIntensityCorr(const HKLIntensityCorr & old):
ScatteringCorr(*old.mpData)
{
  mReflStore.clear();
  for(int i=0;i<old.mReflStore.size();i++) {
    IntensityCorrData *pNewData = new IntensityCorrData;
    const ReflData &d = old.mReflStore.at(i);
    *pNewData = *((IntensityCorrData*)d.data);
    mReflStore.add(d.H,d.K,d.L,d.x,(void*)pNewData);
  }
}

HKLIntensityCorr::~HKLIntensityCorr()
{
  /*for(int i=0;i<mReflStore.size();i++) {
    IntensityCorrData *pData = (IntensityCorrData*) mReflStore.at(i).data;
    ObjCryst::RefinablePar *pPar = &GetPar(&(pData->val));
    RemovePar(pPar);
    delete pPar;
  }*/
  Reset();
}

const string & HKLIntensityCorr::GetName() const
{
   //So far, we do not need a personalized name...
   const static string mName="HKLIntensityCorr";
   return mName;
}

const string & HKLIntensityCorr::GetClassName() const
{
   const static string className="HKLIntensityCorr";
   return className;
}

void HKLIntensityCorr::SetHKLIntensityCorr(int h,int k,int l,REAL val,bool fixed)
{
  // try to find (hkl) reflection in the "store"
  int ind = mReflStore.find(h,k,l,0.);

  // generate parameter name
  string name;
  {
    stringstream s;
    s << mpData->GetName() << "_Ihkl_";
    s << h << "_" << k << "_" << l;
    name = s.str();
  }
   
  // if refl. found use it else add new refl. into the "store"
  IntensityCorrData *pData = 0;
  if(ind>=0) {
    pData = (IntensityCorrData*) mReflStore.at(ind).data;
  }
  else {
    pData = new IntensityCorrData;
    mReflStore.add(h,k,l,0.,(void*)pData);

    RefinablePar tmp(name,&(pData->val),0.,1.e5,
		     gpRefParTypeScattDataCorrHKLIntensity,
		     REFPAR_DERIV_STEP_ABSOLUTE,
		     true,false,false);
    tmp.AssignClock(mClockMaster);
    tmp.SetDerivStep(0.001);
    this->AddPar(tmp);
  }

  // set params
  RefinablePar &par = GetPar(name);
  par.SetValue(val);
  par.SetIsFixed(fixed);
  
  // Is this refl currently used by ScatteringData object?
  pData->i = -1;
  for(int i=0;i<mpData->GetNbRefl();i++)
    if(int(mpData->GetH()(i))==h && int(mpData->GetK()(i))==k &&
       int(mpData->GetL()(i))==l) {
      pData->i = i;
      par.SetIsUsed(true);
      break;
    }
  
  mClockCorrCalc.Reset();
}

void HKLIntensityCorr::BeginOptimization(const bool allowApproximations,
					 const bool enableRestraints)
{
  RefinableObj::BeginOptimization(allowApproximations,enableRestraints);

  for(int i=0;i<mReflStore.size();i++) {
    IntensityCorrData *pData=(IntensityCorrData*)mReflStore.at(i).data;
    int H=mReflStore.at(i).H,K=mReflStore.at(i).K,L=mReflStore.at(i).L;
    pData->i = -1;
    GetPar(&(pData->val)).SetIsUsed(false);
    for(long j=0;j<mpData->GetNbRefl();j++)
      if(int(mpData->GetH()(j))==H && int(mpData->GetK()(j))==K &&
	 int(mpData->GetL()(j))==L) {
	pData->i = j;
	GetPar(&(pData->val)).SetIsUsed(true);
	break;
      }
  }
}

void HKLIntensityCorr::CalcCorr() const
{
  //const CrystVector_REAL &theta=mpData->GetTheta();
  if(mpData->GetClockTheta()<mClockCorrCalc &&
     mClockMaster<mClockCorrCalc) return;
  VFN_DEBUG_MESSAGE("HKLIntensityCorr::CalcCorr()",10)
  TAU_PROFILE("HKLIntensityCorr::CalcCorr()","void ()",TAU_DEFAULT);
  if(mpData->GetClockTheta()>mClockIndexesCalc) {
    for(int i=0;i<mReflStore.size();i++) {
      IntensityCorrData *pData=(IntensityCorrData*)mReflStore.at(i).data;
      int H=mReflStore.at(i).H,K=mReflStore.at(i).K,L=mReflStore.at(i).L;
      pData->i = -1;
      //GetPar(&(pData->val)).SetIsUsed(false);
      for(long j=0;j<mpData->GetNbRefl();j++)
	if(int(mpData->GetH()(j))==H && int(mpData->GetK()(j))==K &&
	   int(mpData->GetL()(j))==L) {
	  pData->i = j;
	  //GetPar(&(pData->val)).SetIsUsed(true);
	  break;
	}
    }
    mClockIndexesCalc.Click();
  }
  mCorr.resize(mpData->GetNbRefl());
  mCorr = 1.;
  for(int i=0;i<mReflStore.size();i++) {
    IntensityCorrData *pData=(IntensityCorrData*)mReflStore.at(i).data;
    if(pData->i>=0) mCorr(pData->i) = pData->val;
  }
  #ifdef __DEBUG__
  stringstream s;
  s << "   HKLIntensityCorr: " << endl;
  s<<setw(10)<<"2theta"<<setw(6)<<"h"<<setw(4)<<"k"<<setw(4)<<"l";
  s<<setw(12)<<"corr"<<endl;
  for(long i=0;i<mpData->GetNbRefl();i++) {
    s<<setw(10)<<2*mpData->GetTheta()(i)*RAD2DEG;
    s<<setw(6)<<(int)mpData->GetH()(i);
    s<<setw(4)<<(int)mpData->GetK()(i);
    s<<setw(4)<<(int)mpData->GetL()(i);
    s<<setw(12)<<mCorr(i)<<endl;
  }
  VFN_DEBUG_MESSAGE("HKLIntensityCorr::CalcCorr(): "<<endl<<s.str(),11)
  #endif
  mClockCorrCalc.Click();
}

void HKLIntensityCorr::Reset()
{
  for(int i=0;i<mReflStore.size();i++) {
    IntensityCorrData *pData=(IntensityCorrData*)mReflStore.at(i).data;
    this->RemovePar(&(this->GetPar(&(pData->val))));
    delete pData;
  }
  mReflStore.clear();
  mClockMaster.Click();
}

const ReflStore& HKLIntensityCorr::GetReflStore()const
{
  return mReflStore;
}

// PowderPatternDiffraction
PowderPatternDiffraction::PowderPatternDiffraction()
:mOmega(-1.),mCorrAbsorption(*this),mCorrTexture(*this),
 mCorrHKLIntensity(*this)
{
  mClockMaster.AddChild(mCorrTexture.GetClockMaster());
  this->AddSubRefObj(mCorrTexture);
  this->AddSubRefObj(mCorrHKLIntensity);
  //mCorrHKLIntensity.RegisterClient(*this);

  mReflProfFact = 2.0;
  mReflProfMinRelIntensity = 0.001;
}

PowderPatternDiffraction::PowderPatternDiffraction(const PowderPatternDiffraction &old):
ObjCryst::PowderPatternDiffraction(old),
mOmega(old.mOmega),mCorrAbsorption(old.mCorrAbsorption),
mCorrTexture(old.mCorrTexture),
mCorrHKLIntensity(old.mCorrHKLIntensity)
{
  //this->AddSubRefObj(mCorrTextureMarchDollase);
  //this->SetProfile(old.mpReflectionProfile->CreateCopy());
  //mClockMaster.AddChild(mClockProfilePar);
  //mClockMaster.AddChild(mClockLorentzPolarSlitCorrPar);
  //mClockMaster.AddChild(mpReflectionProfile->GetClockMaster());
  mClockMaster.AddChild(mCorrTexture.GetClockMaster());
  this->AddSubRefObj(mCorrTexture);
  this->AddSubRefObj(mCorrHKLIntensity);
  //mCorrHKLIntensity.RegisterClient(*this);
}

PowderPatternDiffraction* PowderPatternDiffraction::CreateCopy()const
{
   return new PowderPatternDiffraction(*this);
}

const string& PowderPatternDiffraction::GetClassName() const
{
  const static string className="MStruct::PowderPatternDiffraction";
  return className;
} 

void PowderPatternDiffraction::SetCrystal(ObjCryst::Crystal &crystal)
{
  ObjCryst::PowderPatternDiffraction::SetCrystal(crystal);
  mCorrTexture.SetCrystal(crystal);
}

void PowderPatternDiffraction::SetAbsorptionCorrParams(REAL thickness,
					 REAL depth,REAL absfactor,REAL omega)
{
  mOmega = omega;
  mCorrAbsorption.SetAbsorptionCorrParams(thickness,depth,absfactor,omega);
  mClockIntensityCorr.Reset();
  mClockIhklCalc.Reset();
}

void PowderPatternDiffraction::SetTextureCorrParams(REAL omega)
{
  mOmega = omega;
  mCorrTexture.SetTextureCorrParams(omega);
  mClockIntensityCorr.Reset();
  mClockIhklCalc.Reset();
}

void PowderPatternDiffraction::AddTextureCorrPhase(REAL fraction,
					   const CrystVector_REAL& params, bool bForceTextureSymmetry)
{
  mCorrTexture.AddPhase(fraction,params,bForceTextureSymmetry);
}

void PowderPatternDiffraction::SetHKLIntensityCorrParams(int h, int k, int l,
							 REAL val,bool fixed)
{
  mCorrHKLIntensity.SetHKLIntensityCorr(h,k,l,val,fixed);
  mClockIntensityCorr.Reset();
  mClockIhklCalc.Reset();
}

void PowderPatternDiffraction::GenerateHKLIntensityCorrForAllReflections(const REAL relIntensity)
{
  Prepare();
  this->CalcIhkl();
  CrystVector_REAL intensity = mIhklCalc;
  for(int i=0;i<GetNbRefl();i++) intensity(i) *= mMultiplicity(i);
  REAL averageIntensity = intensity.sum()/GetNbRefl();
  for(int i=0;i<GetNbRefl();i++) {
    bool fixed = relIntensity<0. || intensity(i)<relIntensity*averageIntensity;
    mCorrHKLIntensity.SetHKLIntensityCorr((int)mH(i),(int)mK(i),(int)mL(i),
					  1.,fixed);
  }
  mClockIntensityCorr.Reset();
  mClockIhklCalc.Reset();
}

void PowderPatternDiffraction::PrintHKLIntensityCorr(ostream& s)const
{
  this->GetFhklCalcSq();

	s<<"# "<<setw(4)<<"h"<<setw(4)<<"k"<<setw(4)<<"l";
  s<<setw(14)<<"2Theta(deg)"<<setw(12)<<"|Fhkl|^2";
  s<<setw(12)<<"Icor"<<setw(8)<<"fixed"<<"\n";
  const ReflStore &reflStore=mCorrHKLIntensity.GetReflStore();
  for(int i=0;i<reflStore.size();i++) {
    s<<setw(6)<<reflStore.at(i).H<<setw(4)<<reflStore.at(i).K;
    s<<setw(4)<<reflStore.at(i).L;
    const HKLIntensityCorr::IntensityCorrData *data = 
      (const HKLIntensityCorr::IntensityCorrData*)reflStore.at(i).data;
    if((data->i<0) || (data->i)>=GetNbRefl()) {
      s<<setw(14)<<data->i<<"\n";
      continue;
    }
    s<<setw(14)<<setprecision(3)<<fixed<<2*mTheta(data->i)*RAD2DEG;
    s<<setw(12)<<setprecision(2)<<scientific<<mFhklCalcSq(data->i)*mMultiplicity(data->i);
    s<<setw(12)<<setprecision(2)<<fixed<<data->val;
    s<<setw(8)<<mCorrHKLIntensity.GetPar(&(data->val)).IsFixed()<<"\n";
  }
}

void PowderPatternDiffraction::WriteHKLIntensityCorrToFile(const char* filename) const
{
  // copy content of the file
  ostringstream s;
  {
    char buffer[1024];
    ifstream f(filename);
    while(f.getline(buffer,1022)) s << buffer << '\n';
    f.close();
  }

  // WriteHKLIntensityCorrToFile
  ofstream f(filename);
  this->PrintHKLIntensityCorr(f);
  
  // save original file content at the end of the file
  f<<'\n'<<s.str();

  f.close();
}

void PowderPatternDiffraction::ReadHKLIntensityCorrFromFile(
						  const char* filename)
{
  mCorrHKLIntensity.Reset();

  ifstream f(filename);
  string line;
  while (getline(f,line)) {
    if (line.empty()) break;
    if (line.at(0)=='#') continue;
    int h, k, l, fixed;
    REAL theta2, fhkl2, corr;
    istringstream s(line);
    s >> h >> k >> l >> theta2 >> fhkl2 >> corr >> fixed;
    if (s.fail()==true) break;
    this->SetHKLIntensityCorrParams(h,k,l,corr,fixed==1);
  }
  f.close();
}
		       
  //CrystVector_REAL PowderPatternDiffraction::GetIncidenceIngle
  //                                   (const CrystVector_REAL& theta)
  //{
  //CrystVector_REAL omega(theta.numElements());
  //CrystVector_REAL omega = (mOmega>0.) ? mOmega : theta/2.;
  //return omega;
  //}

const CrystVector_REAL& PowderPatternDiffraction::GetLSQDeriv(
			      const unsigned int n, RefinablePar &par)
{
  string str = string(par.GetName());
  string::size_type loc;

  loc = str.find("_Ihkl_", 0); // CompName_Ihkl_h_k_l
  if(loc == string::npos) {
    return ObjCryst::PowderPatternDiffraction::GetLSQDeriv(n,par);
  }
  else {
    // get (hkl)
    istringstream s(str.data()+loc+6);
    int h, k, l;
    s >> h; s.get(); s >> k; s.get(); s >> l;
    
    // find reflection index
    int ind = -1;
    for(int i=0;i<GetNbRefl();i++)
      if(int(mH(i))==h && int(mK(i))==k && int(mL(i))==l) { ind = i; break; }
    if(ind<0) {
      mLSQDeriv.resize(mPowderPatternCalc.numElements());
      mLSQDeriv = 0.;
      return mLSQDeriv;
    }
  
    // calculation (trick)
    REAL val = par.GetValue();
    if(fabs(val)>FLT_MIN*1.e4) {
      CrystVector_REAL corr = mIntensityCorr;
      mIntensityCorr = 0.; mIntensityCorr(ind) = corr(ind)/val;
      mClockIntensityCorr.Click(); mClockMaster.Click();
      mLSQDeriv = GetPowderPatternCalc();
      mIntensityCorr = corr;
      mClockIntensityCorr.Click(); mClockMaster.Click();
    }
    else {
      SetHKLIntensityCorrParams(h,k,l,1.); CalcIntensityCorr();
      REAL tmp = mIntensityCorr(ind); 
      mIntensityCorr = 0.; mIntensityCorr(ind) = tmp;
      mClockIntensityCorr.Click(); mClockMaster.Click();
      mLSQDeriv = GetPowderPatternCalc();
      SetHKLIntensityCorrParams(h,k,l,val); CalcIntensityCorr();
      mClockMaster.Click();
    }
    return mLSQDeriv;
  }
}

// TODO: many routines shoul be const
void PowderPatternDiffraction::PrintHKLInfo (ostream &s)
{
  this->CalcIhkl();

  // TODO:: doesn't work for TOF data
  
  // get MStruct::ReflectionProfile Object
  if(mpReflectionProfile->GetClassName().compare("MStruct::ReflectionProfile")!=0)
    return;

  MStruct::ReflectionProfile &reflProf =
    *dynamic_cast<MStruct::ReflectionProfile*>(mpReflectionProfile);
  
  // print header
  s<<"#"<<setw(7)<<"h"<<setw(3)<<"k"<<setw(3)<<"l";
  s<<setw(12)<<"2Theta"<<setw(12)<<"|Fhkl|^2"<<setw(12)<<"Ihkl";
  s<<setw(12)<<"FWHM(deg)"<<setw(12)<<"B(deg)"<<"\n";
  // print data
  for(int irefl=0; irefl<this->GetNbRefl(); irefl++) {
    int h = (int)mH(irefl), k = (int)mK(irefl), l = (int)mL(irefl);
    s<<setw(8)<<h<<setw(3)<<k<<setw(3)<<l;

    REAL x0=mpParentPowderPattern->STOL2X(mSinThetaLambda(irefl));
    REAL center = mpParentPowderPattern->X2XCorr(x0);
    center += reflProf.GetPositionCorr(center,h,k,l);
    s<<fixed<<setprecision(3)<<setw(12)<<center*RAD2DEG;
    
    s<<setw(12)<<setprecision(2)<<scientific<<mFhklCalcSq(irefl)*mMultiplicity(irefl);
    
    s<<setw(12)<<setprecision(2)<<scientific<<mIhklCalc(irefl);
    
    s<<setw(12)<<setprecision(3)<<fixed;
    s<<reflProf.GetFullProfileWidth(0.5,center,h,k,l)*RAD2DEG;
    
    s<<setw(12)<<setprecision(3)<<fixed;
    s<<reflProf.GetIntegralWidth(center,h,k,l)*RAD2DEG;

    s<<"\n";
  }
}

// TODO: many routines shoul be const
void PowderPatternDiffraction::PrintHKLInfo2 (ostream &s, const REAL accur) const
{
  this->CalcIhkl();
  
  // get ReflectionProfile Object

  ObjCryst::ReflectionProfile &reflProf = *mpReflectionProfile;
  
  // print header
  s<<"#"<<setw(7)<<"h"<<setw(3)<<"k"<<setw(3)<<"l";
  s<<setw(12)<<"center"<<setw(12)<<"intensity"<<setw(12);
  s<<setw(12)<<"fwhm"<<setw(12)<<"beta"<<setw(12)<<"asym"<<"\n";
  
  // profile calculation parametrs
  int n = 1001;
	std::vector<REAL> width(this->GetNbRefl());
  
  // calculate data for all diffractions
  std::vector<PeakParams> params0(this->GetNbRefl());
	std::vector<PeakParams> params1(this->GetNbRefl());
  
  for(int irefl=0; irefl<this->GetNbRefl(); irefl++) {
  	
  	int h = (int)mH(irefl), k = (int)mK(irefl), l = (int)mL(irefl);
  	
  	// get peak center
    REAL x0=mpParentPowderPattern->STOL2X(mSinThetaLambda(irefl));
    REAL center = mpParentPowderPattern->X2XCorr(x0);
    
  	// range of calculated data
  	width[irefl] = reflProf.GetFullProfileWidth(0.001,center,h,k,l);
    
    // calculate profile
		CrystVector_REAL x(n);
		const REAL step = width[irefl]/(n-1);
		REAL *p = x.data();
		REAL xx = -int((n-1)/2) * step + center;
		for(int i=0; i<n; i++) { *p = xx; p++; xx += step; }
		
		CrystVector_REAL y = reflProf.GetProfile(x,center,h,k,l);
		
		// calculate peak parameters
		params1[irefl] = CalcPeakParams (x,y);
		
  } // irefl
  
  const REAL eps = numeric_limits<REAL>::epsilon();
  
  bool accuracyReached = accur<=0.0;
  int ncalc = 0;
  
  // if calculation accuracy is presribed than multiple calculation must be done
  while (!accuracyReached && ncalc<=10) {
  	
		// increase the number of calculation points and data range :TODO: optimise
  	n = 2*(n-1)+1;
  	
  	// calculate again data for all diffractions
  	for(int irefl=0; irefl<this->GetNbRefl(); irefl++) {
  		
  		width[irefl] *= M_SQRT2;
  		
  		// save the old calculated values
  		params0[irefl] = params1[irefl];
  		
  		int h = (int)mH(irefl), k = (int)mK(irefl), l = (int)mL(irefl);
    
	  	// get peak center
	    REAL x0=mpParentPowderPattern->STOL2X(mSinThetaLambda(irefl));
	    REAL center = mpParentPowderPattern->X2XCorr(x0);
	    
	    // calculate profile
			CrystVector_REAL x(n);
			const REAL step = width[irefl]/(n-1);
			REAL *p = x.data();
			REAL xx = -int((n-1)/2) * step + center;
			for(int i=0; i<n; i++) { *p = xx; p++; xx += step; }
			
			CrystVector_REAL y = reflProf.GetProfile(x,center,h,k,l);
			
			// calculate peak parameters
			params1[irefl] = CalcPeakParams (x,y);
			
  	} // irefl
  
	  // check if prescribed accuracy was reached
	  accuracyReached = true;
	  
	  for(int irefl=0; irefl<this->GetNbRefl(); irefl++) {
	  	
	  	REAL val0, val1;
	  	
	  	val0 = params0[irefl].intensity;
	  	val1 = params1[irefl].intensity;
	  	if(abs(val1*accur)>eps && abs(val0-val1)>abs(accur*val1)) accuracyReached = false;
	  	
	  	val0 = params0[irefl].xmax;
	  	val1 = params1[irefl].xmax;
	  	if(abs(val1*accur)>eps && abs(val0-val1)>abs(accur*val1)) accuracyReached = false;
	  	
	  	val0 = params0[irefl].fwhm;
	  	val1 = params1[irefl].fwhm;
	  	if(abs(val1*accur)>eps && abs(val0-val1)>abs(accur*val1)) accuracyReached = false;
	  	
	  	val0 = params0[irefl].intensity/params0[irefl].ymax;
	  	val1 = params1[irefl].intensity/params1[irefl].ymax;
	  	if(abs(val1*accur)>eps && abs(val0-val1)>abs(accur*val1)) accuracyReached = false;
	  	
	  	val0 = params0[irefl].asym;
	  	val1 = params1[irefl].asym;
	  	if(abs(val1*accur)>eps && abs(val0-val1)>abs(accur*val1)) accuracyReached = false;
	  	
	  } // irefl
	  
	  ncalc++;
  	
  } // if (accur>=0.)
  
  if(!accuracyReached)
  	cerr << "Warning (PowderPatternDiffraction::PrintHKLInfo(...)): accuracy " << accur << " not reached!" << "\n";
  	 
  // print data
  for(int irefl=0; irefl<this->GetNbRefl(); irefl++) {
    // print (hkl)
    int h = (int)mH(irefl), k = (int)mK(irefl), l = (int)mL(irefl);
    s<<setw(8)<<h<<setw(3)<<k<<setw(3)<<l;
    
		// print peak parametrs
    s<<fixed<<setprecision(3)<<setw(12)<<params1[irefl].xmax*RAD2DEG; // :TODO: is here RAD2DEG correctc for not x-ray data?
    
    s<<setw(12)<<setprecision(2)<<scientific<<params1[irefl].intensity;
    
    s<<fixed<<setprecision(3)<<setw(12)<<params1[irefl].fwhm*RAD2DEG;
    
    s<<fixed<<setprecision(3)<<setw(12)<<params1[irefl].intensity/params1[irefl].ymax*RAD2DEG;

	  s<<fixed<<setprecision(2)<<setw(12)<<params1[irefl].asym;
	  
    s<<"\n";
  }
}

void PowderPatternDiffraction::CalcIntensityCorr () const
{
  VFN_DEBUG_MESSAGE("MStruct::PowderPatternDiffraction::CalcIntensityCorr()",
		    10)
  bool needRecalc=false;
  bool thetaModified=false;

  this->CalcSinThetaLambda();
  if(mClockIntensityCorr<mClockTheta) { needRecalc=true; thetaModified=true; }

  // HKL list not modified and peak position not moved significantly
  if(needRecalc==true && mClockIntensityCorr>mClockHKL) {
    REAL ds = 0.;
    for(int i=0;i<mTheta.numElements();i++) {
      REAL ds0 = fabs((mTheta(i)-mIntensityCorrTheta(i))/tan(mTheta(i)));
      if(ds<ds0) ds = ds0;
    }
    if(ds<0.04) { needRecalc=false; thetaModified=true; }
  }

  if(mClockIntensityCorr<mCorrAbsorption.GetClockCorr()) needRecalc=true;
  if(mClockIntensityCorr<mCorrTexture.GetClockCorr()) needRecalc=true;
  if(mClockIntensityCorr<mCorrHKLIntensity.GetClockCorr()) needRecalc=true;

  // HKLIntensityCorr params
  bool HKLIntensityCorrNeedRecalc=false;
  if(mClockIntensityCorr<mCorrHKLIntensity.GetClockMaster()) {
    needRecalc=true;
    HKLIntensityCorrNeedRecalc=true;
  }

  // TextureCorr params
  bool TextureCorrNeedRecalc=false;
  if(mClockIntensityCorr<mCorrTexture.GetClockMaster()) {
    needRecalc=true;
    TextureCorrNeedRecalc=true;
  }

  if(needRecalc==false) {

    if(this->GetRadiation().GetWavelengthType()==WAVELENGTH_TOF) {
      if(mClockIntensityCorr<mCorrTOF.GetClockCorr()) needRecalc=true;
    }
    else
    {
      if(mClockIntensityCorr<mCorrLorentz.GetClockCorr()) needRecalc=true;
      if(this->GetRadiation().GetRadiationType()==RAD_XRAY)
	if(mClockIntensityCorr<mCorrPolar.GetClockCorr()) needRecalc=true;
      if(mClockIntensityCorr<mCorrSlitAperture.GetClockCorr()) needRecalc=true;
    }
   
    if(mCorrTextureMarchDollase.GetNbPhase()>0)
      if(mClockIntensityCorr<mCorrTextureMarchDollase.GetClockCorr()) needRecalc=true;
   
    // only reason to recalculate correction is that theta positions
    // (maybe) were slightly changed but we ignore this matter
    if(needRecalc==false) {
      mClockIntensityCorr.Click();
      mClock2IntensityCorr.Click();
      return;
    }
  }
  
  if(needRecalc==true) mClockIntensityCorr.Reset(); // force recalc
  
  ObjCryst::PowderPatternDiffraction::CalcIntensityCorr();
  
  // if he didn't recalc and I don't want to recalc than we go away
  if(mClockIntensityCorr<mClock2IntensityCorr && needRecalc==false) return;

  // I should add my corrections (absorption, texture)
  const CrystVector_REAL *mpCorr[3];
  
  mpCorr[0]=&(mCorrAbsorption.GetCorr(thetaModified));
  mpCorr[1]=&(mCorrTexture.GetCorr(thetaModified || TextureCorrNeedRecalc));
  mpCorr[2]=&(mCorrHKLIntensity.GetCorr(thetaModified || HKLIntensityCorrNeedRecalc));
  
  mIntensityCorr *= *(mpCorr[0]);
  mIntensityCorr *= *(mpCorr[1]);
  mIntensityCorr *= *(mpCorr[2]);
	
  mIntensityCorrTheta = mTheta;
  mClock2IntensityCorr.Click();
  VFN_DEBUG_MESSAGE("PowderPatternDiffraction::CalcIntensityCorr():finished",10)

}

// ReflectionProfileComponent
ReflectionProfileComponent::ReflectionProfileComponent():
mpParentReflectionProfile(0)
{}

void ReflectionProfileComponent::SetParentReflectionProfile
                    (const ReflectionProfile& s) {
  //if(mpParentPowderPatternDiffraction!=0) 
  //    mClockMaster.RemoveChild(mpParentPowderPattern->GetIntegratedProfileLimitsClock());
   mpParentReflectionProfile = &s;
   //mClockMaster.AddChild(mpParentPowderPattern->GetIntegratedProfileLimitsClock());
   //mClockMaster.AddChild(mpParentPowderPattern->GetClockPowderPatternPar());
   //mClockMaster.AddChild(mpParentPowderPattern->GetClockPowderPatternXCorr());
   //mClockMaster.AddChild(mpParentPowderPattern->GetClockPowderPatternRadiation());
}

const ReflectionProfile& ReflectionProfileComponent::GetParentReflectionProfile()const {
  return *mpParentReflectionProfile;
}

bool ReflectionProfileComponent::IsRealSpaceType()const {
  return true;
}

bool ReflectionProfileComponent::IsAnisotropic()const {
  return false;
}

REAL ReflectionProfileComponent::GetPositionCorr(const REAL xcenter,
					     const REAL h, const REAL k, const REAL l)const
{
  return 0.;
}

// SizeBroadeningEffect
SizeBroadeningEffect::SizeBroadeningEffect():
mM(1000.), mSigma(0.3)
{
  InitParameters();
}

CrystVector_REAL SizeBroadeningEffect::GetProfile(const CrystVector_REAL &x,
						  const REAL xcenter,
						  const REAL h, const REAL k, const REAL l)
{
  int nbPoints = x.numElements(); 
  CrystVector_REAL profile(nbPoints);

  // calc Four. coefs
  // ref: G.Ribarik,T.Ungar,J.Gubicza,J.Appl.Cryst.(2001).34,669-676:MWP-fit
  // size effect

  const REAL *p1 = x.data();
  REAL *p2 = profile.data();
  for(int i=0;i<nbPoints;i++) {
    double L = *p1++; L = abs(L);
    if (L<1.e-4)
      *p2++ = 1.0;
    else {
      double ln = M_SQRT1_2*log(L/mM)/mSigma;
      double A = L*L*L*erfc(ln);
      A -= 3*exp(2*mSigma*mSigma)*L*mM*mM*erfc(ln-M_SQRT2*mSigma);
      A += 2*exp(4.5*mSigma*mSigma)*mM*mM*mM*erfc(ln-3*M_SQRT1_2*mSigma);
      A *= exp(-4.5*mSigma*mSigma)/(4*mM*mM*mM);
      *p2++ = (REAL) A;
    } 
  }

  if (bsavecalc && xcenter>=xcenterlimits[0] && xcenter<=xcenterlimits[1]) {
    ofstream F("profileAS.dat");
    F<<"# m="<<mM<<",sigma="<<mSigma<<",xcenter="<<xcenter*RAD2DEG<<endl;
    for(int i=0;i<nbPoints;i++)
      F<<setw(18)<<x(i)<<setw(18)<<profile(i)<<endl;
    F.close();
  }

  return profile;
}

REAL SizeBroadeningEffect::GetApproxFWHM(const REAL xcenter,
					 const REAL h, const REAL k, const REAL l)const
{
  const Radiation &r = GetParentReflectionProfile().
    GetParentPowderPatternDiffraction().GetRadiation();
  
  return 1.3*r.GetWavelength()(0)/mM/cos(0.5*xcenter);
}

bool SizeBroadeningEffect::IsRealSpaceType()const {
  return true;
}

void SizeBroadeningEffect::SetProfilePar(const REAL m, const REAL sigma)
{
  mM = 10.*m;
  mSigma = sigma;
  mClockMaster.Click();  
}

void SizeBroadeningEffect::InitParameters()
{
  {
    RefinablePar tmp("M", &mM, 5., 3.e3,
                     gpRefParTypeScattDataProfileWidth,
                     REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false,0.1);
    tmp.AssignClock(mClockMaster);
    tmp.SetDerivStep(1.);
    this->AddPar(tmp);
  }
  {
    RefinablePar tmp("Sigma", &mSigma, 0.05, 1.5,
                     gpRefParTypeScattDataProfileWidth,
                     REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false,1.);
    tmp.AssignClock(mClockMaster);
    tmp.SetDerivStep(0.03);
    this->AddPar(tmp);
  }
}

////////////////////////////////////////////////////////////////////////
//
//    SizeDistribPowderPatternDiffraction
//
////////////////////////////////////////////////////////////////////////

SizeDistribPowderPatternDiffraction::SizeDistribPowderPatternDiffraction ()
  : mpSizeDistribReflProf(NULL)
{}

SizeDistribPowderPatternDiffraction::SizeDistribPowderPatternDiffraction (const SizeDistribPowderPatternDiffraction & old)
  : MStruct::PowderPatternDiffraction(old), mpSizeDistribReflProf(NULL)
{}

SizeDistribPowderPatternDiffraction * SizeDistribPowderPatternDiffraction::CreateCopy () const
{
  return new SizeDistribPowderPatternDiffraction(*this);
}

const string & SizeDistribPowderPatternDiffraction::GetClassName () const
{
  static const string className = "MStruct::SizeDistribPowderPatternDiffraction";
  return className;
}

const CrystVector_REAL & SizeDistribPowderPatternDiffraction::GetLSQDeriv (const unsigned int n, ObjCryst::RefinablePar & par)
{
  
  mLSQDeriv.resize( this->GetParentPowderPattern().GetPowderPatternX().numElements() );
  mLSQDeriv = 0.;

  string str = par.GetName();
  string::size_type loc = str.find("pD_", 0); // pD_ProfileName_nb
 
  if( loc != string::npos && mpSizeDistribReflProf!=NULL ) {
    VFN_DEBUG_MESSAGE("MStruct::SizeDistribPowderPatternDiffraction::GetLSQDeriv: "<<
		      "Derivative with respect to a SizeDistrib (histogram) value detected: "<< str,11)
    
    // Confirm from the parameter address that this is really the SizeDistrib parameter
    // and identify the (histogram) value related.
    const CrystVector_REAL & dist = mpSizeDistribReflProf->GetDistribution();
    long indPar = -1;
    if( par.GetPointer()>=dist.data() && par.GetPointer()<=(dist.data()+dist.numElements()-1) )
      indPar = long(par.GetPointer()-dist.data());
    if( indPar<0 ) {
      cerr << "Error: MStruct::SizeDistribPowderPatternDiffraction::GetLSQDeriv(...) asked to calculate "
	   << "derivative with respect to variable name=" << par.GetName()
	   << ". Can not find the related SizeDistrib parameter.\n";
      throw ObjCrystException("MStruct::SizeDistribPowderPatternDiffraction::GetLSQDeriv(...): \
             Wrong input argument.");
    }
    // Check clocks
    if( mvClockLSQDerivCalculated[indPar]>mOtherParamsClock )
      // the stored value can be used
      mLSQDeriv = mvLSQDerivSizeDistrib[indPar];
    else {
      /*cout << par.GetName() << " derivative will be recalculated. " << '\n';*/
      // the derivative has to be calculated
      
      // reference to the SizeDistribution
      const CrystVector_REAL &dist = mpSizeDistribReflProf->GetDistribution();
      // copy distribution (histogram) values
      const CrystVector_REAL dist_current(dist);
      // set distribution to P(*) = 0; P(indPar) = 1
      for(int i=0; i<dist.numElements(); i++)
	mpSizeDistribReflProf->GetPar( dist.data()+i ).SetValue(0.);
      par.SetValue(1.);
      // calculate pattern
      mLSQDeriv = this->GetPowderPatternCalc();
      // restore the distribution (histogram) values
      for(int i=0; i<dist.numElements(); i++)
	mpSizeDistribReflProf->GetPar( dist.data()+i ).SetValue(dist_current(i));

      // save the calculated derivative
      mvLSQDerivSizeDistrib[indPar] = mLSQDeriv;
      mvClockLSQDerivCalculated[indPar].Click();
    }

    // finish the calculation and normalise the derivative
    mLSQDeriv -= this->GetPowderPatternCalc();;
    CrystVector_REAL t = mSizeDistribA0;
    t *= mpSizeDistribReflProf->GetDistribution();
    mLSQDeriv *= mSizeDistribA0(indPar)/t.sum();
    
    return mLSQDeriv;
  }/* else {
    cerr << "Error: MStruct::SizeDistribPowderPatternDiffraction::GetLSQDeriv(...) asked to calculate "
	 << "derivative with respect to variable name=" << par.GetName()
	 << ". This derivative not supported.\n";
    throw ObjCrystException("MStruct::SizeDistribPowderPatternDiffraction::GetLSQDeriv(...): \
             Wrong input argument.");
	     } I_h_k_l ??? */

  return MStruct::PowderPatternDiffraction::GetLSQDeriv(n,par);
}

void SizeDistribPowderPatternDiffraction::Prepare ()
{
  MStruct::PowderPatternDiffraction::Prepare();
  mpSizeDistribReflProf = NULL;

  // Try to find a SizeDistribBroadeningEffect in the list of ReflectionProfileComponents
  if(this->GetProfile().GetClassName()==string("MStruct::ReflectionProfile")) {
    MStruct::ReflectionProfile & reflProf = dynamic_cast<MStruct::ReflectionProfile &>(this->GetProfile());
    vector<long> ind;
    for(long i=0; i<reflProf.GeReflectionProfileComponentNb(); i++)
      if(reflProf.GetReflectionProfileComponent(i).GetClassName()==string("MStruct::SizeDistribBroadeningEffect"))
	ind.push_back(i);
    if(ind.size()>0) {
      mpSizeDistribReflProf =
	dynamic_cast<MStruct::SizeDistribBroadeningEffect*>(&reflProf.GetReflectionProfileComponent(ind[0]));
      VFN_DEBUG_MESSAGE("MStruct::SizeDistribPowderPatternDiffraction::Prepare: "<<
			"SizeDistribBroadeningEffect object found whose clacluation will be optimised : "
			<< mpSizeDistribReflProf->GetName(),11)
    } // if(ind.size()>0)
    if(ind.size()>1) {
      cout << "Warning: (MStruct::SizeDistribPowderPatternDiffraction::Prepare) "
	   << "Multiple ("<<ind.size()<<") MStruct::SizeDistribBroadeningEffect(s) found "
	   << "in the (" << this->GetName() << ") object. Calculation for only the first one ("
	   << mpSizeDistribReflProf->GetName() << ") will be optimised.";
    } // if(ind.size()>1)
  } // if(...)

  // Create a Clock that triggers chages in all possibly related refinable parameters except
  // that of the SizeDistrib model
  {
    // Remove all Child Clocks
    // TODO:: Do this more clearly. Calling destructors explicitly looks ugly.
    mOtherParamsClock.~RefinableObjClock(); mOtherParamsClock = ObjCryst::RefinableObjClock();// we need new clocks
  }

  // Connect "mOtherParamsClock" with appropriate parametrs.
  if(mpSizeDistribReflProf!=NULL) {
    // This is a little bit tricky to save code here
    // 1) build an ObjRegistry. 2) go through the ObjRegistry and register Children Clocks.

    // It is assumed that the SizeDistrib values are stored continuosly in the memory.
    // Hence the first and last element adresses are used later for their identification.
    const CrystVector_REAL & dist = mpSizeDistribReflProf->GetDistribution();
    const REAL *pSizeDistribBegin = dist.data();
    const REAL *pSizeDistribEnd = dist.data()+(dist.numElements()-1);

    ObjCryst::ObjRegistry< ObjCryst::RefinableObj > tmpReg;
    ObjCryst::RefObjRegisterRecursive ((ObjCryst::RefinableObj &)*this, tmpReg);
    for(long iObj=0; iObj<tmpReg.GetNb(); iObj++) {
      ObjCryst::RefinableObj & refObj = tmpReg.GetObj(iObj);
      for(long iPar=0; iPar<refObj.GetNbPar(); iPar++) {
	ObjCryst::RefinablePar & par = refObj.GetPar(iPar);
	if(par.GetPointer()<pSizeDistribBegin || par.GetPointer()>pSizeDistribEnd)
	  mOtherParamsClock.AddChild( par.GetClock() ); // this is not a SizeDistrib parameter
      } // iPar
    } // iObj 
  }

  // Crete Clocks when the LSQDerivatives were calculated last time
  if(mpSizeDistribReflProf!=NULL) {
    mvClockLSQDerivCalculated.clear(); // remove old clocks
    // Create new clocks and Reset them 
    mvClockLSQDerivCalculated.resize( mpSizeDistribReflProf->GetDistribution().numElements() );
    for_each( mvClockLSQDerivCalculated.begin(), mvClockLSQDerivCalculated.end(),
	      mem_fun_ref(&ObjCryst::RefinableObjClock::Reset) );
  }
  
  // Prepare mpSizeDistribReflProf to store calculated derivatives
  if(mpSizeDistribReflProf!=NULL)
    mvLSQDerivSizeDistrib.resize( mpSizeDistribReflProf->GetDistribution().numElements() );

  // Prepare the A(0) values used later for derivative calculation/normalisation
  if(mpSizeDistribReflProf!=NULL) {
    mSizeDistribA0.resize( mpSizeDistribReflProf->GetDistribution().numElements() );
    CrystVector_REAL D1, D2;
    mpSizeDistribReflProf->GetDistributionBins(D1, D2);
    for(int i=0; i<D1.numElements(); i++)
      mSizeDistribA0(i) = (pow(D1(i),2)+pow(D2(i),2))/pow(D1(i)+D2(i),2);
  }
}

//----------------------------------
//----------------------------------
//----------------------------------

// SizeDistribBroadeningEffect
SizeDistribBroadeningEffect::SizeDistribBroadeningEffect()
  :mD1(0), mD2(0), mDistrib(0), mLSQAlpha(0.), mLSQConstraintScale(0.),
   mDmax(0.), mDistIntegral(0.), mVolumeDistIntegral(0.), mCurvIntStep(0.),
   mBeginEndOptimizationCalled(0)
{
  //SetLSQRegularizationOption(LSQRegOpt_None); // Deprecated
  mLSQConstraintScale = ObjCryst::RefinableObj::mDefaultLSQConstraintScale;
}

const string& SizeDistribBroadeningEffect::GetClassName()const
{
  static const string className = "MStruct::SizeDistribBroadeningEffect";
  return className;
}

CrystVector_REAL SizeDistribBroadeningEffect::GetProfile(const CrystVector_REAL &x,
						  const REAL xcenter,
						  const REAL h, const REAL k, const REAL l)
{
  int nbPoints = x.numElements(); 
  CrystVector_REAL profile(nbPoints);

  // size effect

	if(mDistrib.numElements()==0) {
		profile = 1.;
	} else {
	  const REAL *p1 = x.data();
	  REAL *p2 = profile.data();
	  double t = 0;
	  for(int j=0; j<mDistrib.numElements(); j++)
	  	//t += (mD1(j)+mD2(j))*(pow(mD1(j),2)+pow(mD2(j),2)) * mDistrib(j)/pow(mD1(j)+mD2(j),3);
	  	t += (pow(mD1(j),2)+pow(mD2(j),2))/pow(mD1(j)+mD2(j),2) * mDistrib(j);
	  for(int i=0;i<nbPoints;i++) {
	    double L = *p1++; L = abs(L);
	    double A = 0.;
	    if(L<1.e-4)
				A = 1.0;
			else
	    	for(int j=0; j<mDistrib.numElements(); j++) {
	    		/*double D = (L<mD1(j)) ? mD1(j) : L;
						if(L<mD2(j))
	      			A += 2.0*(0.5*(pow(mD2(j),4)-pow(D,4))-L*(pow(mD2(j),3)-pow(D,3))+pow(L,3)*(mD2(j)-D))
	      		     / (pow(mD1(j),3) + pow(mD1(j),2)*mD2(j) + mD1(j)*pow(mD2(j),2) + pow(mD2(j),3))
	      		     * mDistrib(j) / t;
	      	*/ 
	      	/*if(L<mD1(j))
	      		A += ( 1. - (2.*L*(pow(mD1(j),2)+mD1(j)*mD2(j)+pow(mD2(j),2))-2*pow(L,3))/(mD1(j)+mD2(j))
	      		      /  (pow(mD1(j),2)+pow(mD2(j),2)) ) * mDistrib(j)/t;
	        else if(L<mD2(j))
	        	A += pow(mD2(j)-L,3)*(mD2(j)+L)/(mD2(j)-mD1(j)) *mDistrib(j)/t /(mD1(j)+mD2(j))/(pow(mD1(j),2)+pow(mD2(j),2));
	        */
	        /*if(L<mD1(j))
	      		A += ( (mD1(j)+mD2(j))*(pow(mD1(j),2)+pow(mD2(j),2)) 
	      		       - 2.*L*(pow(mD1(j),2)+mD1(j)*mD2(j)+pow(mD2(j),2)) + 2*pow(L,3) ) * mDistrib(j)/pow(mD1(j)+mD2(j),3)/t;
	      	else if(L<mD2(j))
	        	A += pow(mD2(j)-L,3)*(mD2(j)+L)/(mD2(j)-mD1(j)) *mDistrib(j)/pow(mD1(j)+mD2(j),3)/t;
	       	*/
	       	if(L<mD1(j))
	      		A += ( (pow(mD1(j),2)+pow(mD2(j),2))/pow(mD1(j)+mD2(j),2) 
	      		       - (2.*L*(pow(mD1(j),2)+mD1(j)*mD2(j)+pow(mD2(j),2)) - 2*pow(L,3))/pow(mD1(j)+mD2(j),3) ) * mDistrib(j)/t;
	      	else if(L<mD2(j))
	        	A += pow(mD2(j)-L,3)*(mD2(j)+L)/(mD2(j)-mD1(j))/pow(mD1(j)+mD2(j),3) *mDistrib(j)/t;
	      }
	    *p2++ = (REAL) A;
	  } 
	}
	
  if (bsavecalc && xcenter>=xcenterlimits[0] && xcenter<=xcenterlimits[1]) {
    ofstream F("profileASD.dat");
    F<<"# xcenter="<<xcenter*RAD2DEG<<", ApproxFWHM="<<GetApproxFWHM(xcenter,h,k,l)*RAD2DEG<<endl;
    for(int i=0; i<mDistrib.numElements(); i++)
    	F<<"# "<<setw(8)<<mD1(i)/10.<<setw(8)<<mD2(i)/10.<<setw(8)<<mDistrib(i)<<"\n";
    for(int i=0;i<nbPoints;i++)
      F<<setw(18)<<x(i)<<setw(18)<<profile(i)<<endl;
    F.close();
  }

  return profile;
}

REAL SizeDistribBroadeningEffect::GetApproxFWHM(const REAL xcenter,
					 const REAL h, const REAL k, const REAL l)const
{
  const Radiation &r = GetParentReflectionProfile().
    GetParentPowderPatternDiffraction().GetRadiation();
  
  REAL t1 = 0., t2 = 0.;
  for(int i=0; i<mDistrib.numElements(); i++) {
  	t1 += 0.25*pow(mD1(i)+mD2(i),2) * (mD2(i)-mD1(i)) * mDistrib(i);
  	t2 += (mD2(i)-mD1(i))*mDistrib(i)*pow(mD1(i)+mD2(i),3)/8.;
  }
  
  return (mDistrib.numElements()>0) ? 1.3*r.GetWavelength()(0)/cos(0.5*xcenter)*t1/t2 : 0.0;
}

bool SizeDistribBroadeningEffect::IsRealSpaceType()const {
  return true;
}

void SizeDistribBroadeningEffect::SetDistribution(const CrystVector_REAL d1, const CrystVector_REAL d2,
						  const CrystVector_REAL distrib, const CrystVector_int fixed)
{
  if(d1.numElements()!=d2.numElements() ||
     d2.numElements()!=distrib.numElements() ||
     distrib.numElements()!=fixed.numElements()) {
    cerr << "Error: Vectors defining crystallite size distribution must be of same size.\n";
    cerr << "Size distribution not set!" << endl;
    return;
  }
	
  // remove refinable parameters of the old distribution
  for(int i=0; i<mDistrib.numElements(); i++) {
    long n = this->FindPar(&mDistrib(i));
    if (n>=0 && n<this->GetNbPar()) this->RemovePar(&(this->GetPar(n)));
  }
	
  // set the new distribution
  mD1 = d1;
  mD2 = d2;
  mDistrib = distrib;
  
  // find Dmax and min(widthD)
  REAL Dmax = -1., Dw = 1e6;
  for(int i=0; i<distrib.numElements(); i++) {
    if ( mD2(i)>Dmax ) Dmax = mD2(i);
    if ( (mD2(i)-mD1(i))<Dw ) Dw = mD2(i)-mD1(i);
  }
  mDmax = Dmax;
  mCurvIntStep = Dw/2;

  // set new refinable parameters
  for(int i=0; i<distrib.numElements(); i++) {
    // generate parameter name
    string name;
    {
      stringstream s;
      s << "pD_" << this->GetName() << "_" << i;
      name = s.str();
    }
    // add new refinable paramer
    {
      RefinablePar tmp(name, &mDistrib(i), 0., 1.e7,
		       gpRefParTypeScattDataProfileSizeDistrib,
		       REFPAR_DERIV_STEP_ABSOLUTE,true,fixed(i),true,false,1.0);
      tmp.AssignClock(mClockMaster);
      tmp.SetDerivStep(0.01);
      this->AddPar(tmp);
      this->GetPar(name).SetIsFixed(fixed(i));
      this->GetPar(name).Print();
    }  
  }
  
  // rebuild the list of LSQ Regularization Operators
  this->RebuildLSQRegOpList();
  
  // reset clocks
  mClockDistIntegralCalc.Reset();
}

const CrystVector_REAL & SizeDistribBroadeningEffect::GetDistribution () const
{
  return mDistrib;
}

void SizeDistribBroadeningEffect::GetDistributionBins (CrystVector_REAL & D1, CrystVector_REAL & D2) const
{
  D1 = mD1;
  D2 = mD2;
}

REAL SizeDistribBroadeningEffect::GetDistribIntegral() const
{
  this->CalcDistIntegral();
  return mDistIntegral;
}

REAL SizeDistribBroadeningEffect::GetVolumeDistribIntegral() const
{
  this->CalcDistIntegral();
  return mVolumeDistIntegral;
}

void SizeDistribBroadeningEffect::CalcDistIntegral() const
{
  // distribution is connected with the master clock
  if(mClockDistIntegralCalc<this->GetClockMaster()) {
    // recalculate integrals
    mDistIntegral = 0.;
    mVolumeDistIntegral = 0.;

    for(int i=0; i<mDistrib.numElements(); i++) {
      const REAL w = mD2(i)-mD1(i);
      const REAL d = 0.5*(mD1(i)+mD2(i));
      mDistIntegral += mDistrib(i)/pow(d,3); //*w
      mVolumeDistIntegral += mDistrib(i); //*w
    }
    mClockDistIntegralCalc.Click();
  }
}

void SizeDistribBroadeningEffect::ReadDistributionFromFile(
						  const char* filename)
{
  CrystVector_REAL vD1(10), vD2(10), vDistrib(10);
  CrystVector_int vFixed(10);
	
  ifstream f(filename);
  string line;
  int n = 0;
  while (getline(f,line)) {
    if (line.empty()) break;
    if (line.at(0)=='#') continue;
    int fixed = 1;
    REAL d1, d2, v;
    istringstream s(line);
    s >> d1 >> d2 >> v >> fixed; // >> sigma;
    if (s.fail()==true) break;
    if(n>=vD1.numElements()) {
      vD1.resizeAndPreserve(n+10);
      vD2.resizeAndPreserve(n+10);
      vDistrib.resizeAndPreserve(n+10);
      vFixed.resizeAndPreserve(n+10);
    }
    vD1(n) = d1*10.; vD2(n) = d2*10.; vDistrib(n) = v; vFixed(n) = fixed;
    n++;
  }
  f.close();
  vD1.resizeAndPreserve(n);
  vD2.resizeAndPreserve(n);
  vDistrib.resizeAndPreserve(n);
  vFixed.resizeAndPreserve(n);
  
  this->SetDistribution(vD1,vD2,vDistrib,vFixed);
}

void SizeDistribBroadeningEffect::WriteDistributionToFile(const char* filename) const
{
  // copy content of the file
  ostringstream s;
  {
    char buffer[1024];
    ifstream f(filename);
    while(f.getline(buffer,1022)) s << buffer << '\n';
    f.close();
  }

  // write current distribution into the file
  ofstream f(filename);
  
  // write header
  f << "#     D1(nm)" << setw(10) << "D2(nm)" << setw(12) << "distrib" << setw(8) << "fixed" << "\n";
  f << showpoint << fixed;
  const REAL *p = mDistrib.data();
  for(int i=0; i<mDistrib.numElements(); i++) {
    const RefinablePar &par = this->GetPar(p+i);
    f << setprecision(2) << setw(12) << mD1(i)/10. << setw(10) << mD2(i)/10.;
    f << setprecision(3) << setw(12) << mDistrib(i) << setw(8) << (int) par.IsFixed();
    f << setprecision(3) << setw(12) << par.GetHumanSigma();
    f << scientific << setprecision(2) << setw(12) << mDistrib(i)/pow(0.5*(mD2(i)+mD1(i))/10.,3) << fixed;
    f << "\n";
  }
  
  // save original file content at the end of the file
  f<<'\n'<<s.str();

  f.close();
}

void SizeDistribBroadeningEffect::SetLSQConstraintScale(const REAL scale)
{
  mLSQConstraintScale = scale;
}

void SizeDistribBroadeningEffect::BuildDistribution(const REAL Dmin, const REAL Dmax,
						    const int NbIntervals, const string spacing)
{
  CrystVector_REAL D1(NbIntervals);
  CrystVector_REAL D2(NbIntervals);
  CrystVector_REAL  G(NbIntervals); // f(D)*D^3*(D2-D1)
  CrystVector_int Fixed(NbIntervals);
  
  Fixed = 0;
  
  if( (Dmin<=1e-7) || (Dmax<Dmin) || (NbIntervals<1) ) {
    cerr << "< MStruct::SizeDistribBroadeningEffect::GenerateDistribution(...)\n";
    cerr << "\t"<<" Dmin="<<Dmin<<", Dmax="<<Dmax<<", NbIntervals="<<NbIntervals<<"\n";
    cerr << "\t Wrong input argument. >"<<endl;
    throw ObjCrystException("MStruct::SizeDistribBroadeningEffect::GenerateDistribution(...): \
             Wrong input argument.");
  }
  
  // genarete x-axis values
  if ( spacing == string("linear") ) {
    // linear
    for(int i=0; i<NbIntervals; i++) {
      D1(i) = Dmin + i     * (Dmax-Dmin)/NbIntervals;
      D2(i) = Dmin + (i+1) * (Dmax-Dmin)/NbIntervals;
    }

  } else if ( spacing == string("log") ) {
    // logarithmic
    for(int i=0; i<NbIntervals; i++) {
      D1(i) = exp( log(Dmin) + i     * ( log(Dmax/Dmin)/NbIntervals ) );
      D2(i) = exp( log(Dmin) + (i+1) * ( log(Dmax/Dmin)/NbIntervals ) );
    }
  } else if ( spacing == string("sqrt") ) {
    // square root
    for(int i=0; i<NbIntervals; i++) {
      D1(i) = pow( sqrt(Dmin) + i     * ( (sqrt(Dmax)-sqrt(Dmin))/NbIntervals ) ,2);
      D2(i) = pow( sqrt(Dmin) + (i+1) * ( (sqrt(Dmax)-sqrt(Dmin))/NbIntervals ) ,2);
    }
  } else {
    // unknown option
    cerr << "< MStruct::SizeDistribBroadeningEffect::GenerateDistribution(...)\n";
    cerr << "\t"<<"spacing="<<spacing<<" \t"<<"Unknown x spacing type. >"<<endl;
    throw ObjCrystException("MStruct::SizeDistribBroadeningEffect::GenerateDistribution(...): \
             Wrong input argument.");
  }
  
  REAL t = 0.;
  for(int i=0; i<NbIntervals; i++) {
    REAL D = 0.5*(D1(i)+D2(i));
    //G(i) = (1.-cos(2.*M_PI*(D-D1(0))/(D2(NbIntervals-1)-D1(0)))) / (D2(NbIntervals-1)-D1(0)) * pow(D,3);
    G(i) = (1.-cos(2.*M_PI*(D-D1(0))/(D2(NbIntervals-1)-D1(0)))) * pow(D,3) / (D2(i)-D1(i));
    t += G(i);
  }
  G *= 1.e5/t;
  
  this->SetDistribution(D1,D2,G,Fixed);
}

unsigned int SizeDistribBroadeningEffect::GetNbLSQConstraints() const
{
  return 1;
}

void SizeDistribBroadeningEffect::GetLSQConstraint(const unsigned int n,
						   std::vector< const ObjCryst::RefinablePar* > &parList,
						   CrystVector_REAL &coef) const
{
  parList.clear();
  coef.resize(mDistrib.numElements());
  
  for(int i=0; i<mDistrib.numElements(); i++) {
    parList.push_back(&(this->GetPar(mDistrib.data()+i)));
    coef(i) = 1.; // abs(mD2(i)-mD1(i));
  }

  //coef *= mLSQConstraintScale/coef.max();
}

unsigned int SizeDistribBroadeningEffect::GetNbLSQRegularizationOperator(const unsigned int LSQfunc) const
{
  return mLSQRegOpList.size();
}

const LSQRegularizationOperator &
            SizeDistribBroadeningEffect::GetLSQRegularizationOperator(const unsigned int nOp,
								      const unsigned int LSQfunc) const
{
  if( !(nOp<mLSQRegOpList.size()) ) {
    cerr<<"MStruct::SizeDistribBroadeningEffect::GetLSQRegularizationOperator(nOp="<<nOp;
    cerr<<", LSQfunc="<<LSQfunc<<"): Only "<<mLSQRegOpList.size()<<" operator(s) aviable.\n";
    throw ObjCrystException("MStruct::SizeDistribBroadeningEffect::GetLSQRegularizationOperator(...): \
Wrong input argument.");
  }
  
  // Normalise the LSQ regularization matrix

  CrystMatrix_REAL H = mLSQRegOpList[nOp].GetRegularizationOperatorMatrix();
  REAL scale = 1.;

  switch ( mLSQRegTypeList[nOp] ) {
    
  case LSQRegOpt_DistribDeriv: // Integral of the Arithmetic Distribution Derivative
    scale = pow(1./(mDmax*this->GetDistribIntegral()),2)/mDmax;
    break;

  case LSQRegOpt_VolumeDistribDeriv: // Integral of the Volume weighted Distribution Derivative
    scale = mDmax*pow(mDmax/this->GetVolumeDistribIntegral(),2);
    break;

  } // switch

  // rescale matrix
  H *= scale/mLSQRegOpNorm[nOp];
  // save the current norm
  mLSQRegOpNorm[nOp] = scale;

  mLSQRegOpList[nOp].SetRegularizationOperatorMatrix(H);

  return mLSQRegOpList[nOp];
}

void SizeDistribBroadeningEffect::GlobalOptRandomMove(const REAL mutationAmplitude,
							    const RefParType *type)
{
  VFN_DEBUG_MESSAGE("MStruct::SizeDistribBroadeningEffect::GlobalOptRandomMove:Begin",11);

  GenerateRandomDistrib();

  VFN_DEBUG_EXIT("MStruct::SizeDistribBroadeningEffect::GlobalOptRandomMove:End",11);
}

LSQRegularizationOperator SizeDistribBroadeningEffect::CreateLSQRegOpVolumeDistribDeriv(const REAL weight) const
{
  if( mDistrib.numElements()==0 )
    return EmptyLSQRegularizationOperatorObj;

  // prepare regularization operator
  const int nbBins = mDistrib.numElements();
  CrystMatrix_REAL H(nbBins,nbBins);

  /*  | 2  -1   0   0   0    0 |
      |-1   2  -1   0   0    0 |
      | 0  -1   2  -1   0    0 |
      | 0   0  -1   2  -1    0 |
      | 0   0   0  -1   2   -1 |
      | 0   0   0   0  -1    2 | */

  // points at bins centers, derivates are calculated between them
  CrystVector_REAL D(nbBins);
  for(int i=0; i<nbBins; i++)
    D(i) = 0.5*(mD1(i)+mD2(i));

  H = 0.;
  //const REAL scale = mDmax*pow(mDmax/this->GetVolumeDistribIntegral(),2);
  const REAL scale = 1.; // scaled dynamically in the GetLSQRegularizationOperator(...)

  // diagonal
  H(0,0) = scale*(1./(D(1)-D(0))+1./(D(0)-0)) /pow(mD2(0)-mD1(0),2);
  for(int i=1; i<nbBins-1; i++)
    H(i,i) = scale*(1./(D(i)-D(i-1))+1./(D(i+1)-D(i))) /pow(mD2(i)-mD1(i),2);
  H(nbBins-1,nbBins-1) = scale*(1./(D(nbBins-1)-D(nbBins-2))+
				1./(mDmax-D(nbBins-1))) /pow(mD2(nbBins-1)-mD1(nbBins-1),2);
  // upper diagonal
  for(int i=0; i<nbBins-1; i++)
    H(i,i+1) = -1.*scale/(D(i+1)-D(i)) /(mD2(i)-mD1(i))/(mD2(i+1)-mD1(i+1));
  // lower diagonal
  for(int i=0; i<nbBins-1; i++)
    H(i+1,i) = -1.*scale/(D(i+1)-D(i)) /(mD2(i)-mD1(i))/(mD2(i+1)-mD1(i+1));

  // create parameters list
  std::vector< const RefinablePar * > parList;
  for(int i=0; i<nbBins; i++) {
    parList.push_back( &(this->GetPar(mDistrib.data()+i)) );
  }

  // add new LSQRegularizationOperator
  LSQRegularizationOperator regOp("SizeDistrib "+this->GetName()+" LSQRegOp - volume distrib. derivative");
  regOp.SetRegularizationOperatorMatrix(H);
  regOp.SetRegularizationOperatorWeight(weight);
  regOp.SetParamList(parList);
  
  return regOp;
}

LSQRegularizationOperator SizeDistribBroadeningEffect::CreateLSQRegOpDistribDeriv(const REAL weight) const
{
  if( mDistrib.numElements()==0 )
    return EmptyLSQRegularizationOperatorObj;

  // prepare regularization operator
  const int nbBins = mDistrib.numElements();
  CrystMatrix_REAL H(nbBins,nbBins);

  /*  | 2  -1   0   0   0    0 |
      |-1   2  -1   0   0    0 |
      | 0  -1   2  -1   0    0 |
      | 0   0  -1   2  -1    0 |
      | 0   0   0  -1   2   -1 |
      | 0   0   0   0  -1    2 | */

  // points at bins centers, derivates are calculated between them
  CrystVector_REAL D(nbBins);
  for(int i=0; i<nbBins; i++)
    D(i) = 0.5*(mD1(i)+mD2(i));

  H = 0.;
  //const REAL scale = pow(1./(mDmax*this->GetDistribIntegral()),2)/mDmax;
  const REAL scale = 1.; // scaled dynamically in the GetLSQRegularizationOperator(...)

  // diagonal
  H(0,0) = scale*(1./(D(1)-D(0))+1./(D(0)-0)) /pow(D(0)/mDmax,6) /pow(mD2(0)-mD1(0),2);
  for(int i=1; i<nbBins-1; i++)
    H(i,i) = scale*(1./(D(i)-D(i-1))+1./(D(i+1)-D(i))) /pow(D(i)/mDmax,6) /pow(mD2(i)-mD1(i),2);
  H(nbBins-1,nbBins-1) = scale*(1./(D(nbBins-1)-D(nbBins-2))+
				1./(mDmax-D(nbBins-1))) /pow(D(nbBins-1)/mDmax,6) /pow(mD2(nbBins-1)-mD1(nbBins-1),2);
  // upper diagonal
  for(int i=0; i<nbBins-1; i++)
    H(i,i+1) = -1.*scale/(D(i+1)-D(i)) /pow(D(i)/mDmax,3)/pow(D(i+1)/mDmax,3) /(mD2(i)-mD1(i))/(mD2(i+1)-mD1(i+1));
  // lower diagonal
  for(int i=0; i<nbBins-1; i++)
    H(i+1,i) = -1.*scale/(D(i+1)-D(i)) /pow(D(i+1)/mDmax,3)/pow(D(i)/mDmax,3) /(mD2(i)-mD1(i))/(mD2(i+1)-mD1(i+1));

  // create parameters list
  std::vector< const RefinablePar * > parList;
  for(int i=0; i<nbBins; i++)
    parList.push_back( &(this->GetPar(mDistrib.data()+i)) );

  // add new LSQRegularizationOperator
  LSQRegularizationOperator regOp("SizeDistrib "+this->GetName()+" LSQRegOp - distrib. derivative");
  regOp.SetRegularizationOperatorMatrix(H);
  regOp.SetRegularizationOperatorWeight(weight);
  regOp.SetParamList(parList);
  
  return regOp;
}

void SizeDistribBroadeningEffect::GenerateRandomDistrib()
{
  VFN_DEBUG_MESSAGE("MStruct::SizeDistribBroadeningEffect::GenerateRandomDistrib:Begin",11);
  
  // for each bin select a random value from uniform distribution in the interval [0,1)
  REAL *p = mDistrib.data();
  //for(int i=0; i<<mDistrub.numElements(); i++) { *p = real_rand()/pow(i+0.5,3); p++; } // v2
  for(int i=0; i<mDistrib.numElements(); i++) {
    if (!(this->GetPar(p).IsFixed())) *p = real_rand()*pow(0.5*(mD1(i)+mD2(i))/10.,3); // v1
    p++;
  }

  // normalise the distrution so sum(i=0..Nbins-1, P(i)*BinWidth)==1 holds
  REAL sum = 0.;
  for(int i=0; i<mDistrib.numElements(); i++)
    sum += mDistrib(i)/pow(0.5*(mD1(i)+mD2(i))/10.,3) * fabs(mD2(i)-mD1(i));

  if( sum>1.e-7 ) mDistrib *= 1./sum; // NOTE: Here also fixed params. are affected

  // Distribution has changed, reflection profile should be recomputed
  //mClockReflProfCalc.Reset();

  // Object has changed
  mClockMaster.Click();

  VFN_DEBUG_EXIT("MStruct::SizeDistribBroadeningEffect::GenerateRandomDistrib:End",11);
}


void SizeDistribBroadeningEffect::AddLSQRegularizationMethod(const int option, const REAL weight)
{
  VFN_DEBUG_MESSAGE("MStruct::SizeDistribBroadeningEffect::AddLSQRegularizationMethod("
		    <<option<<","<<weight<<"):Begin",11);
  
  // add new in the list
  mLSQRegTypeList.push_back(option);
  mLSQRegWeightList.push_back(weight);
  
  // rebuild the RegOpList
  this->RebuildLSQRegOpList();

  VFN_DEBUG_EXIT("MStruct::SizeDistribBroadeningEffect::AddLSQRegularizationMethod:End",11);
}

REAL SizeDistribBroadeningEffect::CalcCurvInt() const
{
  VFN_DEBUG_MESSAGE("MStruct::SizeDistribBroadeningEffect::CalcCurvInt():Begin",11);

  // calcualte unsigned curvature integral - approximate data with cubic spline

  // create cubic spline approximation of the size distribution
  // TODO:: be carefull ! - noncontinous distributions with "holes" are not handled correctly

  const int nb = mD1.numElements();
  CrystVector_REAL x(nb+2);
  CrystVector_REAL y(nb+2); // zero values will be added to the begin and end 
    
  // assign data, Distrib(D=0)==0, Distrib(Dmax=0)==0
  x(0) = 0.; y(0) = 0.;
  x(nb+1) = mDmax; y(nb+1) = 0.;
    
  REAL *px = x.data()+1; REAL *py = y.data()+1;
  const REAL *pd1 = mD1.data(); const REAL *pd2 = mD2.data(); const REAL *pD = mDistrib.data();
  for(int i=0; i<nb; i++) {
    *px = (*pd1+*pd2)/2.;
    *py = *pD/pow(*px,3);
    px++; py++; pd1++; pd2++; pD++;
  }
    
  // cubic spline approximation with zero derivatives at D=0 and D=Dmax points
  CubicSpline spline(x,y,0.,0.);

  // calculate integral of the size distribution unsigned curvature over
  // the whole distribution interval

  // first derivative
  const int nint = int(mDmax/mCurvIntStep)+1;
  const CrystVector_REAL yy1 = spline.Derivative(0.,mCurvIntStep,nint);
  // second derivative
  const CrystVector_REAL yy2 = spline.Derivative(0.,mCurvIntStep,nint);
  
  // integrate using trapezoidal rule
  REAL sum = 0.;
  const REAL *p1 = yy1.data()+1;
  const REAL *p2 = yy2.data()+1;
  for(int k=1; k<nint-1; k++)
    sum += abs(*p2)/(1. + (*p1)*(*p1));
  sum += 0.5 * ( abs(yy2(0))/(1. + yy1(0)*yy1(0)) + abs(yy2(nint-1))/(1. + yy1(nint-1)*yy1(nint-1)) );
  sum *= mCurvIntStep;
  
  VFN_DEBUG_EXIT("MStruct::SizeDistribBroadeningEffect::CalcCurvInt():End:"<<sum,11);

  return sum;
}

void SizeDistribBroadeningEffect::RebuildLSQRegOpList ()
{
  // mDmax and mCurvIntStep are set at same time as the distribution is defined
  
  // clear the list of regularization operators
  mLSQRegOpList.clear();
  mLSQRegOpNorm.clear();

  // build new operators
  for(int iOp=0; iOp<mLSQRegTypeList.size(); iOp++) {
    
    // prepare attributes for the given option
    switch ( mLSQRegTypeList[iOp] ) {

    case LSQRegOpt_None: // no regularization - nothing to do

      break;

    case LSQRegOpt_DistribDeriv: // simple arithmetic distribution derivative optimization

      mLSQRegOpList.push_back( CreateLSQRegOpDistribDeriv(mLSQRegWeightList[iOp]) );
      mLSQRegOpNorm.push_back(1.);
      break;

    case LSQRegOpt_VolumeDistribDeriv: // simple (volume weighted) distribution derivative optimization
    
      mLSQRegOpList.push_back( CreateLSQRegOpVolumeDistribDeriv(mLSQRegWeightList[iOp]) );
      mLSQRegOpNorm.push_back(1.);
      break;
    
    case LSQRegOpt_BothDistribDeriv: // both arithetic distribution and volume weighted
                                     // distribution derivatives optimization
 
      mLSQRegOpList.push_back( CreateLSQRegOpVolumeDistribDeriv(mLSQRegWeightList[iOp]) );
      mLSQRegOpNorm.push_back(1.);
      mLSQRegOpList.push_back( CreateLSQRegOpDistribDeriv(mLSQRegWeightList[iOp]) );
      mLSQRegOpNorm.push_back(1.);
      break;

    case LSQRegOpt_CurvIntegral: // curvature integral minimization
    
      // unsupportd
      break;
  
    default: // error

      cerr << "< MStruct::SizeDistribBroadeningEffect::RebuildLSQRegOpList()\n";
      cerr << "\t"<<"iOp:"<<iOp<<", option="<<mLSQRegTypeList[iOp]<<" \t"<<"Unknown option. >"<<endl;
      throw ObjCrystException("MStruct::SizeDistribBroadeningEffect::RebuildLSQRegOpList(): \
             Logical error.");
      break;

    } // switch

  } // iOp
}

void SizeDistribBroadeningEffect::BeginOptimization (const bool allowApproximations,
						     const bool enableRestraints)
{
  mBeginEndOptimizationCalled++;

  // When called for the first time, print Constraints and Regularization Info
  if(mBeginEndOptimizationCalled==1) {
    this->PrintConstraintsStatistics();
    this->PrintRegularizationStatistics();
  }
}

void SizeDistribBroadeningEffect::EndOptimization ()
{
  mBeginEndOptimizationCalled--;

  // When called for the last time, print Constraints and Regularization Info
  if(mBeginEndOptimizationCalled==0) {
    this->PrintConstraintsStatistics();
    this->PrintRegularizationStatistics();
  }
}

void SizeDistribBroadeningEffect::PrintRegularizationStatistics () const
{
  if(this->GetNbLSQRegularizationOperator(0)==0) return;

  cout << "SizeDistribBroadeningEffect: " << this->GetName() << " Regularization Statistics\n";
  cout << " ------------------------------------------------------------------------------- \n";
  for(int iOp=0; iOp<this->GetNbLSQRegularizationOperator(0); iOp++) {
    const ObjCryst::LSQRegularizationOperator & regOp = this->GetLSQRegularizationOperator(iOp,0);
    REAL ChiSq = regOp.GetValue();
    REAL Lambda = regOp.GetRegularizationOperatorWeight();
    cout << "(non-weighted)RegOp-ChiSq" << iOp << ": " << ChiSq;
    cout << "\t\tWeight: " << Lambda << "\t\t(" << regOp.GetName() << ")\n";
  }
}

void SizeDistribBroadeningEffect::PrintConstraintsStatistics () const
{
  if(this->GetNbLSQConstraints()==0) return;

  cout << "SizeDistribBroadeningEffect: " << this->GetName() << " Constraints Statistics\n";
  cout << " ------------------------------------------------------------------------------- \n";
  for(int iCon=0; iCon<this->GetNbLSQConstraints(); iCon++) {

    std::vector< const ObjCryst::RefinablePar* > parList;
    CrystVector_REAL coef;
    
    this->GetLSQConstraint(iCon, parList, coef);

    REAL val = 0.;
    for(int i=0; i<parList.size(); i++)
      val += coef(i) * parList[i]->GetValue();

    cout << "Constraint" << iCon << ": " << val;
    cout << "\t\t(VolumeDistribInteg: " << this->GetVolumeDistribIntegral();
    cout << ", DistribInteg: " << this->GetDistribIntegral() << ")\n";
  }
}

//----------------------------------
//----------------------------------
//----------------------------------

// RandomSizeDistribBroadeningEffect
RandomSizeDistribBroadeningEffect::RandomSizeDistribBroadeningEffect()
  :mNbins(0),mBinWidth(0.),mDistrib(0)
{}

const string& RandomSizeDistribBroadeningEffect::GetClassName()const
{
  static const string className = "MStruct::RandomSizeDistribBroadeningEffect";
  return className;
}

void RandomSizeDistribBroadeningEffect::SetDistribution(const int Nbins, const REAL Dmax,
							const CrystVector_REAL &distrib)
{
  VFN_DEBUG_MESSAGE("MStruct::RandomSizeDistribBroadeningEffect::SetDistribution(...):Begin",11);
	
  bool rebuild = ( mNbins != Nbins );
	
  if (rebuild) {
    // Number of params. has changed. The old params. and params. sets have to be removed.
    
    // Remove all possible own RefinablePar sets
    EraseAllParamSet();

    // Remove all registred parameters connected with the old size distribution
    try {
      for(int i=0; i<mNbins; i++) {
	ObjCryst::RefinablePar &par = this->GetPar(&mDistrib(i));
	this->RemovePar(&par);
      }
    }
    catch (std::exception &e) {
      cerr << "< MStruct::RandomSizeDistribBroadeningEffect::SetDistribution()\n";
      cerr << "Unexpected exception: " << e.what() << "\n";
      cerr << "Unexpected exception thrown during removing old parameters from the object.\n >" << endl; 
      throw ObjCrystException("MStruct::RandomSizeDistribBroadeningEffect::SetDistribution(): Program error.");
    }
  } // rebuild==true

  // set new distribution
  mBinWidth = (Nbins!=0) ? Dmax*10./Nbins : 0.; // nm->A
  mDistrib.resize(Nbins);
  mNbins = Nbins;
  
  if( distrib.numElements()==mNbins && fabs(mBinWidth)>1e-6 ) {
    mDistrib = distrib;
    // check if distribution is properly normalised
    mDistrib *= 1./(mDistrib.sum()*mBinWidth);
  }
  else // generate new random set
    GenerateRandomDistrib(); // be careful here - maybe the parameters are not initialised

  if (rebuild) {
    // Number of params. has changed. Object has to be reinitialised.
    Init();
  }

  mClockMaster.Click(); // Object has been changed

  mClockReflProfCalc.Reset(); // Distribution has changed, reflection profile should be recomputed
	
  VFN_DEBUG_EXIT("MStruct::RandomSizeDistribBroadeningEffect::SetDistribution:End",11);
}

void RandomSizeDistribBroadeningEffect::GenerateRandomDistrib()
{
  VFN_DEBUG_MESSAGE("MStruct::RandomSizeDistribBroadeningEffect::GenerateRandomDistrib:Begin",11);
  
  // for each bin select a random value from uniform distribution in the interval [0,1)
  REAL *p = mDistrib.data();
  //for(int i=0; i<mNbins; i++) { *p = real_rand()/pow(i+0.5,3); p++; } // v2
  for(int i=0; i<mNbins; i++) { *p = real_rand(); p++; } // v1

  // normalise the distrution so sum(i=0..Nbins-1, P(i)*BinWidth)==1 holds
  if ( mNbins>0 && fabs(mBinWidth)>0. ) mDistrib *= 1./(mDistrib.sum()*mBinWidth);

  // Distribution has changed, reflection profile should be recomputed
  mClockReflProfCalc.Reset();

  // Object has changed
  mClockMaster.Click();

  VFN_DEBUG_EXIT("MStruct::RandomSizeDistribBroadeningEffect::GenerateRandomDistrib:End",11);
}

CrystVector_REAL RandomSizeDistribBroadeningEffect::GetProfile(const CrystVector_REAL &x,
							       const REAL xcenter,
							       const REAL h, const REAL k, const REAL l)
{
  VFN_DEBUG_MESSAGE("MStruct::RandomSizeDistribBroadeningEffect::GetProfile():Begin",11);

  int nbPoints = x.numElements(); 
  CrystVector_REAL profile(nbPoints);

  // size effect

  if(mNbins==0 || fabs(mBinWidth)<1e-6) {
    profile = 1.;
  } else {
    const REAL *p1 = x.data();
    REAL *p2 = profile.data();
    double t = 0;
    for(int j=0; j<mDistrib.numElements(); j++)
      t += (2*j+1)*( pow(double(j+1),2) + pow(double(j),2) ) * mDistrib(j);//pow(int,int)-not working with gcc4.1

    for(int i=0;i<nbPoints;i++) {
      double x = *p1++; x = fabs(x)/mBinWidth;
      double A = 0.;
      if(x<1.e-4)
	A = 1.0;
      else {
	for(int j=0; j<mDistrib.numElements(); j++) {
	
	  if(x<j)
	    A += ( (j+1+x)*pow(j+1-x,3) - (j+x)*pow(j-x,3) ) / t * mDistrib(j);
	  else if(x<j+1)
	    A += (j+1+x)*pow(j+1-x,3) / t * mDistrib(j);
	} // for(int j=0; j<mDistrib.numElements(); j++)
      }
      *p2++ = (REAL) A;
    } // for(int i=0;i<nbPoints;i++)
  } // if(mNbins==0 || fabs(mBinWidth)<1e-6)

  if (bsavecalc && xcenter>=xcenterlimits[0] && xcenter<=xcenterlimits[1]) {
    ofstream F("profileARSD.dat");
    F<<"# xcenter="<<xcenter*RAD2DEG<<", ApproxFWHM="<<GetApproxFWHM(xcenter,h,k,l)*RAD2DEG<<" (deg)"<<endl;
    for(int i=0; i<mDistrib.numElements(); i++)
      F<<"# "<<setw(8)<<i*mBinWidth/10.<<setw(8)<<(i+1)*mBinWidth/10.<<"  "<<setw(16)<<mDistrib(i)<<"\n";
    for(int i=0;i<nbPoints;i++)
      F<<setw(18)<<x(i)<<setw(18)<<profile(i)<<endl;
    F.close();
  }

  VFN_DEBUG_EXIT("MStruct::RandomSizeDistribBroadeningEffect::GetProfile():End",11);
  return profile;
}

REAL RandomSizeDistribBroadeningEffect::GetApproxFWHM(const REAL xcenter,
					 const REAL h, const REAL k, const REAL l)const
{
  VFN_DEBUG_MESSAGE("MStruct::RandomSizeDistribBroadeningEffect::GetApproxFWHM:Begin",11);

  const Radiation &r =
    GetParentReflectionProfile().GetParentPowderPatternDiffraction().GetRadiation();
  
  REAL t1 = 0., t2 = 0.;
  for(int i=0; i<mDistrib.numElements(); i++) {
    t1 += (2*i+1)*( pow(double(i+1),2) + pow(double(i),2) ) * mDistrib(i);//pow(int,int)-not working with gcc4.1
    t2 += 0.3*(5*i*(i*(pow(double(i),2)+2*i+2)+1)+1) * mDistrib(i);//pow(int,int)-not working with gcc4.1
  }
  t2 *= mBinWidth;
  
  const REAL fwhm = (mNbins>0 && fabs(mBinWidth)>0.) ? 1.3*r.GetWavelength()(0)/cos(0.5*xcenter)*t1/t2 : 0.0;

  VFN_DEBUG_EXIT("MStruct::RandomSizeDistribBroadeningEffect::GetApproxFWHM:End:"<<fwhm*RAD2DEG,11);
  return fwhm;
}

bool RandomSizeDistribBroadeningEffect::IsRealSpaceType()const {
  return true;
}

void RandomSizeDistribBroadeningEffect::GlobalOptRandomMove(const REAL mutationAmplitude,
							    const RefParType *type)
{
  VFN_DEBUG_MESSAGE("MStruct::RandomSizeDistribBroadeningEffect::GlobalOptRandomMove:Begin",11);

  GenerateRandomDistrib();

  VFN_DEBUG_EXIT("MStruct::RandomSizeDistribBroadeningEffect::GlobalOptRandomMove:End",11);
}


void RandomSizeDistribBroadeningEffect::Init()
{
  VFN_DEBUG_MESSAGE("MStruct::RandomSizeDistribBroadeningEffect::Init:Begin",11);
  // add parameters (distribution values for each bin)
  for(int i=0; i<mDistrib.numElements(); i++) {
    // generate parameter name
    string name;
    {
      stringstream s;
      s << "rsD_";
      if(!(this->GetName().empty())) s << "_";
      s << i;
      name = s.str();
    }
    // create a new parameter
    RefinablePar tmp(name,&mDistrib(i),0.,1.e4,
		     gpRefParTypeScattDataProfileWidth,
		     REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false,1.);
    tmp.AssignClock(mClockMaster);
    tmp.SetDerivStep(0.01);
    //tmp.SetGlobalOptimStep(.005);
    this->AddPar(tmp);
  }
  VFN_DEBUG_EXIT("MStruct::RandomSizeDistribBroadeningEffect::Init:End",11);
}

//----------------------------------
//----------------------------------
//----------------------------------

#ifdef __DEPRECATED_DoubleComponentBroadeningEffect__

// DoubleComponentBroadeningEffect
DoubleComponentBroadeningEffect::DoubleComponentBroadeningEffect():
mWeigth(0.)
{}

DoubleComponentBroadeningEffect::SetComponents(const ReflectionProfileComponent* effect1, const ReflectionProfileComponent* effect2)
{
	// check if the input pointers are not NULL
	if(effect1==NULL || effect1==NULL) {
		cerr << "< MStruct::DoubleComponentBroadeningEffect::SetComponents(...)\n";
		cerr << "\t"<<"effect1="<<effect1<<", effect2="<<effect2<<"\n";
		cerr << "\t"<<"NULL pointer to a given effect. >" << endl;
		throw ObjCrystException("MStruct::DoubleComponentBroadeningEffect::SetComponents(...): \
             Wrong input argument.");
	}
	
	// check if the given effects are of same type (real/reciprocal space type)
	try {
		if(effect1->IsRealSpaceType()!=effect2->IsRealSpaceType()) {
			cerr << "< MStruct::DoubleComponentBroadeningEffect::SetComponents(...)\n";
			cerr << "\t"<<"effect1="<<effect1<<", effect2="<<effect2<<"\n";
			cerr << "\t"<<"The given effects are not of same (real/reciprocal space) type. >" << endl;
			throw ObjCrystException("MStruct::DoubleComponentBroadeningEffect::SetComponents(...): \
             	Wrong input argument.");
		}
	}
	catch(std::exception &e) {
		err << "< MStruct::DoubleComponentBroadeningEffect::SetComponents(...)\n";
		cerr << "\t"<<"effect1="<<effect1<<", effect2="<<effect2<<"\n";
		cout << "Unexpected exception: " << e.what() << "\n";
		cout << "Unexpected exception thrown during setting the DoubleComponentBroadeningEffect Object.\n >" << endl; 
		throw;
	}
	
	mEffect1 = effect1;
	mEffect2 = effect2;
}

void DoubleComponentBroadeningEffect::SetProfilePar(const REAL weight)
{
	mWeigtht = weigth;
  mClockMaster.Click();  
}


REAL DoubleComponentBroadeningEffect::GetApproxFWHM(const REAL xcenter,
		     																						const REAL h, const REAL k, const REAL l)const
{
	// check if broadening effect components are set
	if(mEffect1==NULL || mEffect1==NULL) {
		cerr << "< MStruct::DoubleComponentBroadeningEffect::GetApproxFWHM(...)\n";
		cerr << "\t"<<"Broadening effect components are not set. >" << endl;
		throw ObjCrystException("MStruct::DoubleComponentBroadeningEffect::GetApproxFWHM(...): \
             Object not properly initialised.");
	}
	
	return = sqrt( (1.-mWeight)*pow(mEffect1->GetApproxFWHM(xcenter,h,k,l),2)
								 + mWeight*pow(mEffect2->GetApproxFWHM(xcenter,h,k,l),2) );
}

bool IsRealSpaceType()const
{
	// check if broadening effect components are set
	if(mEffect1==NULL || mEffect1==NULL) {
		cerr << "< MStruct::DoubleComponentBroadeningEffect::IsRealSpaceType()\n";
		cerr << "\t"<<"Broadening effect components are not set. >" << endl;
		throw ObjCrystException("MStruct::DoubleComponentBroadeningEffect::IsRealSpaceType(): \
             Object not properly initialised.");
	}
	
	// both components has to be of same type (already checked)
	return mEffect1->IsRealSpaceType();
}

bool IsAnisotropic()const
{
	// check if broadening effect components are set
	if(mEffect1==NULL || mEffect1==NULL) {
		cerr << "< MStruct::DoubleComponentBroadeningEffect::IsAnisotropic()\n";
		cerr << "\t"<<"Broadening effect components are not set. >" << endl;
		throw ObjCrystException("MStruct::DoubleComponentBroadeningEffect::IsAnisotropic(): \
             Object not properly initialised.");
	}
	
	// if any of two effects is anisotropic than is anisotropic
	return (mEffect1->IsAnisotropic() || mEffect2->IsAnisotropic());
} 

#endif // __DEPRECATED_DoubleComponentBroadeningEffect__

//----------------------------------
//----------------------------------
//----------------------------------

// DislocationBroadeningEffectSvB
DislocationBroadeningEffectSvB::DislocationBroadeningEffectSvB()
  :mReOrMWilk(100.), mRou(0.0001), mCg0(1.), mQ1(0.), mQ2(0.), mUseMWilk(false), mFormula(0),
   mArgument(0), mKaganerEta0(2.2), mKaganerEta1(0.), mCellType(MSTRUCT_NONE), mACr(1.)
{
  InitParameters(false);
}

void DislocationBroadeningEffectSvB::SetParentReflectionProfile(const ReflectionProfile &s)
{
  // Call the superclass method to ensure functionality.
  ReflectionProfileComponent::SetParentReflectionProfile(s);
  
  // Init auxiliary parameters (cell type, lenfth of Burgers vector)
  SetAuxParameters();
}

CrystVector_REAL DislocationBroadeningEffectSvB::GetProfile(const CrystVector_REAL &x,
							    const REAL xcenter,
							    const REAL h, const REAL k, const REAL l)
{
  // coeficients for "Full" Wilkens approx formula - series expansion up to x^15
  static const double WilkCoef15[6] = { 32./11025    /M_PI, 2./19845.   /M_PI, 4./343035.   /M_PI,
					 5./2208492. /M_PI, 1./1673100. /M_PI, 7./36067200. /M_PI};

  // nb of points and vector for calc. result
  int nbPoints = x.numElements(); 
  CrystVector_REAL profile(nbPoints);
  
  // First of all we need to check if the object is properly initialized:
  // Auxiliary variables are properly set.
  
  // Generally we need access to the parent ReflectionProfile and ParentPowderPatternDiffraction
  // objects and their GetRadiation and UnitCell objects.

  try {
    // Init auxiliary parameters (cell type, lenfth of Burgers vector) if necessary,
    // calculate length of the difraction vector and dislocation contrast factor
    double s0 = 0., Chkl = 0., t = 0.;
		
    // - these all done by the PrepareCalcAuxParams method
    PrepareCalcAuxParams(xcenter,h,k,l,s0,Chkl,t);
    
    // dislocation model params
    REAL rou = 0.0001, re = 100., MWilk = 1.;
    if (mUseMWilk==true) {
      // MWilk and rou used as parameters
      MWilk = mReOrMWilk; rou = mRou; re = (rou<1.e-7) ? 1.e4 : MWilk/sqrt(rou);
    } else {
      // rou and re used as parameters
      rou = mRou; re = mReOrMWilk; MWilk = re*sqrt(rou);
    }

    double eta;
    // ---- T = Y * (Z*x/Re)^2 * f(Z*x/Re) ----
    
    double Z, Y;
    switch (mArgument) {
    case 0: // x/Re
      Z = 1.;
      Y = 0.5*M_PI*Chkl*(s0*s0)*(mb*mb)*(MWilk*MWilk);
      break;
    case 1: // (b*g)*x/Re
      Z = mb*s0;
      Y = 0.5*M_PI*Chkl*(MWilk*MWilk);
      break;
    case 2: // (sqrt(Chkl)*b*g)*x/Re
      Z = sqrt(Chkl)*mb*s0;
      Y = 0.5*M_PI*(MWilk*MWilk);
      break;
    default: // unknown option
      profile = 1.;
      throw ObjCrystException("DislocationBroadeningEffectSvB::GetProfile: Wrong model-argument type!");
      break;
    } // switch (mArgument)

    if( rou>=1.e-7 ) { // calculation
      const REAL *p1 = x.data();
      REAL *p2 = profile.data();

      switch (mFormula) {
      case 0: // van Berkum --------------------

	// Using van Berkum's formula for x<=1
	// (the 3rd order approx of the full Wilkens's formula)
	for(int i=0;i<nbPoints;i++) {
	  double x = Z*fabs(*p1)/re;
	  if (x<1.e-4)
	    *p2 = 1.;
	  else {
	    double a = (x<=1) ? (-log(x) + 1.75 - M_LN2 + (1./6. - 32./225./M_PI*x)*x*x)
	      : (256./45./M_PI - (11./24. + log(2.*x)/4.)/x)/x;
	    //: (256./45./M_PI - (11./24. + M_LN2*x/4.)/x)/x;
	    a *= Y*(x*x);
	    *p2 = (REAL) exp(-a);
	  }
	  p1++; p2++;
	}
	break;

      case 1: // full Wilkens --------------------

	// Full Wilkens formula for x<=1
	// approx: same as van Berkum's series expansion, but done up to the 15th order,
	// which should ensure accuracy better than 1e-7
	for(int i=0;i<nbPoints;i++) {
	  double x = Z*fabs(*p1)/re;
	  if (x<1.e-4)
	    *p2 = 1.;
	  else {
	    double a = 0.;
	    if(x<=1) {
	      double x2 = x*x;
	      a += (-log(x) + 1.75 - M_LN2 + (1./6. - 32./225./M_PI*x)*x2);
	      double xt = x*x2*x2;
	      for(int iw=0; iw<6; iw++) {
		double t = WilkCoef15[iw]*xt;
		if( t<1e-8 ) break;
		a += t;
		xt *= x2;
	      } // iw
	    } else
	      a += (256./45./M_PI - (11./24. + log(2.*x)/4.)/x)/x;
	    // (256./45./M_PI - (11./24. + M_LN2*x/4.)/x)/x;

	    a *= Y*(x*x);
	    *p2 = (REAL) exp(-a);
	  }
	  p1++; p2++;
	} // i
	break;

      case 2: // Kaganer --------------------
	
	eta = mKaganerEta0 + mKaganerEta1/(Z*re);

	for(int i=0;i<nbPoints;i++) {
	  double x = Z*fabs(*p1)/re;
	  if (x<1.e-4)
	    *p2 = 1.;
	  else {
	    double a = -log( x/(eta+x) );
	    
	    a *= Y*(x*x);
	    *p2 = (REAL) exp(-a);
	  }
	  p1++; p2++;
	} // i
	break;
	
      default: // default - error --------------------
	
	profile = 1.;
	throw ObjCrystException("DislocationBroadeningEffectSvB::GetProfile: Wrong model-formula type!");
	break;

      } // switch
    } else
      profile = 1.;
    
    if (bsavecalc && xcenter>=xcenterlimits[0] && xcenter<=xcenterlimits[1]) {
      ofstream F("profileAD.dat");
      F<<"# re="<<re<<",rou="<<rou<<",MWilk="<<MWilk<<",Cg0="<<mCg0<<",q1="<<mQ1<<",q2="<<mQ2;
      F<<", (hkl)=("<<h<<k<<l<<"), t="<<t<<", Chkl="<<Chkl<<", b="<<mb<<" ,s0="<<s0;
      F<<", approxFWHM="<<GetApproxFWHM(xcenter,h,k,l)*RAD2DEG;
      F<<",xcenter="<<xcenter*RAD2DEG<<endl;
      for(int i=0;i<nbPoints;i++)
      	F<<setw(18)<<x(i)<<setw(18)<<profile(i)<<endl;
      F.close(); }
  }
	
  catch(std::exception &e) {
    cout << "< MStruct::DislocationBroadeningEffectSvB::GetProfiles(...)\n";
    cout << "exception: " << e.what() << "\n";
    cout << "Exception thrown during calculation of the Fourier Coefficients.\n \
				 Maybe some objects are not properly initialized or some parametrs\n \
				 have wrong values. >" << endl; 
    throw;
  }
  
  return profile;
}

REAL DislocationBroadeningEffectSvB::GetApproxFWHM(const REAL xcenter,
					 const REAL h, const REAL k, const REAL l)const
{
  // First of all we need to check if the object is properly initialized:
  // Auxiliary variables are properly set.
  
  // Generally we need access to the parent ReflectionProfile and ParentPowderPatternDiffraction
  // objects and their GetRadiation and UnitCell objects.
	
  double fwhm = 0.;
	
  try {
    // Init auxiliary parameters (cell type, lenfth of Burgers vector) if necessary,
    // calculate length of the difraction vector and dislocation contrast factor
    double s0 = 0., Chkl = 0., t = 0.;
		
    // - these all done by the PrepareCalcAuxParams method
    PrepareCalcAuxParams(xcenter,h,k,l,s0,Chkl,t);
    
    // dislocation model params
    REAL rou = 0.0001, re = 100., MWilk = 1.;
    if (mUseMWilk==true) {
      // MWilk and rou used as parameters
      MWilk = mReOrMWilk; rou = mRou; re = (rou<1e-7) ? 1.e4 : MWilk/sqrt(rou);
    } else {
      // rou and re used as parameters
      rou = mRou; re = mReOrMWilk; MWilk = re*sqrt(rou);
    }

    // Calculation
    const double y = sqrt(rou*Chkl)*re*mb*s0;
    // numerical approximation with absolute precision app. 1e-3 in the range y from the interval (0.,5.)
    // MSVC have no erf(REAl) hence the identity erf(x) = 1 - erfc(x) is used (erfc() is implemented in ObjCryst++) 
    fwhm = 1.8865 * 0.5 * (y*(1.-erfc(sqrt(M_LN2)*y/0.5622)) + 0.5622/sqrt(M_LN2*M_PI)*(exp(-M_LN2*pow(y/0.5622,2))-1.)) +
      1.2343 * M_1_PI * (y*atan(y/2.3364) - 0.5 * 2.3364 * log(1. + pow(y/2.3364,2)));
    fwhm /= re; // (1/A)
  }
	
  catch(std::exception &e) {
    cout << "< MStruct::DislocationBroadeningEffectSvB::GetApproxFWHM(...)\n";
    cout << "exception: " << e.what() << "\n";
    cout << "Exception thrown during calculation of the FWHM guess for the dislocation effect.\n \
				 Maybe some objects are not properly initialized or some parametrs\n \
				 have wrong values. >" << endl; 
    throw;
  }
	
  const Radiation &r = GetParentReflectionProfile().
    GetParentPowderPatternDiffraction().GetRadiation();
  
  return fwhm*r.GetWavelength()(0)/cos(0.5*xcenter); // (rad)
}

void DislocationBroadeningEffectSvB::PrepareCalcAuxParams(const REAL xcenter, const REAL h, const REAL k, const REAL l,
																								 				  double &s0, double &Chkl, double &t)const
{
  // First of all we need to check if the object is properly initialized:
  // Auxiliary variables are properly set.
	
  // Generally we need access to the parent ReflectionProfile and ParentPowderPatternDiffraction
  // objects and their GetRadiation and UnitCell objects.
	
  // Init auxiliary parameters (cell type, lenfth of Burgers vector) if necessary
  if(GetParentReflectionProfile().GetParentPowderPatternDiffraction().
     GetCrystal().GetClockMaster()<mClockAuxParams) SetAuxParameters();

  // Get radiation and wavelength
  const Radiation &rad = GetParentReflectionProfile().
    GetParentPowderPatternDiffraction().GetRadiation();
  const REAL Lambda = rad.GetWavelength()(0);

  // length of the difraction vector
  s0 = 2.*sin(0.5*xcenter)/Lambda;  		
  
  // dislocation contrast factor
  Chkl = 0.;
  t = 0.;
	
  // ref: T.Ungar,J.Gubicza,G.Ribarik,A.Borbely,J.Appl.Cryst.(2001)34,298-310
  // ref: I.C.Dragomir,T.Ungar,J.Appl.Cryst(2002)35,556-564
	
  switch(mCellType) {
  case MSTRUCT_FCC:
  case MSTRUCT_BCC:
  case MSTRUCT_SC: {
    t = REAL((h*h)*(k*k)+(k*k)*(l*l)+(l*l)*(h*h))/pow(h*h+k*k+l*l,2);
    Chkl = mCg0*(1. + mQ1*t); }
    break;
  case MSTRUCT_HCP:
    t = 0.5*(l*l)/((h*h)+(h*k)+(k*k)+(mACr*mACr)*(l*l)); 
    Chkl = mCg0*(1. + mQ1*t + mQ2*(t*t));
    break;
  case MSTRUCT_NONE:
    Chkl = mCg0;
    break;
  default:
    throw ObjCrystException("Invalid cell type.");
    break;
  }
}

void DislocationBroadeningEffectSvB::SetProfilePar(const REAL reOrRou, const REAL rou, const REAL cg0,
						   const REAL q1, const REAL q2)
{
  mReOrMWilk = (mUseMWilk==true) ? reOrRou : 10.*reOrRou;
  mRou = 0.01*rou;
  mCg0 = cg0;
  mQ1 = q1;
  mQ2 = q2;
  mClockMaster.Click();  
}

void DislocationBroadeningEffectSvB::SetUseMWilk(const bool useMWilk)
{
  mUseMWilk = useMWilk;
  InitParameters(true);
  mClockMaster.Click();
}

void DislocationBroadeningEffectSvB::SetFormula(const int formula, const int arg)
{
  mFormula = formula;
  mArgument = arg;
  InitParameters(true);
  mClockMaster.Click();
}

void DislocationBroadeningEffectSvB::InitParameters(const bool reinitialize)
{
  if (reinitialize==true) {
    // remove old parameters

    long ipar;
    ipar = this->FindPar( &mReOrMWilk ); if(ipar>=0) this->RemovePar(&this->GetPar(ipar));
    ipar = this->FindPar( &mRou ); if(ipar>=0) this->RemovePar(&this->GetPar(ipar));
    ipar = this->FindPar( &mCg0 ); if(ipar>=0) this->RemovePar(&this->GetPar(ipar));
    ipar = this->FindPar( &mQ1 ); if(ipar>=0) this->RemovePar(&this->GetPar(ipar));
    ipar = this->FindPar( &mQ2 ); if(ipar>=0) this->RemovePar(&this->GetPar(ipar));
    ipar = this->FindPar( &mKaganerEta0 ); if(ipar>=0) this->RemovePar(&this->GetPar(ipar));
    ipar = this->FindPar( &mKaganerEta1 ); if(ipar>=0) this->RemovePar(&this->GetPar(ipar));
  }

  if (mUseMWilk==true)
  {
    RefinablePar tmp("MWilk", &mReOrMWilk, 0., 100.,
                     gpRefParTypeScattDataProfileWidth,
                     REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false,1.0);
    tmp.AssignClock(mClockMaster);
    tmp.SetDerivStep(0.05);
    this->AddPar(tmp);
  }
  else
  {
    RefinablePar tmp("Re", &mReOrMWilk, 10., 1.e4,
                     gpRefParTypeScattDataProfileWidth,
                     REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false,0.1);
    tmp.AssignClock(mClockMaster);
    tmp.SetDerivStep(5.0);
    this->AddPar(tmp);
  }

  {
    RefinablePar tmp("Rou", &mRou, 0., 0.01,
                     gpRefParTypeScattDataProfileWidth,
                     REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false,1.e2);
    tmp.AssignClock(mClockMaster);
    tmp.SetDerivStep(1.e-3);
    this->AddPar(tmp);
  }
  {
    RefinablePar tmp("Cg0", &mCg0, 0., 100.,
                     gpRefParTypeScattDataProfileWidth,
                     REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false,1.);
    tmp.AssignClock(mClockMaster);
    tmp.SetDerivStep(0.1);
    this->AddPar(tmp);
  }
  {
    RefinablePar tmp("Q1", &mQ1, -10., 10.,
                     gpRefParTypeScattDataProfileWidth,
                     REFPAR_DERIV_STEP_ABSOLUTE,false,true,true,false,1.);
    tmp.AssignClock(mClockMaster);
    tmp.SetDerivStep(0.05);
    this->AddPar(tmp);
  }
  {
    RefinablePar tmp("Q2", &mQ2, -10., 10.,
                     gpRefParTypeScattDataProfileWidth,
                     REFPAR_DERIV_STEP_ABSOLUTE,false,true,true,false,1.);
    tmp.AssignClock(mClockMaster);
    tmp.SetDerivStep(0.05);
    this->AddPar(tmp);
  }
  if(mFormula==2) { // formula = Kaganer
    RefinablePar tmp("KaganerEta0", &mKaganerEta0, 0.1, 5.,
                     gpRefParTypeScattDataProfileWidth,
                     REFPAR_DERIV_STEP_ABSOLUTE,false,true,true,false,1.);
    tmp.AssignClock(mClockMaster);
    tmp.SetDerivStep(0.05);
    this->AddPar(tmp);
  }
  if(mFormula==2) { // formula = Kaganer
    RefinablePar tmp("KaganerEta1", &mKaganerEta1, 0., 10.,
                     gpRefParTypeScattDataProfileWidth,
                     REFPAR_DERIV_STEP_ABSOLUTE,false,true,true,false,1.);
    tmp.AssignClock(mClockMaster);
    tmp.SetDerivStep(0.05);
    this->AddPar(tmp);
  }
}

void DislocationBroadeningEffectSvB::SetAuxParameters()const
{
  // We need to know mainly the cell type, the Burgers vector length and a/c ratio
  // in the case of hexagonal hcp structure
  
  // To calc their proper values we need access to Crystal/UnitCell object.
  
  // We can get the Crystal/UnitCell object trough ParentPowderPatternDiffraction object,
  // but we need the object ParentReflectionProfile and the ParentPowderPatternDiffractio
  // to be accessible.
  try {
    const UnitCell &uc = GetParentReflectionProfile().GetParentPowderPatternDiffraction().GetCrystal();
	                           
    // We want to detect the lattice type (fcc,bcc,hcp) and set the length
    // of the Burgers vector.
	
    // At present time only hcp and cubic structures are supported.
    
    // hcp Mg has space group number 194, cubic groups have numbers 195-230,
    // primitive cubic structures have 1 translation vector, bcc structures have 2 and fcc 4.
	
    const int sgnb = uc.GetSpaceGroup().GetSpaceGroupNumber();
    
    if (sgnb >= 195) {
      // cubic structure
      const REAL aa = uc.GetLatticePar(0);
      switch (uc.GetSpaceGroup().GetNbTranslationVectors()) {
      case 1:
	mCellType = MSTRUCT_SC;
	mb = aa;
	break;
      case 2:
	mCellType = MSTRUCT_BCC;
	mb = sqrt(3.)/2.*aa; 
	break;
      case 4:
	mCellType = MSTRUCT_FCC;
	mb = aa/M_SQRT2;;
	break;
      default:
	mCellType = MSTRUCT_NONE;
	mb = 1.;
	// something strange
	cout << "Warning: MStruct::DislocationBroadeningEffectSvB:" << "Assuming cubic cell, \
				             but can not recognize type (fcc,bcc,sc),\n\t the length of the Burgers vector \
				             set to 1." << endl; 
	break;
      }
      mQ2 = 0.;
    }
    else if (sgnb == 194) {
      // hexagonal cell
      const REAL aa = uc.GetLatticePar(0);
      const REAL cc = uc.GetLatticePar(3);
      // P 6_3 / m m c (typical hcp group)
      mCellType = MSTRUCT_HCP;
      mb = aa;
      // a/c ratio
      mACr = aa/cc;
    }
    else {
      // not supported structure
      mCellType = MSTRUCT_NONE;
      mb = 1.;
      // something strange
      cout << "Warning: MStruct::DislocationBroadeningEffectSvB:" << "Only cubic and hcp structures \
			         supported.\n\t The cell type was not recognised to be any of them, the length of Burgers vector \
			         set to 1.\n\t Width of diffraction lines will not be anisotropic." << endl;
      mQ1 = 0.;
      mQ2 = 0.;
    }
    
    mClockAuxParams.Click();
  }
	
  catch(std::exception &e) {
    cout << "< MStruct::DislocationBroadeningEffectSvB::SetAuxParameters()\n";
    cout << "exception: " << e.what() << "\n";
    cout << "Maybe a parent ReflectionProfile object to this broadenig component\n \
				 or its parent PowderpatterDiffraction object or its Crystal object\n \
				 have not been set yet. Without them a cell type, the length of Burgers vector\n \
				 or other params. can not be calculated. >" << endl; 
    throw;;
  }
}

bool DislocationBroadeningEffectSvB::IsRealSpaceType()const
{
  return true;
}

bool DislocationBroadeningEffectSvB::IsAnisotropic()const
{
  return (abs(mQ1)>1.e-4 || abs(mQ2)>1.e-4);
}

/*
unsigned int DislocationBroadeningEffectSvB::GetNbLSQConstraints() const
{
return 1;
}

void DislocationBroadeningEffectSvB::GetLSQConstraint(const unsigned int n,
						      std::vector< const ObjCryst::RefinablePar* > &parList,
						      CrystVector_REAL &coef) const
{
  parList = std::vector< const ObjCryst::RefinablePar* >(2, NULL);
  coef.resize(2);
  
  // constraint: sqrt(rou) * Re = const
  
  parList[0] = &(this->GetPar(&mRou));
  coef(0) = 0.5*mRe;
  
  parList[1] = &(this->GetPar(&mRe));
  coef(1) = mRou;

  coef *= 1.e6*mDefaultLSQConstraintScale/coef(0);

  cout<<"contraint: "<<"Re: "<<mRe<<", rou: "<<mRou<<", c(0): "<<coef(0)<<", c(1): "<<coef(1)<<endl;
}
*/

// FaultsBroadeningEffectFCC
FaultsBroadeningEffectFCC::FaultsBroadeningEffectFCC():
mAlpha(0.), mBeta(0.), mpUnitCell(0)
{
	InitParameters();
}

void FaultsBroadeningEffectFCC::SetParentReflectionProfile(const ReflectionProfile &s)
{
	// Call the superclass method to ensure functionality.
	ReflectionProfileComponent::SetParentReflectionProfile(s);
	
	// Init auxiliary parameters (check cell type)
	SetAuxParameters();
}

REAL FaultsBroadeningEffectFCC::GetApproxFWHM(const REAL xcenter,
					 																		const REAL h, const REAL k, const REAL l)const
{
	// If (hkl)==(000) is supplied to the method, no broadening is reported.
	if(abs(h)<1e-6 && abs(k)<1e-6 && abs(l)<1e-6) return 0.;
	
	// First of all we need to check if the object is properly initialized:
	// Auxiliary variables are properly set.
	
	// Generally we need access to the parent ReflectionProfile and ParentPowderPatternDiffraction
	// objects and their GetRadiation and UnitCell objects.
	
	double fwhm = 0.;
	
	try {
		// Init auxiliary parameters (cell type) and calculate length of the difraction vector
		double s0 = 0.;
		
		// - these all done by the PrepareCalcAuxParams method (mAbsL0,mSign and mCount vectors also prepared)
		PrepareCalcAuxParams(xcenter,h,k,l,s0);
		
		// Calcualtion
		
		// We use classical Warrens formulas, peak shift is also included into the calculated fwhm
		// becouse {HKL} componets shift are included in the Fourier coefficients and not in the position
		// corrections
		
		int tsum = 0.;
		
		for(int igroup=0; igroup<mSign.numElements(); igroup++)
			if(mSign(igroup)!=0) tsum += mAbsL0(igroup);
		
		// FWHM - Warren formula (in reciprocal space units)
		//fwhm += (1.5*mAlpha+mBeta)/(h*h+k*k+l*l)/s0/mCount.sum()*tsum;
		fwhm += (1.5*mAlpha+mBeta)/s0/mCount.sum()*tsum / M_PI;   // ( FWHM ~ 1/pi * 1/D )

		tsum = 0.;
		
		for(int igroup=0; igroup<mSign.numElements(); igroup++)
			if(mSign(igroup)!=0) tsum += mSign(igroup)*mAbsL0(igroup);
			
		// shift - Warren formula (in reciprocal space units)
		fwhm += abs( mAlpha * sqrt(3.)/4./M_PI * s0/(h*h+k*k+l*l)/mCount.sum() * tsum );
	}
	
	catch(std::exception &e) {
		cout << "< MStruct::FaultsBroadeningEffectFCC::GetApproxFWHM(...)\n";
		cout << "exception: " << e.what() << "\n";
		cout << "Exception thrown during calculation of the FWHM guess for the faulting effects.\n \
				 Maybe some objects are not properly initialized or some parametrs\n \
				 have wrong values. >" << endl; 
		throw;
	}
	
	const Radiation &r = GetParentReflectionProfile().
  	GetParentPowderPatternDiffraction().GetRadiation();
  
  return fwhm*r.GetWavelength()(0)/cos(0.5*xcenter); // (rad)
}

void FaultsBroadeningEffectFCC::PrepareCalcAuxParams(const REAL xcenter, const REAL h, const REAL k, const REAL l,
														     					 				   double &s0)const
{
	// If (hkl)==(000) is supplied to the method, nothing done.
	if(abs(h)<1e-6 && abs(k)<1e-6 && abs(l)<1e-6) return;
	
	// First of all we need to check if the object is properly initialized:
	// Auxiliary variables are properly set.
	
	// Generally we need access to the parent ReflectionProfile and ParentPowderPatternDiffraction
	// objects and their GetRadiation and UnitCell objects.
	
	// Init auxiliary parameters (cell type, pointer to the UnitCell object ... ) if necessary
	if(GetParentReflectionProfile().GetParentPowderPatternDiffraction().
		GetCrystal().GetClockMaster()<mClockAuxParams) SetAuxParameters();

	// Get radiation and wavelength
	const Radiation &rad = GetParentReflectionProfile().
  			GetParentPowderPatternDiffraction().GetRadiation();
	const REAL Lambda = rad.GetWavelength()(0);

	// length of the difraction vector
	s0 = 2.*sin(0.5*xcenter)/Lambda;  		
	
	// ref: L.Velterop,R.Delhez,Th.H.deKeijser,E.J.Mittemeijer,D.Reefman,J.Appl.Cryst.(2000)33,296-306
	// ref: P.Scardi,M.Leoni,Acta.Cryst(2002)A58,190-200
	// ref: B.E.Warren, X-RAY DIFFRACTION, Addison-Wesley 1969
	
	// Equivalent (hkl) reflections should be generated, they need to be divided into
	// the groups according to the sign and value of |L0|=|h+k+l| and this information about
	// the size of |L0|, sign of |L0|=3*J+(0,+/-1) and the number of reflection in each group
	// should be stored (ref: Velterop et al.)
	
	CrystMatrix_REAL equivRefl = mpUnitCell->GetSpaceGroup().GetAllEquivRefl(h,k,l,false,false);
	
	const int maxAbsL0 = int(abs(h)+abs(k)+abs(l));
	CrystMatrix_REAL signMatrix(3,maxAbsL0+1);
	
	signMatrix = 0;
	
	for(int irefl=0; irefl<equivRefl.rows(); irefl++) {
		int absL0 = (int) abs( equivRefl(irefl,0)+equivRefl(irefl,1)+equivRefl(irefl,2) );
		signMatrix(absL0-3*(absL0/3),absL0)++;
	}
	
	int nbgroups = 0;
	for(int iAbsL0=0; iAbsL0<=maxAbsL0; iAbsL0++)
		for(int isign=0; isign<3; isign++) if (signMatrix(isign,iAbsL0)>0) nbgroups++;
	
	mAbsL0 = CrystVector_int(nbgroups);
	mCount = CrystVector_int(nbgroups);
	mSign  = CrystVector_int(nbgroups);
	
	int nn = 0;
	for(int iAbsL0=0; iAbsL0<=maxAbsL0; iAbsL0++) {
		if (signMatrix(2,iAbsL0)>0) {
			mAbsL0(nn) = iAbsL0;
			mCount(nn) = signMatrix(2,iAbsL0);
			mSign(nn) = -1;
			nn++;
		}
		if (signMatrix(0,iAbsL0)>0) {
			mAbsL0(nn) = iAbsL0;
			mCount(nn) = signMatrix(0,iAbsL0);
			mSign(nn) =  0;
			nn++;
		} 
		if (signMatrix(1,iAbsL0)>0) {
			mAbsL0(nn) = iAbsL0;
			mCount(nn) = signMatrix(1,iAbsL0);
			mSign(nn) =  1;
			nn++;
		}
	}
}

void FaultsBroadeningEffectFCC::SetProfilePar(const REAL alpha, const REAL beta)
{
  mAlpha = alpha;
  mBeta = beta;
  mClockMaster.Click();  
}

void FaultsBroadeningEffectFCC::InitParameters()
{
  {
    RefinablePar tmp("Alpha", &mAlpha, 0., 1.,
                     gpRefParTypeScattDataProfileWidth,
                     REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false,1.);
    tmp.AssignClock(mClockMaster);
    tmp.SetDerivStep(0.005);
    this->AddPar(tmp);
  }
  {
    RefinablePar tmp("Beta", &mBeta, 0., 1.,
                     gpRefParTypeScattDataProfileWidth,
                     REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false,1.);
    tmp.AssignClock(mClockMaster);
    tmp.SetDerivStep(0.005);
    this->AddPar(tmp);
  }
}

void FaultsBroadeningEffectFCC::SetAuxParameters()const
{
	// We need to check the cell type only
	
	// To get all information we need access to Crystal/UnitCell object.
	
	// We can get the Crystal/UnitCell object trough ParentPowderPatternDiffraction object,
	// but we need the object ParentReflectionProfile and the ParentPowderPatternDiffraction
	// to be accessible.
	// Pointer to the UnitCell object is stored to allow simple access to this object for eg.
	// equivalent diffractions generation.
	try {
		
		mpUnitCell = &GetParentReflectionProfile().GetParentPowderPatternDiffraction().GetCrystal();
	                           
		// We want to detect the lattice type (fcc,bcc,hcp)
	
		// Only cubic fcc structures are accepted.
	
		// hcp Mg has space group number 194, cubic groups have numbers 195-230,
		// primitive cubic structures have 1 translation vector, bcc structures have 2 and fcc 4.
	
		const int sgnb = mpUnitCell->GetSpaceGroup().GetSpaceGroupNumber();
	
		if (sgnb < 195 || mpUnitCell->GetSpaceGroup().GetNbTranslationVectors() != 4) {
			// not cubic fcc structure
			cout << "Warning: MStruct::FaultsBroadeningEffectFCC:" << "Only cubic fcc cell supported." << endl; 
			mClockAuxParams.Click();
		}
	}
	
	catch(std::exception &e) {
		cout << "< MStruct::FaultsBroadeningEffectFCC::SetAuxParameters()\n";
		cout << "exception: " << e.what() << "\n";
		cout << "Maybe a parent ReflectionProfile object to this broadenig component\n \
				 		 or its parent PowderpatterDiffraction object or its Crystal object\n \
				     have not been set yet. Without them a cell type or other inforamtion \
				 		 can not be checked. >" << endl; 
		throw ;
	}
}

bool FaultsBroadeningEffectFCC::IsRealSpaceType()const
{
	// Default value is true.
	return true;
}

bool FaultsBroadeningEffectFCC::IsAnisotropic()const
{
	// Default value is true if fault probabilities are nonzero.
	return (abs(mAlpha)>1.e-6 || abs(mBeta)>1.e-6);
}

// FaultsBroadeningEffectVelteropFCC
FaultsBroadeningEffectVelteropFCC::FaultsBroadeningEffectVelteropFCC():
mAlpha(0.), mBeta(0.), mpUnitCell(0)
{
	InitParameters();
}

void FaultsBroadeningEffectVelteropFCC::SetParentReflectionProfile(const ReflectionProfile &s)
{
	// Call the superclass method to ensure functionality.
	ReflectionProfileComponent::SetParentReflectionProfile(s);
	
	// Init auxiliary parameters (check cell type)
	SetAuxParameters();
}

CrystVector_REAL FaultsBroadeningEffectVelteropFCC::GetProfile(const CrystVector_REAL &x,
						  									const REAL xcenter,
						  									const REAL h, const REAL k, const REAL l)
{
	// nb of points and vector for calc. result
	int nbPoints = x.numElements(); 
	CrystVector_REAL profile(2*nbPoints);

	profile = 0.;
  
	// First of all we need to check if the object is properly initialized:
	// Auxiliary variables are properly set.
	
	// Generally we need access to the parent ReflectionProfile and ParentPowderPatternDiffraction
	// objects and their GetRadiation and UnitCell objects.
	
	try {
		// Init auxiliary parameters (check cell type, divide equivalent reflections into groups,
		// prepare Sign, |L0| and Count vectors for all groups ) if necessary,
		// calculate length of the difraction vector
		double s0 = 0.;
		
		// - these all done by the PrepareCalcAuxParams method
		PrepareCalcAuxParams(xcenter,h,k,l,s0);
  	
  	// calculation
		const REAL beta0 = sqrt(3.-6.*mBeta-pow(mBeta,2)-12.*mAlpha+12.*pow(mAlpha,2));
  	const REAL gamma = atan(beta0/(1.-mBeta));
  	const REAL logZ = log(sqrt(1.-2*mBeta-3.*mAlpha+3.*pow(mAlpha,2)));
  	const REAL tt = (h*h+k*k+l*l)/s0;

		const double mult = mCount.sum();
  	
 		const REAL *p1 = x.data();
 		REAL *p2 = profile.data();
 		for(int i=0;i<nbPoints;i++) {
   		double x = *p1;
   		double A = 0., B = 0.;
			if (abs(x)<1.e-4) {
     			A = 1.;
     			B = 0.;
			} else {
				for(int igroup=0; igroup<mAbsL0.numElements(); igroup++) {
					if(mSign(igroup)==0) {
						// unaffected component
						A += mCount(igroup)/mult; // mult is double type (real division)
					} else {
						// affected component
						double dx = tt/mAbsL0(igroup);
						double a = mCount(igroup)/mult*exp(logZ*abs(x)/dx);
						double b = mSign(igroup)*mBeta/beta0*x/abs(x)*a;
						// shift
						double ds = -mSign(igroup)*(gamma/2./M_PI-1./6.)/dx;
						double c = cos(2.*M_PI*ds*x); // FFT: Integrate f(x)*exp(-2*pi*q*x) dx
						double s = sin(2.*M_PI*ds*x);
						A += a*c - b*s;
						B += b*c + a*s;
					}
				}
   		}
   		*p2 = A; p2++;
   		*p2 = B; p2++;
   		p1++;
 		}

 		if (bsavecalc && xcenter>=xcenterlimits[0] && xcenter<=xcenterlimits[1]) {
 			CrystMatrix_REAL profiles;
 			CrystVector_REAL shifts;
 			
 			GetAllGroupsProfiles(profiles,shifts,x,xcenter,h,k,l);
 			
	   	ofstream F("profileAF.dat");
			F<<"# alpha="<<mAlpha<<",beta="<<mBeta;
			F<<",(hkl)=("<<h<<k<<l<<"),s0="<<s0;
			F<<",approxFWHM="<<GetApproxFWHM(xcenter,h,k,l)*RAD2DEG;
			F<<",xcenter="<<xcenter*RAD2DEG<<endl;
			for(int igroup=0; igroup<profiles.rows(); igroup++)
			{	F<<"# "<<setw(4)<<mAbsL0(igroup)<<setw(3)<<mSign(igroup)<<setw(6)<<mCount(igroup);
				F<<setw(14)<<shifts(igroup)<<endl;
			}
  		for(int i=0;i<nbPoints;i++)
      	//F<<setw(18)<<x(i)<<setw(18)<<profile(2*i)<<setw(18)<<profile(2*i+1)<<endl;
      {
      	F<<setw(18)<<x(i)<<setw(18)<<profile(2*i)<<setw(18)<<profile(2*i+1);
      	for(int igroup=0; igroup<profiles.rows(); igroup++)
      		 F<<setw(18)<<profiles(igroup,2*i)<<setw(18)<<profiles(igroup,2*i+1);
      	F<<endl;
      }
    	F.close(); }
	}
	
	catch(std::exception &e) {
		cout << "< MStruct:: FaultsBroadeningEffectVelteropFCC::GetProfile(...)\n";
		cout << "exception: " << e.what() << "\n";
		cout << "Exception thrown during calculation of the Fourier Coefficients.\n \
				 Maybe some objects are not properly initialized or some parametrs\n \
				 have wrong values. >" << endl; 
		throw;
	}

  return profile;
}

void FaultsBroadeningEffectVelteropFCC::GetAllGroupsProfiles(CrystMatrix_REAL &profiles,
					 																									 CrystVector_REAL &shifts,
																														 const CrystVector_REAL &x,
			      								  					const REAL xcenter, const REAL h, const REAL k, const REAL l)
{
	// nb of points and vector for calc. result
		int nbPoints = x.numElements(); 
		
	// First of all we need to check if the object is properly initialized:
	// Auxiliary variables are properly set.
	
	// Generally we need access to the parent ReflectionProfile and ParentPowderPatternDiffraction
	// objects and their GetRadiation and UnitCell objects.
	
	try {
		// Init auxiliary parameters (check cell type, divide equivalent reflections into groups,
		// prepare Sign, |L0| and Count vectors for all groups ) if necessary,
		// calculate length of the difraction vector
		double s0 = 0.;
		
		// - these all is done by the PrepareCalcAuxParams method
		PrepareCalcAuxParams(xcenter,h,k,l,s0);
		
		// prapre space for output data
		profiles.resize(mAbsL0.numElements(),2*nbPoints);
		shifts.resize(mAbsL0.numElements());
		
		profiles = 0.;
  	
  	// calculation
		const REAL beta0 = sqrt(3.-6.*mBeta-pow(mBeta,2)-12.*mAlpha+12.*pow(mAlpha,2));
  	const REAL gamma = atan(beta0/(1.-mBeta));
  	const REAL logZ = log(sqrt(1.-2*mBeta-3.*mAlpha+3.*pow(mAlpha,2)));
  	const REAL tt = (h*h+k*k+l*l)/s0;

		const double mult = mCount.sum();

		// for each group
		for(int igroup=0; igroup<mAbsL0.numElements(); igroup++) {
			if(mSign(igroup)==0) {
				// unaffected component
				for(int i=0;i<nbPoints;i++)
					profiles(igroup,2*i) += mCount(igroup)/mult; // mult is double type (real division)
				shifts(igroup) = 0.;
			} else {
				// affected component
				double dx = tt/mAbsL0(igroup);
				
				const REAL *p1 = x.data();
				
				for(int i=0;i<nbPoints;i++) {
					double x = *p1;
					if(abs(x)<1.e-4) {
						profiles(igroup,2*i) = mCount(igroup)/mult;
					} else {
						profiles(igroup,2*i)   = mCount(igroup)/mult*exp(logZ*abs(x)/dx);
						profiles(igroup,2*i+1) = mSign(igroup)*mBeta/beta0*x/abs(x)*profiles(igroup,2*i);
					}
					p1++;
				} // i
				
				// shift
				shifts(igroup) = -mSign(igroup)*(gamma/2./M_PI-1./6.)/dx;
						
			} // affected component	
		} // igroup
	
	} // try
	
	catch(std::exception &e) {
		cout << "< MStruct:: FaultsBroadeningEffectVelteropFCC::GetAllGroupsProfiles(...)\n";
		cout << "exception: " << e.what() << "\n";
		cout << "Exception thrown during calculation of the Fourier Coefficients.\n \
				 Maybe some objects are not properly initialized or some parametrs\n \
				 have wrong values. >" << endl; 
		throw;
	}
}

			      															
REAL FaultsBroadeningEffectVelteropFCC::GetApproxFWHM(const REAL xcenter,
					 const REAL h, const REAL k, const REAL l)const
{
	// First of all we need to check if the object is properly initialized:
	// Auxiliary variables are properly set.
	
	// Generally we need access to the parent ReflectionProfile and ParentPowderPatternDiffraction
	// objects and their GetRadiation and UnitCell objects.
	
	double fwhm = 0.;
	
	try {
		// Init auxiliary parameters (cell type, lenfth of Burgers vector) if necessary,
		// calculate length of the difraction vector and dislocation contrast factor
		double s0 = 0.;
		
		// - these all done by the PrepareCalcAuxParams method (mAbsL0,mSign and mCount vectors also prepared)
		PrepareCalcAuxParams(xcenter,h,k,l,s0);
		
		// Calcualtion
		
		// We use classical Warrens formulas, peak shift is also included into the calculated fwhm
		// becouse {HKL} componets shift are included in the Fourier coefficients and not in the position
		// corrections
		
		int tsum = 0.;
		
		for(int igroup=0; igroup<mSign.numElements(); igroup++)
			if(mSign(igroup)!=0) tsum += mAbsL0(igroup);
		
		// FWHM - Warren formula (in reciprocal space units)
		fwhm += (1.5*mAlpha+mBeta)/(h*h+k*k+l*l)/s0/mCount.sum()*tsum;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
		
		tsum = 0.;
		
		for(int igroup=0; igroup<mSign.numElements(); igroup++)
			if(mSign(igroup)!=0) tsum += mSign(igroup)*mAbsL0(igroup);
			
		// shift - Warren formula (in reciprocal space units)
		fwhm += abs( mAlpha * sqrt(3.)/4./M_PI * s0/(h*h+k*k+l*l)/mCount.sum() * tsum );
	}
	
	catch(std::exception &e) {
		cout << "< MStruct::FaultsBroadeningEffectVelteropFCC::GetApproxFWHM(...)\n";
		cout << "exception: " << e.what() << "\n";
		cout << "Exception thrown during calculation of the FWHM guess for the faulting effects.\n \
				 Maybe some objects are not properly initialized or some parametrs\n \
				 have wrong values. >" << endl; 
		throw;
	}
	
	const Radiation &r = GetParentReflectionProfile().
  	GetParentPowderPatternDiffraction().GetRadiation();
  
  return fwhm*r.GetWavelength()(0)/cos(0.5*xcenter); // (rad)
}

void FaultsBroadeningEffectVelteropFCC::PrepareCalcAuxParams(const REAL xcenter, const REAL h, const REAL k, const REAL l,
																								 				  double &s0)const
{
	// First of all we need to check if the object is properly initialized:
	// Auxiliary variables are properly set.
	
	// Generally we need access to the parent ReflectionProfile and ParentPowderPatternDiffraction
	// objects and their GetRadiation and UnitCell objects.
	
	// Init auxiliary parameters (cell type, pointer to the UnitCell object ... ) if necessary
	if(GetParentReflectionProfile().GetParentPowderPatternDiffraction().
		GetCrystal().GetClockMaster()<mClockAuxParams) SetAuxParameters();

	// Get radiation and wavelength
	const Radiation &rad = GetParentReflectionProfile().
  			GetParentPowderPatternDiffraction().GetRadiation();
	const REAL Lambda = rad.GetWavelength()(0);

	// length of the difraction vector
	s0 = 2.*sin(0.5*xcenter)/Lambda;  		
	
	// ref: L.Velterop,R.Delhez,Th.H.deKeijser,E.J.Mittemeijer,D.Reefman,J.Appl.Cryst.(2000)33,296-306
	// ref: P.Scardi,M.Leoni,Acta.Cryst(2002)A58,190-200
	
	// Equivalent (hkl) reflections should be generated, they need to be divided into
	// the groups according to the sign and value of |L0|=|h+k+l| and this information about
	// the size of |L0|, sign of |L0}=3*J(0,+/-1) and the number of reflection in each group
	// should be stored (ref: Velterop et al.)
	
	CrystMatrix_REAL equivRefl = mpUnitCell->GetSpaceGroup().GetAllEquivRefl(h,k,l,false,false);
	
	const int maxAbsL0 = int(abs(h)+abs(k)+abs(l));
	CrystMatrix_REAL signMatrix(3,maxAbsL0+1);
	
	signMatrix = 0;
	
	for(int irefl=0; irefl<equivRefl.rows(); irefl++) {
		int absL0 = (int) abs( equivRefl(irefl,0)+equivRefl(irefl,1)+equivRefl(irefl,2) );
		signMatrix(absL0-3*(absL0/3),absL0)++;
	}
	
	int nbgroups = 0;
	for(int iAbsL0=0; iAbsL0<=maxAbsL0; iAbsL0++)
		for(int isign=0; isign<3; isign++) if (signMatrix(isign,iAbsL0)>0) nbgroups++;
	
	mAbsL0 = CrystVector_int(nbgroups);
	mCount = CrystVector_int(nbgroups);
	mSign  = CrystVector_int(nbgroups);
	
	int nn = 0;
	for(int iAbsL0=0; iAbsL0<=maxAbsL0; iAbsL0++) {
		if (signMatrix(2,iAbsL0)>0) {
			mAbsL0(nn) = iAbsL0;
			mCount(nn) = signMatrix(2,iAbsL0);
			mSign(nn) = -1;
			nn++;
		}
		if (signMatrix(0,iAbsL0)>0) {
			mAbsL0(nn) = iAbsL0;
			mCount(nn) = signMatrix(0,iAbsL0);
			mSign(nn) =  0;
			nn++;
		} 
		if (signMatrix(1,iAbsL0)>0) {
			mAbsL0(nn) = iAbsL0;
			mCount(nn) = signMatrix(1,iAbsL0);
			mSign(nn) =  1;
			nn++;
		}
	}
}

void FaultsBroadeningEffectVelteropFCC::SetProfilePar(const REAL alpha, const REAL beta)
{
  mAlpha = alpha;
  mBeta = beta;
  mClockMaster.Click();  
}

void FaultsBroadeningEffectVelteropFCC::InitParameters()
{
  {
    RefinablePar tmp("Alpha", &mAlpha, 0., 1.,
                     gpRefParTypeScattDataProfileWidth,
                     REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false,1.);
    tmp.AssignClock(mClockMaster);
    tmp.SetDerivStep(0.005);
    this->AddPar(tmp);
  }
  {
    RefinablePar tmp("Beta", &mBeta, 0., 1.,
                     gpRefParTypeScattDataProfileWidth,
                     REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false,1.);
    tmp.AssignClock(mClockMaster);
    tmp.SetDerivStep(0.005);
    this->AddPar(tmp);
  }
}

void FaultsBroadeningEffectVelteropFCC::SetAuxParameters()const
{
	// We need to check the cell type only
	
	// To get all information we need access to Crystal/UnitCell object.
	
	// We can get the Crystal/UnitCell object trough ParentPowderPatternDiffraction object,
	// but we need the object ParentReflectionProfile and the ParentPowderPatternDiffraction
	// to be accessible.
	// Pointer to the UnitCell object is stored to allow simple access to this object for eg.
	// equivalent diffractions generation.
	try {
		
		mpUnitCell = &GetParentReflectionProfile().GetParentPowderPatternDiffraction().GetCrystal();
	                           
		// We want to detect the lattice type (fcc,bcc,hcp)
	
		// At present time only cubic fcc structures are supported.
	
		// hcp Mg has space group number 194, cubic groups have numbers 195-230,
		// primitive cubic structures have 1 translation vector, bcc structures have 2 and fcc 4.
	
		const int sgnb = mpUnitCell->GetSpaceGroup().GetSpaceGroupNumber();
	
		if (sgnb < 195 || mpUnitCell->GetSpaceGroup().GetNbTranslationVectors() != 4) {
			// not cubic fcc structure
			cout << "Warning: MStruct::FaultsBroadeningEffectVelteropFCC:" << "Only cubic fcc cell supported." << endl; 
			mClockAuxParams.Click();
		}
	}
	
	catch(std::exception &e) {
		cout << "< MStruct::FaultsBroadeningEffectVelteropFCC::SetAuxParameters()\n";
		cout << "exception: " << e.what() << "\n";
		cout << "Maybe a parent ReflectionProfile object to this broadenig component\n \
				 		 or its parent PowderpatterDiffraction object or its Crystal object\n \
				     have not been set yet. Without them a cell type or other inforamtion \
				 		 can not be checked. >" << endl; 
		throw ;
	}
}

bool FaultsBroadeningEffectVelteropFCC::IsRealSpaceType()const
{
	return true;
}

bool FaultsBroadeningEffectVelteropFCC::IsAnisotropic()const
{
	return (abs(mAlpha)>1.e-6 || abs(mBeta)>1.e-6);
}

// FaultsBroadeningEffectFCCBaloghUngar

FaultsBroadeningEffectFCCBaloghUngar::FaultsBroadeningEffectFCCBaloghUngar():
mFaultsType(0)
{}

CrystVector_REAL FaultsBroadeningEffectFCCBaloghUngar::GetProfile(const CrystVector_REAL &x,
						  									const REAL xcenter,
						  									const REAL h, const REAL k, const REAL l)
{
	// nb of points and vector for calc. result
	int nbPoints = x.numElements(); 
	CrystVector_REAL profile(2*nbPoints);

	profile = 0.;
	
	try {
		
		// Access to the parent ReflectionProfile and ParentPowderPatternDiffraction
		// objects and their GetRadiation and UnitCell objects is necessary.
	
		// Get radiation and wavelength
		const Radiation &rad = GetParentReflectionProfile().
  							GetParentPowderPatternDiffraction().GetRadiation();
		const REAL Lambda = rad.GetWavelength()(0);
	
		// calcualtion
		
		// get parametrs
		REAL weight[4], fwhm[4], shift0[4]; // alpha(fraction),fwhm and shift0 (1/A)
		GetSubComponentsPar(mAlpha,mFaultsType,h,k,l,weight,fwhm,shift0);
		
		// add all subcomponent
		for(int icomp=0; icomp<4; icomp++) {
			// get a weight, fwhm and shift
			const REAL w = weight[icomp];
			const REAL f = -M_PI*fwhm[icomp]; // (1/A)
			const REAL s = 2*M_PI*shift0[icomp]; // (1/A)
			// calculate Four.coef.
			const REAL *px = x.data();
			REAL *pA = profile.data();
			if (abs(w)>1e-4) { // 'non-virtual' component
				if (abs(f)<1e-7) // unbroadened comp. (possibly there is a shift)
					if (abs(s)<1e-7) // unbroadened comp. with no shift
						for(int i=0; i<nbPoints; i++) { *pA += w; pA += 2; }
					else // unbroadened comp. with a nonzero shift
						for(int i=0; i<nbPoints; i++) { *pA += w*cos(s*(*px)); pA++; *pA += w*sin(s*(*px)); px++; pA++; }
				else // broadened and possibly shifted comp.
					if (abs(s)<1e-7) // broadened comp. with no shift
						for(int i=0; i<nbPoints; i++) { *pA += w*exp(f*abs(*px)); px++; pA += 2; }
					else // broadened comp. with a nonzero shift
						for(int i=0; i<nbPoints; i++) {
							REAL t = w*exp(f*abs(*px));
							*pA += t*cos(s*(*px));
							pA++;
							*pA += t*sin(s*(*px));
							px++; pA++;
						}
			}
		} // icomp
		
		if (bsavecalc && xcenter>=xcenterlimits[0] && xcenter<=xcenterlimits[1]) {
 			CrystMatrix_REAL profiles;
 			CrystVector_REAL shifts;
 			
 			GetAllGroupsProfiles(profiles,shifts,x,xcenter,h,k,l);
 			
	   	ofstream F("profileAF.dat");
			F<<"# alpha="<<mAlpha;
			F<<",(hkl)=("<<h<<k<<l<<"),s0="<<2*sin(0.5*xcenter)/Lambda;
			F<<",approxFWHM="<<GetApproxFWHM(xcenter,h,k,l)*RAD2DEG;
			F<<",xcenter="<<xcenter*RAD2DEG<<endl;
			F<<"# Waren groups: |L0|,sign,count"<<endl;
			for(int igroup=0; igroup<mAbsL0.numElements(); igroup++)
				F<<"# "<<setw(4)<<mAbsL0(igroup)<<setw(3)<<mSign(igroup)<<setw(6)<<mCount(igroup)<<endl;
			F<<"# Balogh&Ungar groups: shifts"<<endl;
			for(int igroup=0; igroup<profiles.rows(); igroup++)
				F<<"# "<<setw(14)<<shifts(igroup)<<endl;
  		for(int i=0;i<nbPoints;i++)
      	//F<<setw(18)<<x(i)<<setw(18)<<profile(2*i)<<setw(18)<<profile(2*i+1)<<endl;
      {
      	F<<setw(18)<<x(i)<<setw(18)<<profile(2*i)<<setw(18)<<profile(2*i+1);
      	for(int igroup=0; igroup<profiles.rows(); igroup++)
      		 F<<setw(18)<<profiles(igroup,i)<<setw(18)<<0.0;
      	F<<endl;
      }
    	F.close(); } // save profile into file block
 
	} // try block

	catch(std::exception &e) {
		cout << "< MStruct:: FaultsBroadeningEffectFCCBaloghUngar::GetProfile(...)\n";
		cout << "exception: " << e.what() << "\n";
		cout << "Exception thrown during calculation of the Fourier Coefficients.\n \
				 Maybe some objects are not properly initialized or some parametrs\n \
				 have wrong values. >" << endl; 
		throw;
	}
	
  return profile;
}

void FaultsBroadeningEffectFCCBaloghUngar::GetAllGroupsProfiles(CrystMatrix_REAL &profiles,
					 																									 CrystVector_REAL &shifts,
																														 const CrystVector_REAL &x,
			      								  					const REAL xcenter, const REAL h, const REAL k, const REAL l)
{
	// nb of points and vector for calc. result
	int nbPoints = x.numElements(); 
	
	shifts.resize(4);
	profiles.resize(4,nbPoints);
	
	// Access to the parent ReflectionProfile and ParentPowderPatternDiffraction
	// objects and their GetRadiation and UnitCell objects is necessary.
	
	// Get radiation and wavelength
	const Radiation &rad = GetParentReflectionProfile().
  			GetParentPowderPatternDiffraction().GetRadiation();
	const REAL Lambda = rad.GetWavelength()(0);

	// calculatation
	
	// get parametrs
	REAL weight[4], fwhm[4], shift0[4];
	GetSubComponentsPar(mAlpha,mFaultsType,h,k,l,weight,fwhm,shift0); // alpha(fraction),fwhm and shift0 (1/A)
	
	// for each component
	for(int icomp=0; icomp<4; icomp++) {
		// set a shift
		shifts(icomp) = shift0[icomp] / (cos(0.5*xcenter)/Lambda); // 2theta(rad)
		// calculate Four.coef.
		const REAL *px = x.data();
		REAL *pA = profiles.data() + icomp*nbPoints;
		if (abs(weight[icomp])<1e-4) // 'virtual' component
			for(int i=0; i<nbPoints; i++) { *pA = 0.; pA++; }
		else if (abs(fwhm[icomp])<1e-7) { // unbroadened comp. (possibly there is a shift)
			const REAL w = weight[icomp];
			for(int i=0; i<nbPoints; i++) { *pA = w; pA++; }
		} else { // broadened and possibly shifted comp.
			const REAL w = weight[icomp];
			const REAL t = - M_PI * fwhm[icomp];	
			for(int i=0; i<nbPoints; i++) { *pA = w*exp(t*abs(*px)); pA++; px++; }
		}
	}
	
}

void FaultsBroadeningEffectFCCBaloghUngar::SetProfilePar(const int type, const REAL alpha)
{
  FaultsBroadeningEffectFCC::SetProfilePar(alpha,0.);
  mFaultsType = type;
  mClockMaster.Click();  
}

void FaultsBroadeningEffectFCCBaloghUngar::GetSubComponentsPar(const REAL alpha, const int type,
																															 const REAL h, const REAL k, const REAL l,
													 																		 REAL weight[4], REAL fwhm[4], REAL shift[4])const
{
	// If (hkl)==(000) is supplied to the method, nothing done.
	if(abs(h)<1e-6 && abs(k)<1e-6 && abs(l)<1e-6) {
		for(int i=0; i<4; i++) { weight[i]=0.; fwhm[i]=0.; shift[i]=0.; }
		weight[0] = 1.;
	}
	
	// Get parameters table for the given type of defects
	const ParametersTable *ptable = NULL;
	switch (type) {
		case INTRINSIC: ptable = &TableIntrinsic; break;
		case TWINS: ptable = &TableTwins; break;
		case EXTRINSIC: ptable = &TableExtrinsic; break;
		default: ptable = NULL;
	}
	
	if(ptable==NULL) {
		cerr << "< MStruct::FaultsBroadeningEffectFCCBaloghUngar::GetSubComponentsPar(...)\n";
		cerr << "\t"<<"type="<<type<<", alpha="<<alpha<<", h="<<h<<", k="<<k<<", l="<<l<<"\n";
		cerr << "\t"<<"Defect type not identified. >" << endl;
		throw ObjCrystException("MStruct::FaultsBroadeningEffectFCCBaloghUngar::GetSubComponentsPar(...): \
             Wrong input argument.");
	}
	
	// Find the given hkl reflection in the table list of hkl reflections
	int ind = -1;
	{
		int H = int(h), K = int(k), L = int(l);
		// sort the indexes
		if(K<L) { int t = K; K = L; L = t; }   
		if(H<K) { int t = K; K = H; H = t; }	
		if(K<L) { int t = K; K = L; L = t; }
		// search the hkl list
		for(int i=0; i<ptable->Nhkl; i++)
			if ((ptable->hkl[i][0]==H) && (ptable->hkl[i][1]==K) && (ptable->hkl[i][2]==L)) { ind = i; break; }
		if(ind<0) {
			cerr << "< MStruct::FaultsBroadeningEffectFCCBaloghUngar::GetSubComponentsPar(...)\n";
			cerr << "\t"<<"type="<<type<<", alpha="<<alpha<<", h="<<h<<", k="<<k<<", l="<<l<<"\n";
			cerr << "\t"<<"Reflection not found in the \"Levente Balogh\" parameters table. >" << endl;
			throw ObjCrystException("MStruct::FaultsBroadeningEffectFCCBaloghUngar::GetSubComponentsPar(...): \
             Wrong input argument.");
		}
	}
	
	// Get weigths, fwhms and shifts for all(4) reflection subcomponents
	for(int i=0; i<4; i++) {
		// weight
		weight[i] = ptable->weight[ind][i];
		// fwhm
		fwhm[i] = ((((ptable->fwhm[i][ind][0]*alpha + ptable->fwhm[i][ind][1])*alpha + ptable->fwhm[i][ind][2])*alpha
							+ ptable->fwhm[i][ind][3])*alpha + ptable->fwhm[i][ind][4])*alpha + ptable->fwhm[i][ind][5];
		// shift
		shift[i] = ((((ptable->shift[i][ind][0]*alpha + ptable->shift[i][ind][1])*alpha + ptable->shift[i][ind][2])*alpha
							+ ptable->shift[i][ind][3])*alpha + ptable->shift[i][ind][4])*alpha + ptable->shift[i][ind][5];
	}
}

// PseudoVoigtBroadeningEffect

PseudoVoigtBroadeningEffect::PseudoVoigtBroadeningEffect():
mAccuracy(0.005),mAsymXMax(60.*DEG2RAD)
{
  /*
  {
    RefinablePar tmp("param", &mParam, 5., 3.e3,
                     gpRefParTypeScattDataProfileWidth,
                     REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false,0.1);
    tmp.AssignClock(mClockMaster);
    tmp.SetDerivStep(1.);
    this->AddPar(tmp);
  }
  RefinableObj::Print();
  */
}

CrystVector_REAL PseudoVoigtBroadeningEffect::GetProfile(
							 const CrystVector_REAL &x,
							 const REAL xcenter,
							 const REAL h, const REAL k, const REAL l)
{
  // TODO:: should be faster and more simple, integrated intensity
  // TODO:: improve self determination calc. range and optimal step

  // prepare for calculation

  // reflection profile pseudo-Voigt object for internal calculations
  // (to prevent recursive calling of working subroutines)
  ReflectionProfilePseudoVoigt reflprof = ReflectionProfilePseudoVoigt(*this);
	/*if(xcenter>=mAsymXMax) {
		reflprof.GetPar("Asym0").SetValue(1.);
		reflprof.GetPar("Asym1").SetValue(0.);
		reflprof.GetPar("Asym2").SetValue(0.);
	}*/
	
	// maybe there are limits for asymmmetry parametrs, hence I prefer to set it manually
	const REAL asym = (xcenter<mAsymXMax) ? 
         						this->GetPar("Asym0").GetValue()
         					 +this->GetPar("Asym1").GetValue()/sin(xcenter)
         					 +this->GetPar("Asym2").GetValue()/pow(sin(xcenter),(REAL)2.) : 1.;
  
  reflprof.GetPar("Asym0").SetValue(asym);
	reflprof.GetPar("Asym1").SetValue(0.);
	reflprof.GetPar("Asym2").SetValue(0.);
        						
  // get prescribed nb. of points in which the profile should be evaluated
  int nbPoints = x.numElements();

  // get radiation and wavelength
  const Radiation &rad = GetParentReflectionProfile().
    GetParentPowderPatternDiffraction().GetRadiation();
  REAL Lambda = rad.GetWavelength()(0);

  REAL factor = cos(xcenter/2)/Lambda;

  // estimate calculation range and optimal step

  // calc range for the given accuracy 0.04
  REAL tthwidth0 = 2*reflprof.GetFullProfileWidth(0.01,xcenter,h,k,l);
  
  // calc range for the given accuracy (mAccuracy)
  REAL tthwidth = 2*reflprof.GetFullProfileWidth(mAccuracy,xcenter,h,k,l);

  // we calculate at least in the prescribed limits
  tthwidth = max(tthwidth,(x(nbPoints-1)-x(0))/factor); // not a good idea (why?)
      // I need it in the case that this is the only one effect,
      // increasing range to achiave presciribed accuracy do not
      // converge in the other case
  
  // if it is suitable to use more than prescribed number of points than
  // increase nbPoints
  nbPoints = (tthwidth0<tthwidth) ? int(tthwidth/tthwidth0*nbPoints+1) : nbPoints;

	// if there is less than approx 10 points per FWHM than increase nb. of points
	const REAL fwhm = this->GetApproxFWHM(xcenter,h,k,l);
	nbPoints = (fwhm/10<tthwidth/(nbPoints-1)) ? int(tthwidth/(fwhm/10)+1) : nbPoints;

  // create a vector of tth-data where profile will be evaluated
  // data will be regulary spaced in 2theta

  CrystVector_REAL tth(nbPoints);

  REAL tthstep = tthwidth/(nbPoints-1);

  {
    REAL val = xcenter - tthwidth/2;
    REAL *p1 = tth.data();
    for(int i=0;i<nbPoints;i++) { *p1 = val; p1++; val += tthstep; }
  }

  // calc profile
  CrystVector_REAL profile = reflprof.GetProfile(tth,xcenter,h,k,l);

  // result will contain at first x-data and after them y-data (intensities)
  CrystVector_REAL result(2*nbPoints);

  {
    // calc values of x-data points in the s-space
    
    REAL *p1 = result.data();
    const REAL *p2 = tth.data();

    for(int i=0;i<nbPoints;i++) {
      *p1 = factor * (*p2-xcenter); // ds = cos(theta)/Lambda * d(2theta)
      p1++; p2++;
    }

    // calc values of y-data points (intensiies) in the s-space

    p2 = profile.data();

    for(int i=0;i<nbPoints;i++) {
      *p1 = *p2 / factor; // dI(s) = (cos(theta)/Lambda)^-1 * dI(2theta)
      p1++; p2++;
    }
  }

  if (bsavecalc && xcenter>=xcenterlimits[0] && xcenter<=xcenterlimits[1]) {
    ofstream F("profilePpv.dat");
    F<<"# xcenter="<<xcenter*RAD2DEG<<endl;
    for(int i=0;i<nbPoints;i++)
      F<<setw(18)<<result(i)<<setw(18)<<result(nbPoints+i)<<endl;
    F.close();
  }

  return result;
}

REAL PseudoVoigtBroadeningEffect::GetApproxFWHM(const REAL xcenter,
				const REAL h, const REAL k, const REAL l)const
{
  // TODO:: should be faster
  // possible solution: if I can acces PV profile params I can calc
  // fwhm directly
  //const REAL fwhm=sqrt( mCagliotiW
  //                      +mCagliotiV*tan(center/2.0)
  //                      +mCagliotiU*pow(tan(center/2.0),2));

  // reflection profile pseudo-Voigt object for internal calculations
  // (to prevent recursive calling of working subroutines)
  ReflectionProfilePseudoVoigt reflprof = ReflectionProfilePseudoVoigt(*this);

	return reflprof.GetFullProfileWidth(0.5,xcenter,h,k,l);
}

bool PseudoVoigtBroadeningEffect::IsRealSpaceType()const {
  return false;
}

void PseudoVoigtBroadeningEffect::SetAsymXMax(const REAL asymXMax)
{
	mAsymXMax = asymXMax;
	mClockMaster.Click();
}

// PseudoVoigtBroadeningEffectA
PseudoVoigtBroadeningEffectA::PseudoVoigtBroadeningEffectA():
mCagliotiU(0.), mCagliotiV(0.), mCagliotiW(0.09*DEG2RAD*DEG2RAD),
mPseudoVoigtEta0(0.5), mPseudoVoigtEta1(0.), mAsym0(1.), mAsym1(0.),
mAsym2(0.),mAsymXMax(60.*DEG2RAD)
{
  InitParameters();
}

CrystVector_REAL PseudoVoigtBroadeningEffectA::GetProfile(
				     const CrystVector_REAL &x,
				     const REAL xcenter,
				     const REAL h, const REAL k, const REAL l)
{
  // calc asymetry
  const REAL asym = (xcenter<mAsymXMax) ? 
         mAsym0+mAsym1/sin(xcenter)+mAsym2/pow(sin(xcenter),(REAL)2.) : 1.;

  const bool is_asym = fabs(asym-1.)>1.e-4; // ??? unclear !!!

  const int nbPoints = x.numElements(); 
  CrystVector_REAL profile((is_asym==true) ? 2*nbPoints : nbPoints);
 
  // get radiation and wavelength
  const Radiation &rad = GetParentReflectionProfile().
    GetParentPowderPatternDiffraction().GetRadiation();
  const REAL Lambda = rad.GetWavelength()(0);

  const REAL factor = cos(xcenter/2)/Lambda;
  
  // calc Four. coefs
  // ref: P.Scardi,L.Matteo,J.Appl.Cryst.(1999).32,671-682:Fourier modelling...
  // size effect

  const REAL *p1 = x.data();
  REAL *p2 = profile.data();

  const REAL hwhm=sqrt( mCagliotiW
                       +mCagliotiV*tan(xcenter/2.0) // ds = cos(theta)/Lambda * d(2theta)
                       +mCagliotiU*pow(tan(xcenter/2.0),2))/2. * factor;
  const REAL factorPhi = 1./sqrt(M_PI*M_LN2); // ? should I use this ???
  const REAL ni = mPseudoVoigtEta0+mPseudoVoigtEta1*xcenter;
  const REAL w = (abs(ni)<FLT_EPSILON) ? 0. : 1./(1.+factorPhi*(1.-ni)/ni);

  const REAL factorGsq = M_PI*M_PI*hwhm*hwhm/log(2.);
  const REAL factorC   = 2.*M_PI*hwhm;

  if(is_asym==false) {

    for(int i=0;i<nbPoints;i++) {
      double L = *p1++; L = abs(L);
      *p2++ = REAL((1.-w)*exp(-factorGsq*L*L)+w*exp(-factorC*L));
    }

  }
  else { // is asymmetric
    
    // left side
    const REAL c1 = 2.*asym/(1.+asym);
    const REAL factorG1 = sqrt(factorGsq)*c1;
    const REAL factorG1sq = factorG1*factorG1;
    const REAL factorC1 = factorC*c1;

    // right side
    const REAL c2 = 2./(1.+asym);
    const REAL factorG2 = sqrt(factorGsq)*c2;
    const REAL factorG2sq = factorG2*factorG2;
    const REAL factorC2 = factorC*c2;
    
    for(int i=0;i<nbPoints;i++) {
      double L = *p1++;
      // real part
      // left side
      *p2  = REAL((1.-w)*exp(-factorG1sq*L*L)+w*exp(-factorC1*fabs(L)))*c1;
      // right side
      *p2 += REAL((1.-w)*exp(-factorG2sq*L*L)+w*exp(-factorC2*fabs(L)))*c2;

      p2++;

      // imaginary part
      double arg1, arg2,arg3;
      // left side
       //arg1 = factorG1*L*2.;
			arg1 = factorG1*L;
      arg2 = factorC1*L;
      arg3 = -arg2;
      //*p2  = REAL((1.-w)*daw_(&arg1)+w*(expei_(&arg2)-expei_(&arg3))/M_PI)*c1;
       //*p2  = REAL((1.-w)*func_daw(arg1)+w*((fabs(arg2)<1.e-5) ? 0. : (exp(-arg2)*func_ei(arg2)-exp(-arg3)*func_ei(arg3)))/M_PI)*c1;
      *p2  = -REAL((1.-w)*M_2_SQRTPI*func_daw(arg1)+w*((fabs(arg2)<1.e-5) ? 0. : (exp(-arg2)*func_ei(arg2)-exp(-arg3)*func_ei(arg3)))/M_PI)*c1;
			// right side
       //arg1 = factorG2*L*2.;
			arg1 = factorG2*L;
      arg2 = factorC2*L;
      arg3 = -arg2;
      //*p2 -= REAL((1.-w)*daw_(&arg1)+w*(expei_(&arg2)-expei_(&arg3))/M_PI)*c2;
       //*p2 -= REAL((1.-w)*func_daw(arg1)+w*((fabs(arg2)<1.e-5) ? 0. : (exp(-arg2)*func_ei(arg2)-exp(-arg3)*func_ei(arg3)))/M_PI)*c2;
			*p2 += REAL((1.-w)*M_2_SQRTPI*func_daw(arg1)+w*((fabs(arg2)<1.e-5) ? 0. : (exp(-arg2)*func_ei(arg2)-exp(-arg3)*func_ei(arg3)))/M_PI)*c2;

      p2++;
    }   
    
    profile *= 0.5;

  }

  if (bsavecalc && xcenter>=xcenterlimits[0] && xcenter<=xcenterlimits[1]) {
    ofstream F("profileAG.dat");
    const REAL sc = RAD2DEG*RAD2DEG;
    F<<"# W="<<mCagliotiW*sc<<",U="<<mCagliotiU*sc<<",V="<<mCagliotiV*sc;
    F<<",Eta0="<<mPseudoVoigtEta0<<",Eta1="<<mPseudoVoigtEta1;
    F<<",xcenter="<<xcenter*RAD2DEG<<endl;
    if(is_asym==false) {
      for(int i=0;i<nbPoints;i++)
	F<<setw(18)<<x(i)<<setw(18)<<profile(i)<<endl;
    }
    else {
      for(int i=0;i<nbPoints;i++)
	F<<setw(18)<<x(i)
	           <<setw(18)<<profile(2*i)<<setw(18)<<profile(2*i+1)<<endl;
    }
    F.close();
  }

  return profile;
}

REAL PseudoVoigtBroadeningEffectA::GetApproxFWHM(const REAL xcenter,
			 const REAL h, const REAL k, const REAL l)const
{
  return sqrt( mCagliotiW
              +mCagliotiV*tan(xcenter/2.0)
              +mCagliotiU*pow(tan(xcenter/2.0),2));
}

bool PseudoVoigtBroadeningEffectA::IsRealSpaceType()const {
  return true;
}

void PseudoVoigtBroadeningEffectA::SetProfilePar(const REAL fwhmCagliotiW,
						 const REAL fwhmCagliotiU,
						 const REAL fwhmCagliotiV,
						 const REAL eta0,
						 const REAL eta1)
{
   mCagliotiU=fwhmCagliotiU;
   mCagliotiV=fwhmCagliotiV;
   mCagliotiW=fwhmCagliotiW;
   mPseudoVoigtEta0=eta0;
   mPseudoVoigtEta1=eta1;
   mClockMaster.Click();
}

void PseudoVoigtBroadeningEffectA::SetAsymXMax(const REAL asymXMax)
{
	mAsymXMax = asymXMax;
	mClockMaster.Click();
}

void PseudoVoigtBroadeningEffectA::InitParameters()
{
  {
    RefinablePar tmp("U",&mCagliotiU,-1/RAD2DEG/RAD2DEG,1./RAD2DEG/RAD2DEG,
		     gpRefParTypeScattDataProfileWidth,
		     REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false,
		     RAD2DEG*RAD2DEG);
    tmp.AssignClock(mClockMaster);
    tmp.SetDerivStep(1e-7);
    this->AddPar(tmp);
  }
  {
    RefinablePar tmp("V",&mCagliotiV,-1/RAD2DEG/RAD2DEG,1./RAD2DEG/RAD2DEG,
		     gpRefParTypeScattDataProfileWidth,
		     REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false,
		     RAD2DEG*RAD2DEG);
    tmp.AssignClock(mClockMaster);
    tmp.SetDerivStep(1e-7);
    this->AddPar(tmp);
  }
  {
    RefinablePar tmp("W",&mCagliotiW,0,1./RAD2DEG/RAD2DEG,
		     gpRefParTypeScattDataProfileWidth,
		     REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false,
		     RAD2DEG*RAD2DEG);
    tmp.AssignClock(mClockMaster);
    tmp.SetDerivStep(1e-7);
    this->AddPar(tmp);
  }
  {
    RefinablePar tmp("Eta0",&mPseudoVoigtEta0,0,1.,
		     gpRefParTypeScattDataProfileType,
		     REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false);
    tmp.AssignClock(mClockMaster);
    tmp.SetDerivStep(1e-4);
    this->AddPar(tmp);
  }
  {
    RefinablePar tmp("Eta1",&mPseudoVoigtEta1,-1,1.,
		     gpRefParTypeScattDataProfileType,
		     REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false);
    tmp.AssignClock(mClockMaster);
    tmp.SetDerivStep(1e-4);
    this->AddPar(tmp);
  }
  {
    RefinablePar tmp("Asym0",&mAsym0,0.01,10.,
		     gpRefParTypeScattDataProfileAsym,
		     REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false);
    tmp.AssignClock(mClockMaster);
    tmp.SetDerivStep(1e-4);
    this->AddPar(tmp);
  }
  {
    RefinablePar tmp("Asym1",&mAsym1,-1.0,1.0,
		     gpRefParTypeScattDataProfileAsym,
		     REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false);
    tmp.AssignClock(mClockMaster);
    tmp.SetDerivStep(1e-4);
    this->AddPar(tmp);
  }
  {
    RefinablePar tmp("Asym2",&mAsym2,-1.0,1.0,
		     gpRefParTypeScattDataProfileAsym,
		     REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false);
    tmp.AssignClock(mClockMaster);
    tmp.SetDerivStep(1e-4);
    this->AddPar(tmp);
  }
}

// HKLPseudoVoigtBroadeningEffectA
HKLPseudoVoigtBroadeningEffectA::HKLPseudoVoigtBroadeningEffectA()
{}

HKLPseudoVoigtBroadeningEffectA::~HKLPseudoVoigtBroadeningEffectA()
{
  for(int i=0;i<mReflStore.size();i++) {
    HKLProfilePar *pData=(HKLProfilePar*)mReflStore.at(i).data;
    this->RemovePar(&(this->GetPar(&(pData->dx))));
    this->RemovePar(&(this->GetPar(&(pData->fwhm))));
    this->RemovePar(&(this->GetPar(&(pData->eta))));
    delete pData;
  }
  mReflStore.clear();
}

CrystVector_REAL HKLPseudoVoigtBroadeningEffectA::GetProfile(
				     const CrystVector_REAL &x,
				     const REAL xcenter,
				     const REAL h, const REAL k, const REAL l)
{
  int nbPoints = x.numElements(); 
  CrystVector_REAL profile(nbPoints);
 
  const REAL *p1 = x.data();
  REAL *p2 = profile.data();

  // try to find (hkl) reflection in the "store"
  int ind = mReflStore.find((int)h,(int)k,(int)l,xcenter);
   
  // if refl. found use stored params else return delta function
  if(ind>=0) {
  
    const HKLProfilePar *pData= (HKLProfilePar*) mReflStore.at(ind).data;

    // get radiation and wavelength
    const Radiation &rad = GetParentReflectionProfile().
      GetParentPowderPatternDiffraction().GetRadiation();
    const REAL Lambda = rad.GetWavelength()(0);

    const REAL factor = cos(xcenter/2)/Lambda;
  
    // calc Four. coefs
    // ref: P.Scardi,L.Matteo,J.Appl.Cryst.(1999).32,671-682:Fourier modelling...

    const REAL hwhm= pData->fwhm/2.* factor;// ds = cos(theta)/Lambda * d(2theta)
    const REAL factorPhi = 1./sqrt(M_PI*M_LN2);
    const REAL ni = pData->eta;
    const REAL w = (abs(ni)<FLT_EPSILON) ? 0. : 1./(1.+factorPhi*(1.-ni)/ni);

    const REAL factorG = M_PI*M_PI*hwhm*hwhm/log(2.);
    const REAL factorC = 2.*M_PI*hwhm;

    for(int i=0;i<nbPoints;i++) {
      double L = *p1++; L = abs(L);
      *p2++ = REAL((1.-w)*exp(-factorG*L*L)+w*exp(-factorC*L)); 
    }

  } else {
    // delta function
    for(int i=0;i<nbPoints;i++) {*p2++ = 1.;p1++;}
  }


  if (bsavecalc && xcenter>=xcenterlimits[0] && xcenter<=xcenterlimits[1]) {
    ofstream F("profileAGhkl.dat");
    F<<"# dx="<<0.<<",fwhm="<<0.<<",eta="<<0.;
    F<<",xcenter="<<xcenter*RAD2DEG<<",(HKL)=("<<h<<" "<<k<<" "<<l<<")"<<endl;
    for(int i=0;i<nbPoints;i++)
      F<<setw(18)<<x(i)<<setw(18)<<profile(i)<<endl;
    F.close();
  }

  return profile;
}

REAL HKLPseudoVoigtBroadeningEffectA::GetApproxFWHM(const REAL xcenter,
			 const REAL h, const REAL k, const REAL l)const
{
  // try to find (hkl) reflection in the "store"
  int ind = mReflStore.find((int)h,(int)k,(int)l,xcenter);
  
  REAL fwhm = 0.;

  // if refl. found use stored params else return 0. for delta function
  if(ind>=0) {
    const HKLProfilePar *pData= (HKLProfilePar*) mReflStore.at(ind).data;
    fwhm= pData->fwhm;
  }

  return fwhm;
}

bool HKLPseudoVoigtBroadeningEffectA::IsRealSpaceType()const {
  return true;
}

void HKLPseudoVoigtBroadeningEffectA::SetProfilePar(int h,int k,int l,
						    REAL dx,REAL fwhm,REAL eta,
						    bool dx_fixed,
						    bool fwhm_fixed,
						    bool eta_fixed)
{
  // try to find (hkl) reflection in the "store"
  int ind = mReflStore.find(h,k,l,0.);

  const ObjCryst::PowderPatternDiffraction *pDiffData = &(mpParentReflectionProfile->
							 GetParentPowderPatternDiffraction());

  // generate parameter name
  string name;
  {
    stringstream s;
    s << pDiffData->GetName() << "_Phkl_";
    s << h << "_" << k << "_" << l;
    name = s.str();
  }
   
  // if refl. found use it else add new refl. into the "store"
  HKLProfilePar *pData = 0;
  if(ind>=0) {
    pData = (HKLProfilePar*) mReflStore.at(ind).data;
  }
  else {
    pData = new HKLProfilePar;
    mReflStore.add(h,k,l,0.,(void*)pData);

    {
    RefinablePar tmp(string(name+"_dx").c_str(),&(pData->dx),-5.*DEG2RAD,5.*DEG2RAD,
		     gpRefParTypeScattDataProfileWidth, // TODO:: change this
		     REFPAR_DERIV_STEP_ABSOLUTE,
		     true,false,false,false,RAD2DEG);
    tmp.AssignClock(mClockMaster);
    tmp.SetDerivStep(0.001);
    this->AddPar(tmp);
    }

    {
    RefinablePar tmp(string(name+"_fwhm").c_str(),&(pData->fwhm),0.*DEG2RAD,10.*DEG2RAD,
		     gpRefParTypeScattDataProfileWidth, // TODO:: change this
		     REFPAR_DERIV_STEP_ABSOLUTE,
		     true,false,false,false,RAD2DEG);
    tmp.AssignClock(mClockMaster);
    tmp.SetDerivStep(0.001);
    this->AddPar(tmp);
    }

    {
    RefinablePar tmp(string(name+"_eta").c_str(),&(pData->eta),0.,1.,
		     gpRefParTypeScattDataProfileWidth, // TODO:: change this
		     REFPAR_DERIV_STEP_ABSOLUTE,
		     true,false,false,false,1.);
    tmp.AssignClock(mClockMaster);
    tmp.SetDerivStep(0.05);
    this->AddPar(tmp);
    }
  }
  
  // Is this refl currently used by ScatteringData object?
  bool isused = false;
  for(int i=0;i<pDiffData->GetNbRefl();i++)
    if(int(pDiffData->GetH()(i))==h && int(pDiffData->GetK()(i))==k &&
       int(pDiffData->GetL()(i))==l) {
      isused = true;
      break;
    }

  // set params
  {
  RefinablePar &par = GetPar(string(name+"_dx"));
  par.SetValue(dx);
  par.SetIsFixed(dx_fixed);
  par.SetIsUsed(isused);
  }

  {
  RefinablePar &par = GetPar(string(name+"_fwhm"));
  par.SetValue(fwhm);
  par.SetIsFixed(fwhm_fixed);
  par.SetIsUsed(isused);
  }

  {
  RefinablePar &par = GetPar(string(name+"_eta"));
  par.SetValue(eta);
  par.SetIsFixed(eta_fixed);
  par.SetIsUsed(isused);
  }

  mClockMaster.Click();
}

REAL HKLPseudoVoigtBroadeningEffectA::GetPositionCorr(const REAL xcenter,
				     const REAL h, const REAL k, const REAL l)const
{ 
  // try to find (hkl) reflection in the "store"
  int ind = mReflStore.find((int)h,(int)k,(int)l,xcenter);
   
  REAL xcorr = 0.;

  // if refl. found use stored params else return 0.
  if(ind>=0) {
    const HKLProfilePar *pData = (HKLProfilePar*) mReflStore.at(ind).data;
    xcorr += pData->dx;
  }

  return xcorr;
}

// ReflectionProfile
ReflectionProfile::ReflectionProfile (const UnitCell& cell,
				      const Radiation& radiation):
mpUnitCell(&cell), mpRadiation(&radiation), mOmega(-1.), mFactor(1.), // 9
mN(1024),mNbPoints(0),mvL(1024),mvs(1024),mpParentPowderPatternDiffraction(0)
{
  this->mClockMaster.AddChild( mClockStressCorrQ );
  mReflectionProfileComponentRegistry.SetName("ReflectionProfile Components");
  mStressCoeff.resize(5);
  mStressCoeff = 0.; mStressCoeff(1) = -1.;
  mvIn.reserve(1024);
  mvOut.reserve(1024);
  mvProfile.reserve(1024);
  InitParameters();
}

ReflectionProfile::ReflectionProfile(const ReflectionProfile &old):
mpUnitCell(old.mpUnitCell),
mpRadiation(old.mpRadiation), mOmega(old.mOmega), mStressCoeff(old.mStressCoeff),
mFactor(old.mFactor), mN(old.mN),
mvIn(old.mvIn), mvOut(old.mvOut), mdL(old.mdL), mds(old.mds), mQ(old.mQ),
ms0(old.ms0), ms1(old.ms1), mLambda(old.mLambda), mNbPoints(old.mNbPoints),
mvL(old.mvL),
mvs(old.mvs),
mReflStore(old.mReflStore),
mReflectionProfileComponentRegistry(old.mReflectionProfileComponentRegistry),
mpParentPowderPatternDiffraction(old.mpParentPowderPatternDiffraction)
{
  this->mClockMaster.AddChild( mClockStressCorrQ );
  InitParameters();
}
  
ReflectionProfile* ReflectionProfile::CreateCopy()const
{
   return new ReflectionProfile(*this);
}

ReflectionProfile::~ReflectionProfile()
{
  for(int i=0;i<mReflectionProfileComponentRegistry.GetNb();i++)
  {
    mReflectionProfileComponentRegistry.GetObj(i).DeRegisterClient(*this);
    this->RemoveSubRefObj(mReflectionProfileComponentRegistry.GetObj(i));
    delete &(mReflectionProfileComponentRegistry.GetObj(i)); // !!!! BUG
          // ----> no delete - it should be deleted by a creater of the component object
  }
  for(int i=0;i<mReflStore.size();i++) {
    ReflCalcData *pData = (ReflCalcData*) mReflStore.at(i).data;
    delete pData;
  }
  mReflStore.clear();
}

const string& ReflectionProfile::GetClassName () const
{
  const static string className="MStruct::ReflectionProfile";
  return className;
}

void ReflectionProfile::AddReflectionProfileComponent(ReflectionProfileComponent &comp)
{
  VFN_DEBUG_ENTRY("ReflectionProfile::AddReflectionProfileComponent():"<<comp.GetName(),11)
  //comp.SetParentReflectionProfile(*this);
  this->AddSubRefObj(comp);
  comp.RegisterClient(*this);
  mClockMaster.AddChild(comp.GetClockMaster());
  mClockReflectionProfileCalc.Reset();
  mReflectionProfileComponentRegistry.Register(comp);

  this->UpdateDisplay();
  VFN_DEBUG_EXIT("ReflectionProfile::AddReflectionProfileComponent():"<<comp.GetName(),11)
}

long ReflectionProfile::GeReflectionProfileComponentNb() const
{ return mReflectionProfileComponentRegistry.GetNb(); }
const ReflectionProfileComponent &  ReflectionProfile::GetReflectionProfileComponent (long n) const
{ return mReflectionProfileComponentRegistry.GetObj(n); }
ReflectionProfileComponent &  ReflectionProfile::GetReflectionProfileComponent (long n)
{ return mReflectionProfileComponentRegistry.GetObj(n); }
const ReflectionProfileComponent &  ReflectionProfile::GetReflectionProfileComponent (const string &objName) const
{ return mReflectionProfileComponentRegistry.GetObj(objName); }
ReflectionProfileComponent &  ReflectionProfile::GetReflectionProfileComponent (const string &objName)
{ return mReflectionProfileComponentRegistry.GetObj(objName); }

bool ReflectionProfile::PrepareForCalc(const CrystVector_REAL &x,
				       const REAL xcenter, 
				       REAL h, REAL k, REAL l)
{
  VFN_DEBUG_ENTRY("ReflectionProfile::PrepareForCalc(): xcenter:"
		  <<xcenter*RAD2DEG<<", h,k,l:"<<h<<","<<k<<","<<l,11)
  mNbPoints = x.numElements();

  // Wavelength of the experiment (in Angstroems)
  //mLambda = mpRadiation->GetWavelength()(0); // ???
  mLambda = mpParentPowderPatternDiffraction->GetRadiation().GetWavelength()(0); // ???
	
  ms0 = 2*sin(xcenter/2)/mLambda;

  // range of calc profile in the reciprocal units
  ms1 = max(abs(2*sin(x(0)/2)/mLambda-ms0),
	    abs(2*sin(x(mNbPoints-1)/2)/mLambda-ms0));
  
  // step in the direct space
  mdL = 0.5/ms1/sqrt(mFactor);

  // nb of points used for FT
  mN = 2*int(mFactor*(mNbPoints-1)/2.+1);

  // using complex DFT - from FFTW lib
  // step in the direct space ds=1/dL/n; s=k*ds; L=j*dL; k,j=-n/2+1,..,0,..n/2
  
  // step in the reciprocal space
  mds = 1.0/mdL/mN;

	// now it will be checked if there are at least 10 points per line fwhm,
	// if not, the calculation grid will be changed
	REAL fwhm = 0.;
	for(int icomp=0;icomp<mReflectionProfileComponentRegistry.GetNb();icomp++) {
    const ReflectionProfileComponent &comp = mReflectionProfileComponentRegistry.GetObj(icomp);
		fwhm += comp.GetApproxFWHM(xcenter,h,k,l);
	}
	fwhm *= cos(xcenter/2)/mLambda;
	if(fwhm>std::numeric_limits<REAL>::epsilon() && fwhm/mds<10) {
		// new value of mN and ds
		mN = 2*int(mds/(fwhm/10)*mN/2);
		mds = 1.0/mdL/mN;
	}

  // if something has been modified clear the "store"
  if(mClockReflectionProfileCalc<mClockMaster) {
    for(int i=0;i<mReflStore.size();i++) {
      ReflCalcData *pData= (ReflCalcData*)mReflStore.at(i).data;
      delete pData;
    }
    mReflStore.clear();
    VFN_DEBUG_MESSAGE("MStruct::ReflectionProfile::PrepareForCalc(): "<<
		      "ReflStore CLEARED - ALL reflections will be "<<
		      "RECALCULATED.",11)
  }
  else {
    VFN_DEBUG_MESSAGE("MStruct::ReflectionProfile::PrepareForCalc(): "<<
		      "Nothing modified - reflection will loaded "<<
		      "or calculated.",11)
  }

  // try to find reflection in the "store"
  int ind = mReflStore.find((int)h,(int)k,(int)l,xcenter);
  if(ind>=0) { // reflection found
    // get data for the found reflection
    ReflCalcData *pData = (ReflCalcData*) mReflStore.at(ind).data;
    // does the center,range or step changed significantly?
    if(fabs((ms0-pData->s0)/ms0)<0.04 && fabs((ms1-pData->s1)/ms1)<0.04 && 
       fabs((mds-pData->ds)/mds)<0.04) {
      // use already calc data as an aproximation (prepare data)
      mds = pData->ds;
      mN  = pData->N;
      mvProfile = pData->vProfile;
      //if ((int)mvProfile.capacity()<mN) mvProfile.reserve(mN);
	  if ((int)mvProfile.size()<mN) mvProfile.resize(mN);
      double *p1=&(pData->vProfile[0]), *p2=&(mvProfile[0]);
      for(int i=0;i<mN;i++,p1++,p2++) *p2 = *p1;
      VFN_DEBUG_MESSAGE("ReflectionProfile::PrepareForCalc()"<<
			"Diffraction found-will be loaded.",11)
      VFN_DEBUG_EXIT("ReflectionProfile::PrepareForCalc()():false",11)
      // we dont't need to calc profile (just approx)
      return false;
    } else { // reflection found but we don't want to use it
      // remove data from the "store"
      // one exception acceptable - if x of the found reflection 
      // is close to the first or the last point of the pattern
      // than maybe it is a testing reflection and such reflections
      // are not removed
      REAL xmin = mpParentPowderPatternDiffraction->
	              GetParentPowderPattern().GetPowderPatternXMin();
      REAL xmax = mpParentPowderPatternDiffraction->
	              GetParentPowderPattern().GetPowderPatternXMax();
      if(fabs((mReflStore.at(ind).x-xmin)/tan(xmin))>=0.04 &&
	 fabs((mReflStore.at(ind).x-xmax)/tan(xmax))>=0.04) {
	VFN_DEBUG_MESSAGE("ReflectionProfile::PrepareForCalc()"<<
			  "Diffraction found-will be removed.",11)
	mReflStore.erase(ind);
	delete pData;
      } else
	VFN_DEBUG_MESSAGE("ReflectionProfile::PrepareForCalc()"<<
			  "Diffraction found-identified as an exception.",11);
    }
  }
  
  // continue normal computation

  // prepare memory
  //if ((int)mvIn.capacity()<mN) mvIn.resize(mN);
  //if ((int)mvOut.capacity()<mN) mvOut.reserve(mN);
  //if ((int)mvProfile.capacity()<mN) mvProfile.reserve(mN);
  if ((int)mvIn.size()<mN) mvIn.resize(mN);
  if ((int)mvOut.size()<mN) mvOut.resize(mN);
  if ((int)mvProfile.size()<mN) mvProfile.resize(mN);

  // prepare real space grid (sometimes not necessary)
  mvL.resize(mN);
  REAL *p1=mvL.data();
  for(int i=0;i<=mN/2;i++) *p1++ = i*mdL;
  for(int i=mN/2-1;i>=1;i--) *p1++ = -i*mdL;
  // prepare reciprocal space grid (sometimes not necessary)
  mvs.resize(mN);
  p1=mvs.data();
  for(int i=0;i<=mN/2;i++) *p1++ = i*mds;
  for(int i=mN/2-1;i>=1;i--) *p1++ = -i*mds;

  // some usefull constants
  mpUnitCell->MillerToOrthonormalCoords(h,k,l);
  mQ=sqrt(h*h+k*k+l*l);

  VFN_DEBUG_EXIT("ReflectionProfile::PrepareForCalc()():true",11)
  return true;
}

void ReflectionProfile::SetParentPowderPatternDiffraction
                     (const ObjCryst::PowderPatternDiffraction& s) {
  //if(mpParentPowderPatternDiffraction!=0) 
  //    mClockMaster.RemoveChild(mpParentPowderPattern->GetIntegratedProfileLimitsClock());
   mpParentPowderPatternDiffraction = &s;
   if ( GetName().empty() )
   	 SetName(string("reflProf_")+s.GetCrystal().GetName());
   mReflStore.clear();
   //mClockMaster.AddChild(mpParentPowderPattern->GetIntegratedProfileLimitsClock());
   //mClockMaster.AddChild(mpParentPowderPattern->GetClockPowderPatternPar());
   //mClockMaster.AddChild(mpParentPowderPattern->GetClockPowderPatternXCorr());
   //mClockMaster.AddChild(mpParentPowderPattern->GetClockPowderPatternRadiation());
}

const ObjCryst::PowderPatternDiffraction& ReflectionProfile::GetParentPowderPatternDiffraction()const {
  return *mpParentPowderPatternDiffraction;
}


//void ReflectionProfile::SetProfilePar(const REAL m, const REAL sigma)
//{}

bool ReflectionProfile::IsAnisotropic()const { return false; }

CrystVector_REAL ReflectionProfile::GetProfile (const CrystVector_REAL &x,
							const REAL xcenter,
							const REAL H,
							const REAL K,
							const REAL L)
{ 
  bool needRecalc = PrepareForCalc(x,xcenter,H,K,L);

  REAL factor = cos(xcenter/2)/mLambda;

  if(needRecalc) {
	
  /* calculation */

  // no efect (TODO:: not necessary)
  for(int i=0;i<mN;i++) { mvIn[i] = 1.; mvOut[i] = 1.; }

  // appprox. FWHM of all included effects
  REAL fwhm = 0.;

  // add effects
  //TODO::can be faster for symmetric profiles
  for(int k=0;k<mReflectionProfileComponentRegistry.GetNb();k++) {
    ReflectionProfileComponent &comp = mReflectionProfileComponentRegistry.GetObj(k);
    if (!comp.IsRealSpaceType()) continue;
    fwhm += comp.GetApproxFWHM(xcenter,H,K,L);
    CrystVector_REAL profile=comp.GetProfile(mvL,xcenter,H,K,L);
    const REAL *p=profile.data();
    if (profile.numElements()==mN) // Four. coef. are real
      for(int i=0;i<mN;i++) mvIn[i] *= (double) *p++;
     else if (profile.numElements()==2*mN) // Four. coef are complex
       for(int i=0;i<mN;i++) {
	 mvIn[i] *= complex<double>((double)p[0],(double)p[1]);
	 p += 2;
       }
    else // program error - throw exception
      throw ObjCrystException("MStruct::ReflectionProfile::GetProfile(...)\
             Wrong nb. of calculated profile points.");
  }
  
/// ------ 
  if (bsavecalc && xcenter>=xcenterlimits[0] && xcenter<=xcenterlimits[1]) {
    ofstream F("profileA.dat");
    F<<"# xcenter="<<xcenter*RAD2DEG<<endl;
    for(int i=0;i<mN;i++)
      F<<setw(18)<<mvL(i)<<setw(18)<<mvIn[i].real()<<setw(18)<<mvIn[i].imag()<<endl;
    F.close();
  }

  fftw_plan p = fftw_plan_dft_1d (mN,
				  reinterpret_cast<fftw_complex*>(&mvIn[0]),
				  reinterpret_cast<fftw_complex*>(&mvOut[0]),
				  FFTW_FORWARD, FFTW_ESTIMATE);

  fftw_execute(p);

  if (bsavecalc && xcenter>=xcenterlimits[0] && xcenter<=xcenterlimits[1]) {
    ofstream F("profileF.dat");
    F<<"# xcenter="<<xcenter*RAD2DEG<<endl;
    for(int i=0;i<mN;i++)
      F<<setw(18)<<mvs(i)<<setw(18)<<mvOut[i].real()*mdL<<setw(18)<<mvOut[i].imag()*mdL<<endl;
    F.close();
  }

  // copy real part of the calc. profile and do a fft-shift
  for(int i=0;i<=mN/2;i++) mvProfile[mN/2-1+i] = mvOut[i].real()*mdL;
  for(int i=mN/2+1;i<=mN-1;i++) mvProfile[i-mN/2-1] = mvOut[i].real()*mdL;
 // this part I need 
  { // common order of x-data
    REAL *p=mvs.data();
    for(int i=-mN/2+1;i<=mN/2;i++) *p++=i*mds;
  }

  if (bsavecalc && xcenter>=xcenterlimits[0] && xcenter<=xcenterlimits[1]) {
    ofstream F("profileF1.dat");
    F<<"# xcenter="<<xcenter*RAD2DEG<<endl;
    for(int i=0;i<mN;i++)
      F<<setw(18)<<mvs(i)<<setw(18)<<mvProfile[i]<<endl;
    F.close();
  }
/// ------ if (fwhm*factor<mds) unneccessary
  // will be used later for direct convolution
  vector<double> fx(mN); // TODO:: sometimes not necessary
  {
    const REAL *p=mvs.data();
    for(int i=0;i<mN;i++) fx[i]=*p++;
  }

  // direct convolution in recip. space - TODO:: integ. intensity. ???
  for(int k=0;k<mReflectionProfileComponentRegistry.GetNb();k++) {
    ReflectionProfileComponent &comp = mReflectionProfileComponentRegistry.GetObj(k);
    if (comp.IsRealSpaceType()) continue;
    REAL fwhm1 = comp.GetApproxFWHM(xcenter,H,K,L);
    CrystVector_REAL profile=( (fwhm1*factor>=mds) || (mReflectionProfileComponentRegistry.GetNb()==1) ) //nonsence ???!!!						
    							? comp.GetProfile(mvs,xcenter,H,K,L)
    							: CrystVector_REAL(0);
    int n = profile.numElements();
    if (1 & ((unsigned int) n)) // is odd (error)
       throw ObjCrystException("MStruct::ReflectionProfile::GetProfile(...)\
             Wrong nb. of calculated profile points.");
    else
      n = n/2;
    // copy data
    vector<double> gx(n), gy(n);
    const REAL *p=profile.data();
    for(int i=0;i<n;i++) gx[i]=*p++;
    for(int i=0;i<n;i++) gy[i]=*p++;
    // if the direct space profile (f) is a "delta-function"
    if (fwhm*factor<mds) {
      for(int i=0;i<mN;i++) mvProfile[i]=interp1(gx,gy,mvs(i));
    } 
    else { // convolution
      if (fwhm1*factor<mds) continue; // g is a "delta-function"
      // copy physical profile
      // I not use vector properly !!! (cheat) (mvProfile.size() == 0 !!!)
      vector<double> fy(mN);
      for(int i=0;i<mN;i++) fy[i]=mvProfile[i];
      // convolution (always the 1st profile in conv2 routine is more
      // broad than the 2nd one)
      if (fwhm>=fwhm1)
	for(int i=0;i<mN;i++) mvProfile[i]=conv2(fx,fy,gx,gy,mvs(i));
      else
	for(int i=0;i<mN;i++) mvProfile[i]=conv2(gx,gy,fx,fy,mvs(i));
    }
    fwhm += fwhm1;
  }
 
  if (bsavecalc && xcenter>=xcenterlimits[0] && xcenter<=xcenterlimits[1]) {
    ofstream F("profileP1.dat");
    F<<"# xcenter="<<xcenter*RAD2DEG<<endl;
    for(int i=0;i<mN;i++)
      F<<setw(18)<<mvs(i)<<setw(18)<<mvProfile[i]<<endl;
    F.close();
  }

  /* free memory */
  fftw_destroy_plan (p);

  // save calc profile into the "store"
  ReflCalcData *pData = new ReflCalcData;
  pData->s1=ms1; pData->s0=ms0; pData->ds=mds;
  pData->N=mN; pData->vProfile=mvProfile;
  //if ((int)pData->vProfile.capacity()<mN) pData->vProfile.reserve(mN);
  if ((int)pData->vProfile.size()<mN) pData->vProfile.resize(mN);
  double *p1=&(mvProfile[0]), *p2=&(pData->vProfile[0]);
  for(int i=0;i<mN;i++,p1++,p2++) *p2 = *p1;
  mReflStore.add((int)H,(int)K,(int)L,xcenter,(void*)pData);

  } // if(needRecalc)

  CrystVector_REAL profile(mNbPoints);

  
  // position corr
  REAL xcorr = GetPositionCorr(xcenter,H,K,L);

  // add postion corr of all components - not necessary to recalc always
  for(int k=0;k<mReflectionProfileComponentRegistry.GetNb();k++) {
    ReflectionProfileComponent &comp = mReflectionProfileComponentRegistry.GetObj(k);
    xcorr += comp.GetPositionCorr(xcenter,H,K,L);
  }
  
  // linear interpolation 
  for(int i=0;i<mNbPoints;i++) {
    double s = ((x(i)-xcorr)>=0.0) && ((x(i)-xcorr)<=M_PI) ? 2.*sin((x(i)-xcorr)/2)/mLambda - ms0
    	: factor * (x(i)-xcorr-xcenter);
    int m = int(floor(s/mds));
    REAL y1 = (-mN/2<=(m-1) && (m-1)<mN/2) ? mvProfile[mN/2+m-1]
                                           : (-mN/2>(m-1) ? mvProfile[0] : mvProfile[mN-1]);
    REAL y2 = (-mN/2<=m     &&     m<mN/2) ? mvProfile[mN/2+m]
                                           : (-mN/2>m     ? mvProfile[0] : mvProfile[mN-1]);
    profile(i) = y1 + (s-m*mds)/mds*(y2-y1);
  }

  profile *= factor;

  mClockReflectionProfileCalc.Click();

  if (bsavecalc && xcenter>=xcenterlimits[0] && xcenter<=xcenterlimits[1]) {
    ofstream F("profileP.dat");
    F<<"# xcenter="<<xcenter*RAD2DEG<<endl;
    for(int i=0;i<mNbPoints;i++)
      F<<setw(18)<<(x(i)-xcorr)*RAD2DEG<<setw(18)<<profile(i)<<endl;
    F.close();
  }

  return profile;
}

REAL ReflectionProfile::GetApproxFWHM(const REAL xcenter,
		     const REAL h, const REAL k, const REAL l)const
{
  REAL fwhm = 0.;
  // sum appprox. FWHMs from all effects
  for(int k=0;k<mReflectionProfileComponentRegistry.GetNb();k++) {
    const ReflectionProfileComponent &comp = mReflectionProfileComponentRegistry.GetObj(k);
    fwhm += comp.GetApproxFWHM(xcenter,h,k,l);
  }

  return fwhm;
}

REAL ReflectionProfile::GetFullProfileWidth (const REAL relativeIntensity,
					     const REAL xcenter,
					     const REAL h,const REAL k,
					     const REAL l)
{
  // copy of code in ObjCryst files (ReflectionProfile.cpp)
  VFN_DEBUG_ENTRY("MStruct::ReflectionProfile::GetFullProfileWidth()",10)
  const int nb=100; // 1000
  const int halfnb=nb/2;
  CrystVector_REAL x(nb);
  REAL n=10.0; // 5.0
  const REAL fwhm=GetApproxFWHM(xcenter,h,k,l);
  CrystVector_REAL prof;
  while(true)
  {
    REAL *p=x.data();
    const REAL tmp=fwhm*n/nb;
    for(int i=0;i<nb;i++) *p++ = tmp*(i-halfnb);
    x+=xcenter;
    prof=this->GetProfile(x,xcenter,h,k,l);
    const REAL max=prof.max();
    const REAL test=max*relativeIntensity;
    int n1=0,n2=0;
    if((prof(0)<test)&&(prof(nb-1)<test))
    {
      p=prof.data();
      while(*p<test){ p++; n1++;n2++;}
      n1--;
      n2=nb-1; p=prof.data()+n2;
      while(*p<test){ p--; n2--;}
      n2++;
      VFN_DEBUG_EXIT("MStruct::ReflectionProfile::GetFullProfileWidth():"<<x(n2)-x(n1),10)
      return x(n2)-x(n1);
      // DUMMY !!!
	/*
      REAL profilewidth = x(n2)-x(n1);
      // position corr
      REAL xcorr = GetPositionCorr(xcenter,h,k,l);
      // add postion corr of all components
      for(int k=0;k<mReflectionProfileComponentRegistry.GetNb();k++) {
	ReflectionProfileComponent &comp = mReflectionProfileComponentRegistry.GetObj(k);
	xcorr += comp.GetPositionCorr(xcenter,h,k,l);
      }
      return profilewidth+2.*fabs(xcorr);
      */
    }
    VFN_DEBUG_MESSAGE("MStruct::ReflectionProfile::GetFullProfileWidth():"<<max<<","<<test
		      <<endl<<FormatVertVector<REAL>(x,prof),10)
    n*=2.0;
    //if(n>200) exit(0);
  }
}

REAL ReflectionProfile::GetIntegralWidth (const REAL xcenter,
					  const REAL h,const REAL k,const REAL l) 
{
  // profile is normalised to unity - just find a value in the maximum
  // it is not possible to simply calculate the value in one point
  // hence it is necessary calc the whole profile
  
  const REAL fullwidth = this->GetFullProfileWidth(0.001,xcenter,h,k,l);
  const int n = 511; // 1023
  const REAL dx = fullwidth/(n-1);
  CrystVector_REAL x(n);
  REAL *p = x.data();
  for(int i=-n/2;i<=n/2;i++) { *p = i*dx; p++; }
  x += xcenter;
  CrystVector_REAL prof = this->GetProfile(x,xcenter,h,k,l);
  
  /*
  if((abs(h-1.)+abs(k-0.)+abs(l-1.))<1e-7) {
    ofstream s("prof_101.dat");
    for(int i=0;i<x.numElements();i++) {
      s<<setw(12)<<scientific<<setprecision(3)<<x(i);
      s<<setw(12)<<scientific<<setprecision(3)<<prof(i)<<"\n";
    }
    s.close();
  }
  */
  return prof.sum()*dx/prof.max();
}

void ReflectionProfile::SetIncidenceAngle(const REAL omega)
{
	mOmega = omega;
	mClockMaster.Click();
}

REAL ReflectionProfile::GetIncidenceAngle(const REAL xcenter) const
{
	return (mOmega>0.) ? mOmega : xcenter/2;
}

REAL ReflectionProfile::GetPositionCorr(const REAL xcenter,
				        const REAL h, const REAL k, const REAL l) const
{
  // TODO: TOF data ?, tetragonal only
  /* Where is Psi here?
  REAL corr = 0.; //d2Theta

  const REAL ac = mpUnitCell->GetLatticePar(0)/mpUnitCell->GetLatticePar(2);  
  const REAL G0 = h*h+k*k*+ac*ac*l*l;
  if(G0>FLT_EPSILON) {
    const REAL GH = (h*h*k*k+ac*ac*(h*h+k*k)*l*l)/G0;
    const REAL HH = (ac*ac*l*l)/G0;

    corr += mStressCoeff(0);
    corr += mStressCoeff(1)*GH + (mStressCoeff(2)+mStressCoeff(3)*HH)*HH;
    corr *= -2.*tan(xcenter/2);
  }
  */
  
  REAL corr = 0.; //d2Theta
  //REAL tt = xcenter/2 - this->GetParentPowderPatternDiffraction().
  //							GetParentPowderPattern().GetIncidenceAngle(); ???
  REAL omega = (mOmega>0.) ? mOmega : xcenter/2;
  REAL tt = xcenter/2 - omega;
  tt = sin(tt);
  
  // mStressCoeff(0) - stress
  // (mStressCoeff(1) >= 0.) ? Reuss/Voigt model weigth : isotropic case
  // mStressCoeff(2) S11  /  E
  // mStressCoeff(3) S12  /  ni
  // mStressCoeff(4) S44

  REAL S1, S2;

  if (mStressCoeff(1)>=0.) {
	  // Reuss/Voigt
	  const REAL GG = (h*h*k*k + k*k*l*l + l*l*h*h)/pow(h*h + k*k + l*l,2);
	  const REAL s11 = mStressCoeff(2);
	  const REAL s12 = mStressCoeff(3);
	  const REAL s44 = mStressCoeff(4);

	  S1 = (1.-mStressCoeff(1))*(s11*(2*s11+2*s12-s44)+s12*(3*s44-4*s12))/
			                    (2*s44+6*(s11-s12))  +
			   mStressCoeff(1) *(s12+GG*(s11-s12-0.5*s44));
	  S2 = (1.-mStressCoeff(1))*2.*(5*(s11-s12)*s44)/(2*s44+6*(s11-s12))  +
			   mStressCoeff(1) *2.*(s11-s12-3*GG*(s11-s12-0.5*s44));
  }
  else {
	  // isotropic case
	  const REAL E = mStressCoeff(2);
	  const REAL ni = mStressCoeff(3);
	  S1 = - ni/E;
	  S2 = 2.*(1+ni)/E;
  }

  REAL ee = (fabs(mStressCoeff(0))>FLT_MIN) ? (0.5*S2*tt*tt + 2.*S1)*mStressCoeff(0) : 0.;
  corr = -2.*tan(xcenter/2)*ee;
  
  return corr;
}

void ReflectionProfile::InitParameters()
{
  {
    RefinablePar tmp("q0",&mStressCoeff(0),-1,1.,
		     gpRefParTypeScattDataProfileType,
		     REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false);
    tmp.AssignClock(mClockStressCorrQ);
    tmp.SetDerivStep(1e-4);
    this->AddPar(tmp);
  }
  {
    RefinablePar tmp("q1",&mStressCoeff(1),-1,1.,
		     gpRefParTypeScattDataProfileType,
		     REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false);
    tmp.AssignClock(mClockStressCorrQ);
    tmp.SetDerivStep(1e-4);
    this->AddPar(tmp);
  }
  {
    RefinablePar tmp("q2",&mStressCoeff(2),-1,1.,
		     gpRefParTypeScattDataProfileType,
		     REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false);
    tmp.AssignClock(mClockStressCorrQ);
    tmp.SetDerivStep(1e-4);
    this->AddPar(tmp);
  }
  {
    RefinablePar tmp("q3",&mStressCoeff(3),-1,1.,
		     gpRefParTypeScattDataProfileType,
		     REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false);
    tmp.AssignClock(mClockStressCorrQ);
    tmp.SetDerivStep(1e-4);
    this->AddPar(tmp);
  }
  {
    RefinablePar tmp("q4",&mStressCoeff(4),-1,1.,
		     gpRefParTypeScattDataProfileType,
		     REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false);
    tmp.AssignClock(mClockStressCorrQ);
    tmp.SetDerivStep(1e-4);
    this->AddPar(tmp);
  }
}

////////////////////////////////////////////////////////////////////////
//
//    ReflectionPositionCorrBase
//
////////////////////////////////////////////////////////////////////////

ReflectionPositionCorrBase::ReflectionPositionCorrBase()
{}

ReflectionPositionCorrBase::ReflectionPositionCorrBase(const ReflectionPositionCorrBase & old)
//:ReflectionProfileComponent(old) // ReflectionProfileComponent has no copy constructor and the RefinableObj copy constructor is not implemented
{}

const string& ReflectionPositionCorrBase::GetClassName()const
{
	const static string name = "MStruct::ReflectionPositionCorrBase";
	
	return name;
}

CrystVector_REAL ReflectionPositionCorrBase::GetProfile(const CrystVector_REAL &x,
			      										const REAL xcenter,
			      										const REAL h, const REAL k, const REAL l)
{
	/* This method should be never called. If called the method prints a warning
	 * message and return a profile with only one nonzero point - the one most closed
	 * to the xcenter. This can cause an unpredicable small reflection shifts. 
	 */
	
	cerr << "Warning: MStruct::ReflectionPositionCorrBase::GetProfile(...) method should be never called.\n";
	cerr << "\tThis can cause an unpredicable small reflection shifts." << endl;
	
	unsigned long nbPoints = x.numElements(); 
	CrystVector_REAL profile(nbPoints);
 
	const REAL *p1 = x.data();
	
	unsigned long icenter = (unsigned long) nbPoints/2;
	REAL min_dist = x(nbPoints-1)-x(0);
	
	for(unsigned long i=0; i<nbPoints; i++) {
		if (abs((*p1)-xcenter)<min_dist) { icenter = i; min_dist = abs((*p1)-xcenter); }
		p1++;
	}
	
	profile = 0.;
	profile(icenter) = 1.;
	
	return profile;
}

REAL ReflectionPositionCorrBase::GetApproxFWHM(const REAL xcenter,
								     		   const REAL h,
									   		   const REAL k,
									   		   const REAL l)const
{
	return 0.;
}

bool ReflectionPositionCorrBase::IsRealSpaceType()const
{
	return false;
}

REAL ReflectionPositionCorrBase::GetPositionCorr(const REAL xcenter,
		       					 				 										 const REAL h, const REAL k, const REAL l)const
{
	return 0.; // by default
}

////////////////////////////////////////////////////////////////////////
//
//    RefractionPositionCorr
//
////////////////////////////////////////////////////////////////////////

RefractionPositionCorr::RefractionPositionCorr()
:mDensity(1.),mpCrystal(NULL),mAbsDensity(0.),mChi0ValueChoice(CHI0_VALUE),mConsiderCrystalFixed(true),mChi0(0.,0.),
mInitializationFlag(0)
{
	InitParameters();
	mClockMaster.AddChild(mClockChi0);
}

RefractionPositionCorr::RefractionPositionCorr(const RefractionPositionCorr & old)
:ReflectionPositionCorrBase(old),mDensity(old.mDensity),mpCrystal(NULL),mAbsDensity(old.mAbsDensity),
mChi0ValueChoice(old.mChi0ValueChoice),mConsiderCrystalFixed(true),mChi0(old.mChi0),mClockChi0(old.mClockChi0),
mInitializationFlag(old.mInitializationFlag)
{
	InitParameters();
	if ( old.mChi0ValueChoice==CHI0_CRYSTAL && old.mpCrystal!=NULL )
		SetCrystal(*old.mpCrystal,old.mConsiderCrystalFixed);
	mClockMaster.AddChild(mClockChi0);
}

/*RefractionPositionCorr::~RefractionPositionCorr()
{
		mClockMaster.RemoveChild(mClockChi0);
}*/

const string& RefractionPositionCorr::GetClassName()const
{
	const static string name = "MStruct::RefractionPositionCorr";
	
	return name;
}

REAL RefractionPositionCorr::GetPositionCorr(const REAL xcenter,
							       					 							 const REAL h, const REAL k, const REAL l)const
{
	// calculate (update) chi0 value (without forcing recalculation)  
	GetChi0(false);
	
	// immediatelly return zero in case of a zero chi0
	if (abs(mChi0)<=1e-7) return 0.;
	
	// get incidence angle for the current 2Theta position
	const REAL alpha = this->GetParentReflectionProfile().GetIncidenceAngle(xcenter);
	
	const REAL delta = -mChi0.real()/2 * mDensity;
	const REAL beta  = -mChi0.imag()/2 * mDensity;
 
 	//if (alpha < sqrt(2*delta)) return alpha; // corr( 2Theta ) (rad)
 	
	const REAL dalpha2 = pow(alpha,2)-2*delta; 
	
	// corr( 2Theta ) (rad)
	return alpha - M_SQRT1_2 * sqrt( dalpha2 + sqrt( pow(dalpha2,2)+4*pow(beta,2) ) );
}

void RefractionPositionCorr::SetCrystal(const ObjCryst::Crystal & crystal, const bool fixed)
{
	// Set new Crystal
	try {
		mChi0 = complex<REAL>(0.,0.);
		mFormula.SetFormula("",cout);
		mpCrystal = &crystal;
		mClockChi0.Reset();
  	mConsiderCrystalFixed = fixed;
  	mChi0ValueChoice = CHI0_CRYSTAL;
  	mInitializationFlag = 0;
	}
	catch (std::exception &e) {
		cerr << "< MStruct::RefractionPositionCorr::SetCrystal(...)\n";
		cerr << "  Unexpected exception when setting Crystal: " << e.what() << "\n";
		cerr << ">" << endl;
		throw ObjCrystException("MStruct::RefractionPositionCorr::SetCrystal(...): Unexpected error.");
  }
	
	// Force chi0 value recalculation (also absolute density is recalculated)
	///GetChi0(true); // object could here still not be connected with diffraction data, radiation, etc.
	
	// Print density and chi0 value

	/*ostringstream os;
	os << "MStruct::RefractionPositionCorr::SetCrystal(...): Chi0 and absolute density computed for Crystal: " << mpCrystal->GetName() << "\n";
	os << "\t" << "chi0: " << scientific << setprecision(4) << mChi0 << " (n=1-delata-ii*beta~=1+chi0/2)\n";
	os << "\t" << "critical angle: " << FormatFloat(sqrt(-mChi0.real() * mDensity)*RAD2DEG,4,2) << "(deg)" << "\n";
	os << "\t" << "density: " << FormatFloat(mAbsDensity,6,3) << " (g/cm3)\n";
	cout << os.str() << flush;*/
		
		/*const ObjCryst::ScatteringPowerAtom * pNewScattPowAtom =
							new ObjCryst::ScatteringPowerAtom("Mg","Mg_chi0"); 
		
		const REAL F_forward = pNewScattPowAtom->GetForwardScatteringFactor(radiationType);
		const REAL F_resonantReal = pNewScattPowAtom->GetResonantScattFactReal(scattData)(0); //:TODO: More than one wavelength
		const REAL F_resonantImag = pNewScattPowAtom->GetResonantScattFactImag(scattData)(0); //:TODO: More than one wavelength
		
		cout << " F_f: " << F_forward << ", F_rR: " << F_resonantReal << ", F_rI: " << F_resonantImag << endl;  
		*/
}

void RefractionPositionCorr::SetChemicalFormula(const string & formula, const REAL density)
{
	try {
		mpCrystal = NULL;
		mConsiderCrystalFixed = true;
		mInitializationFlag = 0;
		mChi0 = complex<REAL>(0.,0.);
		
		mChi0ValueChoice = CHI0_CHEM_FORMULA;
		
		mFormula.SetFormula(formula,cout);
		
		mAbsDensity = density;
		
		mClockChi0.Reset();
	}
	catch (std::exception &e) {
		cerr << "< MStruct::RefractionPositionCorr::SetChemicalFormula(...)\n";
		cerr << "  Unexpected exception when setting Chemical formula: " << e.what() << "\n";
		cerr << ">" << endl;
		throw ObjCrystException("MStruct::RefractionPositionCorr::SetChemicalFormula(...): Unexpected error.");
  }
}

REAL CalcUnitCellMass(const ObjCryst::Crystal& crystal)
{
  // get Scattering Component List to obtain Atom scatters
  const ObjCryst::ScatteringComponentList & sc_list = crystal.GetScatteringComponentList();
					
  REAL mass = 0.;
					
  for (long icomp=0; icomp<sc_list.GetNbComponent(); icomp++) {
    
    const ObjCryst::ScatteringPower * pScattPow = sc_list(icomp).mpScattPow;
								
    if ( pScattPow->GetClassName()=="ScatteringPowerAtom" ) {
      const ObjCryst::ScatteringPowerAtom * pScattPowAtom =
	dynamic_cast<const ObjCryst::ScatteringPowerAtom*>(pScattPow);
      // TODO:: Atom scatter population could be calculated also using Dynamical populacy correction without
      //        calling SpaceGroup::GetAllSymmetrics(...) method but DynPopCorr looks working wrongly
      //        in the case of an one atom crystal (e.g. Al, Cu, Mg, etc.)
      ///const REAL popu = sc_list(icomp).mOccupancy*sc_list(icomp).mDynPopCorr;
							
      // find atom in the periodic table to get its atomic weight 
      cctbx::eltbx::tiny_pse::table tpse(pScattPowAtom->GetSymbol());
      // number of equivalent atom site symmetry positions
      const int nb_sym_pos = crystal.GetSpaceGroup().GetAllSymmetrics(sc_list(icomp).mX,sc_list(icomp).mY,
								      sc_list(icomp).mZ,false,false,true).rows(); 
      //cout << "Scatterer nb. " << icomp << ", symbol: " << pScattPowAtom->GetSymbol();
      //cout << ", x, y, z: " << sc_list(icomp).mX << ", " << sc_list(icomp).mY << ", " << sc_list(icomp).mZ << ", sym: " << nb_sym_pos;
      //cout << ", occ: " << sc_list(icomp).mOccupancy << ", DynPopCorr: " << sc_list(icomp).mDynPopCorr << ", popu: " << popu << ", weight: " << tpse.weight() << endl;
					
      // atomic mass constant (ref: http://physics.nist.gov/cuu)
      // m_u = 1.660 538 782 x 10^-27 kg
      // 
      //  density = Sum_atoms { nb_sym_pos*occ*weight } * m_u / UnitCell.volume
      //  density(g/cm3)   exponents   (27-3) - (3*10-3*2) = 24 - 24 = 0
      
      REAL t = nb_sym_pos * sc_list(icomp).mOccupancy;
      mass += t * tpse.weight();
    }
  } // for icomp
  
  mass *= 1.660538782; // (1e-24 g)

  return mass; // (1e-24 g)
}

void RefractionPositionCorr::SetChi0(const complex<REAL> & chi0)
{
	mpCrystal = NULL;
	mFormula.SetFormula("",cout);
	mAbsDensity = 0.;
	mConsiderCrystalFixed = true;
	mInitializationFlag = 0;
	
	mChi0ValueChoice = CHI0_VALUE;
	
	mChi0 = chi0;
	mClockChi0.Click();
}

const complex< REAL > & RefractionPositionCorr::GetChi0(const bool forceReCalc)const
{
	try {
		
		switch (mChi0ValueChoice) {
			
			case CHI0_VALUE:
				// chi0 value set directly - nothing to do
				break;
				
			case CHI0_CRYSTAL:
				// chi0 value calculated from crystal structure
				if (mpCrystal==NULL) {
					cerr << "< MStruct::RefractionPositionCorr::GetChi0(...)\n";
					cerr << "  No crystal set. Can not calculate Chi0 value. \n";
					cerr << ">" << endl;
					throw ObjCrystException("MStruct::RefractionPositionCorr::GetChi0(...): Unexpected error.");
				}
				// check if (re)calculation is necessary
				if ( forceReCalc==true || mInitializationFlag<1 ||
						 ( mConsiderCrystalFixed==false && mClockChi0<mpCrystal->GetClockLatticePar() 
						                                && mClockChi0<mpCrystal->GetMasterClockScatteringPower()
						 															  && mClockChi0<mpCrystal->GetClockScattCompList() ) ) {
					// calculate chi0 value and crystal density
					
					// TODO:: Scattering factor could be calculated also directly using "FOX-build-in" methods, unfortunately
					//    calling calculation of resonant factors for (000) reflection is not straightforward
					//    and it looks that there are some additional problems with dynamical populacy correction
					//    in the case of one-atomic simple structures
					
					// forward scattering factor
					REAL S_forward = 0.;
					// resonant scattering factor
					REAL S_resonantReal = 0., S_resonantImag = 0.;
		
					// get scattering data (for calculation of the resonant scattering factor)
					const ObjCryst::PowderPatternDiffraction & scattData = GetParentReflectionProfile().GetParentPowderPatternDiffraction();
					// get radiation type (for calculation of the forward scattering factor)
					const ObjCryst::RadiationType & radiationType = GetParentReflectionProfile().GetParentPowderPatternDiffraction().GetRadiationType();
		
					// get Scattering Component List to obtain Atom scatters
					const ObjCryst::ScatteringComponentList & sc_list = mpCrystal->GetScatteringComponentList();
					
					mAbsDensity = 0.;
					
					for (long icomp=0; icomp<sc_list.GetNbComponent(); icomp++) {
		
						const ObjCryst::ScatteringPower * pScattPow = sc_list(icomp).mpScattPow;
								
						if ( pScattPow->GetClassName()=="ScatteringPowerAtom" ) {
							const ObjCryst::ScatteringPowerAtom * pScattPowAtom = dynamic_cast<const ObjCryst::ScatteringPowerAtom*>(pScattPow);
							// TODO:: Atom scatter population could be calculated also using Dynamical populacy correction without calling
							//          SpaceGroup::GetAllSymmetrics(...) method but DynPopCorr looks working wrongly in the case of an one atom crystal (e.g. Al, Cu, Mg, etc.)
							const REAL popu = sc_list(icomp).mOccupancy*sc_list(icomp).mDynPopCorr;
							
							// find atom in the periodic table to get its atomic weight 
							cctbx::eltbx::tiny_pse::table tpse(pScattPowAtom->GetSymbol());
							// number of equivalent atom site symmetry positions
							const int nb_sym_pos = mpCrystal->GetSpaceGroup().GetAllSymmetrics(sc_list(icomp).mX,sc_list(icomp).mY,sc_list(icomp).mZ,false,false,true).rows(); 
							//cout << "Scatterer nb. " << icomp << ", symbol: " << pScattPowAtom->GetSymbol();
							//cout << ", x, y, z: " << sc_list(icomp).mX << ", " << sc_list(icomp).mY << ", " << sc_list(icomp).mZ << ", sym: " << nb_sym_pos;
							//cout << ", occ: " << sc_list(icomp).mOccupancy << ", DynPopCorr: " << sc_list(icomp).mDynPopCorr << ", popu: " << popu << ", weight: " << tpse.weight() << endl;
					
							// atomic mass constant (ref: http://physics.nist.gov/cuu)
							// m_u = 1.660 538 782 x 10^-27 kg
							// 
							//  density = Sum_atoms { nb_sym_pos*occ*weight } * m_u / UnitCell.volume
							//  density(g/cm3)   exponents   (27-3) - (3*10-3*2) = 24 - 24 = 0
							
							REAL t = nb_sym_pos * sc_list(icomp).mOccupancy;
							mAbsDensity += t * tpse.weight();
					
							// scattering factor calculation
							S_forward += t * pScattPowAtom->GetForwardScatteringFactor(radiationType);
							S_resonantReal += t * pScattPowAtom->GetResonantScattFactReal(scattData)(0); //:TODO: More than one wavelength
							S_resonantImag += t * pScattPowAtom->GetResonantScattFactImag(scattData)(0); //:TODO: More than one wavelength
						} else {
							cerr << "< MStruct::RefractionPositionCorr::GetChi0(...)\n";
							cerr << "  Warning: During calculation of the absolute density of the Crystal: " << mpCrystal->GetName() << "\n";
							cerr << "    a scattering component: " << pScattPow->GetName() << " of not ScatteringPowerAtom type: " << pScattPow->GetClassName() << " found." << "\n";
							cerr << "    The scatterrer not included in the density calculation." << "\n";
							cerr << ">" << endl;
						}
					} // for icomp
					
					//cout << "Nb. of symm.: " << mpCrystal->GetSpaceGroup().GetNbSymmetrics() << endl;
					//cout << "Volume: " << mpCrystal->GetVolume() << endl;
					mAbsDensity *= 1.660538782 / mpCrystal->GetVolume(); // (g/cm3)
					//cout << "density: " << density << " (g/cm3)" << endl;
					//cout << "ForwardScatteringFactor: " << S_forward << endl;
					//cout << "ResonantScattFactReal: " << S_resonantReal << endl;
					//cout << "ResonantScattFactImag " << S_resonantImag << endl;
		
					// chi0 calculation (classical electron radius rel = 2.8179e-5 A)

					mChi0 = complex<REAL>(S_forward+S_resonantReal,S_resonantImag);
	
					mChi0 *= - 2.8179e-5 * pow(scattData.GetRadiation().GetWavelength()(0),2) / (M_PI * mpCrystal->GetVolume());
					
					mClockChi0.Click();
					
					if (mInitializationFlag==0) {
						// Print density and chi0 value
						ostringstream os;
						os << "MStruct::RefractionPositionCorr::GetChi0(...): Chi0 and absolute density computed for Crystal: " << mpCrystal->GetName() << "\n";
						os << "\t" << "chi0: " << scientific << setprecision(4) << mChi0 << " (n=1-delata-ii*beta~=1+chi0/2)\n";
						os << "\t" << "critical angle: " << FormatFloat(sqrt(-mChi0.real() * mDensity)*RAD2DEG,4,2) << "(deg)" << "\n";
						os << "\t" << "density: " << FormatFloat(mAbsDensity,6,3) << " (g/cm3)\n";
						cout << os.str() << flush;
					} // if (mInitializationFlag==0)
					
					mInitializationFlag = 1;
					
				} // if - (re)calculation was necessary
				
				break;
			
			case CHI0_CHEM_FORMULA:
				// chi0 value calculated from chemical formula and absolute material density
				
				// check if (re)calculation is necessary
				if ( forceReCalc==true || mInitializationFlag<1 || mClockChi0<mFormula.GetClock() ) {
					
					mChi0 = complex<REAL>(0.,0.);
					
					// The complex dielectric permittivity value chi0 can be calculated from the total complex
					// scattering power S(formula) of the molecule represented by the chemical formula as:
					//
					//            chi0(formula) = -rel*Lambda^2 / ( PI*V(formula) ) * S(formula),
					//
					// where rel ~ 2.8179e-5 A is the classical electron radius and V(formula) represents volume
					// per one molecule. This can be calculated as:
					// 
					//            V(fomula) = m_u * Sum_atoms { atomic_weight } / density,
					//
					// where m_u ~ 1.66 10^-27 kg is the atomic mass constant. If density value (stored here
					// in g/cm^3) is written in kg/A^3 the exponent factors (-27) cancel themselves.  
					
					// forward scattering factor
					REAL S_forward = 0.;
					// resonant scattering factor
					REAL S_resonantReal = 0., S_resonantImag = 0.;
					// sum of atomic weights
					REAL t = 0.;
		
					// get scattering data (for calculation of the resonant scattering factor)
					const ObjCryst::PowderPatternDiffraction & scattData = GetParentReflectionProfile().GetParentPowderPatternDiffraction();
					// get radiation type (for calculation of the forward scattering factor)
					const ObjCryst::RadiationType & radiationType = GetParentReflectionProfile().GetParentPowderPatternDiffraction().GetRadiationType();
					
					// loop over all molecule (chemical formula) atom scatterers
					vector< ObjCryst::ScatteringPowerAtom >::const_iterator it1 = mFormula.GetScatteringPowerAtomList().begin();
					const vector< ObjCryst::ScatteringPowerAtom >::const_iterator it1_end = mFormula.GetScatteringPowerAtomList().end();
					vector< int >::const_iterator it2 = mFormula.GetScatteringPowerAtomCountList().begin();
					for( ; it1 != it1_end; it1++, it2++) {
						S_forward += *it2 * it1->GetForwardScatteringFactor(radiationType);
						S_resonantReal += *it2 * it1->GetResonantScattFactReal(scattData)(0); //:TODO: More than one wavelength
						S_resonantImag += *it2 * it1->GetResonantScattFactImag(scattData)(0); //:TODO: More than one wavelength
						// find the atom in the periodic table to get its atomic weight 
						cctbx::eltbx::tiny_pse::table tpse(it1->GetSymbol());
						t += *it2 * tpse.weight();
					}
					
					// volume per one molecule (formula) - A^3
					t = 1.660538782 * t / mAbsDensity;
					
					// permittivity value 
					mChi0 = complex<REAL>(S_forward+S_resonantReal,S_resonantImag);
	
					mChi0 *= - 2.8179e-5 * pow(scattData.GetRadiation().GetWavelength()(0),2) / (M_PI * t);
					
					mClockChi0.Click();
					
					if (mInitializationFlag==0) {
						// Print density and chi0 value
						ostringstream os;
						os << "MStruct::RefractionPositionCorr::GetChi0(...): Chi0 computed for chemical formula: " << mFormula.GetFormula() << "\n";
						os << "\t" << "chi0: " << scientific << setprecision(4) << mChi0 << " (n=1-delata-ii*beta~=1+chi0/2)\n";
						os << "\t" << "critical angle: " << FormatFloat(sqrt(-mChi0.real() * mDensity)*RAD2DEG,4,2) << "(deg)" << "\n";
						os << "\t" << "density: " << FormatFloat(mAbsDensity,6,3) << " (g/cm3)\n";
						cout << os.str() << flush;
					} // if (mInitializationFlag==0)
								
					mInitializationFlag = 1;
					
				} // if - (re)calculation was necessary 
				
				break;
					
			default:
				cerr << "< MStruct::RefractionPositionCorr::GetChi0(...)\n";
				cerr << "  Can not recognise method how Chi0 value should be calculated. Choice: " << mChi0ValueChoice << "\n";
				cerr << ">" << endl;
				throw ObjCrystException("MStruct::RefractionPositionCorr::GetChi0(...): Unexpected error.");
				break;
		}
		
	} // try
	catch (std::exception &e) {
		cerr << "< MStruct::RefractionPositionCorr::GetChi0(...)\n";
		cerr << "  Unexpected exception when retrieving chi0 value: " << e.what() << "\n";
		cerr << ">" << endl;
		throw ObjCrystException("MStruct::RefractionPositionCorr::GetChi0(...): Unexpected error.");
  }
  
	return mChi0;
}

void RefractionPositionCorr::SetParams(const REAL density)
{
	mDensity = density;
	mClockMaster.Click();
}

void RefractionPositionCorr::InitParameters()
{
	{
		string name(GetName());
		if (!name.empty()) name = name + "_";
    RefinablePar tmp(name+"Density",&mDensity,1e-2,1e2,
		     gpRefParTypeScattDataCorrPos,
		     REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false);
    tmp.AssignClock(mClockMaster);
    tmp.SetDerivStep(1e-3);
    this->AddPar(tmp);
  }
}

////////////////////////////////////////////////////////////////////////
//
//   RefractionPositionCorr::ChemicalFormula 
//
////////////////////////////////////////////////////////////////////////

RefractionPositionCorr::ChemicalFormula::ChemicalFormula()
:mName("formula")
{}

RefractionPositionCorr::ChemicalFormula::ChemicalFormula(const ChemicalFormula & old)
:mName(old.mName), mFormula(old.mFormula), mvScatt(old.mvScatt), mvScattCount(old.mvScattCount),
mClock(old.mClock)
{}

const string & RefractionPositionCorr::ChemicalFormula::GetClassName()const
{
	const static string name = "MStruct::RefractionPositionCorr::ChemicalFormula";
	
	return name;
}

void RefractionPositionCorr::ChemicalFormula::SetName(const string & name)
{
	mName = name;
}

void RefractionPositionCorr::ChemicalFormula::SetFormula(const string & formula, ostream & ss)
{
	// clear data
	mFormula.clear();
	mvScatt.clear();
	mvScattCount.clear();
	
	// set formula string
	mFormula = formula;
	
	// analyse the formula string and build atom and atoms count lists
	try {
		// remove leading and terminating blank spaces
		string::size_type ind;
		
		ind = mFormula.find_first_not_of(" \t\n");
		if(ind!=string::npos) mFormula = mFormula.substr(ind);
		
		ind = mFormula.find_last_not_of(" \t\n");
		if(ind!=string::npos) mFormula = mFormula.substr(0,ind+1);
		
		map<const string, int> elements;
		
		// parse the formula string
		while (!mFormula.empty()) {
			// current found part of the formula
			string part = "";
			// remove leading blank spaces
			ind = mFormula.find_first_not_of(" \t\n");
			if(ind!=string::npos) mFormula = mFormula.substr(ind);
			// capital letters or blank spaces are possible chemical elements separators
			ind = mFormula.find_first_of("ABCDEFGHIJKLMNOPQRSTUVWXYZ ",1);
			// if search stoped at the 2nd position and it is not a blank space, try search once more
			if(ind!=string::npos && ind==1 && mFormula[ind]!=' ') ind = mFormula.find_first_of("ABCDEFGHIJKLMNOPQRSTUVWXYZ ",2);
			// if nothing found - this is the last part
			if(ind==string::npos) {
				part = mFormula.substr(0);
				mFormula = "";
			}
			// if something found, try to extract a formula part
			else {
				part = mFormula.substr(0,ind);
				mFormula = mFormula.substr(ind);
			}
			
			// process the extracted part of the formula
			
			// check last chars of the string - if it is a number, this means an element count
			int count = 1;
			ind = part.length()-1;
			for(string::reverse_iterator it=part.rbegin(); it!=part.rend(); ++it, ind--)
				if(*it=='+' || *it=='-' || isalpha(*it)) break;
			if(ind<part.length()-1) {
				count = atoi(part.substr(ind+1).c_str());
				// cut the string
				part = part.substr(0,ind+1);
			}
			
			// save possible scatterer name
			if(elements.find(part)==elements.end())
				elements[part] = count;
			else
				elements[part] += count;
		
		} // while (!mFormula.empty())
		
		// recover formula string
		mFormula = formula;
		
		// remove leading and terminating blank spaces
		ind = mFormula.find_first_not_of(" \t\n");
		if(ind!=string::npos) mFormula = mFormula.substr(ind);
		
		ind = mFormula.find_last_not_of(" \t\n");
		if(ind!=string::npos) mFormula = mFormula.substr(0,ind+1);
		
		// initialise all attributes by a formula represented by elements-map
		if(!mFormula.empty()) {
			ss << "\tChemical formula: " << mFormula << "\n";
			ss << "\t\tel.name ... el.count ... scattering power symbol" << "\n";
			// using preferably reverse iterator
			// if map-object is inserting items to its end, gives almost same order as in the fomula
			for( map<const string, int>::reverse_iterator it = elements.rbegin(); it != elements.rend(); ++it ) {
				string name = (*it).first + (mName.empty() ? "" : "_"+mName);
				mvScatt.push_back(ObjCryst::ScatteringPowerAtom(name,(*it).first));
				mvScattCount.push_back((*it).second);
	      ss << "\t\t" << (*it).first << " ... " << (*it).second << " ... " << mvScatt.back().GetSymbol() << "\n";
	    }
		} // if(!mFormula.empty())
	}
	catch (std::exception &e) {
		cerr << "< MStruct::RefractionPositionCorr::ChemicalFormula::SetFormula(...)\n";
		cerr << "  Unexpected exception when parsing chemical formula '" << formula << "' (name: " << mName << "): " << e.what() << "\n";
		cerr << ">" << endl;
		throw ObjCrystException("MStruct::RefractionPositionCorr::ChemicalFormula::SetFormula(...): Unexpected error.");
  }
  
	// set clock
	mClock.Click();
}

const string & RefractionPositionCorr::ChemicalFormula::GetFormula()const
{
	return mFormula;
}

const vector< ObjCryst::ScatteringPowerAtom > & RefractionPositionCorr::ChemicalFormula::GetScatteringPowerAtomList()const
{
	return mvScatt;
}

const vector< int > & RefractionPositionCorr::ChemicalFormula::GetScatteringPowerAtomCountList()const
{
	return mvScattCount;
}

const ObjCryst::RefinableObjClock & RefractionPositionCorr::ChemicalFormula::GetClock()const
{
	return mClock;
}

////////////////////////////////////////////////////////////////////////
//
//    ResidualStressPositionCorrection
//
////////////////////////////////////////////////////////////////////////

ResidualStressPositionCorrection::ResidualStressPositionCorrection()
:mStress(0.),pXECsObj(NULL)
{
	InitParameters();
}

ResidualStressPositionCorrection::ResidualStressPositionCorrection(const ResidualStressPositionCorrection & old)
:ReflectionPositionCorrBase(old),mStress(old.mStress),pXECsObj(old.pXECsObj)
{
	InitParameters();
	if( pXECsObj != NULL ) {
		this->AddSubRefObj(*pXECsObj);
		pXECsObj->RegisterClient(*this);
	} 
}

ResidualStressPositionCorrection::~ResidualStressPositionCorrection()
{
	if( pXECsObj != NULL ) {
		pXECsObj->DeRegisterClient(*this);
		this->RemoveSubRefObj(*pXECsObj);
	}
}

const string& ResidualStressPositionCorrection::GetClassName()const
{
	const static string name = "MStruct::ResidualStressPositionCorrection";
	
	return name;
}

void ResidualStressPositionCorrection::SetXECsObj(XECsObj & obj)
{
	if( pXECsObj != NULL ) {
		pXECsObj->DeRegisterClient(*this);
		this->RemoveSubRefObj(*pXECsObj);
	}
	
	pXECsObj = &obj;
	
	this->AddSubRefObj(obj);
	obj.RegisterClient(*this);
	
	mClockMaster.Reset();
}

void ResidualStressPositionCorrection::SetParams(const REAL stress)
{
	mStress = stress;
	mClockMaster.Click();
}

REAL ResidualStressPositionCorrection::GetPositionCorr(const REAL xcenter,
							       					 												 const REAL h, const REAL k, const REAL l)const
{
	// immediatelly return zero in case of a zero stress
	if (abs(mStress)<=1e-6) return 0.;
	
	// get incidence angle for the current 2Theta position
	const REAL omega = this->GetParentReflectionProfile().GetIncidenceAngle(xcenter);
	
	// calculate sin(psi)^2
  REAL tt = xcenter/2 - omega;
  tt = sin(tt);
 
  // calc XECs
  if( pXECsObj == NULL ) {
  	cerr << "< MStruct::ResidualStressPositionCorrection::GetPositionCorr(...)\n";
		cerr << "  xcenter=" << xcenter*RAD2DEG << "h=" << h << ",k=" << k << ",l=" << l << "\n";
		cerr << "  No XECs Object was set. >" << endl; 
		throw ObjCrystException("MStruct::ResidualStressPositionCorrection::GetPositionCorr(...): Object not properly initialised.");
  }
  
  REAL s1 = 0., s2 = 0.; // (1/GPa)
  try {
  	pXECsObj->GetXECs(s1,s2,h,k,l);
  }
  catch (std::exception &e) {
		cerr << "< MStruct::ResidualStressPositionCorrection::GetPositionCorr(...)\n";
		cerr << "  xcenter=" << xcenter*RAD2DEG << "h=" << h << ",k=" << k << ",l=" << l << "\n";
		cerr << "  Unexpected exception during XECs computation: " << e.what() << "\n";
		cerr << ">" << endl;
		throw ObjCrystException("MStruct::ResidualStressPositionCorrection::GetPositionCorr(...): Unexpected error.");
  }
	
	REAL ee = (0.5*s2*pow(tt,2) + 2*s1) * mStress; 
	
	return -2.*tan(xcenter/2)*ee; // corr( 2Theta(rad) )
}

void ResidualStressPositionCorrection::InitParameters()
{
	{
		string name(GetName());
		if (!name.empty()) name = name + "_";
    RefinablePar tmp(name+"Stress",&mStress,-1.e6,1.e6,
		     gpRefParTypeScattDataCorrPos,
		     REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false);
    tmp.AssignClock(mClockMaster);
    tmp.SetDerivStep(1e-4);
    this->AddPar(tmp);
  }
}

////////////////////////////////////////////////////////////////////////
//
//    XECsObj
//
////////////////////////////////////////////////////////////////////////

ObjCryst::ObjRegistry<XECsObj> gXECsObjRegistry("List of all XECs objects");

XECsObj::XECsObj()
{
	gXECsObjRegistry.Register(*this);
}

XECsObj::XECsObj(const XECsObj & old)
//:ObjCryst::RefinableObj(old) // Defined not implemented. 
{
	gXECsObjRegistry.Register(*this);
}

XECsObj::~XECsObj()
{
	gXECsObjRegistry.DeRegister(*this);
}

const string& XECsObj::GetClassName()const
{
	const static string name = "MStruct::XECsObj";
	
	return name;
}

////////////////////////////////////////////////////////////////////////
//
//    XECsIsotropic
//
////////////////////////////////////////////////////////////////////////

XECsIsotropic::XECsIsotropic()
:mE(211.),mNi(0.33)
{
	InitParameters();
}

XECsIsotropic::XECsIsotropic(const XECsIsotropic & old)
:XECsObj(old),mE(old.mE),mNi(old.mNi)
{
	InitParameters();
}

const string& XECsIsotropic::GetClassName()const
{
	const static string name = "MStruct::XECsIsotropic";
	
	return name;
}

void XECsIsotropic::GetXECs(REAL & s1, REAL & s2, const REAL h, const REAL k, const REAL l)const
{
	s1 = -mNi/mE;
	s2 = 2.*(1.+mNi)/mE;
}

void XECsIsotropic::SetElasticParameters(const REAL E, const REAL ni)
{
	mE = E;
	mNi = ni;
	mClockMaster.Click();
}
	
void XECsIsotropic::InitParameters()
{
	{
    RefinablePar tmp("E",&mE,0.,1.e6,
		     gpRefParTypeMaterialData,
		     REFPAR_DERIV_STEP_ABSOLUTE,true,true,false,false);
    tmp.AssignClock(mClockMaster);
    tmp.SetDerivStep(1.);
    this->AddPar(tmp);
  }
  
  {
    RefinablePar tmp("ni",&mNi,-10.,10.,
		     gpRefParTypeMaterialData,
		     REFPAR_DERIV_STEP_ABSOLUTE,true,true,false,false);
    tmp.AssignClock(mClockMaster);
    tmp.SetDerivStep(1e-3);
    this->AddPar(tmp);
  }
}

////////////////////////////////////////////////////////////////////////
//
//    XECsReussVoigt
//
////////////////////////////////////////////////////////////////////////
 
const int XECsReussVoigt::MCij_a_triclinic[][2] = {{1,1},{1,2},{1,3},{1,4},{1,5},{1,6},{2,2},{2,3},{2,4},{2,5},{2,6},{3,3},{3,4},{3,5},{3,6},{4,4},{4,5},{4,6},{5,5},{5,6},{6,6}};
const int XECsReussVoigt::MCij_a_monoclinic_axis_c[][2] = {{1,1},{1,2},{1,3},{1,6},{2,2},{2,3},{2,6},{3,3},{3,6},{4,4},{4,5},{5,5},{6,6}};
const int XECsReussVoigt::MCij_a_monoclinic_axis_b[][2] = {{1,1},{1,2},{1,3},{1,5},{2,2},{2,3},{2,5},{3,3},{3,5},{4,4},{4,6},{5,5},{6,6}};
const int XECsReussVoigt::MCij_a_orthorhombic[][2] = {{1,1},{1,2},{1,3},{2,2},{2,3},{3,3},{4,4},{5,5},{6,6}};
const int XECsReussVoigt::MCij_a_tetragonal_low[][2] = {{1,1},{1,2},{1,3},{1,6},{3,3},{4,4},{6,6}};
const int XECsReussVoigt::MCij_a_tetragonal_high[][2] = {{1,1},{1,2},{1,3},{3,3},{4,4},{6,6}};
const int XECsReussVoigt::MCij_a_trigonal_low[][2] = {{1,1},{1,2},{1,3},{1,4},{1,5},{3,3},{4,4}};
const int XECsReussVoigt::MCij_a_trigonal_high[][2] = {{1,1},{1,2},{1,3},{1,4},{3,3},{4,4}};
const int XECsReussVoigt::MCij_a_hexagonal[][2] = {{1,1},{1,2},{1,3},{3,3},{4,4}};
const int XECsReussVoigt::MCij_a_cubic[][2] = {{1,1},{1,2},{4,4}};

const XECsReussVoigt::MC_indices XECsReussVoigt::MCij_undefined(0, NULL);
const XECsReussVoigt::MC_indices XECsReussVoigt::MCij_triclinic(21, XECsReussVoigt::MCij_a_triclinic);
const XECsReussVoigt::MC_indices XECsReussVoigt::MCij_monoclinic_axis_c(13, XECsReussVoigt::MCij_a_monoclinic_axis_c);
const XECsReussVoigt::MC_indices XECsReussVoigt::MCij_monoclinic_axis_b(13, XECsReussVoigt::MCij_a_monoclinic_axis_b);
const XECsReussVoigt::MC_indices XECsReussVoigt::MCij_orthorhombic(9, XECsReussVoigt::MCij_a_orthorhombic);
const XECsReussVoigt::MC_indices XECsReussVoigt::MCij_tetragonal_low(7, XECsReussVoigt::MCij_a_tetragonal_low);
const XECsReussVoigt::MC_indices XECsReussVoigt::MCij_tetragonal_high(6, XECsReussVoigt::MCij_a_tetragonal_high);
const XECsReussVoigt::MC_indices XECsReussVoigt::MCij_trigonal_low(7, XECsReussVoigt::MCij_a_trigonal_low);
const XECsReussVoigt::MC_indices XECsReussVoigt::MCij_trigonal_high(6, XECsReussVoigt::MCij_a_trigonal_high);
const XECsReussVoigt::MC_indices XECsReussVoigt::MCij_hexagonal(5, XECsReussVoigt::MCij_a_hexagonal);
const XECsReussVoigt::MC_indices XECsReussVoigt::MCij_cubic(3, XECsReussVoigt::MCij_a_cubic);

XECsReussVoigt::XECsReussVoigt()
:mWeight(0.),mpUnitCell(NULL),mConsiderUnitCellFixed(true),mGroupType(LAUE_UNDEFINED),
mvCijConstantsNames(0),mCijValues(0),mCijMatrix(0,6),mSijMatrix(0,6),
mS1Voigt(-0.001),mS2Voigt(0.01),mMCij(MCij_undefined)
{
	InitParameters();
	mClockMaster.AddChild(mClockCijValues);
}

XECsReussVoigt::XECsReussVoigt(const XECsReussVoigt & old)
:mWeight(old.mWeight),mpUnitCell(NULL),mConsiderUnitCellFixed(true),mGroupType(LAUE_UNDEFINED),
mvCijConstantsNames(0),mCijValues(0),mCijMatrix(0,6),mSijMatrix(0,6),
mS1Voigt(-0.001),mS2Voigt(0.01),mMCij(MCij_undefined),mReflStore(old.mReflStore)
{
	InitParameters();
	mClockMaster.AddChild(mClockCijValues);
	SetUnitCell(*old.mpUnitCell, old.mConsiderUnitCellFixed);
}

XECsReussVoigt::~XECsReussVoigt()
{
	// clear Reflection Store
	for(int i=0;i<mReflStore.size();i++) {
		REAL *pData = (REAL*) mReflStore.at(i).data;
		delete[] pData;
	}
	mReflStore.clear();

	//mClockMaster.RemoveChild(mClockCijValues);
}

const string& XECsReussVoigt::GetClassName()const
{
	const static string name = "MStruct::XECsReussVoigt";
	
	return name;
}

void XECsReussVoigt::GetXECs(REAL & s1, REAL & s2, const REAL h, const REAL k, const REAL l)const
{
	VFN_DEBUG_MESSAGE("MStruct::XECsReussVoigt::GetXECs(...),h:"<<h<<",k:"<<k<<",l:"<<l<<":Begin",11);

	// Check if UnitCell was set 
	if (mpUnitCell==NULL) {
		cerr << "< MStruct::XECsReussVoigt::GetXECs()\n";
		cerr << "  No UnitCell Object. UnitCell Object should be set first.\n";
		cerr << ">" << endl;
		throw ObjCrystException("MStruct::XECsReussVoigt::GetXECs(): Object not initialised.");
	}
	
	// Calculate matrices and Voigt XECS
	CalcMatrices();
	
	// Check if the calculated matrices are still older than the oldest data in the store,
	// or if UnitCell lattice parameters were not changed if not mConsiderUnitCellFixed.
	// If something was changed clear the store.
	if ( mClockReflStore<mClockMatricesCalc
	     || (!mConsiderUnitCellFixed && mClockReflStore<mpUnitCell->GetClockLatticePar()) )
	{
		VFN_DEBUG_MESSAGE("MStruct::XECsReussVoigt::GetXECs(...):Clearing Stored data",11);

		// clear Reflection Store
		for(int i=0; i<mReflStore.size(); i++) {
			REAL *pData = (REAL*) mReflStore.at(i).data;
			delete[] pData;
		}
		mReflStore.clear();
		mClockReflStore.Reset();
	}

	// Try to find the given reflection in the Reflection Store
	int ind = mReflStore.find(h,k,l,0.);

	if ( ind >= 0 ) {
		
		VFN_DEBUG_MESSAGE("MStruct::XECsReussVoigt::GetXECs(...):Loading data from the Store",11);

		// Get data from the Reflection store
		REAL *pData = (REAL*) mReflStore.at(ind).data;
		s1 = pData[0];
		s2 = pData[1];

	} else {

		VFN_DEBUG_MESSAGE("MStruct::XECsReussVoigt::GetXECs(...):Calculating Reuss XECs",11);

		// Calculate Reuss XECs
		try {
			// Using simplified equations from [1], assuming only s1 = s2 = sigma (stress) nonzero
			//  Ref: [1] J.C.Popa, J.Appl.Cryst. 33 (2000) 103-107
				
			static const REAL rou[6] = { 1, 1, 1, 2, 2, 2 };
			static const REAL delta[6] = { 1, 1, 1, 0, 0, 0 };
			
			// Calculate hkl directional cosines (in crystallite coordinates)
			REAL A[3] = { h , k, l };
			mpUnitCell->MillerToOrthonormalCoords(A[0],A[1],A[2]);
			const REAL H = sqrt( pow(A[0],2) + pow(A[1],2) + pow(A[2],2) );
			A[0] /= H; A[1] /= H; A[2] /= H;
		
			REAL E[6] = { pow(A[0],2), pow(A[1],2), pow(A[2],2), A[1]*A[2], A[2]*A[0], A[0]*A[1] };
			
			s1 = 0.;
			s2 = 0.;
			
			// (eq. (7)<-(3a)<-(2a)<-(13), using (8) and consider s1=s2=sigma, all others si=0
			//     and B = ( sin(psi), 0, cos(psi) )
			// s1 = 1/2 * sum(i,j) (delta(j)-E(j))*E(i)*rou(i)*Sij*rou(j)
			// s2 = sum(i,j) (3*E(j)-delta(j))*E(i)*rou(i)*Sij*rou(j)
			for(int i=0; i<6; i++)
				for(int j=0; j<6; j++) {
					s1 += (delta[j]-E[j])*E[i]*rou[i]*mSijMatrix(i,j)*rou[j];
					s2 += (3.*E[j]-delta[j])*E[i]*rou[i]*mSijMatrix(i,j)*rou[j];
				}
				
			s1 /= 2.;
		}
		catch (std::exception &e) {
			cerr << "< MStruct::XECsReussVoigt::SGetXECs(...)\n";
			cerr << "  h: " << h <<  ", k: " << k << ", l: " << l << endl;
			cerr << "  Unexpected exception while computing Reuss XECs: " << e.what() << "\n";
			cerr << ">" << endl;
			throw ObjCrystException("MStruct::XECsReussVoigt::GetXECs(...): Unexpected error.");
		}

		VFN_DEBUG_MESSAGE("MStruct::XECsReussVoigt::GetXECs(...):Saving data in the Store",11);

		// Save data into the Reflection Store
		REAL * pData = new REAL[2];
		pData[0] = s1;
		pData[1] = s2;
		mReflStore.add((long)h,(long)k,(long)l,0.,(void*)pData);
		if ( mReflStore.size()==1 ) mClockReflStore.Click(); // the first data stored

	}
		
	s1 = (1.-mWeight) * s1 + mWeight * mS1Voigt;
	s2 = (1.-mWeight) * s2 + mWeight * mS2Voigt;

	VFN_DEBUG_MESSAGE("MStruct::XECsReussVoigt::GetXECs(...):End",11);
}

int XECsReussVoigt::SetUnitCell(const ObjCryst::UnitCell & uc, const bool fixed)
{
	VFN_DEBUG_MESSAGE("MStruct::XECsReussVoigt::SetUnitCell(...):Begin",11);
	
	if (mpUnitCell != NULL ) {
		// Old UnitCell settings has to be removed
		try {
			// Remove Refinable parameters connected with the old UnitCell
			for(int i=0; i<mCijValues.numElements(); i++) RemovePar(&GetPar(&mCijValues(i)));
			// Clear internal variables
			mpUnitCell = NULL;
			mConsiderUnitCellFixed = true;
			mGroupType = LAUE_UNDEFINED;
			mvCijConstantsNames.clear();
			mCijValues.resize(0);
			mClockCijValues.Click();
			mCijMatrix.resize(0,6);
			mSijMatrix.resize(0,6);
			mClockMatricesCalc.Reset();
			mMCij = MCij_undefined;
		}
		catch (std::exception &e) {
			cerr << "< MStruct::XECsReussVoigt::SetUnitCell(...)\n";
			cerr << "  Unexpected exception while clearing old UnitCell settings: " << e.what() << "\n";
			cerr << ">" << endl;
			throw ObjCrystException("MStruct::XECsReussVoigt::SetUnitCell(...): Unexpected error.");
  	}
	}
	
	// Set new UnitCell
	try {
		
		mpUnitCell = &uc;
		mClockMaster.Reset();
		
		// Analyse new UnitCell symmetry
		
			// Get cctbx::sgtbx::space_group
			const cctbx::sgtbx::space_group & sg = mpUnitCell->GetSpaceGroup().GetCCTbxSpg();
		
			// Get Laue group code
			const cctbx::sgtbx::matrix_group::code & LaueCode = sg.laue_group_type();
		
			using namespace cctbx::sgtbx::matrix_group;
		
			// Determine the Laue Symmetry from the Laue Code
			if      ( LaueCode == code_1b )
				mGroupType = LAUE_TRICLINIC;
			else if ( LaueCode == code_2_m ) {
				// Identify Unique axis
				
				// Note:: Simple method does't work now properly.
				// mpUnitCell->GetSpaceGroup().GetUniqueAxis();
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
					case 0 : { cerr << "Warning: Unique axis 'a' for monoclinic system not suppoerted." << endl;
										 mGroupType = LAUE_UNDEFINED;
									 } break;
					case 1 : mGroupType = LAUE_MONOCLINIC_AXIS_B; break;
					case 2 : mGroupType = LAUE_MONOCLINIC_AXIS_C; break;
					default : mGroupType = LAUE_UNDEFINED; break;
				}
				
			}
			else if ( LaueCode == code_mmm )
				mGroupType = LAUE_ORTHORHOMBIC;
			else if ( LaueCode == code_4_m )
				mGroupType = LAUE_TETRAGONAL_LOW;
			else if ( LaueCode == code_4_mmm )
				mGroupType = LAUE_TETRAGONAL_HIGH;
			else if ( LaueCode == code_3b )
				mGroupType = LAUE_TRIGONAL_LOW;
			else if ( LaueCode == code_3bm )
				mGroupType = LAUE_TRIGONAL_HIGH;
			else if ( LaueCode == code_6_m || LaueCode == code_6_mmm )
				mGroupType = LAUE_HEXAGONAL;
			else if ( LaueCode == code_m3b || LaueCode == code_m3bm )
				mGroupType = LAUE_CUBIC;
			else
				mGroupType = LAUE_UNDEFINED;

		if ( mGroupType == LAUE_UNDEFINED ) {
			cerr << "< MStruct::XECsReussVoigt::SetUnitCell(...)\n";
			cerr << " Can not identify (Laue) stiffness/complience tensor symmetry from the space group.\n";
			mpUnitCell->GetSpaceGroup().Print();
			cerr << ">" << endl;
			throw ObjCrystException("MStruct::XECsReussVoigt::SetUnitCell(...): Logical error.");
  	}
  
  	mConsiderUnitCellFixed = fixed;
  	
  	// Prepare vector of stiffness constants names
  	
  		// Prepare material constants tensor indices object
	  	switch (mGroupType) {
	  		
	  		case LAUE_TRICLINIC         : mMCij = MCij_triclinic; break;
	  		case LAUE_MONOCLINIC_AXIS_C : mMCij = MCij_monoclinic_axis_c; break;
	  		case LAUE_MONOCLINIC_AXIS_B : mMCij = MCij_monoclinic_axis_b; break;
	  		case LAUE_ORTHORHOMBIC      : mMCij = MCij_orthorhombic; break;
	  		case LAUE_TETRAGONAL_LOW    : mMCij = MCij_tetragonal_low; break;
	  		case LAUE_TETRAGONAL_HIGH   : mMCij = MCij_tetragonal_high; break;
	  		case LAUE_TRIGONAL_LOW      : mMCij = MCij_trigonal_low; break;
	  		case LAUE_TRIGONAL_HIGH     : mMCij = MCij_trigonal_high; break;
	  		case LAUE_HEXAGONAL         : mMCij = MCij_hexagonal; break;
	  		case LAUE_CUBIC             : mMCij = MCij_cubic; break;
	  		
	  		default :	mMCij = MCij_undefined;
	  							throw ObjCrystException("MStruct::XECsReussVoigt::SetUnitCell(...): Logical error.");
	  							break;
	  	}
  	
  		// Prepare vector of constants names
  		const vector< string > & vStrIJ = mMCij.indices_str();
  		const int n = vStrIJ.size();
  		for(int i=0; i<n; i++)
  			mvCijConstantsNames.push_back( "C" + vStrIJ[i] ); // "C" + "ij"
  	
  		// Resize vector of material Stiffness constants
  		mCijValues.resize(n);
  		
  		// Assign default Stiffness constants values
  		{
  			mCijValues = 0.;
  			// Isotropic (E = 200 GPa, ni = 0.33)
  			const REAL E = 200.;
  			const REAL ni = 0.33;
  			const REAL c12 = E*ni/(1+ni)/(1-2*ni);
  			const REAL c44 = E/2./(1+ni);
  			const REAL c11 = c12 + 2.*c44; // E*(1-ni)/(1+ni)/(1-2*ni);
  			  
  			const vector< pair<int,int> > & vInd = mMCij.indices();
  			vector< pair<int,int> >::const_iterator it;
  			
   			// c11
   			it = find( vInd.begin(), vInd.end(), make_pair(1,1) );
   			if ( it != vInd.end() ) mCijValues( it-vInd.begin() ) = c11;
   				else throw ObjCrystException("MStruct::XECsReussVoigt::SetUnitCell(...): Logical error.");
   			it = find( vInd.begin(), vInd.end(), make_pair(2,2) );
   			if ( it != vInd.end() ) mCijValues( it-vInd.begin() ) = c11;
   			it = find( vInd.begin(), vInd.end(), make_pair(3,3) );
   			if ( it != vInd.end() ) mCijValues( it-vInd.begin() ) = c11;
   			// c12
   			it = find( vInd.begin(), vInd.end(), make_pair(1,2) );
   			if ( it != vInd.end() ) mCijValues( it-vInd.begin() ) = c12;
   				else throw ObjCrystException("MStruct::XECsReussVoigt::SetUnitCell(...): Logical error.");
   			it = find( vInd.begin(), vInd.end(), make_pair(1,3) );
   			if ( it != vInd.end() ) mCijValues( it-vInd.begin() ) = c12;
   			it = find( vInd.begin(), vInd.end(), make_pair(2,3) );
   			if ( it != vInd.end() ) mCijValues( it-vInd.begin() ) = c12;
   			// c44
   			it = find( vInd.begin(), vInd.end(), make_pair(4,4) );
   			if ( it != vInd.end() ) mCijValues( it-vInd.begin() ) = c44;
   				else throw ObjCrystException("MStruct::XECsReussVoigt::SetUnitCell(...): Logical error.");
   			it = find( vInd.begin(), vInd.end(), make_pair(5,5) );
   			if ( it != vInd.end() ) mCijValues( it-vInd.begin() ) = c44;
   			it = find( vInd.begin(), vInd.end(), make_pair(6,6) );
   			if ( it != vInd.end() ) mCijValues( it-vInd.begin() ) = c44;
  		}
  		
  		// Initialise new Refinable parameters (Cij values)
  		for(int i=0; i<n; i++) {
  			RefinablePar tmp(mvCijConstantsNames[i],&mCijValues(i),0.,1.e6,
		  	gpRefParTypeMaterialData,
		  	REFPAR_DERIV_STEP_ABSOLUTE,true,true,false,false);
    		tmp.AssignClock(mClockCijValues);
    		tmp.SetDerivStep(1.);
    		this->AddPar(tmp);
  		}
  		
	}
	catch (std::exception &e) {
		cerr << "< MStruct::XECsReussVoigt::SetUnitCell(...)\n";
		cerr << "  Unexpected exception while setting new UnitCell: " << e.what() << "\n";
		cerr << ">" << endl;
		throw ObjCrystException("MStruct::XECsReussVoigt::SetUnitCell(...): Unexpected error.");
  }
	
	VFN_DEBUG_MESSAGE("MStruct::XECsReussVoigt::SetUnitCell(...):End, return " << mGroupType,11);
	
	return mGroupType;
}

void XECsReussVoigt::SetStiffnessConstants(const map< string, REAL > & CijValues)
{
	VFN_DEBUG_MESSAGE("MStruct::XECsReussVoigt::SetStiffnessConstants(...):Begin",11);
	
	// Check if UnitCell was set 
	if (mpUnitCell==NULL) {
		cerr << "< MStruct::XECsReussVoigt::SetStiffnessConstants()\n";
		cerr << "  No UnitCell Object. UnitCell Object should be set first.\n";
		cerr << ">" << endl;
		throw ObjCrystException("MStruct::XECsReussVoigt::SetStiffnessConstants(): Object not initialised.");
	}
	
	try {
		// Set all material Stiffness constants values to zero
		mCijValues = 0.;
		
		// Number of requred constants
		const int n = mCijValues.numElements();
		
		for(int i=0; i<n; i++) {
			// find constant value in the map
			map< string, REAL >::const_iterator it = CijValues.find( mvCijConstantsNames[i] );
			if ( it != CijValues.end() )
				mCijValues(i) = it->second; // set value
			else
				cerr << "Warning: " << mvCijConstantsNames[i] << " value not found, set to zero." << endl;
		}
		
		mClockCijValues.Click();
	}
	catch (std::exception &e) {
		cerr << "< MStruct::XECsReussVoigt::SetStiffnessConstants(...)\n";
		cerr << "  Unexpected exception while setting Stiffness constants values: " << e.what() << "\n";
		cerr << ">" << endl;
		throw ObjCrystException("MStruct::XECsReussVoigt::SetStiffnessConstants(...): Unexpected error.");
  }
  
  // Update/Calculate Matrices
  CalcMatrices();
  
  VFN_DEBUG_MESSAGE("MStruct::XECsReussVoigt::SetStiffnessConstants(...):End",11);
}

const vector< string > & XECsReussVoigt::GetStiffnessConstantsNames()const
{
	// Check if UnitCell was set 
	if (mpUnitCell==NULL) {
		cerr << "< MStruct::XECsReussVoigt::GetStiffnessConstantsNames()\n";
		cerr << "  No UnitCell Object. UnitCell Object should be set first.\n";
		cerr << ">" << endl;
		throw ObjCrystException("MStruct::XECsReussVoigt::GetStiffnessConstantsNames(): Object not initialised.");
	}
	
	return mvCijConstantsNames;
}

const CrystMatrix_REAL & XECsReussVoigt::GetCijMatrix()const
{
	// Check if UnitCell was set 
	if (mpUnitCell==NULL) {
		cerr << "< MStruct::XECsReussVoigt::GetCijMatrix()\n";
		cerr << "  No UnitCell Object. UnitCell Object should be set first.\n";
		cerr << ">" << endl;
		throw ObjCrystException("MStruct::XECsReussVoigt::GetCijMatrix(): Object not initialised.");
	}
	
	// Check if matrices were calculated and are up-to-date, if not, calculate them
	if ( mClockMatricesCalc<mClockCijValues  ) CalcMatrices();
	
	return mCijMatrix;
}

void XECsReussVoigt::SetParams(const REAL weight)
{
	mWeight = weight;
	mClockMaster.Click();
}

void XECsReussVoigt::CalcMatrices()const
{
	VFN_DEBUG_MESSAGE("MStruct::XECsReussVoigt::CalcMatrices():Begin",11);
	
	// If Matrices are newer than material Cij values than no calculation done
	if ( mClockMatricesCalc>mClockCijValues  ) {
		VFN_DEBUG_MESSAGE("MStruct::XECsReussVoigt::CalcMatrices():End",11);
		return;
	}
	
	// Check if UnitCell was set 
	if (mpUnitCell==NULL) {
		cerr << "< MStruct::XECsReussVoigt::CalcMatrices()\n";
		cerr << "  No UnitCell Object. UnitCell Object should be set first.\n";
		cerr << ">" << endl;
		throw ObjCrystException("MStruct::XECsReussVoigt::CalcMatrices(): Object not initialised.");
	}
	
	// Assign Cij matrix by material Cij values, make inversion and create Sij matrix
	
		// Assign 'basic' values of the symmetric tensor
		NEWMAT::SymmetricMatrix nCij(6);
		
		nCij = 0.;
		
		const vector< pair<int,int> > & vInd = mMCij.indices();
		
		for (unsigned int i=0; i<vInd.size(); i++)
			nCij(vInd[i].first,vInd[i].second) = mCijValues(i);
		
		try {
			// Assign symmetry equivalent/presribed values of the symmetric tensor
			switch (mGroupType) {
		  		
		  		case LAUE_TRICLINIC         : break; // nothing to do
		  		case LAUE_MONOCLINIC_AXIS_C : break; // nothing to do
		  		case LAUE_MONOCLINIC_AXIS_B : break; // nothing to do
		  		case LAUE_ORTHORHOMBIC      : break; // nothing to do
		  		case LAUE_TETRAGONAL_LOW    : 
		  			nCij(2,2) =  nCij(1,1);	nCij(2,3) =  nCij(1,3);	nCij(2,6) = -nCij(1,6);
		  			nCij(5,5) =  nCij(4,4);
		  			break;
		  		case LAUE_TETRAGONAL_HIGH   :
		  			nCij(2,2) =  nCij(1,1);	nCij(2,3) =  nCij(1,3);
		  			nCij(5,5) =  nCij(4,4);
		  			break;
		  		case LAUE_TRIGONAL_LOW      : 
		  			nCij(2,2) =  nCij(1,1);	nCij(2,3) =  nCij(1,3); nCij(2,4) = -nCij(1,4); nCij(2,5) = -nCij(1,5);
		  			nCij(4,6) = -nCij(1,5);
		  			nCij(5,5) =  nCij(4,4); nCij(5,6) =  nCij(1,4);
		  			nCij(6,6) = 0.5*(nCij(1,1)-nCij(1,2));
		  			break;
		  		case LAUE_TRIGONAL_HIGH     :
		  			nCij(2,2) =  nCij(1,1);	nCij(2,3) =  nCij(1,3); nCij(2,4) = -nCij(1,4);
		  			nCij(5,5) =  nCij(4,4); nCij(5,6) =  nCij(1,4);
		  			nCij(6,6) = 0.5*(nCij(1,1)-nCij(1,2));
		  			break;
		  		case LAUE_HEXAGONAL         :
		  			nCij(2,2) =  nCij(1,1);	nCij(2,3) =  nCij(1,3);
		  			nCij(5,5) =  nCij(4,4);
		  			nCij(6,6) = 0.5*(nCij(1,1)-nCij(1,2));
		  			break;
		  		case LAUE_CUBIC             :
		  			nCij(1,3) =  nCij(1,2);
		  			nCij(2,2) =  nCij(1,1);	nCij(2,3) =  nCij(1,2);
		  			nCij(3,3) =  nCij(1,1);
		  			nCij(5,5) =  nCij(4,4);
		  			nCij(6,6) =  nCij(4,4);
		  			break;
		  			
		  		default :
		  			throw ObjCrystException("MStruct::XECsReussVoigt::CalcMatrices(): Logical error.");
		  			break;
		  	}
			
			// Save Cij Matrix
			mCijMatrix.resize(6,6);
			for(int i=0; i<6; i++) for(int j=0; j<6; j++) mCijMatrix(i,j) = nCij(i+1,j+1);
			
			// Prepare Newmat Cij matrix for inversion (according to Popa's definition) 
			//      ==>   partitione matrix in four blocks and multiply elements by
			//            _   _
			//           | 1 2 |
		  //           |_2 4_|
		  //
		  for(int i=1; i<=3; i++) for(int j=4; j<=6; j++) nCij(i,j) *= 2.;
		  for(int i=4; i<=6; i++) for(int j=i; j<=6; j++) nCij(i,j) *= 4.;
		  
		  /*cout << "Cij' Matrix (GPa):" << endl;
		  for(int i=0; i<6; i++) {
				for(int j=0; j<6; j++)
					cout << setw(10) << nCij(i+1,j+1);
					cout << endl;
			}*/
			
		  // Calculate Sij (Popa) by Cij matrix inversion
		  nCij = nCij.i();
		  
		  // Save Sij Matrix
			mSijMatrix.resize(6,6);
			for(int i=0; i<6; i++) for(int j=0; j<6; j++) mSijMatrix(i,j) = nCij(i+1,j+1);
		  
		  /*cout << "Sij Matrix (1/TPa):" << endl;
		  for(int i=0; i<6; i++) {
				for(int j=0; j<6; j++)
					cout << setw(10) << mSijMatrix(i,j)*1e3;
				cout << endl;
			}*/
	
		// Calculate Voigt XECs
		//     Ref: [1] J.C.Popa, J.Appl.Cryst. 33 (2000) 103-107
		
		
		// c11v = (c11 + c22 + c33)/5 + 2*(c12+c13+c23+2*c44+2*c55+2*c66)/15    ([1] eq. 12a)
		const REAL c11v = (mCijMatrix(0,0) + mCijMatrix(1,1) + mCijMatrix(2,2))/5. +
		  2.*(mCijMatrix(0,1)+mCijMatrix(0,2)+mCijMatrix(1,2)+2.*mCijMatrix(3,3)+2.*mCijMatrix(4,4)+2.*mCijMatrix(5,5))/15.;
		// c12v = (c11+c22+c33-2*c44-2*c55-2*c66)/15 + 4*(c12+c13+c23)          ([1] eq. 12b)
		const REAL c12v = (mCijMatrix(0,0)+mCijMatrix(1,1)+mCijMatrix(2,2)-2.*mCijMatrix(3,3)-2.*mCijMatrix(4,4)-2.*mCijMatrix(5,5))/15. +
		  + 4.*(mCijMatrix(0,1)+mCijMatrix(0,2)+mCijMatrix(1,2))/15.;
		
		// s11v = (c11v+c12v)/(c11v*(c11v+c12v)-2*c12v^2)   ([1] eq. 30a)
		const REAL s11v = (c11v+c12v)/(c11v*(c11v+c12v)-2*pow(c12v,2));
		// s12v = -c12v/(c11v*(c11v+c12v)-2*c12v^2)   ([1] eq. 30b)
		const REAL s12v = -c12v/(c11v*(c11v+c12v)-2*pow(c12v,2));
		
		// Voigt XECs  ([1] simplified eq. 29)
		mS1Voigt = s12v;
		mS2Voigt = 2*(s11v-s12v);
		
		// calculation done
		mClockMatricesCalc.Click();
		
	}
	catch (std::exception &e) {
		cerr << "< MStruct::XECsReussVoigt::CalcMatrices()\n";
		cerr << "  Unexpected exception while computing Stiffness and Compliance Tensors: " << e.what() << "\n";
		cerr << ">" << endl;
		throw ObjCrystException("MStruct::XECsReussVoigt::CalcMatrices(): Unexpected error.");
  }
  
	VFN_DEBUG_MESSAGE("MStruct::XECsReussVoigt::CalcMatrices():End",11);
}

void XECsReussVoigt::InitParameters()
{
	{
    RefinablePar tmp("RV_weight",&mWeight,0.,1.,
		     gpRefParTypeScattDataCorrPos,
		     REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false);
    tmp.AssignClock(mClockMaster);
    tmp.SetDerivStep(1.e-3);
    this->AddPar(tmp);
  }
}

// an interpolation routine
double interp1(const vector<double> &vx, const vector<double> &vy,
	       const double x) {
  if (vx.size()<1) return 0.;

  if (vx.size()<2 && abs(vx.front()-x)<=numeric_limits<double>::epsilon())
    return vy.front();

  //if (x<vx.front() || x>=vx.back()) return 0.;
  if (x<=vx.front()) return vy.front();
  if (x>=vx.back()) return vy.back();

  // find interval
  vector<double>::const_iterator px=vx.begin(), py=vy.begin();
  double x0, x1, y0, y1;
  x1 = *px; y1 = *py;
  while (x>=x1) {
    px++; py++;
    x0 = x1; y0 = y1;
    x1 = *px; y1 = *py;
  }
  // linear interpolation
  return y0+(x-x0)/(x1-x0)*(y1-y0);
}

// a convolution routine
double conv2(const vector<double> &fx,const vector<double> &fy,
	     const vector<double> &gx,const vector<double> &gy,
	     const double x)
{
  // vectors of x have to be increasing
  if (gx.size()<1) return 0.;
  
  vector<double>::const_iterator pgx = gx.begin(), pgy = gy.begin();

  if (gx.size()<2) return interp1(fx,fy,x-(*pgx))*(*pgy); // only multiply

  // convolution
  
  double x0, y0, x1, y1, f0, f1, xi0, yi0, xi1, yi1, r = 0.;
  
  vector<double>::const_reverse_iterator pfx=fx.rbegin(), pfy=fy.rbegin();
  
  if (gx.back()<=x-*pfx || *pgx>=x-fx.front()) return 0.;

  // interpolation of f - find interval and interp. linearly
  while(x-*pgx>=*pfx) { pgx++; pgy++; }
  if (pgx==gx.end()) return 0.;
  x1  = *pgx++; y1  = *pgy++;
  xi1 = *pfx++; yi1 = *pfy++;
  xi0 = *pfx++; yi0 = *pfy++;
  while(x-x1<=xi0 && pfx!=fx.rend())
    { xi1 = xi0; yi1 = yi0; xi0 = *pfx++; yi0 = *pfy++; }
  if (pfx==fx.rend()) return 0.5*r; // TODO:: ???
  f1 = yi0+(x-x1-xi0)/(xi1-xi0)*(yi1-yi0);  // interp1(fx,fy,x-x1)
  while(pgx!=gx.end()) {
    x0 = x1; y0 = y1; f0 = f1;
    x1 = *pgx++; y1 = *pgy++;
    while(x-x1<=xi0 && pfx!=fx.rend())
      { xi1 = xi0; yi1 = yi0; xi0 = *pfx++; yi0 = *pfy++; }
    if (pfx==fx.rend()) return 1./6.*r; /// TODO:: ???
    f1 = yi0+(x-x1-xi0)/(xi1-xi0)*(yi1-yi0); // interp1(fx,fy,x-x1)
    //r += (f0*y0+f1*y1)*(x1-x0);
    r += (6.*f0*y0+3.*y0*(f1-f0)+3.*f0*(y1-y0)+2.*(f1-f0)*(y1-y0))*(x1-x0);
  }
  return 1./6.*r; // 0.5*r
}

// DoubleComponentReflectionProfile
DoubleComponentReflectionProfile::DoubleComponentReflectionProfile():
mWeight(0.0), mComponent1(NULL), mComponent2(NULL)
{
	InitParameters();
}

DoubleComponentReflectionProfile::DoubleComponentReflectionProfile(const DoubleComponentReflectionProfile &old):
mWeight(old.mWeight), mComponent1(old.mComponent1), mComponent2(old.mComponent2)
{
	InitParameters();
}

DoubleComponentReflectionProfile::~DoubleComponentReflectionProfile()
{
	if ( mComponent1 != NULL ) {
		mComponent1->DeRegisterClient(*this);
		this->RemoveSubRefObj(*mComponent1);
		// no delete - it should be deleted by a creater of the component object 
	}
	
	if ( mComponent2 != NULL ) {
		mComponent2->DeRegisterClient(*this);
		this->RemoveSubRefObj(*mComponent2);
		// no delete - it should be deleted by a creater of the component object 
	}
}

const string& DoubleComponentReflectionProfile::GetClassName () const
{
  const static string className="MStruct::DoubleComponentReflectionProfile";
  return className;
}

DoubleComponentReflectionProfile* DoubleComponentReflectionProfile::CreateCopy()const
{
	return new DoubleComponentReflectionProfile(*this);
}

CrystVector_REAL DoubleComponentReflectionProfile::GetProfile(const CrystVector_REAL &x, const REAL xcenter,
																												      const REAL h, const REAL k, const REAL l)
{
	// check if reflection profile components are set
	if(mComponent1==NULL || mComponent1==NULL) {
		cerr << "< DoubleComponentReflectionProfile::GetProfile(...)\n";
		cerr << "\t"<<"x="<<FormatHorizVector<REAL>(x,8,3)<<"\n";
		cerr << "\t"<<"xcenter="<<FormatFloat(xcenter,10,4)<<", h k l="<<h<<" "<<k<<" "<<l<<"\n";
		cerr << "\t"<<"Object ObjCryst::Name="<<this->GetName()<<"\n";
		cerr << "\t"<<"Reflection profile components are not set. >" << endl;
		throw ObjCrystException("DoubleComponentReflectionProfile::GetProfile(...): \
             Object not properly initialised.");
	}
	
	CrystVector_REAL result(x.numElements());
	CrystVector_REAL tmp(x.numElements());
	
	result = 0.;
	
	if (abs(1.-mWeight)>1e-6) { tmp = mComponent1->GetProfile(x,xcenter,h,k,l); tmp *= (1.-mWeight); result += tmp; }
	if (abs(mWeight)>1e-6) { tmp = mComponent2->GetProfile(x,xcenter,h,k,l); tmp *= mWeight; result += tmp; }
		
	return result;
}

REAL DoubleComponentReflectionProfile::GetFullProfileWidth(const REAL relativeIntensity, const REAL xcenter,
			    								 																 const REAL h, const REAL k, const REAL l)
{
	// check if reflection profile components are set
	if(mComponent1==NULL || mComponent1==NULL) {
		cerr << "< DoubleComponentReflectionProfile::GetFullProfileWidth(...)\n";
		cerr << "\t"<<"relativeIntensity="<<relativeIntensity;
		cerr << ", "<<"xcenter="<<FormatFloat(xcenter,10,4)<<", h k l="<<h<<" "<<k<<" "<<l<<"\n";
		cerr << "\t"<<"Object ObjCryst::Name="<<this->GetName()<<"\n";
		cerr << "\t"<<"Reflection profile components are not set. >" << endl;
		throw ObjCrystException("DoubleComponentReflectionProfile::GetFullProfileWidth(...): \
             Object not properly initialised.");
	}
	
	// if it is a "pure" function than this task is simple
	if (abs(1.-mWeight)<1e-6)
		return mComponent2->GetFullProfileWidth(relativeIntensity,xcenter,h,k,l);
	else if (abs(mWeight)<1e-6)
		 return mComponent1->GetFullProfileWidth(relativeIntensity,xcenter,h,k,l);
		 
	// mixed type reflection profile (mWeight > 0 && mWeight < 1)
	
	// calculate a guess for relative intensity of components
	REAL r1, r2;
	{
		// evaluate the functions on a finite interval to find a maximum
		// (take into consideration a possibility that the peak maximu doesn't have to be in the xcenter)
		
		REAL width;
		{
			REAL width1 = mComponent1->GetFullProfileWidth(0.01,xcenter,h,k,l);
			REAL width2 = mComponent2->GetFullProfileWidth(0.01,xcenter,h,k,l);
			width = (width1 > width2) ? width1 : width2;
		}
		
		CrystVector_REAL x(100), y1(100), y2(100);
		{ 
			REAL *p = x.data();
			REAL dx = width/(x.numElements()-1);
			for(int i=0; i<100; i++) { *p = -width/2 + i*dx; p++; }
		}
		y1 =  mComponent1->GetProfile(x,xcenter,h,k,l);
		y2 =  mComponent1->GetProfile(x,xcenter,h,k,l);
		
		// find a maximum of each component and of the mixed profile
		r1 = y1.max(); r2 = y2.max();
		y2 *= mWeight; y1 *= (1.-mWeight);
	
		if(abs(r1)<numeric_limits<REAL>::epsilon() || abs(r2)<numeric_limits<REAL>::epsilon()) {
			cerr << "< DoubleComponentReflectionProfile::GetFullProfileWidth(...)\n";
			cerr << "\t"<<"relativeIntensity="<<relativeIntensity;
			cerr << ", "<<"xcenter="<<FormatFloat(xcenter,10,4)<<", h k l="<<h<<" "<<k<<" "<<l<<"\n";
			cerr << "\t"<<"Object ObjCryst::Name="<<this->GetName()<<"\n";
			cerr << "\t"<<"Zero maximum value (a) reflection profile component. >" << endl;
			throw ObjCrystException("DoubleComponentReflectionProfile::GetFullProfileWidth(...): \
      	       Bad numeric value.");
		}
		
		y2 += y1;
		REAL Imax = y2.max();
		
		r1 = relativeIntensity/(1.-mWeight) * Imax/r1;
		r2 = relativeIntensity/mWeight * Imax/r2;
	}
	
	// calculate the first guess for components widths
	REAL width1 = mComponent1->GetFullProfileWidth(r1,xcenter,h,k,l);
	REAL width2 = mComponent1->GetFullProfileWidth(r2,xcenter,h,k,l);
	
	REAL width = (width1 > width2) ? width1 : width2;
	
	int counter = 0;
	const int n = 100;
	CrystVector_REAL x(n), y(n);
	
	while (counter<10) {
		REAL *p = x.data();
		REAL dx = width/(x.numElements()-1);
		for(int i=0; i<n; i++) { *p = -width/2 + i*dx; p++; }
		y = this->GetProfile(x,xcenter,h,k,l);
		REAL Imax = y.max();
		if ((y(0)/Imax <= relativeIntensity) && (y(n-1)/Imax <= relativeIntensity))
			break;
		else
			width *= 2;
	}
	
	if(counter>=10) {
		cerr << "< DoubleComponentReflectionProfile::GetFullProfileWidth(...)\n";
		cerr << "\t"<<"relativeIntensity="<<relativeIntensity;
		cerr << ", "<<"xcenter="<<FormatFloat(xcenter,10,4)<<", h k l="<<h<<" "<<k<<" "<<l<<"\n";
		cerr << "\t"<<"Object ObjCryst::Name="<<this->GetName()<<"\n";
		cerr << "\t"<<"Can not find x-interval large enough to reach requred relative intensity. >" << endl;
		throw ObjCrystException("DoubleComponentReflectionProfile::GetFullProfileWidth(...): \
    	       Can not reach convergence.");
	}
	
	// find the required profile width more precisly
	
	const int new_n = 4*n;
	
	// evaluate function
	x.resize(new_n);
	y.resize(new_n);
	{
		REAL *p = x.data();
		REAL dx = width/(x.numElements()-1);
		for(int i=0; i<x.numElements(); i++) { *p = -width/2 + i*dx; p++; }
	}
	y = this->GetProfile(x,xcenter,h,k,l);
	
	REAL Ilim = relativeIntensity*y.max();
	
	// left point
	REAL *py = y.data();
	int n1 = 0;
	py++; n1++;
	for(int i=1; i<y.numElements(); i++)
		if(*py>Ilim) break; else { n1++; py++; }
	
	REAL xleft;
	{
		REAL x0 = x(n1-1), x1 = x(n1); // TODO:: n1<x.numElements() ?
		REAL y0 = y(n1-1), y1 = y(n1);
		xleft = x0 + (x1-x0)/(y1-y0) * (Ilim-y0);
	}
	
	// right point
	py = y.data() + (y.numElements()-1);
	n1 = y.numElements()-1;
	py--; n1--;
	for(int i=1; i<y.numElements(); i++)
		if(*py>Ilim) break; else { n1--; py--; }
	
	REAL xright;
	{
		REAL x0 = x(n1), x1 = x(n1+1); // TODO:: n1>=0 ?
		REAL y0 = y(n1), y1 = y(n1+1);
		xright = x0 + (x1-x0)/(y1-y0) * (Ilim-y0);
	}
	
	return (xright-xleft);
}
			    								 
bool DoubleComponentReflectionProfile::IsAnisotropic()const
{
	// check if reflection profile components are set
	if(mComponent1==NULL || mComponent1==NULL) {
		cerr << "< DoubleComponentReflectionProfile::IsAnisotropic()\n";
		cerr << "\t"<<"Object ObjCryst::Name="<<this->GetName()<<"\n";
		cerr << "\t"<<"Reflection profile components are not set. >" << endl;
		throw ObjCrystException("DoubleComponentReflectionProfile::IsAnisotropic(): \
             Object not properly initialised.");
	}
	
	// if any of two components is anisotropic than this reflection profile is anisotropic
	return (mComponent1->IsAnisotropic() || mComponent2->IsAnisotropic());
}
	
void DoubleComponentReflectionProfile::SetComponents(ObjCryst::ReflectionProfile &comp1, ObjCryst::ReflectionProfile &comp2)
{
	VFN_DEBUG_ENTRY("DoubleComponentReflectionProfile::SetComponents(...):"<<comp1.GetName()<<", "<<comp2.GetName(),11)

	// remove old components if there are any
	if ( mComponent1 != NULL ) {
		mComponent1->DeRegisterClient(*this);
		this->RemoveSubRefObj(*mComponent1);
		mClockMaster.RemoveChild(mComponent1->GetClockMaster());
		// no delete - it should be deleted by a creater of the component object
	}
	
	if ( mComponent2 != NULL ) {
		mComponent2->DeRegisterClient(*this);
		this->RemoveSubRefObj(*mComponent2);
		mClockMaster.RemoveChild(mComponent2->GetClockMaster());
		// no delete - it should be deleted by a creater of the component object 
	}
	
  this->AddSubRefObj(comp1);
  comp1.RegisterClient(*this);
  mClockMaster.AddChild(comp1.GetClockMaster());
  
  this->AddSubRefObj(comp2);
  comp2.RegisterClient(*this);
  mClockMaster.AddChild(comp2.GetClockMaster());
  
	mComponent1 = &comp1;
	mComponent2 = &comp2;
	
	this->UpdateDisplay();
  VFN_DEBUG_EXIT("DoubleComponentReflectionProfile::SetComponents(...):"<<comp1.GetName()<<", "<<comp2.GetName(),11)
}

void DoubleComponentReflectionProfile::SetProfileParams(const REAL weight)
{
	mWeight = weight;
	mClockMaster.Click();
}

void DoubleComponentReflectionProfile::InitParameters()
{
	{
    RefinablePar tmp("Fraction",&mWeight,0,1.,
		     gpRefParTypeScattDataProfileType,
		     REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false);
    tmp.AssignClock(mClockMaster);
    tmp.SetDerivStep(1e-4);
    this->AddPar(tmp);
  }
}

////////////////////////////////////////////////////////////////////////
//
//    TextureOdfBase
//
////////////////////////////////////////////////////////////////////////

REAL TextureOdfBase::GetOdfValue(const CrystMatrix_REAL& A,
						        	  	  		 REAL phi1, REAL phi, REAL phi2,
																 const bool AnglesInitialised) const
{
	// This method should return the value of the model ODF.
	// The firt method parameter is the rotation matrix (3x3)
	// which defines the orientation of the selected crystallites
	// that fraction shold be returned. The rotation matrix is
	// the Eulerian rotation matrix representing rotations
	// in the (phi1,phi,phi2) notation. Values of the appropriate
	// angles (phi1,phi,phi2) have to be supplied if the 
	// AnglesInitialised parameter is true. Function
	// RotationTB::GetEulerAngles(phi1,phi,phi2,A) can be used
	// to get the Euler angles if they are not supplied.

	// We don't need the Euler angles.
	 //if(!AnglesInitialised) RotationTB::GetEulerAngles(phi1,phi,phi2,A);
	
	// The default model - no texture
	return 1.;
}

bool TextureOdfBase::IsOdfSymmetric() const
{
	// The default ODF model is nontextured powder which is naturally
	// symmetrised.
	return true;
}

bool TextureOdfBase::IsOdfNormalised() const
{
	// The default ODF model is nontextured powder which is properly
	// normalised.
	return true;
}

////////////////////////////////////////////////////////////////////////
//
//    TextureModelFiber
//
////////////////////////////////////////////////////////////////////////
TextureModelFiber::TextureModelFiber():
mHKL(3),mShapeFuncType(0),mWidthPar(20.*DEG2RAD),mShapePar(1.),mTiltParams(3),
mpUnitCell(0),mHKLDirections(0,3),mTiltDirection(3),mFuncNorm(1.)
{
	mClockMaster.AddChild(mClockFuncParams);
	mClockMaster.AddChild(mClockTiltParams);
	mClockMaster.AddChild(mClockHKL);
	
	mTiltParams = 0.;
	mTiltDirection = 0.; mTiltDirection(2) = 1.;

	// Initialise parameters and options
	Init();
}

TextureModelFiber::~TextureModelFiber()
{
	// remove all parameters
	RemovePar(&GetPar(&mWidthPar));
	RemovePar(&GetPar(&mShapePar));
	RemovePar(&GetPar(&mTiltParams(0)));
	RemovePar(&GetPar(&mTiltParams(1)));
	RemovePar(&GetPar(&mTiltParams(2)));
	// remove the current used unit cell if any 
	if (mpUnitCell!=0) mClockMaster.RemoveChild(mpUnitCell->GetClockMaster());	
}

void TextureModelFiber::Init()
{
	{ 
    string name = "Width_"+GetName();
		RefinablePar tmp(name,&mWidthPar,1.*DEG2RAD,360.*DEG2RAD,
				gpRefParTypeScattDataCorrInt,REFPAR_DERIV_STEP_ABSOLUTE,
				true,true,true,false,RAD2DEG);
		tmp.AssignClock(mClockFuncParams);
		tmp.SetDerivStep(0.5*DEG2RAD);
		//tmp.SetGlobalOptimStep(.05);
		this->AddPar(tmp);
  }
  { 
    string name = "Shape_"+GetName();
		RefinablePar tmp(name,&mShapePar,0.01,10.,
				gpRefParTypeScattDataCorrInt,REFPAR_DERIV_STEP_ABSOLUTE,
				true,true,true,false,1.);
		tmp.AssignClock(mClockFuncParams);
		tmp.SetDerivStep(0.02);
		//tmp.SetGlobalOptimStep(.05);
		this->AddPar(tmp);
  }
  { 
    string name = "Tilt_phi1_"+GetName();
		RefinablePar tmp(name,&mTiltParams(0),-180.*DEG2RAD,180.*DEG2RAD,
				gpRefParTypeScattDataCorrInt,REFPAR_DERIV_STEP_ABSOLUTE,
				true,true,true,true,360.*DEG2RAD);
		tmp.AssignClock(mClockTiltParams);
		tmp.SetDerivStep(3.*DEG2RAD);
		//tmp.SetGlobalOptimStep(.05);
		this->AddPar(tmp);
  }
  { 
    string name = "Tilt_phi_"+GetName();
		RefinablePar tmp(name,&mTiltParams(1),0.*DEG2RAD,180.*DEG2RAD,
				gpRefParTypeScattDataCorrInt,REFPAR_DERIV_STEP_ABSOLUTE,
				true,true,true,true,180.*DEG2RAD);
		tmp.AssignClock(mClockTiltParams);
		tmp.SetDerivStep(3.*DEG2RAD);
		//tmp.SetGlobalOptimStep(.05);
		this->AddPar(tmp);
  }
  { 
    string name = "Tilt_phi2_"+GetName();
		RefinablePar tmp(name,&mTiltParams(2),-180.*DEG2RAD,180.*DEG2RAD,
				gpRefParTypeScattDataCorrInt,REFPAR_DERIV_STEP_ABSOLUTE,
				true,true,true,true,360.*DEG2RAD);
		tmp.AssignClock(mClockTiltParams);
		tmp.SetDerivStep(3.*DEG2RAD);
		//tmp.SetGlobalOptimStep(.05);
		this->AddPar(tmp);
  }
}

void TextureModelFiber::SetUnitCell(const ObjCryst::UnitCell* unitcell)
{
	// Remove the current used unit cell if any 
	if (mpUnitCell!=0) mClockMaster.RemoveChild(mpUnitCell->GetClockMaster());
	// Add the new unit cell 
	mpUnitCell = unitcell;
	if (mpUnitCell==0) return;
	mClockMaster.AddChild(mpUnitCell->GetClockMaster());
}

void TextureModelFiber::SetFiberHKL(const int h, const int k, const int l, const bool forceFriedelLaw)
{
	mHKL(0) = h; mHKL(1) = k; mHKL(2) = l;
	mHKLDirections = CrystMatrix_REAL(0,3);
	// Calculate orthonormal coordinates of all equivalent HKL directions
	if(mpUnitCell==0) return;
	mHKLDirections = mpUnitCell->GetSpaceGroup().GetAllEquivRefl(h,k,l,false,forceFriedelLaw);
	// Transform HKL indexes of equvalent directions to orthonormal coordinates
	for(int ihkl=0; ihkl<mHKLDirections.rows(); ihkl++)
		mpUnitCell->MillerToOrthonormalCoords(mHKLDirections(ihkl,0),mHKLDirections(ihkl,1),mHKLDirections(ihkl,2));
	// Normalise them to unit vectors
	mHKLDirections *= 1./sqrt(pow(mHKLDirections(0,0),2) + pow(mHKLDirections(0,1),2) + pow(mHKLDirections(0,2),2));
	mClockHKL.Click();
}

REAL TextureModelFiber::GetOdfValue(const CrystMatrix_REAL& A,
						        	  	  		 		REAL phi1, REAL phi, REAL phi2,
																 		const bool AnglesInitialised) const
{
	// This method should return the value of the model ODF.
	// The firt method parameter is the rotation matrix (3x3)
	// which defines the orientation of the selected crystallites
	// that fraction shold be returned. The rotation matrix is
	// the Eulerian rotation matrix representing rotations
	// in the (phi1,phi,phi2) notation. Values of the appropriate
	// angles (phi1,phi,phi2) have to be supplied if the 
	// AnglesInitialised parameter is true. Function
	// RotationTB::GetEulerAngles(phi1,phi,phi2,A) can be used
	// to get the Euler angles if they are not supplied.

	// We don't need the Euler angles.
	 //if(!AnglesInitialised) RotationTB::GetEulerAngles(phi1,phi,phi2,A);
	
	// If no Unit cell has been set return default value (1.)
	if(mpUnitCell==0) return 1.;
	
	// Check if the reference direction (tilted normal direction)
	// was correctly calculated
	if(mClockTiltDirection<mClockTiltParams) {
		CrystMatrix_REAL a = RotationTB::GetEulerMatrix(mTiltParams(0),mTiltParams(1),mTiltParams(2));
		mTiltDirection(0) = a(0,2); mTiltDirection(1) = a(1,2); mTiltDirection(2) = a(2,2); 
		mClockTiltDirection.Click();
	}
	
	// Check if the fiber axis HKL indexes were set properly
	if(mHKLDirections.rows()==0) return 1.;
	
	// Check if the normalisation factor was correctly calculated
	if(mClockNormCalculated<mClockFuncParams)	CalcShapeFuncNorm();
	
	// Rotate the HKL fiber axis with the given crystalline orientation
	// and calculatete the angle between the tlted reference direction
	// and this rotated fiber direction
	
	REAL odf = 0.;
	const REAL Xn = pow(REAL(1-cos(mWidthPar/2.)),mShapePar)/M_LN2;
	
	CrystMatrix_REAL Q = RotationTB::MatrixTranspose(mHKLDirections);
	
	Q = RotationTB::MatrixMult(A,Q);
	
	for(int ihkl=0; ihkl<Q.cols(); ihkl++) {
		REAL f = 1.-(mTiltDirection(0)*Q(0,ihkl)+mTiltDirection(1)*Q(1,ihkl)+mTiltDirection(2)*Q(2,ihkl));
		f = exp(-pow(f,mShapePar)/Xn)/mFuncNorm;
		odf += f;
	}
	
	odf /= Q.cols();
	
	return odf;
}

REAL TextureModelFiber::CalcShapeFuncNorm() const
{
	// Calculate the normalisation factor for the used shape function
	
	if(mClockNormCalculated>mClockFuncParams) return mFuncNorm;
	
	// D.Simek et al., J.Appl.Cryst.(2006) 39, 487-501
	// func = exp(-(1-cos(om))^n/X, X = (1-cos(om_hw))/ln2^(1/n)
	
	const REAL Xn = pow(REAL(1-cos(mWidthPar/2.)),mShapePar)/M_LN2;
	
	REAL sum = 0.;
	REAL om;
	const int n = 100;
	
	for(int i=1; i<n; i++) {
		om = i*M_PI/n;
		sum += exp(-pow(REAL(1.-cos(om)),mShapePar)/Xn)*sin(om);
	}
	
	om = 0.;
	sum += 0.5*exp(-pow(REAL(1.-cos(om)),mShapePar)/Xn)*sin(om);
	om = M_PI;
	sum += 0.5*exp(-pow(REAL(1.-cos(om)),mShapePar)/Xn)*sin(om);
	
	mFuncNorm = 0.5*sum*M_PI/n; // integral[0->Pi] dpsi sin(psi) = 2
	mClockNormCalculated.Click();
	
	return mFuncNorm;
}

bool TextureModelFiber::IsOdfSymmetric() const
{
	return true;
}

bool TextureModelFiber::IsOdfNormalised() const
{
	return true;
}

////////////////////////////////////////////////////////////////////////
//
//    TextureOdfNumCalculator
//
////////////////////////////////////////////////////////////////////////
TextureOdfNumCalculator::TextureOdfNumCalculator():
mpOdfModel(0),mpUnitCell(0),mOdfGridParams(3,3),
mSetOdfProjectionIntegStep(3.*DEG2RAD)
{}

TextureOdfNumCalculator::~TextureOdfNumCalculator()
{
	//if (mpOdfModel!=0) mClockMaster.RemoveChild(mpOdfModel->GetClockMaster());
	mpOdfModel = 0;
	//if (mpUnitCell!=0) mClockMaster.RemoveChild(mpUnitCell->GetClockMaster());
	mpUnitCell = 0;
}

void TextureOdfNumCalculator::SetUnitCell(const ObjCryst::UnitCell* unitcell)
{
	// Remove the current used unit cell if any 
	if (mpUnitCell!=0) mClockMaster.RemoveChild(mpUnitCell->GetClockMaster());
	// Add the new unit cell 
	mpUnitCell = unitcell;
	if (mpUnitCell==0) return;
	mClockMaster.AddChild(mpUnitCell->GetClockMaster());
}

void TextureOdfNumCalculator::SetOdfModel(const TextureOdfBase* OdfModel)
{
	// Remove the current used ODF model if any
	if (mpOdfModel!=0) mClockMaster.RemoveChild(mpOdfModel->GetClockMaster());
	// Add new model
	mpOdfModel = OdfModel;
	if (mpOdfModel==0) return;
	mClockMaster.AddChild(mpOdfModel->GetClockMaster());
	
	// Set parametrs for the ODF grid representation
	 // phi1_min, phi1_max, phi1_step
	 mOdfGridParams(0,0) = 0.; mOdfGridParams(0,1) = 360.*DEG2RAD; mOdfGridParams(0,2) = 5.*DEG2RAD;
	 // phi_min, phi_max, phi_step
	 mOdfGridParams(1,0) = 0.; mOdfGridParams(1,1) = 180.*DEG2RAD; mOdfGridParams(1,2) = 3.*DEG2RAD;
	 // phi2_min, phi2_max, phi2_step
	 mOdfGridParams(2,0) = 0.; mOdfGridParams(2,1) = 360.*DEG2RAD; mOdfGridParams(2,2) = 5.*DEG2RAD;
}

REAL TextureOdfNumCalculator::GetOdfNorm() const
{
	if(mpOdfModel==0)
		throw ObjCrystException("TextureOdfNumCalculator: No ODF model!");
	if(mpOdfModel->IsOdfNormalised())
		return 1.;
	else {
		// We have to cvalculate norm ourselfs properly
		
		// Calc the ODF normalization factor
  	REAL odfnfactor = 0.;

  	// Integration of the ODF function over the whole Euler angles space
  	{
    	const REAL phi1min = mOdfGridParams(0,0); const REAL phi1max = mOdfGridParams(0,1); const REAL phi1step = mOdfGridParams(0,2);
    	const REAL phimin  = mOdfGridParams(1,0); const REAL phimax  = mOdfGridParams(1,1); const REAL phistep  = mOdfGridParams(1,2);
    	const REAL phi2min = mOdfGridParams(2,0); const REAL phi2max = mOdfGridParams(2,1); const REAL phi2step = mOdfGridParams(2,2);
    	const int nphi1 = int((phi1max - phi1min)/phi1step);
    	const int nphi  = int((phimax  - phimin )/phistep );
    	const int nphi2 = int((phi2max - phi2min)/phi2step);
    	
    	REAL sum3 = 0.;
    	for(int iphi=0; iphi<=nphi; iphi++)
    	{
      	REAL phi = phimin + iphi*phistep;
				REAL w = ((iphi==0) || (iphi==nphi)) ? 0.5 : 1.0;
      	REAL domega = sin(phi)/8./M_PI/M_PI;
      	REAL sum2 = 0.;
      	
      	for(int iphi2=0; iphi2<=nphi2; iphi2++)
      	{
					REAL phi2 = phi2min + iphi2*phi2step;
					REAL w = ((iphi2==0) || (iphi2==nphi2)) ? 0.5 : 1.0;

					REAL sum1 = 0.;
					
					for(int iphi1=0; iphi1<=nphi1; iphi1++)
					{
	  				REAL phi1 = phi1min + iphi1*phi1step;
						REAL w = ((iphi1==0) || (iphi1==nphi1)) ? 0.5 : 1.0; 
	  				CrystMatrix_REAL a = RotationTB::GetEulerMatrix(phi1,phi,phi2);
	  				sum1 += w*mpOdfModel->GetOdfValue(a,phi1,phi,phi2,true);
	  			}
					
					sum2 += w*sum1*phi1step;
				}
				
      	sum3 += w*sum2*phi2step*domega;
			}
    	odfnfactor = sum3*phistep;	
		}

  	return odfnfactor;
	}
}

void TextureOdfNumCalculator::ExportOdfXPert(ostream& os) const
{
	/* Exports the current ODF model in the Panalytical XPert like
	 * text format - not exactly the same used by XPert-Texture,
	 * however with almost similar data structure.
	 */
	
	// Print header
	os<<"File:"<<"\n";
	os<<"Sample:"<<"\n";
	os<<"Created:\t";
	
	try {
		time_t rawtime;
  	struct tm * timeinfo;

  	time (&rawtime);
  	timeinfo = localtime(&rawtime);
  	char* str_time = asctime(timeinfo);
  	char swday[64], smonth[64];
  	int day, hour, min, sec, year;
  	sscanf(str_time,"%3c %3c %d %d:%d:%d %d",swday,smonth,&day,&hour,&min,&sec,&year);
  	os<<setfill('0')<<setw(2)<<day<<setfill(' ')<<"-"<<smonth<<"-"<<year<<" "<<hour<<":"<<min<<"\n";
	}
	catch(std::exception) {
		// just show the error
		cerr<<"MStruct::TextureOdfNumCalculator::ExportOdfXPert: Error during time conversion!"<<endl;
	}
	
	os<<"Type:\t"<<"ODF data"<<"\n";
	os<<"Origin:\t"<<"MStruct::TextureOdfNumCalculator"<<"\n";
	os<<"\n";
	os<<"ODF calculation parameters"<<"\n";
	os<<"Used pole figures:\t"<<0<<"\n";
	os<<"\t#"<<"\th"<<"\tk"<<"\tl"<<"\tPsi min"<<"\tPsi max"<<"\n";
	os<<"Symmetries:"<<"\n";
	os<<"\t"<<"Crystal:\t";
	if(mpUnitCell!=0)
		os<<mpUnitCell->GetSpaceGroup().GetName()<<"\n";
	else
		os<<"none"<<"\n";
	os<<"\t"<<"Sample:"<<"\t"<<"Triclinic"<<"\n";
	os<<"Sample information:"<<"\n";
	os<<"\tGrid (\xb0):"<<"\t"<<"3.00 x 3.00"<<"\n";
	os<<fixed<<showpoint<<setprecision(3);
	os<<"\ta= "<<((mpUnitCell!=0) ? mpUnitCell->GetLatticePar(0) : 1.)<<"\n";
	os<<"\tb= "<<((mpUnitCell!=0) ? mpUnitCell->GetLatticePar(1) : 1.)<<"\n";
	os<<"\tc= "<<((mpUnitCell!=0) ? mpUnitCell->GetLatticePar(2) : 1.)<<"\n";
	os<<setprecision(2);
	os<<"\tAlpha= "<<((mpUnitCell!=0) ? mpUnitCell->GetLatticePar(3)*RAD2DEG : 90.)<<"\n";
	os<<"\tBeta= "<<((mpUnitCell!=0) ? mpUnitCell->GetLatticePar(3)*RAD2DEG : 90.)<<"\n";
	os<<"\tGamma= "<<((mpUnitCell!=0) ? mpUnitCell->GetLatticePar(3)*RAD2DEG : 90.)<<"\n";
	os<<"\n";
	
	// Print data
	os<<"ODF data table"<<endl;
	
	if(mpUnitCell==0 || mpOdfModel==0)
		return;
		
	const REAL odfnfactor = GetOdfNorm();
	
	os<<fixed<<showpoint<<setprecision(2);
	
	CrystVector_REAL PHI(int((mOdfGridParams(1,1)-mOdfGridParams(1,0))/mOdfGridParams(1,2))+1);
	for(int iphi=0; iphi<PHI.numElements(); iphi++)
			PHI(iphi) = mOdfGridParams(1,0) + iphi*mOdfGridParams(1,2);
			
	CrystVector_REAL Phi2(int((mOdfGridParams(2,1)-mOdfGridParams(2,0))/mOdfGridParams(2,2))+1);
	for(int iphi2=0; iphi2<Phi2.numElements(); iphi2++)
			Phi2(iphi2) = mOdfGridParams(2,0) + iphi2*mOdfGridParams(2,2);
		
	for(int iphi1=0; iphi1<=int((mOdfGridParams(0,1)-mOdfGridParams(0,0))/mOdfGridParams(0,2)); iphi1++) {
		const REAL phi1 = mOdfGridParams(0,0)+iphi1*mOdfGridParams(0,2);
		os<<"\n"<<"Phi1="<<setw(6)<<phi1*RAD2DEG<<"\n";
		os<<"Phi2\\PHI";
		for(int iphi=0; iphi<PHI.numElements(); iphi++)
			os<<setw(9)<<PHI(iphi)*RAD2DEG;
		os<<"\n";
		for(int iphi2=0; iphi2<Phi2.numElements(); iphi2++) {
			os<<setw(8)<<Phi2(iphi2)*RAD2DEG;
			for(int iphi=0; iphi<PHI.numElements(); iphi++) {
				CrystMatrix_REAL a = RotationTB::GetEulerMatrix(phi1,PHI(iphi),Phi2(iphi2)); 
				os<<setw(9)<<mpOdfModel->GetOdfValue(a,phi1,PHI(iphi),Phi2(iphi2),true)/odfnfactor;
			}
			os<<"\n";
		}
	}
}

void TextureOdfNumCalculator::PrepareForOdfProjectionCalc(const REAL h, const REAL k, const REAL l,
																							 						const bool forceFriedelLaw) const
{
	// Generate the rotation matrix transforming a vector in the given (hkl) direction
	// into the sample normal direction (001) for all equivalent orientations (if it
	// is necessary)
	
	if(mpUnitCell==0 || mpOdfModel==0) return;
	
	// chceck if we are already prepared
	if(mpUnitCell->GetClockMaster()<mClockOdfProjectionCalcPrepared &&
		 mpOdfModel->GetClockMaster()<mClockOdfProjectionCalcPrepared &&
		 mLastOdfProjectionHKL.numElements()==3 &&
		 fabs(mLastOdfProjectionHKL(0)-h)<1.e-4 && fabs(mLastOdfProjectionHKL(1)-k)<1.e-4 && fabs(mLastOdfProjectionHKL(2)-l)<1.e-4)
		return;
		
	// clear data
	mAuxInvRotMatHKL.clear();
	
	// set new HKL direction
	mLastOdfProjectionHKL.resize(3);
	mLastOdfProjectionHKL(0) = h; mLastOdfProjectionHKL(1) = k; mLastOdfProjectionHKL(2) = l;
	
	// the highest n-fold order of the symmetry of the (hkl) rotation axis
	// used to truncate the integration period
	int n_fold = 1;
	// only if the special direction (hkl indexes are integer values) something is done
	if(fabs(int(h)-h)<1.e-4 && fabs(int(k)-k)<1.e-4 && fabs(int(l)-l)<1.e-4) {
		// a reference to the CCTBX spacegroup
		const cctbx::sgtbx::space_group& cctbxsg = mpUnitCell->GetSpaceGroup().GetCCTbxSpg();
		
		// cctbx::sgtbx (hkl) vector
		const cctbx::sg_vec3 v = cctbx::sg_vec3(int(h),int(k),int(l));

		// test all symmetry operations
		for (std::size_t i_op = 0; i_op < cctbxsg.order_z(); i_op++) {
			// We are interested only in the rotational/inversion part of the symmetry matrix
			// A list of reciprocal space symmetry group rotations can be obtained as a list
			// of direct space matrices transposed. (see e.g. Giacovazzo - Fund. of Cryst.) 
			cctbx::sgtbx::rot_mx rm = cctbxsg(i_op).r().transpose();
			// test if the symmetry operation conserves the hkl axis
			 //cctbx::sgtbx::rot_mx const & arg1 = rm;
			 //cctbx::sg_vec3 const & arg2 = v; 
			cctbx::sg_vec3 u = rm*v - v;
			if (u.is_zero()==true && abs(rm.type())>n_fold) n_fold = rm.type();
			if (forceFriedelLaw || cctbxsg.is_centric()) {
				u = rm*v + v; // test in case of inversion center
				if (u.is_zero()==true && abs(rm.type())>n_fold) n_fold = rm.type();
			}
		}
	}
	
	// set the integration step and the number of points
	mOdfProjectionIntegPointsNb = int(2.*M_PI/n_fold/mSetOdfProjectionIntegStep);
	mOdfProjectionIntegStep = 2.*M_PI/n_fold/mOdfProjectionIntegPointsNb;

	// equivalent (hkl) directions necessary for calculation
	CrystMatrix_REAL equivHKLs;
	
	if(mpOdfModel->IsOdfSymmetric() &&
		 fabs(int(h)-h)<1.e-4 && fabs(int(k)-k)<1.e-4 && fabs(int(l)-l)<1.e-4 ) {
		// we need (or will use) just a one hkl for the calculation
		equivHKLs = CrystMatrix_REAL(1,3);
		equivHKLs(0,0) = h; equivHKLs(0,1) = k; equivHKLs(0,2) = l; 
	}
	else {
		// we need all equivalent diffractions (note: works correctly only for integer (hkl) indexes) 
		equivHKLs = mpUnitCell->GetSpaceGroup().GetAllEquivRefl(h,k,l,false,forceFriedelLaw);
	}
	
	const int nhkl = equivHKLs.rows();
	
	// transform (HKL)s directions into the orthogonal space and normalise them
	for(int ihkl=0; ihkl<nhkl; ihkl++)
		mpUnitCell->MillerToOrthonormalCoords(equivHKLs(ihkl,0),equivHKLs(ihkl,1),equivHKLs(ihkl,2));
	equivHKLs *= 1./sqrt(pow(equivHKLs(0,0),2)+pow(equivHKLs(0,1),2)+pow(equivHKLs(0,2),2));
	
	// clear old data
	mAuxInvRotMatHKL.clear();
	
	// create the list of rotation matrices transforming the (hkl) directions
	// into the sample normal direction (rotation around the axis hklx(001))
	for(int ihkl=0; ihkl<nhkl; ihkl++) {
		CrystMatrix_REAL a(3,3);
		
		const REAL c = equivHKLs(ihkl,2);
		const REAL s = sqrt(pow(equivHKLs(ihkl,0),2)+pow(equivHKLs(ihkl,1),2));
		
		// are vectors (001) and (hkl) parallel?
		if(s<1.e-4) {
			a = 0.;
			a(0,0) = c;
			a(1,1) = c;
			a(2,2) = c;
		} else {
			const REAL l1 =  equivHKLs(ihkl,1)/s;
			const REAL l2 = -equivHKLs(ihkl,0)/s;
			const REAL l3 = 0.;
			
			a(0,0) = c+(1-c)*pow(l1,2);
			a(1,1) = c+(1-c)*pow(l2,2);
			a(2,2) = c+(1-c)*pow(l3,2);
		
			a(0,1) = -s*l3 + (1-c)*l1*l2;
			a(1,0) =  s*l3 + (1-c)*l1*l2;
		
			a(0,2) =  s*l2 + (1-c)*l1*l3;
			a(2,0) = -s*l2 + (1-c)*l1*l3;
		
			a(1,2) = -s*l1 + (1-c)*l2*l3;
			a(2,1) =  s*l1 + (1-c)*l2*l3;
		}
		
		mAuxInvRotMatHKL.push_back(a);
	}

	mClockOdfProjectionCalcPrepared.Click();
}

REAL TextureOdfNumCalculator::CalcOdfProjection(const REAL psi, const REAL phi,
											   												const REAL h, const REAL k, const REAL l) const
{
	// Calculate the integral (projection) of the ODF function over the crystallites
	// having their crystal direction (hkl) parallel to the sample direction y =
	// ( cos(phi)*sin(psi), sin(phi)*sin(phi), cos(psi) ).
	//
	// If a special crystal direction (hkl indexes are integer values) is supplied
	// the calculation speed is optimised using the rotation symmetry of the (hkl)
	// crystal axis. To the contrary, the sample symmetry has not been yet implemented
	// and used here for the calcualtion speed optimization.
	//
	// Futhermore, if a special crystal direction (hkl indexes are integer values) is supplied
	// and the current used ODF model is not properly symmetrised (IsOdfSymmetric returns false)
	// the projection is averaged over all symmetry equivalent (hkl) directions and the Friedel
	// Symmetry is forced. If the non-integer (hkl) indexes are supplied, the equivalent directions
	// are not generated and if the ODF model is not correctly symmetric the calculated
	// projection doesn't respect crystal symmetry in any way.
	
	// if no Crystal or ODF model supplied return a primitive value
	if(mpUnitCell==0 || mpOdfModel==0) return 1.;

	// prepare object for the calculation	
	PrepareForOdfProjectionCalc(h,k,l,true);
	
	// if the (hkl) reflection makes no sence than return a primitive value 
	if(mAuxInvRotMatHKL.size()==0) return 1.;
	
	// integrate the supplied odf around the given (hkl) crystal direction
	// tilted from the sample normal direction by the angle psi and rotated
	// by the angle phi. Do this integration for all equivalent (hkl) directions
	// that were considered by PrepareForPFCalc method as necessary
	
	// We are using the convention for pole figures measurements
	
	// The (hkl) axis of the diffracting crystallites is parallel with the (001) z-axis
	// of the laboratory reference system. - The matrix denoted as inv(R) which rotates the (hkl)
	// vetor in the orthoromal coordianates into the vector (001) is prepared and stored
	// in the mAuxInvRotMatHKL structure for all necessary symmetry equivalent diffractions.
	//          (001)' = inv(R) * (hkl) 
	
	// The diffracting crystallites axis perpendicular to the (hkl) axis (we can call
	// the (hkl) axis now the z-axis and the perpendicular axis the x-axis) can be
	// rotated around the (hkl) z-axis in an arbitrary way. This is the projection
	// integration angle chi and this rotation is desribed later by the rotation matrix
	//             Dz(chi) - rotation around the z-axis by the angle chi
	
	// The sample normal direction (001) is rotated by angle psi and phi from the laboratory
	// (001) direction which is parallel with the difractiong crystallines (hkl) direction.
	// This is desribed by the Euler rotation matrix
	//                   B(0,psi,phi)
  // The angles psi and phi are the angles used for conventional pole figures measurements.
  
  // To conclude: The series of rotations transforming the axes of diffracting crystalline
  // starts with the rotation of the actual (hkl) direction into the laboratory z-axis (001). This
  // rotation is represented by matrix inv(R). Second step is the rotation of the crystalline aroud
  // this new z-axis by the integration angle chi represented by matrix Dz(chi) and at the end,
  // to transform the reference system of the actual diffracting crystalline into the sample reference
  // system, it si necesary to apply the rotation desribed by matrix B(0,psi,phi). Hence the orientation
  // of the diffracting crystallite is desribed by rotation matrix:
  //            A = B(0,psi,phi) * Dz(chi) * inv(R)
	
	// It is not essential if in this argumentation for simplification the notation of
	// "diffracting crystallites" was used. This was only . Derived matrix "A" desribes
	// a rotational transformation matrix from the crystal reference system into the sample
	// reference system if the choosen (hkl) direction should be parallel to the given sample
	// direction (psi,phi).
	
	const CrystMatrix_REAL B = RotationTB::GetEulerMatrix(0.,psi,phi);
	
	// Dz(chi) matrix
	CrystMatrix_REAL Dz(3,3);
	Dz = 0.; Dz(2,2) = 1.;
		
	REAL sum = 0.;
	
	// for all equivalent (hkl) directions
	for(int ihkl=0; ihkl<mAuxInvRotMatHKL.size(); ihkl++) {
		// integration (periodic function - 0..2Pi/n_fold)
	 	for(int ichi=0; ichi<mOdfProjectionIntegPointsNb; ichi++) {
	 		// create the Dz(chi) matrix
	 		const REAL c = cos(ichi*mOdfProjectionIntegStep);
	 		const REAL s = sin(ichi*mOdfProjectionIntegStep);
	 		Dz(0,0) = c; Dz(1,1) = c; Dz(0,1) = s; Dz(1,0) = -s;   
	 		// the orientation matrix
	 		CrystMatrix_REAL A = RotationTB::MatrixMult(B,Dz);
			A = RotationTB::MatrixMult(A,mAuxInvRotMatHKL[ihkl]);
	 		// value of the ODF
	 		sum += mpOdfModel->GetOdfValue(A);
	 	}
	}
	
	sum *= 1./mOdfProjectionIntegPointsNb/mAuxInvRotMatHKL.size(); // 1/(2Pi) * integral
	
	return sum;
}

CrystMatrix_REAL TextureOdfNumCalculator::CalcPFProjection(const CrystVector_REAL& psi,
																													 const CrystVector_REAL& phi,
												 																	 const REAL h, const REAL k, const REAL l) const
{
	// Calculate the Pole Figure projection - output matrix has indexing: result(psi,phi)
	CrystMatrix_REAL result(psi.numElements(),phi.numElements());
	
	// default value
	result = 1.;
	
	// the calculation idea is same as in CalcOdfProjection method. This method
	// is only slightly optimised version of that more general one - hence
	// see there if you need more info
	
	// if no Crystal or ODF model supplied return a primitive value
	if(mpUnitCell==0 || mpOdfModel==0) return result;

	// prepare object for the calculation	
	PrepareForOdfProjectionCalc(h,k,l,true);
	
	// if the (hkl) reflection makes no sence than return a primitive value 
	if(mAuxInvRotMatHKL.size()==0) return result;
	
	// prepare auxilliary data structures (matrices) to accelerate the calculation
	std::vector< CrystMatrix_REAL > vDziR; // Dz(chi)*inv(R_hkl)
	std::vector< CrystMatrix_REAL > vCx; // Cx(psi)
	std::vector< CrystMatrix_REAL > vBz; // Bz(phi)
	
	// Dz(chi) matrix
	CrystMatrix_REAL Dz(3,3);
	Dz = 0.; Dz(2,2) = 1.;
	
	// for all equivalent (hkl) directions
	for(int ihkl=0; ihkl<mAuxInvRotMatHKL.size(); ihkl++) {
		// for all integration points
	 	for(int ichi=0; ichi<mOdfProjectionIntegPointsNb; ichi++) {
	 		// create the Dz(chi) matrix
	 		const REAL c = cos(ichi*mOdfProjectionIntegStep);
	 		const REAL s = sin(ichi*mOdfProjectionIntegStep);
	 		Dz(0,0) = c; Dz(1,1) = c; Dz(0,1) = s; Dz(1,0) = -s;
	 		// a part of the orentation matrix Dz(chi)*inv(R_hkl)
	 		vDziR.push_back(RotationTB::MatrixMult(Dz,mAuxInvRotMatHKL[ihkl]));
	 	}
	}
	
	// Cx(psi) matrix
	CrystMatrix_REAL Cx(3,3);
	Cx = 0.; Cx(0,0) = 1.;
	
	// for all required psi angles
	for(int ipsi=0; ipsi<psi.numElements(); ipsi++) {
		// Cx(psi) matrix
		const REAL c = cos(psi(ipsi));
	 	const REAL s = sin(psi(ipsi));
	 	Cx(1,1) = c; Cx(2,2) = c; Cx(1,2) = s; Cx(2,1) = -s;
	 	vCx.push_back(Cx);
	}
	
	// Bz(phi) matrix
	CrystMatrix_REAL Bz(3,3);
	Bz = 0.; Bz(2,2) = 1.;
	
	// for all required phi angles
	for(int iphi=0; iphi<phi.numElements(); iphi++) {
		// Bz(phi) matrix
		const REAL c = cos(phi(iphi));
	 	const REAL s = sin(phi(iphi));
	 	Bz(0,0) = c; Bz(1,1) = c; Bz(0,1) = s; Bz(1,0) = -s;
	 	vBz.push_back(Bz);
	}
	
	// PF claculation
	
	result = 0.;
	
	// for all required psi angles
	for(int ipsi=0; ipsi<psi.numElements(); ipsi++) {
	 	// for all equivalent (hkl) directions
		for(int ihkl=0; ihkl<mAuxInvRotMatHKL.size(); ihkl++) {
			// in each integration point
			for(int ichi=0; ichi<mOdfProjectionIntegPointsNb; ichi++) {
	 			// a part of the orientation matrix
	 			CrystMatrix_REAL A1 = RotationTB::MatrixMult(vCx[ipsi],vDziR[ihkl*mOdfProjectionIntegPointsNb+ichi]);
	 			// for all required phi angles
	 			for(int iphi=0; iphi<phi.numElements(); iphi++) {
	 				// full orientation matrix
	 				CrystMatrix_REAL A = RotationTB::MatrixMult(vBz[iphi],A1);
	 				// value of the ODF
	 				result(ipsi,iphi) += mpOdfModel->GetOdfValue(A);
	 			} // iphi
			} // ichi
		} // ihkl
	} // ipsi
	
	result *= 1./mOdfProjectionIntegPointsNb/mAuxInvRotMatHKL.size(); // 1/(2Pi) * integral
	
	return result;
}

void TextureOdfNumCalculator::ExportPFProjectionXPert(ostream& os,
																											const CrystVector_REAL& psi,
																											const CrystVector_REAL& phi,
																											const CrystMatrix_REAL& pf_data,
																											const REAL scale) const
{
	/* Exports the supplied Pole Figure data in the Panalytical XPert like
	 * text format - not exactly the same used by XPert-Texture,
	 * however with almost similar data structure. (It is suitable
	 * to use this method immediately after using this object to
	 * calculate the Pole Figure data by the CalcPFProjection or
	 * CalcOdfProjection method. In this way more correct information
	 * is written in the output stream but still in most cases information
	 * in the header part of the exported stream is only ballast.)
	 */
	
	// Print header
	os<<"File:"<<"\n";
	os<<"Sample:"<<"\n";
	os<<"Created:\t";
	
	try {
		time_t rawtime;
  	struct tm * timeinfo;

  	time (&rawtime);
  	timeinfo = localtime(&rawtime);
  	char* str_time = asctime(timeinfo);
  	char swday[4], smonth[4];
  	int day, hour, min, sec, year;
  	sscanf(str_time,"%3c %3c %d %d:%d:%d %d",swday,smonth,&day,&hour,&min,&sec,&year);
  	os<<setfill('0')<<setw(2)<<day<<setfill(' ')<<"-"<<smonth<<"-"<<year<<" "<<hour<<":"<<min<<"\n";
	}
	catch(std::exception) {
		// just show the error
		cerr<<"MStruct::TextureOdfNumCalculator::ExportPFProjectionXPert: Error during time conversion!"<<endl;
	}

	if(psi.numElements()==0 || phi.numElements()==0 ||
		 psi.numElements()!=pf_data.rows() || phi.numElements()!=pf_data.cols()) {
		os<<"error:"<<"\t"<<"no data or bad data format"<<"\n";
		return;		
	}

	os<<"Type:\t"<<"Calculated pole figure"<<"\n";
	os<<"Origin:\t"<<"MStruct::TextureOdfNumCalculator"<<"\n";
	os<<"\n";
	os<<"Goniometer radius (mm):"<<"\t"<<"320"<<"\n";
	os<<"Sample stage:"<<"\t\t\t"<<"Other"<<"\n";
	os<<"Receiving slit (mm):"<<"\t\t"<<"0.00"<<"\n";
	os<<"Divergence slit (mm):"<<"\t"<<"0.30"<<"\n";
	os<<"Distance focus mask (mm):"<<"\t"<<"100"<<"\n";
	os<<"X-ray tube anode:"<<"\t"<<"Cu"<<"\n";
	os<<"Tube focus:"<<"\t\t"<<"LFF"<<"\n";
	os<<"Generator (kV):"<<"\t"<<"40"<<"\n";
	os<<"Generator (mA):"<<"\t"<<"35"<<"\n";
	os<<"Wavelength (A):"<<"\t"<<"1.5419"<<"\n";
	os<<"\n";
	os<<fixed<<showpoint<<setprecision(2);
	os<<"\t"<<setw(6)<<"Start"<<setw(7)<<"End"<<setw(7)<<"Step"<<"\n";
	os<<"Psi:"<<"\t"<<setw(6)<<psi(0)*RAD2DEG<<setw(7)<<psi(psi.numElements()-1)*RAD2DEG;
	os<<setw(7)<<((psi.numElements()>1) ? psi(1)-psi(0) : 0.)*RAD2DEG<<"\n";
	os<<"Phi:"<<"\t"<<setw(6)<<phi(0)*RAD2DEG<<setw(7)<<phi(phi.numElements()-1)*RAD2DEG;
	os<<setw(7)<<((phi.numElements()>1) ? phi(1)-phi(0) : 0.)*RAD2DEG<<"\n";
	os<<"\n";
	int H = 0, K = 0, L = 0; 
	if(mLastOdfProjectionHKL.numElements()==3) {
		H = int(mLastOdfProjectionHKL(0));
		K = int(mLastOdfProjectionHKL(1));
		L = int(mLastOdfProjectionHKL(2));
	}
	os<<"h k l:"<<"\t\t\t"<<H<<" "<<K<<" "<<L<<"\n";
	os<<"\n";

	os<<"2Theta (\xb0):"<<"\t\t"<<"10.0000"<<"\n";
	os<<"Time per step (s):"<<"\t"<<scale<<"\n";
	os<<"Sample oscillation (mm):"<<"\t"<<"0"<<"\n";
	
	os<<"\n";
	
	// Print data
	
	os<<fixed<<showpoint<<setprecision(2);
	
	// the first line with psi values
	os<<"Phi\\Psi";
	for(int ipsi=0; ipsi<psi.numElements(); ipsi++)	os<<setw(9)<<psi(ipsi)*RAD2DEG;
	os<<"\n";
	// data table
	for(int iphi=0; iphi<phi.numElements(); iphi++) {
		os<<setw(7)<<phi(iphi)*RAD2DEG;
		for(int ipsi=0; ipsi<psi.numElements(); ipsi++)	os<<setw(9)<<pf_data(ipsi,iphi)*scale;
		os<<"\n";
	}
	
}
							 						
////////////////////////////////////////////////////////////////////////
//
//    TextureModelHKL
//
////////////////////////////////////////////////////////////////////////

TextureModelHKL::TextureModelHKL()
:mNbComponents(0),mGlobalPhi1(0.),mGlobalPhi(0.),mGlobalPhi2(0.),
mpCrystal(0),mbforceFriedelSymmetry(false)
{
	this->InitParameters();
	mClockMaster.AddChild(mClockParams);
}

TextureModelHKL::TextureModelHKL(const TextureModelHKL& old)
:mNbComponents(0),
mTextureCalculator(old.mTextureCalculator),
mGlobalPhi1(old.mGlobalPhi1),mGlobalPhi(old.mGlobalPhi),mGlobalPhi2(old.mGlobalPhi2),
mpCrystal(old.mpCrystal),mbforceFriedelSymmetry(old.mbforceFriedelSymmetry)
{
	mClockMaster.AddChild(mClockParams);
	this->InitParameters();
	// add components
	for(int icomp=0; icomp<old.mNbComponents; icomp++)
		AddComponent(*old.mComponentTextureModels_phi[icomp],
								 *old.mComponentTextureModels_phi2[icomp],
								 *old.mComponentTextureModels_phi1[icomp],
								 *old.mComponentFractions[icomp],
								 *old.mComponentTextureParams[icomp]);
}

TextureModelHKL::~TextureModelHKL()
{
	// remove all components
	for(int icomp=0; icomp<mNbComponents; icomp++)
		RemoveComponent(icomp);
}

const string & TextureModelHKL::GetClassName() const
{
	const static string className="MStruct::TextureModelHKL";
   return className;
}

void TextureModelHKL::SetCrystal(const Crystal& crystal, const bool forceFriedelSymmetry)
{
  mpCrystal = &crystal;
	mbforceFriedelSymmetry = forceFriedelSymmetry;
}

/** Prepare the object for the odf function calculation - generete all equivalent
	 * crystal orientations for all HKL components, transform their coordinates to
	 * the normal orthogonal system; create texture reference systems (tilted and rotated
	 * base) for all components; calculate the odf normalization factor.
	 */
void TextureModelHKL::PrepareForOdfCalc()
{
	// Generete all equivalent crystal orientations for all HKL components

	// Clear the lists of all equivalent crystal orientations (HKLz x HKLx)
	mComponentEquivHKLz.clear();
	mComponentEquivHKLx.clear();

	// We need a reference to the CCTBX spacegroup - we will use it to get 
	// symmetry operations to generate all equivalent crystal orientations
	const cctbx::sgtbx::space_group cctbxsg = mpCrystal->GetSpaceGroup().GetCCTbxSpg();

	for(int icomp=0; icomp<mNbComponents; icomp++) {
		// Generate list of all crystal (HKLz x HKLx) equivalent orientations for this component
		std::vector< CrystVector_REAL > equivHKLz;
		std::vector< CrystVector_REAL > equivHKLx;
	
		// Basic orientation
		CrystVector_REAL hklz = *mComponentHKLz[icomp];
		CrystVector_REAL hklx = *mComponentHKLx[icomp];
	
		// Add basic orientation system to the list of equivalent systems
		equivHKLz.push_back(hklz);
		equivHKLx.push_back(hklx);

		// If the algorithm has to force the Friedel symmetry and the crystal
		// has no center of symmetry we should add an inversion to basic system manually
		if(mbforceFriedelSymmetry && !mpCrystal->GetSpaceGroup().IsCentrosymmetric()) {
			CrystVector_REAL t = hklz; t *= -1.; equivHKLz.push_back(t);
			                 t = hklx; t *= -1.; equivHKLx.push_back(t); }

		// Loop over all space group symmetry elements and generate from the basic system
		// new equivalent systems.
		// Vectors of starting orientation
		cctbx::sgtbx::tr_vec uz(int((*mComponentHKLz[icomp])(0)),
														int((*mComponentHKLz[icomp])(1)),
														int((*mComponentHKLz[icomp])(2)));
		cctbx::sgtbx::tr_vec ux(int((*mComponentHKLx[icomp])(0)),
														int((*mComponentHKLx[icomp])(1)),
														int((*mComponentHKLx[icomp])(2)));
		for (std::size_t i_op = 0; i_op < cctbxsg.order_z(); i_op++) {
			// We are interested only in the rotational/inversion part of the symmetry matrix
			// A list of reciprocal space symmetry group rotations can be obtained as a list
			// of direct space matrices transposed. (see e.g. Giacovazzo - Fund. of Cryst.) 
			const cctbx::sgtbx::rot_mx rm = cctbxsg(i_op).r().transpose();
			// Vectors of created orientation
			cctbx::sgtbx::tr_vec vz = uz;
			cctbx::sgtbx::tr_vec vx = ux;

			// Sequence of applaying the symmetry operation
			while (1) {
				// apply rotation
				vz = rm.multiply(vz);
				vx = rm.multiply(vx);
				// test for end of symmetry operation sequence
				if (vz[0]==uz[0] && vz[1]==uz[1] && vz[2]==uz[2] &&
					  vx[0]==ux[0] && vx[1]==ux[1] && vx[2]==ux[2]) break; 
				// test if this orietation is equivalent to some already created orientations
				bool notequiv = true;
				for(size_t ieq=0; ieq<equivHKLz.size(); ieq++) {
					notequiv = int(equivHKLz[ieq](0))!=vz[0] || int(equivHKLz[ieq](1))!=vz[1] ||
						int(equivHKLz[ieq](2))!=vz[2] || int(equivHKLx[ieq](0))!=vx[0] ||
						int(equivHKLx[ieq](1))!=vx[1] || int(equivHKLx[ieq](2))!=vx[2];
					if(!notequiv) break;
				} // ieq
				// test if a new nonequivalent symmetrical configuration found
				if (notequiv) {
					CrystVector_REAL t(3);
					t(0) = vz[0]; t(1) = vz[1]; t(2) = vz[2];
					equivHKLz.push_back(t);
					t(0) = vx[0]; t(1) = vx[1]; t(2) = vx[2];
					equivHKLx.push_back(t);
				}
				// if the algorithm has to force the Friedel symmetry and the crystal
				// has no center of symmetry we should test for inversion symmetry manually
				if(mbforceFriedelSymmetry && !mpCrystal->GetSpaceGroup().IsCentrosymmetric()) {
					// test if this orietation is equivalent to some already created orientations
					bool notequiv = true;
					for(size_t ieq=0; ieq<equivHKLz.size(); ieq++) {
						notequiv = int(equivHKLz[ieq](0))!=-vz[0] || int(equivHKLz[ieq](1))!=-vz[1] ||
							int(equivHKLz[ieq](2))!=-vz[2] || int(equivHKLx[ieq](0))!=-vx[0] ||
							int(equivHKLx[ieq](1))!=-vx[1] || int(equivHKLx[ieq](2))!=-vx[2];
						if(!notequiv) break;
					} // ieq
					// test if a new nonequivalent symmetrical configuration found
					if (notequiv) {
						// new symmetrical configuration found
						CrystVector_REAL t(3);
						t(0) = -vz[0]; t(1) = -vz[1]; t(2) = -vz[2];
						equivHKLz.push_back(t);
						t(0) = -vx[0]; t(1) = -vx[1]; t(2) = -vx[2];
						equivHKLx.push_back(t);
					}
				} // test for inversion symmetry
			} // while (1)
		} // i_op
		
		// Transform axes coordinates from a generally rectilinear reciprocal space 
    // system (A*) to classical orthogonal system (E1* (=E1))
		for(size_t isys=0; isys<equivHKLz.size(); isys++) {
			mpCrystal->MillerToOrthonormalCoords(equivHKLz[isys](0),
																					 equivHKLz[isys](1),
																					 equivHKLz[isys](2));
			mpCrystal->MillerToOrthonormalCoords(equivHKLx[isys](0),
																					 equivHKLx[isys](1),
																					 equivHKLx[isys](2));
		}
		
		// Normalise orietations vectors to unity and orthogonalise the x-vector
		for(size_t isys=0; isys<equivHKLz.size(); isys++) {
			REAL t = sqrt(equivHKLz[isys](0)*equivHKLz[isys](0) +
							 			equivHKLz[isys](1)*equivHKLz[isys](1) +
							 			equivHKLz[isys](2)*equivHKLz[isys](2));
			if (fabs(t)<1.e-7) throw ObjCrystException("TextureModelHKL: Wrong crystal orientation!");
			
			equivHKLz[isys] *= 1/t;
			
			     t = sqrt(equivHKLx[isys](0)*equivHKLx[isys](0) +
							 			equivHKLx[isys](1)*equivHKLx[isys](1) +
							 			equivHKLx[isys](2)*equivHKLx[isys](2));
							 			
			if (fabs(t)<1.e-7) throw ObjCrystException("TextureModelHKL: Wrong crystal orientation!");
			equivHKLx[isys] *= 1/t;
			
			// generally the input HKLz and HKLx vectors don't have to be perpendicular,
			// hence we would like to ensure their orthogonality here
			REAL nt =  equivHKLz[isys](0)*equivHKLx[isys](0) +
							 	 equivHKLz[isys](1)*equivHKLx[isys](1) +
							 	 equivHKLz[isys](2)*equivHKLx[isys](2);
			if ((1.-fabs(nt))<1.e-4) throw ObjCrystException("TextureModelHKL: Wrong crystal orientation!");
			CrystVector_REAL w = equivHKLz[isys]; // HKLz
			w *= -nt; // - dot(HKLz,HKLx) * HKLz
			equivHKLx[isys] += w; // - dot(HKLz,HKLx) * HKLz + HKLx  
			equivHKLx[isys] *= 1./sqrt(1-nt);
		}

		/*cout<<"TextureModelHKL: "<<GetName()<<", symmetry systems ("<<equivHKLz.size()<<"):\n";
		for(size_t i=0; i<equivHKLz.size(); i++) {
			cout<<"\t vz:";
			cout<<equivHKLz[i](0)<<" "<<equivHKLz[i](1)<<" "<<equivHKLz[i](2)<<"\n";
			cout<<"\t vx:";
			cout<<equivHKLx[i](0)<<" "<<equivHKLx[i](1)<<" "<<equivHKLx[i](2)<<"\n";
			cout << "\n";
		}*/
		
		// Store calc results
		mComponentEquivHKLz.push_back(equivHKLz);
		mComponentEquivHKLx.push_back(equivHKLx);
		
	} // icomp
	
	// create texture reference systems (properly rotated base) for all components
	mComponentReferenceSystems.clear();
	
	CrystMatrix_REAL identity(3,3);
	identity = 0.; identity(0,0) = 1.; identity(1,1) = 1.; identity(2,2) = 1.;
	
	for(int icomp=0; icomp<mNbComponents; icomp++) {
		CrystMatrix_REAL referenceSystem = identity;
		// get the Euler angles from component params
		
		
		// save the reference system
		mComponentReferenceSystems.push_back(referenceSystem);
	} // icomp
	
	// calculate the odf normalization factor
}

void TextureModelHKL::AddComponent(const int texture_model_phi,
																	 const int texture_model_phi2,
																	 const int texture_model_phi1,
																	 const REAL fraction,
																	 const CrystVector_REAL& HKLz,
																	 const CrystVector_REAL& HKLx,
																	 const CrystVector_REAL &params)
{
	// number of parameters used for this component (phi1 part)
	int nb_params_phi1 = 0;
		
	switch (texture_model_phi1) {
		case TEXTURE_MODEL_HKL_GAUSS:
			// phi1 func. weight, tilt, fwhm
			nb_params_phi1 += 3;
			break;
		case TEXTURE_MODEL_HKL_SIMEK:
			// phi1 func. weight, tilt, fwhm, shape factor
			nb_params_phi1 += 4;
			break;
		case TEXTURE_MODEL_HKL_BIRKHOLZ:
			// phi1 func. weight, tilt, texture order
			nb_params_phi1 += 3;
			break;
		case TEXTURE_MODEL_HKL_AXIAL:
			break;
		default:
			throw ObjCrystException("Wrong Texture model identificator!");
			break;
	}
	
	// number of parameters used for this component (phi part)
	int nb_params_phi = 0;
		
	switch (texture_model_phi) {
		case TEXTURE_MODEL_HKL_GAUSS:
			// phi func. weight, tilt, fwhm
			nb_params_phi += 3;
			break;
		case TEXTURE_MODEL_HKL_SIMEK:
			// phi func. weight, tilt, fwhm, shape factor
			nb_params_phi += 4;
			break;
		case TEXTURE_MODEL_HKL_BIRKHOLZ:
			// phi func. weight, tilt, texture order
			nb_params_phi += 3;
			break;
		case TEXTURE_MODEL_HKL_AXIAL:
			break;
		default:
			throw ObjCrystException("Wrong Texture model identificator!");
			break;
	}

	// number of parameters used for this component (all params: phi1 + phi + phi2)
	int nb_params = nb_params_phi1 + nb_params_phi; 

	switch (texture_model_phi2) {
		case TEXTURE_MODEL_HKL_GAUSS:
			// phi2 func. weight, rotation, fwhm
			nb_params += 3;
			break;
		case TEXTURE_MODEL_HKL_SIMEK:
			// phi2 func. weight, rotation, fwhm, shape factor
			nb_params += 4;
			break;
		case TEXTURE_MODEL_HKL_BIRKHOLZ:
			// phi2 func. weight, rotation, texture order
			nb_params += 3;
			break;
		case TEXTURE_MODEL_HKL_AXIAL:
			break;
		default:
			throw ObjCrystException("Wrong Texture model identificator!");
			break;
	}
	
	// set new component model type
	mComponentTextureModels_phi1.push_back(new int(texture_model_phi1));
	mComponentTextureModels_phi.push_back(new int(texture_model_phi));
	mComponentTextureModels_phi2.push_back(new int(texture_model_phi2));

	// set component HKL axes 
	mComponentHKLz.push_back(new CrystVector_REAL(3));
	// if HKLz suppled use the supplied value
	if( HKLz.numElements()==3 && abs(HKLz(0))+abs(HKLz(1))+abs(HKLz(2))>0. )
		(*mComponentHKLz[mNbComponents]) = HKLz;
	else {
		(*mComponentHKLz[mNbComponents])(0) = 0.;
		(*mComponentHKLz[mNbComponents])(1) = 0.;
		(*mComponentHKLz[mNbComponents])(2) = 1.;
	}
	// if HKLx suppled use the supplied value
	mComponentHKLx.push_back(new CrystVector_REAL(3));
	if( HKLx.numElements()==3 && abs(HKLx(0))+abs(HKLx(1))+abs(HKLx(2))>0.)
		(*mComponentHKLx[mNbComponents]) = HKLx;
	else {
		(*mComponentHKLx[mNbComponents])(0) = 1.;
		(*mComponentHKLx[mNbComponents])(1) = 0.;
		(*mComponentHKLx[mNbComponents])(2) = 0.;
	}

	// set new component fraction
	mComponentFractions.push_back(new REAL(fraction));
	
	// set component prameters (use default values if params. not supplied)
	mComponentTextureParams.push_back(new CrystVector_REAL(nb_params));
	
	// deafault values - rewriten later if supplied
	
	// phi1 weight - for all models (with exception of the axial model)
	if( texture_model_phi1 != TEXTURE_MODEL_HKL_AXIAL )
		(*mComponentTextureParams[mNbComponents])(0) = 1.;
	// texture function rotation - also for all type of models (with exception of the axial model)
	if( texture_model_phi1 != TEXTURE_MODEL_HKL_AXIAL )
		(*mComponentTextureParams[mNbComponents])(1) = 0.;
	
	// other phi1 func. params - specific for each model
	
	// phi1 fwhm - TEXTURE_MODEL_HKL_GAUSS, TEXTURE_MODEL_HKL_SIMEK
	if( texture_model_phi1 == TEXTURE_MODEL_HKL_GAUSS ||
			texture_model_phi1 == TEXTURE_MODEL_HKL_SIMEK)
		(*mComponentTextureParams[mNbComponents])(2) = 20.*DEG2RAD;
	
	// phi1 - texture shape par. - n - TEXTURE_MODEL_HKL_SIMEK, TEXTURE_MODEL_BIRKHOLZ
	if( texture_model_phi1 == TEXTURE_MODEL_HKL_SIMEK ||
			texture_model_phi1 == TEXTURE_MODEL_HKL_BIRKHOLZ) {
		int tn = (texture_model_phi1==TEXTURE_MODEL_HKL_SIMEK) ? 3 : 2;
		(*mComponentTextureParams[mNbComponents])(tn) = 2.; 
	}
	
	// phi weight - for all models (with exception of the axial model)
	if( texture_model_phi != TEXTURE_MODEL_HKL_AXIAL )
		(*mComponentTextureParams[mNbComponents])(nb_params_phi1) = 1.;
	// texture function rotation - also for all type of models (with exception of the axial model)
	if( texture_model_phi != TEXTURE_MODEL_HKL_AXIAL )
		(*mComponentTextureParams[mNbComponents])(nb_params_phi1+1) = 0.;

	// other rotation func. params - specific for each model
	
	// phi func. fwhm - TEXTURE_MODEL_HKL_GAUSS, TEXTURE_MODEL_HKL_SIMEK
	if( texture_model_phi == TEXTURE_MODEL_HKL_GAUSS ||
			texture_model_phi == TEXTURE_MODEL_HKL_SIMEK) 
		(*mComponentTextureParams[mNbComponents])(nb_params_phi1+2) = 20.*DEG2RAD;
	
	// phi func. - texture shape par. - n - TEXTURE_MODEL_HKL_SIMEK, TEXTURE_MODEL_BIRKHOLZ
	if( texture_model_phi == TEXTURE_MODEL_HKL_SIMEK ||
			texture_model_phi == TEXTURE_MODEL_HKL_BIRKHOLZ) {
		int tn = (texture_model_phi==TEXTURE_MODEL_HKL_SIMEK) ? 3 : 2;
		(*mComponentTextureParams[mNbComponents])(nb_params_phi1+tn) = 2.;
	}

	// phi2 weight - for all models (with exception of the axial model)
	if( texture_model_phi2 != TEXTURE_MODEL_HKL_AXIAL )
		(*mComponentTextureParams[mNbComponents])(nb_params_phi1+nb_params_phi) = 1.;
	// texture function rotation - also for all type of models (with exception of the axial model)
	if( texture_model_phi2 != TEXTURE_MODEL_HKL_AXIAL )
		(*mComponentTextureParams[mNbComponents])(nb_params_phi1+nb_params_phi+1) = 0.;

	// other rotation func. params - specific for each model
	
	// phi2 func. fwhm - TEXTURE_MODEL_HKL_GAUSS, TEXTURE_MODEL_HKL_SIMEK
	if( texture_model_phi2 == TEXTURE_MODEL_HKL_GAUSS ||
			texture_model_phi2 == TEXTURE_MODEL_HKL_SIMEK) 
		(*mComponentTextureParams[mNbComponents])(nb_params_phi1+nb_params_phi+2) = 20.*DEG2RAD;
	
	// phi2 func. - texture shape par. - n - TEXTURE_MODEL_HKL_SIMEK, TEXTURE_MODEL_BIRKHOLZ
	if( texture_model_phi2 == TEXTURE_MODEL_HKL_SIMEK ||
			texture_model_phi2 == TEXTURE_MODEL_HKL_BIRKHOLZ) {
		int tn = (texture_model_phi2==TEXTURE_MODEL_HKL_SIMEK) ? 3 : 2;
		(*mComponentTextureParams[mNbComponents])(nb_params_phi1+nb_params_phi+tn) = 2.;
	}
	
	// if component parameters supplied use them
	if(params.numElements()>0) {
		if((*mComponentTextureParams[mNbComponents]).numElements()!=params.numElements())
			throw ObjCrystException("Wrong number of parameters supplied!");
		(*mComponentTextureParams[mNbComponents]) = params;
	}

	// increase the number of components used
	mNbComponents++;
	
	// Initialize texture component model parameters
	InitComponentParameters(mNbComponents-1);
	
	// set parametrs to be used if it is essential
  if ( texture_model_phi1 != TEXTURE_MODEL_HKL_AXIAL ) this->GetPar(&mGlobalPhi1).SetIsUsed(true);
	if ( texture_model_phi != TEXTURE_MODEL_HKL_AXIAL ) this->GetPar(&mGlobalPhi).SetIsUsed(true);
	if ( texture_model_phi2 != TEXTURE_MODEL_HKL_AXIAL ) this->GetPar(&mGlobalPhi2).SetIsUsed(true);
	
	// let the params. clocks know that something was changed
  mClockParams.Click();
}

void TextureModelHKL::RemoveComponent(const int comp_nb)
{
	// remove the component - remove all parameters connected with this component
	// clear memory allocated for them and remove their cobntainers from vector structures
	
	// test component number
	if ( (comp_nb < 0) || (comp_nb >= mNbComponents) )
		throw ObjCrystException("Wrong texture component number!");
		
	// remove all parameters connected with this component
	{
		RefinablePar *ppar = 0;
		
		// fraction
		ppar = 0;
		try {
			ppar = &GetPar(mComponentFractions[comp_nb]);
			if ( ppar != 0) // it looks like really existing parameter
			RemovePar(ppar);
		}
		catch (std::exception e) {;} // no problem if failed
		
		// texture parameters
		for(int i=0; i<(*mComponentTextureParams[comp_nb]).numElements(); i++) {
			ppar = 0;
			try {
				ppar = &GetPar(&((*mComponentTextureParams[comp_nb])(i)));
				if ( ppar != 0) // it looks like really existing parameter
				RemovePar(ppar);
			}
			catch (std::exception e) {;} // no problem if failed
		}
	}

	// clear memory allocated for texture parameters
	{
		// phi1 type
		delete mComponentTextureModels_phi1[comp_nb];
		// phi type
		delete mComponentTextureModels_phi[comp_nb];
		// phi2 type
		delete mComponentTextureModels_phi2[comp_nb];
		// fraction
		delete mComponentFractions[comp_nb];
		// HKL axes
		delete mComponentHKLz[comp_nb];
		delete mComponentHKLx[comp_nb];
		// texture parameters
		delete mComponentTextureParams[comp_nb];
	}

	// remove parameters containers from vector structures
	{
		// phi1 type
		mComponentTextureModels_phi1.erase(mComponentTextureModels_phi1.begin()+comp_nb);
		// phi type
		mComponentTextureModels_phi.erase(mComponentTextureModels_phi.begin()+comp_nb);
		// phi2 type
		mComponentTextureModels_phi2.erase(mComponentTextureModels_phi2.begin()+comp_nb);
		// fraction
		mComponentFractions.erase(mComponentFractions.begin()+comp_nb);
		// HKL axes
		mComponentHKLz.erase(mComponentHKLz.begin()+comp_nb);
		mComponentHKLx.erase(mComponentHKLx.begin()+comp_nb);
		// texture parameters
		mComponentTextureParams.erase(mComponentTextureParams.begin()+comp_nb);
	}

	// decrease number of parameters
	mNbComponents--;
}

void TextureModelHKL::InitComponentParameters(const int comp_nb)
{
	// test component number
	if ( (comp_nb < 0) || (comp_nb >= mNbComponents) )
		throw ObjCrystException("Wrong texture component number!");
	
	// initialise all parrameters connected with this component
	
	// get string containing component number
	string str_comp_nb;
	{
		ostringstream ss;
		ss << comp_nb;
		str_comp_nb = ss.str();
	}
	
	{ // initialise component fraction
		string name = "Fraction_"+str_comp_nb;
    RefinablePar tmp(name,mComponentFractions[comp_nb],0.,1.,
		     gpRefParTypeScattDataCorrInt,REFPAR_DERIV_STEP_ABSOLUTE,
		     true,true,true,false,1.);
    tmp.AssignClock(mClockParams);
    tmp.SetDerivStep(1e-2);
    //tmp.SetGlobalOptimStep(.05);
    this->AddPar(tmp);
  }

	// number of parameters used for this fraction (phi1 part)
	int nb_params_phi1 = 0;
	
	switch (*mComponentTextureModels_phi1[comp_nb]) {
		case TEXTURE_MODEL_HKL_GAUSS:
			// tilt func. weight, tilt, fwhm
			nb_params_phi1 += 3;
			break;
		case TEXTURE_MODEL_HKL_SIMEK:
			// tilt func. weight, tilt, fwhm, shape factor
			nb_params_phi1 += 4;
			break;
		case TEXTURE_MODEL_HKL_BIRKHOLZ:
			// tilt func. weight, tilt, texture order
			nb_params_phi1 += 3;
			break;
		case TEXTURE_MODEL_HKL_AXIAL:
			break;
		default:
			throw ObjCrystException("Wrong Texture model identificator!");
			break;
	}

	// number of parameters used for this fraction (phi part)
	int nb_params_phi = 0;
	
	switch (*mComponentTextureModels_phi[comp_nb]) {
		case TEXTURE_MODEL_HKL_GAUSS:
			// tilt func. weight, tilt, fwhm
			nb_params_phi += 3;
			break;
		case TEXTURE_MODEL_HKL_SIMEK:
			// tilt func. weight, tilt, fwhm, shape factor
			nb_params_phi += 4;
			break;
		case TEXTURE_MODEL_HKL_BIRKHOLZ:
			// tilt func. weight, tilt, texture order
			nb_params_phi += 3;
			break;
		case TEXTURE_MODEL_HKL_AXIAL:
			break;
		default:
			throw ObjCrystException("Wrong Texture model identificator!");
			break;
	}


	// phi1 weight - for all models (with exception of the axial model)
	if( *mComponentTextureModels_phi1[comp_nb] != TEXTURE_MODEL_HKL_AXIAL ) {
		string name = "Phi1_Weight_"+str_comp_nb;
    RefinablePar tmp(name,&((*mComponentTextureParams[comp_nb])(0)),0.,1.,
		     gpRefParTypeScattDataCorrInt,REFPAR_DERIV_STEP_ABSOLUTE,
		     true,true,true,false,1.);
    tmp.AssignClock(mClockParams);
    tmp.SetDerivStep(1e-2);
    //tmp.SetGlobalOptimStep(.05);
    this->AddPar(tmp);
  }

	// texture function phi1 - also for all type of models (with exception of the axial model)
	if( *mComponentTextureModels_phi1[comp_nb] != TEXTURE_MODEL_HKL_AXIAL ) {
		string name = "Phi1_"+str_comp_nb;
    RefinablePar tmp(name,&((*mComponentTextureParams[comp_nb])(1)),-180.*DEG2RAD,180.*DEG2RAD,
		     gpRefParTypeScattDataCorrInt,REFPAR_DERIV_STEP_ABSOLUTE,
		     true,true,true,true,RAD2DEG,360.*DEG2RAD);
    tmp.AssignClock(mClockParams);
    tmp.SetDerivStep(1.*DEG2RAD);
    //tmp.SetGlobalOptimStep(.05);
    this->AddPar(tmp);
  }

	// other phi1 func. params - specific for each model
	
	// phi1 fwhm - TEXTURE_MODEL_HKL_GAUSS, TEXTURE_MODEL_HKL_SIMEK
	if( *mComponentTextureModels_phi1[comp_nb] == TEXTURE_MODEL_HKL_GAUSS ||
			*mComponentTextureModels_phi1[comp_nb] == TEXTURE_MODEL_HKL_SIMEK) 
	{
		string name = "Phi1_fwhm_"+str_comp_nb;
		RefinablePar tmp(name,&((*mComponentTextureParams[comp_nb])(2)),1.*DEG2RAD,360.*DEG2RAD,
				gpRefParTypeScattDataCorrInt,REFPAR_DERIV_STEP_ABSOLUTE,
				true,true,true,false,RAD2DEG);
		tmp.AssignClock(mClockParams);
		tmp.SetDerivStep(0.5*DEG2RAD);
		//tmp.SetGlobalOptimStep(.05);
		this->AddPar(tmp);
	}
	
	// phi1 - texture shape par. - n - TEXTURE_MODEL_HKL_SIMEK, TEXTURE_MODEL_BIRKHOLZ
	if( *mComponentTextureModels_phi1[comp_nb] == TEXTURE_MODEL_HKL_SIMEK ||
			*mComponentTextureModels_phi1[comp_nb] == TEXTURE_MODEL_HKL_BIRKHOLZ) 
	{
		string name = "Phi1_n_"+str_comp_nb;
		int tn = (*mComponentTextureModels_phi1[comp_nb]==TEXTURE_MODEL_HKL_SIMEK) ? 3 : 2;
		RefinablePar tmp(name,&((*mComponentTextureParams[comp_nb])(tn)),0.1,10.,
				gpRefParTypeScattDataCorrInt,REFPAR_DERIV_STEP_ABSOLUTE,
				true,true,true,false,1.);
		tmp.AssignClock(mClockParams);
		tmp.SetDerivStep(0.1);
		//tmp.SetGlobalOptimStep(.05);
		this->AddPar(tmp);
	}
	
	// phi weight - for all models (with exception of the axial model)
	if( *mComponentTextureModels_phi[comp_nb] != TEXTURE_MODEL_HKL_AXIAL ) {
		string name = "Phi_Weight_"+str_comp_nb;
    RefinablePar tmp(name,&((*mComponentTextureParams[comp_nb])(nb_params_phi1)),0.,1.,
		     gpRefParTypeScattDataCorrInt,REFPAR_DERIV_STEP_ABSOLUTE,
		     true,true,true,false,1.);
    tmp.AssignClock(mClockParams);
    tmp.SetDerivStep(1e-2);
    //tmp.SetGlobalOptimStep(.05);
    this->AddPar(tmp);
  }

	// texture function phi - also for all type of models (with exception of the axial model)
	if( *mComponentTextureModels_phi[comp_nb] != TEXTURE_MODEL_HKL_AXIAL ) {
		string name = "Phi_"+str_comp_nb;
    RefinablePar tmp(name,&((*mComponentTextureParams[comp_nb])(nb_params_phi1+1)),0.,180.*DEG2RAD,
		     gpRefParTypeScattDataCorrInt,REFPAR_DERIV_STEP_ABSOLUTE,
		     true,true,true,true,RAD2DEG,180.*DEG2RAD);
    tmp.AssignClock(mClockParams);
    tmp.SetDerivStep(1.*DEG2RAD);
    //tmp.SetGlobalOptimStep(.05);
    this->AddPar(tmp);
  }

	// other rotation func. params - specific for each model
	
	// phi func. fwhm - TEXTURE_MODEL_HKL_GAUSS, TEXTURE_MODEL_HKL_SIMEK
	if( *mComponentTextureModels_phi[comp_nb] == TEXTURE_MODEL_HKL_GAUSS ||
			*mComponentTextureModels_phi[comp_nb] == TEXTURE_MODEL_HKL_SIMEK) 
	{
		string name = "Phi_fwhm_"+str_comp_nb;
		RefinablePar tmp(name,&((*mComponentTextureParams[comp_nb])(nb_params_phi1+2)),1.*DEG2RAD,360.*DEG2RAD,
				gpRefParTypeScattDataCorrInt,REFPAR_DERIV_STEP_ABSOLUTE,
				true,true,true,false,RAD2DEG);
		tmp.AssignClock(mClockParams);
		tmp.SetDerivStep(0.5*DEG2RAD);
		//tmp.SetGlobalOptimStep(.05);
		this->AddPar(tmp);
	}
	
	// phi func. - texture shape par. - n - TEXTURE_MODEL_HKL_SIMEK, TEXTURE_MODEL_BIRKHOLZ
	if( *mComponentTextureModels_phi[comp_nb] == TEXTURE_MODEL_HKL_SIMEK ||
			*mComponentTextureModels_phi[comp_nb] == TEXTURE_MODEL_HKL_BIRKHOLZ) 
	{
		string name = "Phi_n_"+str_comp_nb;
		int tn = (*mComponentTextureModels_phi[comp_nb]==TEXTURE_MODEL_HKL_SIMEK) ? 3 : 2; 
		RefinablePar tmp(name,&((*mComponentTextureParams[comp_nb])(nb_params_phi1+tn)),0.1,10.,
				gpRefParTypeScattDataCorrInt,REFPAR_DERIV_STEP_ABSOLUTE,
				true,true,true,false,1.);
		tmp.AssignClock(mClockParams);
		tmp.SetDerivStep(0.1);
		//tmp.SetGlobalOptimStep(.05);
		this->AddPar(tmp);
	}

	// phi2 weight - for all models (with exception of the axial model)
	if( *mComponentTextureModels_phi2[comp_nb] != TEXTURE_MODEL_HKL_AXIAL ) {
		string name = "Phi2_Weight_"+str_comp_nb;
    RefinablePar tmp(name,&((*mComponentTextureParams[comp_nb])(nb_params_phi1+nb_params_phi)),0.,1.,
		     gpRefParTypeScattDataCorrInt,REFPAR_DERIV_STEP_ABSOLUTE,
		     true,true,true,false,1.);
    tmp.AssignClock(mClockParams);
    tmp.SetDerivStep(1e-2);
    //tmp.SetGlobalOptimStep(.05);
    this->AddPar(tmp);
  }

	// texture function phi2 - also for all type of models (with exception of the axial model)
	if( *mComponentTextureModels_phi2[comp_nb] != TEXTURE_MODEL_HKL_AXIAL ) {
		string name = "Phi2_"+str_comp_nb;
    RefinablePar tmp(name,&((*mComponentTextureParams[comp_nb])(nb_params_phi1+nb_params_phi+1)),-180.*DEG2RAD,180.*DEG2RAD,
		     gpRefParTypeScattDataCorrInt,REFPAR_DERIV_STEP_ABSOLUTE,
		     true,true,true,true,RAD2DEG,360.*DEG2RAD);
    tmp.AssignClock(mClockParams);
    tmp.SetDerivStep(1.*DEG2RAD);
    //tmp.SetGlobalOptimStep(.05);
    this->AddPar(tmp);
  }

	// other phi2 func. params - specific for each model
	
	// phi2 fwhm - TEXTURE_MODEL_HKL_GAUSS, TEXTURE_MODEL_HKL_SIMEK
	if( *mComponentTextureModels_phi2[comp_nb] == TEXTURE_MODEL_HKL_GAUSS ||
			*mComponentTextureModels_phi2[comp_nb] == TEXTURE_MODEL_HKL_SIMEK) 
	{
		string name = "Phi2_fwhm_"+str_comp_nb;
		RefinablePar tmp(name,&((*mComponentTextureParams[comp_nb])(nb_params_phi1+nb_params_phi+2)),1.*DEG2RAD,360.*DEG2RAD,
				gpRefParTypeScattDataCorrInt,REFPAR_DERIV_STEP_ABSOLUTE,
				true,true,true,false,RAD2DEG);
		tmp.AssignClock(mClockParams);
		tmp.SetDerivStep(0.5*DEG2RAD);
		//tmp.SetGlobalOptimStep(.05);
		this->AddPar(tmp);
	}
	
	// phi2 - texture shape par. - n - TEXTURE_MODEL_HKL_SIMEK, TEXTURE_MODEL_BIRKHOLZ
	if( *mComponentTextureModels_phi2[comp_nb] == TEXTURE_MODEL_HKL_SIMEK ||
			*mComponentTextureModels_phi2[comp_nb] == TEXTURE_MODEL_HKL_BIRKHOLZ) 
	{
		string name = "Phi2_n_"+str_comp_nb;
		int tn = (*mComponentTextureModels_phi2[comp_nb]==TEXTURE_MODEL_HKL_SIMEK) ? 3 : 2;
		RefinablePar tmp(name,&((*mComponentTextureParams[comp_nb])(nb_params_phi1+nb_params_phi+tn)),0.1,10.,
				gpRefParTypeScattDataCorrInt,REFPAR_DERIV_STEP_ABSOLUTE,
				true,true,true,false,1.);
		tmp.AssignClock(mClockParams);
		tmp.SetDerivStep(0.1);
		//tmp.SetGlobalOptimStep(.05);
		this->AddPar(tmp);
	}
	
	// let the params. clocks know that something was changed
  mClockParams.Click();
}

void TextureModelHKL::InitParameters()
{
  {
    RefinablePar tmp("GlobalPhi1",&mGlobalPhi1,-180.*DEG2RAD,180.*DEG2RAD,
		     gpRefParTypeScattDataCorrInt,REFPAR_DERIV_STEP_ABSOLUTE,
		     true,true,false,true,RAD2DEG,360.*DEG2RAD);
		tmp.SetValue(0.);
    tmp.AssignClock(mClockParams);
    tmp.SetDerivStep(1.*DEG2RAD);
    this->AddPar(tmp);
  }
  {
    RefinablePar tmp("GlobalPhi",&mGlobalPhi,0.*DEG2RAD,180.*DEG2RAD,
		     gpRefParTypeScattDataCorrInt,REFPAR_DERIV_STEP_ABSOLUTE,
		     true,true,false,true,RAD2DEG,180.*DEG2RAD);
		tmp.SetValue(0.);
    tmp.AssignClock(mClockParams);
    tmp.SetDerivStep(1.*DEG2RAD);
    this->AddPar(tmp);
  }
  {
    RefinablePar tmp("GlobalPhi2",&mGlobalPhi2,-180.*DEG2RAD,180.*DEG2RAD,
		     gpRefParTypeScattDataCorrInt,REFPAR_DERIV_STEP_ABSOLUTE,
		     true,true,false,true,RAD2DEG,360.*DEG2RAD);
		tmp.SetValue(0.);
    tmp.AssignClock(mClockParams);
    tmp.SetDerivStep(1.*DEG2RAD);
    this->AddPar(tmp);
  }
}

std::vector< TextureModelHKL::ComponentPartParams > TextureModelHKL::SplitHKLComponentParams(
																								CrystVector_REAL& params,
																								const int texture_model_phi,
																								const int texture_model_phi2,
																							 	const int texture_model_phi1)
{
	std::vector< ComponentPartParams > result;
	ComponentPartParams cp;
	
	// number of parameters used for the given component (phi1 part)
	int nb_params_phi1 = 0;
		
	switch (texture_model_phi1) {
		case TEXTURE_MODEL_HKL_GAUSS:
			// phi1 func. weight, tilt, fwhm
			nb_params_phi1 += 3;
			break;
		case TEXTURE_MODEL_HKL_SIMEK:
			// phi1 func. weight, tilt, fwhm, shape factor
			nb_params_phi1 += 4;
			break;
		case TEXTURE_MODEL_HKL_BIRKHOLZ:
			// phi1 func. weight, tilt, texture order
			nb_params_phi1 += 3;
			break;
		case TEXTURE_MODEL_HKL_AXIAL:
			break;
		default:
			throw ObjCrystException("Wrong Texture model identificator!");
			break;
	}
	
	// number of parameters used for the given component (phi part)
	int nb_params_phi = 0;
		
	switch (texture_model_phi) {
		case TEXTURE_MODEL_HKL_GAUSS:
			// phi func. weight, tilt, fwhm
			nb_params_phi += 3;
			break;
		case TEXTURE_MODEL_HKL_SIMEK:
			// phi func. weight, tilt, fwhm, shape factor
			nb_params_phi += 4;
			break;
		case TEXTURE_MODEL_HKL_BIRKHOLZ:
			// phi func. weight, tilt, texture order
			nb_params_phi += 3;
			break;
		case TEXTURE_MODEL_HKL_AXIAL:
			break;
		default:
			throw ObjCrystException("Wrong Texture model identificator!");
			break;
	}

	// number of parameters used for the given component (phi2 part)
	int nb_params_phi2 = 0; 

	switch (texture_model_phi2) {
		case TEXTURE_MODEL_HKL_GAUSS:
			// phi2 func. weight, rotation, fwhm
			nb_params_phi2 += 3;
			break;
		case TEXTURE_MODEL_HKL_SIMEK:
			// phi2 func. weight, rotation, fwhm, shape factor
			nb_params_phi2 += 4;
			break;
		case TEXTURE_MODEL_HKL_BIRKHOLZ:
			// phi2 func. weight, rotation, texture order
			nb_params_phi2 += 3;
			break;
		case TEXTURE_MODEL_HKL_AXIAL:
			break;
		default:
			throw ObjCrystException("Wrong Texture model identificator!");
			break;
	}
	
	if((nb_params_phi1+nb_params_phi+nb_params_phi2) != params.numElements())
		throw ObjCrystException("Wrong number of component params!");
		
	// phi1
	cp.params.clear();
	cp.model = texture_model_phi1;
	
	// phi1 weight - for all models (with exception of the axial model)
	if( texture_model_phi1 != TEXTURE_MODEL_HKL_AXIAL )
		cp.params.push_back(&params(0));
	// texture function rotation - also for all type of models (with exception of the axial model)
	if( texture_model_phi1 != TEXTURE_MODEL_HKL_AXIAL )
		cp.params.push_back(&params(1));
	
	// other phi1 func. params - specific for each model
	
	// phi1 fwhm - TEXTURE_MODEL_HKL_GAUSS, TEXTURE_MODEL_HKL_SIMEK
	if( texture_model_phi1 == TEXTURE_MODEL_HKL_GAUSS ||
			texture_model_phi1 == TEXTURE_MODEL_HKL_SIMEK)
		cp.params.push_back(&params(2));
	
	// phi1 - texture shape par. - n - TEXTURE_MODEL_HKL_SIMEK, TEXTURE_MODEL_BIRKHOLZ
	if( texture_model_phi1 == TEXTURE_MODEL_HKL_SIMEK ||
			texture_model_phi1 == TEXTURE_MODEL_HKL_BIRKHOLZ) {
		int tn = (texture_model_phi1==TEXTURE_MODEL_HKL_SIMEK) ? 3 : 2;
		cp.params.push_back(&params(tn)); 
	}
	
	result.push_back(cp);
	
	// phi
	cp.params.clear();
	cp.model = texture_model_phi;
	
	// phi weight - for all models (with exception of the axial model)
	if( texture_model_phi != TEXTURE_MODEL_HKL_AXIAL )
		cp.params.push_back(&params(nb_params_phi1));
	// texture function rotation - also for all type of models (with exception of the axial model)
	if( texture_model_phi != TEXTURE_MODEL_HKL_AXIAL )
		cp.params.push_back(&params(nb_params_phi1+1));

	// other rotation func. params - specific for each model
	
	// phi func. fwhm - TEXTURE_MODEL_HKL_GAUSS, TEXTURE_MODEL_HKL_SIMEK
	if( texture_model_phi == TEXTURE_MODEL_HKL_GAUSS ||
			texture_model_phi == TEXTURE_MODEL_HKL_SIMEK) 
		cp.params.push_back(&params(nb_params_phi1+2));
	
	// phi func. - texture shape par. - n - TEXTURE_MODEL_HKL_SIMEK, TEXTURE_MODEL_BIRKHOLZ
	if( texture_model_phi == TEXTURE_MODEL_HKL_SIMEK ||
			texture_model_phi == TEXTURE_MODEL_HKL_BIRKHOLZ) {
		int tn = (texture_model_phi==TEXTURE_MODEL_HKL_SIMEK) ? 3 : 2;
		cp.params.push_back(&params(nb_params_phi1+tn));
	}

	result.push_back(cp);
	
	// phi2
	cp.params.clear();
	cp.model = texture_model_phi2;
	
	// phi2 weight - for all models (with exception of the axial model)
	if( texture_model_phi2 != TEXTURE_MODEL_HKL_AXIAL )
		cp.params.push_back(&params(nb_params_phi1+nb_params_phi));
	// texture function rotation - also for all type of models (with exception of the axial model)
	if( texture_model_phi2 != TEXTURE_MODEL_HKL_AXIAL )
		cp.params.push_back(&params(nb_params_phi1+nb_params_phi+1));

	// other rotation func. params - specific for each model
	
	// phi2 func. fwhm - TEXTURE_MODEL_HKL_GAUSS, TEXTURE_MODEL_HKL_SIMEK
	if( texture_model_phi2 == TEXTURE_MODEL_HKL_GAUSS ||
			texture_model_phi2 == TEXTURE_MODEL_HKL_SIMEK) 
		cp.params.push_back(&params(nb_params_phi1+nb_params_phi+2));
	
	// phi2 func. - texture shape par. - n - TEXTURE_MODEL_HKL_SIMEK, TEXTURE_MODEL_BIRKHOLZ
	if( texture_model_phi2 == TEXTURE_MODEL_HKL_SIMEK ||
			texture_model_phi2 == TEXTURE_MODEL_HKL_BIRKHOLZ) {
		int tn = (texture_model_phi2==TEXTURE_MODEL_HKL_SIMEK) ? 3 : 2;
		cp.params.push_back(&params(nb_params_phi1+nb_params_phi+tn));
	}
	
	result.push_back(cp);
	
	return result;
}

// ----------------------------------

// Initialisation of some const static data

// Levente Balogh parameters tables for fallting broadening effects in fcc crystals
// ref.: Levente Balogh, Gabor Ribarik, Tamas Ungar, Stacking faults and twins boundaries in fcc crystals
//       determined by x-ray diffraction profile analysis, J.Appl.Phys. (2006) 100, 023512

const FaultsBroadeningEffectFCCBaloghUngar::ParametersTable FaultsBroadeningEffectFCCBaloghUngar::TableIntrinsic =
{ 
  //13, // field Nhkl
  { // field hkl
    { 1, 1, 1 },
    { 2, 0, 0 },
    { 2, 2, 0 },
    { 3, 1, 1 },
    { 2, 2, 2 },
    { 4, 0, 0 },
    { 3, 3, 1 },
    { 4, 2, 0 },
    { 4, 2, 2 },
    { 4, 4, 0 },
    { 5, 3, 1 },
    { 6, 2, 0 },
    { 5, 3, 3 },
  },
  { // field weight
    { 0.25, 0.75, 0.00, 0.00 },
    { 0.00, 1.00, 0.00, 0.00 },
    { 0.50, 0.50, 0.00, 0.00 },
    { 0.50, 0.25, 0.25, 0.00 },
    { 0.25, 0.75, 0.00, 0.00 },
    { 0.00, 1.00, 0.00, 0.00 },
    { 0.00, 0.25, 0.50, 0.25 },
    { 0.50, 0.50, 0.00, 0.00 },
    { 0.25, 0.50, 0.25, 0.00 },
    { 0.50, 0.50, 0.00, 0.00 },
    { 0.50, 0.25, 0.25, 0.00 },
    { 0.00, 0.50, 0.50, 0.00 },
    { 0.00, 0.25, 0.50, 0.25 },
  },
  { // field fwhm  (1/A), alpha(fraction)
    { // field fwhm(1)
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
    },
    { // field fwhm(2)
      { +2.03513000e+01, -9.03079000e+00, +1.78081000e+00, -6.95345000e-02, +8.10219900e-02, +0.00000000e+00 },
      { +3.42893000e+01, -1.58436000e+01, +2.96642000e+00, -1.88284000e-01, +1.40343540e-01, +0.00000000e+00 },
      { +2.59494000e+01, -1.49005000e+01, +3.12098000e+00, -1.80226000e-01, +1.97590850e-01, +0.00000000e+00 },
      { +1.22867000e+01, -5.26286000e+00, +1.02901000e+00, -3.98497000e-02, +4.24542200e-02, +0.00000000e+00 },
      { +1.88532000e+01, -8.81948000e+00, +1.66969000e+00, -1.10012000e-01, +8.09866700e-02, +0.00000000e+00 },
      { +1.84035000e+01, -1.05516000e+01, +2.22317000e+00, -1.23315000e-01, +1.39730090e-01, +0.00000000e+00 },
      { +2.19697000e+03, -5.04755000e+02, +4.29458000e+01, -1.53573700e+00, +2.39385490e-01, +0.00000000e+00 },
      { +1.44747000e+01, -6.79429000e+00, +1.29504000e+00, -8.59916000e-02, +6.27462800e-02, +0.00000000e+00 },
      { +1.49619000e+01, -8.57738000e+00, +1.81581000e+00, -9.90487000e-02, +1.14078340e-01, +0.00000000e+00 },
      { +2.60562000e+01, -1.49592000e+01, +3.11804000e+00, -1.90198000e-01, +1.97585910e-01, +0.00000000e+00 },
      { +6.94159000e+00, -3.01017000e+00, +5.87894000e-01, -2.30794000e-02, +2.38053300e-02, +0.00000000e+00 },
      { +1.16393000e+01, -6.66820000e+00, +1.41932000e+00, -7.64934000e-02, +8.83799900e-02, +0.00000000e+00 },
      { +6.96817000e+00, -2.94947000e+00, +5.70431000e-01, -2.29212000e-02, +2.15428400e-02, +0.00000000e+00 },
    },
    { // field fwhm(3)
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +2.79860000e+01, -1.60769000e+01, +3.33919000e+00, -2.04947000e-01, +2.10669260e-01, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +1.00914000e+01, -4.31240000e+00, +8.35914000e-01, -3.34921000e-02, +3.23806700e-02, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +3.02472000e+01, -1.73292000e+01, +3.61770000e+00, -2.17811000e-01, +2.28175730e-01, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { -1.23150000e+05, +1.17582400e+04, -3.76633600e+02, +4.40254000e+00, +1.51595900e-01, +0.00000000e+00 },
      { +3.02473000e+01, -1.63997000e+01, +3.23750000e+00, -1.97970000e-01, +1.77226580e-01, +0.00000000e+00 },
      { +1.83455000e+01, -9.94389000e+00, +1.95509000e+00, -1.24175000e-01, +1.06838310e-01, +0.00000000e+00 },
    },
    { // field fwhm(4)
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +2.75973000e+01, -1.49823000e+01, +2.94173000e+00, -1.84084000e-01, +1.60738890e-01, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +3.98319000e+01, -2.15694000e+01, +4.27489000e+00, -2.55329000e-01, +2.34980740e-01, +0.00000000e+00 },
    },
  },
  { // field shift (1/A), alpha(fraction)
    { // field shift(1)
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
    },
    { // field shift(2)
      { +5.23762000e+00, -2.46753000e+00, +5.47458000e-01, -9.73368000e-03, +2.38248300e-02, +0.00000000e+00 },
      { +1.05841000e+00, -1.50509000e+00, +3.86143000e-01, -1.18973000e-01, -3.51583300e-02, +0.00000000e+00 },
      { -1.02602000e-01, -5.03705000e-01, +4.17185000e-01, +2.11592000e-02, +5.66692300e-02, +0.00000000e+00 },
      { +3.04070000e+00, -1.36825000e+00, +2.58289000e-01, +5.48767000e-03, +1.16333000e-02, +0.00000000e+00 },
      { -1.26645000e+00, +4.56282000e-01, -1.11276000e-01, -3.02131000e-02, -2.18281000e-02, +0.00000000e+00 },
      { +1.96119000e+01, -9.59126000e+00, +1.74965000e+00, -6.82585000e-02, +4.16005300e-02, +0.00000000e+00 },
      { +5.59563000e+03, -1.14012000e+03, +8.17572000e+01, -2.35381600e+00, +8.81735700e-02, +0.00000000e+00 },
      { +1.33015000e-01, -1.42181000e-01, -1.79960000e-02, -2.51967000e-02, -1.70237300e-02, +0.00000000e+00 },
      { +7.92576000e+00, -3.85236000e+00, +7.84621000e-01, -1.64797000e-02, +3.33150400e-02, +0.00000000e+00 },
      { +1.44314000e+01, -7.85189000e+00, +1.36452000e+00, -2.19488000e-01, -4.90442500e-02, +0.00000000e+00 },
      { +4.54175000e+00, -2.24194000e+00, +4.23464000e-01, -2.05176000e-02, +7.44306000e-03, +0.00000000e+00 },
      { -6.69149000e+00, +3.25809000e+00, -4.39870000e-01, +5.63658000e-02, +2.46061700e-02, +0.00000000e+00 },
      { -3.17642000e+00, +1.54881000e+00, -2.41826000e-01, +2.79779000e-02, +5.34884000e-03, +0.00000000e+00 },
    },
    { // field shift(3)
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +1.63492000e+01, -7.88138000e+00, +1.16935000e+00, -1.95934000e-01, -5.39516800e-02, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { -2.26041000e+00, +1.29238000e+00, -2.30218000e-01, +3.77035000e-02, +7.85792000e-03, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { -2.01389000e+01, +9.28538000e+00, -1.59297000e+00, -3.88196000e-02, -6.14535600e-02, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { -4.73573000e+04, +4.56906000e+03, -1.39083700e+02, +1.28453000e+00, +4.89675900e-02, +0.00000000e+00 },
      { -2.67762000e+00, -2.91447000e-02, +1.46173000e-01, -1.26184000e-01, -4.55926000e-02, +0.00000000e+00 },
      { -3.83127000e+00, +1.34308000e+00, -2.31874000e-01, -3.89303000e-02, -2.88579700e-02, +0.00000000e+00 },
    },
    { // field shift(4)
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { -4.65714000e+00, +1.18039000e+00, -6.46495000e-02, -1.04397000e-01, -4.12468800e-02, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +4.13119000e+00, -1.36594000e+00, +3.44026000e-02, -1.45090000e-01, -6.07575200e-02, +0.00000000e+00 },
    },
  },
};

const FaultsBroadeningEffectFCCBaloghUngar::ParametersTable FaultsBroadeningEffectFCCBaloghUngar::TableTwins =
{ 
  //13, // field Nhkl
  { // field hkl
    { 1, 1, 1 },
    { 2, 0, 0 },
    { 2, 2, 0 },
    { 3, 1, 1 },
    { 2, 2, 2 },
    { 4, 0, 0 },
    { 3, 3, 1 },
    { 4, 2, 0 },
    { 4, 2, 2 },
    { 4, 4, 0 },
    { 5, 3, 1 },
    { 6, 2, 0 },
    { 5, 3, 3 },
  },
  { // field weight
    { 0.25, 0.75, 0.00, 0.00 },
    { 0.00, 1.00, 0.00, 0.00 },
    { 0.50, 0.50, 0.00, 0.00 },
    { 0.50, 0.25, 0.25, 0.00 },
    { 0.25, 0.75, 0.00, 0.00 },
    { 0.00, 1.00, 0.00, 0.00 },
    { 0.00, 0.25, 0.50, 0.25 },
    { 0.50, 0.50, 0.00, 0.00 },
    { 0.25, 0.50, 0.25, 0.00 },
    { 0.50, 0.50, 0.00, 0.00 },
    { 0.50, 0.25, 0.25, 0.00 },
    { 0.00, 0.50, 0.50, 0.00 },
    { 0.00, 0.25, 0.50, 0.25 },
  },
  { // field fwhm (1/A), alpha(fraction)
    { // field fwhm(1)
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
    },
    { // field fwhm(2)
      { +1.62936000e+01, -7.40305000e+00, +1.57710000e+00, -7.43350000e-02, +5.54438300e-02, +0.00000000e+00 },
      { +3.01676000e+01, -1.30966000e+01, +2.78409000e+00, -1.26340000e-01, +9.58772200e-02, +0.00000000e+00 },
      { +3.37887000e+01, -1.60818000e+01, +3.49823000e+00, -1.60788000e-01, +1.35287670e-01, +0.00000000e+00 },
      { +1.21453000e+01, -4.87727000e+00, +1.00149000e+00, -4.63229000e-02, +2.90917100e-02, +0.00000000e+00 },
      { +1.91976000e+01, -8.04327000e+00, +1.70202000e+00, -7.69075000e-02, +5.54320300e-02, +0.00000000e+00 },
      { +2.47609000e+01, -1.15627000e+01, +2.51116000e+00, -1.14687000e-01, +9.56713200e-02, +0.00000000e+00 },
      { +1.64034000e+02, -6.31225000e+01, +9.51175000e+00, -4.55760000e-01, +1.57921920e-01, +0.00000000e+00 },
      { +1.58910000e+01, -6.53804000e+00, +1.37524000e+00, -6.23322000e-02, +4.29952300e-02, +0.00000000e+00 },
      { +2.07080000e+01, -9.56413000e+00, +2.07403000e+00, -9.44960000e-02, +7.81292600e-02, +0.00000000e+00 },
      { +4.99766000e+01, -2.07844000e+01, +4.26735000e+00, -1.93556000e-01, +1.35835570e-01, +0.00000000e+00 },
      { +5.16180000e+00, -2.28916000e+00, +4.84775000e-01, -2.30693000e-02, +1.62691900e-02, +0.00000000e+00 },
      { +1.67229000e+01, -7.64106000e+00, +1.65076000e+00, -7.57535000e-02, +6.05810300e-02, +0.00000000e+00 },
      { +7.07083000e+00, -2.86902000e+00, +5.73086000e-01, -2.78423000e-02, +1.48201800e-02, +0.00000000e+00 },
    },
    { // field fwhm(3)
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +5.06722000e+01, -2.14211000e+01, +4.42210000e+00, -2.01212000e-01, +1.44718070e-01, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +9.95809000e+00, -4.00909000e+00, +8.08215000e-01, -3.80309000e-02, +2.21942200e-02, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +5.76066000e+01, -2.40939000e+01, +4.92337000e+00, -2.24521000e-01, +1.56877200e-01, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +1.46891000e+03, -3.38938000e+02, +2.93236000e+01, -1.02460300e+00, +1.23287610e-01, +0.00000000e+00 },
      { +4.79581000e+01, -1.99117000e+01, +4.02802000e+00, -1.85072000e-01, +1.21724500e-01, +0.00000000e+00 },
      { +2.97033000e+01, -1.21099000e+01, +2.46790000e+00, -1.12283000e-01, +7.33859600e-02, +0.00000000e+00 },
    },
    { // field fwhm(4)
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +4.24099000e+01, -1.77767000e+01, +3.61112000e+00, -1.66398000e-01, +1.10367410e-01, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +6.33539000e+01, -2.64494000e+01, +5.31055000e+00, -2.44743000e-01, +1.61402250e-01, +0.00000000e+00 },
    },
  },
  { // field shift (1/A), alpha(fraction)
    { // field shift(1)
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
    },
    { // field shift(2)
      { +7.23948000e+00, -3.41120000e+00, +6.32089000e-01, -2.96264000e-02, +1.15938000e-03, +0.00000000e+00 },
      { -3.98106000e+00, +1.68449000e+00, -3.78391000e-01, -3.96181000e-03, -4.76613000e-04, +0.00000000e+00 },
      { +1.45469000e+01, -7.63632000e+00, +1.55997000e+00, -8.77133000e-02, +3.36975000e-03, +0.00000000e+00 },
      { -2.91168000e+00, +1.60709000e+00, -2.74118000e-01, +2.87067000e-02, -4.73286000e-04, +0.00000000e+00 },
      { -4.76993000e+00, +2.00709000e+00, -3.43966000e-01, +1.70127000e-03, -2.01500000e-04, +0.00000000e+00 },
      { +1.28405000e+01, -6.19358000e+00, +1.17462000e+00, -6.36353000e-02, +2.53967000e-03, +0.00000000e+00 },
      { +9.03153000e+01, -3.53602000e+01, +5.07303000e+00, -2.63579000e-01, +6.53281000e-03, +0.00000000e+00 },
      { -5.71832000e-01, +3.32453000e-01, -1.24961000e-01, -1.27853000e-03, -4.85763000e-04, +0.00000000e+00 },
      { +1.13017000e+01, -5.27570000e+00, +9.50949000e-01, -4.42757000e-02, +1.43466000e-03, +0.00000000e+00 },
      { -3.36612000e+00, +1.66191000e+00, -4.69021000e-01, -1.35050000e-02, -2.47514000e-04, +0.00000000e+00 },
      { +2.00549000e+00, -9.12150000e-01, +1.48181000e-01, -1.89729000e-03, -6.30655000e-05, +0.00000000e+00 },
      { +2.29341000e+00, -1.01949000e+00, +2.21941000e-01, +2.15516000e-03, +1.88177000e-04, +0.00000000e+00 },
      { +3.38928000e+00, -1.45551000e+00, +2.20175000e-01, -6.42614000e-03, -1.46785000e-05, +0.00000000e+00 },
    },
    { // field shift(3)
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { -9.88129000e+00, +3.82094000e+00, -6.46763000e-01, -2.17040000e-02, +4.88740000e-04, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +1.28750000e+00, -3.84142000e-01, +4.33839000e-02, +6.57769000e-03, -5.50572000e-05, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { -1.99690000e+00, +1.19333000e-01, -6.59046000e-02, -6.36204000e-02, +1.21462000e-03, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { -4.56561000e+02, +7.05828000e+01, -2.48459000e+00, -2.16177000e-02, +2.69606000e-03, +0.00000000e+00 },
      { -4.66220000e+00, +1.74930000e+00, -3.46663000e-01, -2.63196000e-02, +4.78499000e-04, +0.00000000e+00 },
      { +1.82218000e+00, -1.13198000e+00, +1.46681000e-01, -3.75473000e-02, +8.03234000e-04, +0.00000000e+00 },
    },
    { // field shift(4)
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { -2.42764000e+00, +8.16366000e-01, -2.17437000e-01, -2.58941000e-02, +2.45820000e-04, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { -1.50368000e+01, +5.73605000e+00, -8.87333000e-01, -2.09924000e-02, +7.46912000e-04, +0.00000000e+00 },
    },
  },
};

const FaultsBroadeningEffectFCCBaloghUngar::ParametersTable FaultsBroadeningEffectFCCBaloghUngar::TableExtrinsic =
{ 
  //13, // field Nhkl
  { // field hkl
    { 1, 1, 1 },
    { 2, 0, 0 },
    { 2, 2, 0 },
    { 3, 1, 1 },
    { 2, 2, 2 },
    { 4, 0, 0 },
    { 3, 3, 1 },
    { 4, 2, 0 },
    { 4, 2, 2 },
    { 4, 4, 0 },
    { 5, 3, 1 },
    { 6, 2, 0 },
    { 5, 3, 3 },
  },
  { // field weight
    { 0.25, 0.75, 0.00, 0.00 },
    { 0.00, 1.00, 0.00, 0.00 },
    { 0.50, 0.50, 0.00, 0.00 },
    { 0.50, 0.25, 0.25, 0.00 },
    { 0.25, 0.75, 0.00, 0.00 },
    { 0.00, 1.00, 0.00, 0.00 },
    { 0.00, 0.25, 0.50, 0.25 },
    { 0.50, 0.50, 0.00, 0.00 },
    { 0.25, 0.50, 0.25, 0.00 },
    { 0.50, 0.50, 0.00, 0.00 },
    { 0.50, 0.25, 0.25, 0.00 },
    { 0.00, 0.50, 0.50, 0.00 },
    { 0.00, 0.25, 0.50, 0.25 },
  },
  { // field fwhm (1/A), alpha(fraction)
    { // field fwhm(1)
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
    },
    { // field fwhm(2)
      { +7.79192000e+02, -2.83305000e+02, +3.68187000e+01, -1.89927500e+00, +1.09656620e-01, +0.00000000e+00 },
      { +4.04177000e+02, -1.59144000e+02, +2.32303000e+01, -1.26204800e+00, +1.60003600e-01, +0.00000000e+00 },
      { +3.71561000e+01, -1.77621000e+01, +3.12905000e+00, -2.00576000e-01, +1.97550690e-01, +0.00000000e+00 },
      { +4.45287000e+02, -1.30168000e+02, +1.43620000e+01, -6.15185000e-01, +4.96624800e-02, +0.00000000e+00 },
      { +3.80615000e+02, -1.51588000e+02, +2.19968000e+01, -1.20718400e+00, +1.01079360e-01, +0.00000000e+00 },
      { +3.07850000e+01, -1.39635000e+01, +2.44288000e+00, -1.56872000e-01, +1.39866310e-01, +0.00000000e+00 },
      { +2.47005000e+02, -7.99191000e+01, +1.00756000e+01, -4.71504000e-01, +1.64841040e-01, +0.00000000e+00 },
      { +3.74880000e+02, -1.50027000e+02, +2.17249000e+01, -1.19837900e+00, +8.30548900e-02, +0.00000000e+00 },
      { +2.79206000e+01, -1.23053000e+01, +2.13347000e+00, -1.35750000e-01, +1.14306030e-01, +0.00000000e+00 },
      { +8.02726000e+02, -3.00761000e+02, +4.16671000e+01, -2.18970800e+00, +2.31369230e-01, +0.00000000e+00 },
      { +1.39791000e+02, -4.69428000e+01, +5.81120000e+00, -2.87210000e-01, +2.75080400e-02, +0.00000000e+00 },
      { +2.49249000e+01, -1.06888000e+01, +1.82987000e+00, -1.15143000e-01, +8.87134200e-02, +0.00000000e+00 },
      { +2.67235000e+02, -7.86026000e+01, +8.63763000e+00, -3.69775000e-01, +2.59281100e-02, +0.00000000e+00 },
    },
    { // field fwhm(3)
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +2.34393000e+02, -7.47082000e+01, +1.06287000e+01, -5.02728000e-01, +2.15078900e-01, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +6.95033000e+01, -2.01619000e+01, +2.45990000e+00, -1.22602000e-01, +3.31627900e-02, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +5.50542000e+02, -1.91280000e+02, +2.58740000e+01, -1.28334700e+00, +2.44988840e-01, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +2.40773000e+02, -8.12522000e+01, +1.00064000e+01, -4.96114000e-01, +1.70020350e-01, +0.00000000e+00 },
      { +3.97593000e+02, -1.37411000e+02, +1.87178000e+01, -9.21214000e-01, +1.88646210e-01, +0.00000000e+00 },
      { +3.99843000e+02, -1.48328000e+02, +2.06028000e+01, -1.06848000e+00, +1.22879200e-01, +0.00000000e+00 },
    },
    { // field fwhm(4)
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +3.22150000e+02, -1.09315000e+02, +1.33815000e+01, -6.62212000e-01, +2.30661410e-01, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +3.08935000e+02, -9.75487000e+01, +1.35217000e+01, -6.36725000e-01, +2.40395000e-01, +0.00000000e+00 },
    },
  },
  { // field shift (1/A), alpha(fraction)
    { // field shift(1)
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
    },
    { // field shift(2)
      { +6.83445000e+00, -2.59934000e+00, +3.25813000e-01, -4.31359000e-02, -2.18135700e-02, +0.00000000e+00 },
      { -4.08307000e+00, +1.61099000e+00, -1.32192000e-01, +5.86227000e-02, +3.78333200e-02, +0.00000000e+00 },
      { +7.79096000e+00, -3.05383000e+00, +3.66297000e-01, -1.16600000e-01, -5.23653200e-02, +0.00000000e+00 },
      { +2.33063000e+01, -9.71916000e+00, +1.44241000e+00, -1.04433000e-01, -9.37563000e-03, +0.00000000e+00 },
      { -1.05957000e+00, +4.14803000e-01, -6.29484000e-03, +3.14773000e-02, +2.19299800e-02, +0.00000000e+00 },
      { +4.68920000e+00, -1.77238000e+00, +2.06214000e-01, -7.78824000e-02, -3.69780300e-02, +0.00000000e+00 },
      { +2.67453000e+01, -6.57247000e+00, +6.02558000e-01, +3.33334000e-02, +4.39563500e-02, +0.00000000e+00 },
      { -1.05466000e+00, +6.30494000e-01, -1.07401000e-01, +3.79783000e-02, +1.63358600e-02, +0.00000000e+00 },
      { +3.10459000e+00, -1.45868000e+00, +2.28187000e-01, -6.93777000e-02, -3.01409200e-02, +0.00000000e+00 },
      { +1.46202000e+00, -1.33635000e+00, +4.44979000e-01, +2.71811000e-02, +5.47492100e-02, +0.00000000e+00 },
      { -2.24182000e+00, +1.31713000e+00, -2.67873000e-01, +1.49757000e-02, -7.17721000e-03, +0.00000000e+00 },
      { +2.01031000e+00, -4.81804000e-01, -3.58263000e-02, -2.81750000e-02, -2.43111400e-02, +0.00000000e+00 },
      { +4.06056000e+00, -1.96864000e+00, +3.13920000e-01, -2.60073000e-02, -5.51450000e-03, +0.00000000e+00 },
    },
    { // field shift(3)
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { -4.18964000e+00, -1.62680000e+00, +8.04274000e-01, -8.58715000e-03, +5.96749600e-02, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { -2.76521000e+00, +3.82216000e+00, -8.09601000e-01, +4.61418000e-02, -9.92090000e-03, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +5.87287000e+00, -2.91163000e+00, +6.56120000e-01, +2.01689000e-02, +6.36899600e-02, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { -1.11129000e+02, +3.11490000e+01, -3.07102000e+00, +4.61065000e-02, -4.62481200e-02, +0.00000000e+00 },
      { +1.52975000e+01, -6.96675000e+00, +1.24887000e+00, -3.07316000e-02, +5.06514600e-02, +0.00000000e+00 },
      { +1.00740000e+00, -9.32438000e-01, +3.24644000e-01, +4.97430000e-03, +3.03117600e-02, +0.00000000e+00 },
    },
    { // field shift(4)
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { -3.46639000e+01, +9.45195000e+00, -9.18152000e-01, -8.17188000e-02, -6.02532400e-02, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +0.00000000e+00 },
      { +1.11684000e+01, -5.34330000e+00, +1.05710000e+00, -5.67544000e-03, +6.61915900e-02, +0.00000000e+00 },
    },
  },
};

// ----------------------------------

// polyfit auxiliary function
CrystVector_REAL polyfit(const CrystVector_REAL &x, const CrystVector_REAL &y, const int np)
/** Finds coefficients of a polynomial p(x) of the order np fitting the the data y.
 * The polynomial coefficients are stored in the descending powers:
 * p(x) = p(0)*x^np + p(1)*x^(np-1) + ... + p(np) .
 * 
 * (not optimised yet)
 */
{
	CrystVector_REAL result(np+1);
	
	{
  	// build the fitting matrix and the rigth hand side
  	NEWMAT::SymmetricMatrix A(np+1);
  	NEWMAT::ColumnVector b(np+1);
  	
  	A = 0.0;
  	b = 0.0;
  	
  	for(int i=0; i<x.numElements(); i++) {
  		REAL yi = y(i);
  		REAL xi = x(i);
  	  for(int m=0; m<np+1; m++) {
  		  b(m+1) += yi * pow(xi,np-m);
  			for(int n=0; n<=m; n++) {
  				A(m+1,n+1) += pow(xi,np-m)*pow(xi,np-n); //c:TODO: optimise
  			} // n
  		} // m
  	} // i
  	
  	// solve the matrix equation
  	NEWMAT::ColumnVector c = A.i() * b;
  	
  	// store the result
  	for(int m=0; m<np+1; m++) result(m) = c(m+1);
	}
	
	return result; 
}

// Auxiliary routine to calculate simple peak parmeters of numerical data representing a peak-like function. 
PeakParams CalcPeakParams(const CrystVector_REAL &x, const CrystVector_REAL &y, const int nbApproxPoints)
{
  // find maximum value and minimum value
  const REAL ymin = y.min();
  //const int imax = y.imax(); :TODO: check new version of code for CrystVector::imax( , )
  REAL imax;
  {
  	const REAL *p = y.data();
  	imax = 0;
  	REAL vmax = *p;
  	for(int i=0; i<y.numElements(); i++) { if(*p>vmax) { vmax = *p; imax = i; }; p++; }
  }
  REAL xmax = x(imax);
  REAL ymax = y(imax)-ymin;
  
  // fitting data in the maximum (using nbApproxPoints around the maximum)
  // by quadratic polynomial to get better guess of the maximum value and position
  {
  	// prepare data (pick up the selected points) 
  	int ind1 = imax - (nbApproxPoints-1)/2;
  	ind1 = (ind1 < 0) ? 0 : ind1;
  	int ind2 = imax + (nbApproxPoints-1)/2;
  	ind2 = (ind2 > x.numElements()-1) ? x.numElements()-1 : ind2;
		CrystVector_REAL tx(ind2-ind1+1);
		CrystVector_REAL ty(ind2-ind1+1);
		// x-scaling
		REAL xsc = x(ind2)-x(ind1); 
		// copy and scale data
		for(int i=0; i<tx.numElements(); i++) {	tx(i) = (x(ind1+i)-xmax)/xsc; ty(i) = y(ind1+i); }
		///for(int i=0; i<tx.numElements(); i++) { tx(i) = i; ty(i) = 1000.-i*i; }
		// data fitting with cubic polynomial
		CrystVector_REAL p = polyfit(tx,ty,2);
		// calculation of the better guess of xmax and ymax values
		REAL xmax1 = -p(1)/2./p(0);
		ymax = p(0)*pow(xmax1,2) + p(1)*xmax1 + p(2) - ymin;
		xmax = xmax1 * xsc + xmax;
  }
  
  // calcualtion of the integrated intensity and asymmetry
  REAL intensity = 0.0;
  REAL asym = 1.0;
 	{
 		// calc intensity on the lower side of the peak
  	REAL intensity1 = 0.0;
  	REAL intensity2 = 0.0;
  	REAL xi1, xi2, yi1, yi2;
  	
  	xi2 = x(0);
  	yi2 = y(0);
  	for(int i=1; i<x.numElements(); i++) {
  			xi1 = xi2; yi1 = yi2;
  			xi2 = x(i); yi2 = y(i);
  			if(xi2>=xmax) break;
  			intensity1 += 0.5*(xi2-xi1)*(yi2+yi1);
  	}
  	xi2 = xmax; yi2 = ymax;
  	intensity1 += 0.5*(xi2-xi1)*(yi2+yi1);
  	intensity1 -= (xmax-x(0))*ymin;
  	
  	// calc intensity on the upper side of the peak
  	xi1 = x(x.numElements()-1);
  	yi1 = y(x.numElements()-1);
  	for(int i=x.numElements()-2; i>=0; i--) {
  			xi2 = xi1; yi2 = yi1;
  			xi1 = x(i); yi1 = y(i);
  			if(xi1<=xmax) break;
  			intensity2 += 0.5*(xi2-xi1)*(yi2+yi1);
  	}
  	xi1 = xmax; yi1 = ymax;
  	intensity2 += 0.5*(xi2-xi1)*(yi2+yi1);
  	intensity2 -= (x(x.numElements()-1)-xmax)*ymin;

		// intensity
		intensity = intensity1 + intensity2;
  	
  	// calc asymmetry (from integrated intensity)
  	asym = intensity1/intensity2;
 	}
 	
 	// calculate fwhm
 	REAL fwhm = 0.0;
 	{
 		int i;
 		
 		// find lower (left) hwhm value
 		for(i = 0; i<x.numElements(); i++) if ((y(i)-ymin)>=ymax/2.0) break;
 		// linear interpolation to get better guess
 		REAL x1 = x(i-1) + (x(i)-x(i-1)) * (ymax/2.0+ymin-y(i-1))/(y(i)-y(i-1));
 		
 		// find upper (right) hwhm value
 		for(i = x.numElements()-1; i>=0; i--) if ((y(i)-ymin)>=ymax/2.0) break;
 		// linear interpolation to get better guess
 		REAL x2 = x(i) + (x(i+1)-x(i)) * (ymax/2.0+ymin-y(i))/(y(i+1)-y(i));
 		
 		// fwhm
 		fwhm = x2 - x1;
 	}
 	
  PeakParams p;
  p.xmax = xmax;
	p.ymax = ymax;
	p.intensity = intensity;
	p.fwhm = fwhm;
	p.asym = asym;
  
  return p;
}

} // namespace MStruct

//***********EXPLICIT INSTANTIATION*******************//
#include "RefinableObj/ObjRegistry_T.cpp"
template class ObjCryst::ObjRegistry<MStruct::ReflectionProfileComponent>;
template class ObjCryst::ObjRegistry<MStruct::XECsObj>;

//************ SOME FUNCTIONS ***********************//

// double func_ei(const double x) - calculates value of the Exponential Integral:
//
//   Ei(x) = -integral (from t=-x to t=infinity) (exp(-t)/t)
//
// for any real argument x.
//
// For aguments x < 0. an approximative method from the SPECFUN module ei.f
// written by W. J. Code is used.
// For agument x = 0. function returns -numeric_limits<double>::infinity().
// For aguments x > 0. an approximative method from Numerical Recipes
// is utilized.
//
// September 25, 2007, Zdenek Matej
//------------------------------------------------------------------------
double func_ei(const double x)
{
  //----------------------------------------------------------------------
  // Coefficients  for -1.0 <= x < 0.0
  //----------------------------------------------------------------------
  static const double a[7] =
    { 1.1669552669734461083368e+2, 2.1500672908092918123209e+3,
      1.5924175980637303639884e+4, 8.9904972007457256553251e+4,
      1.5026059476436982420737e+5,-1.4815102102575750838086e+5,
      5.0196785185439843791020e+0 };
  static const double b[6] =
    { 4.0205465640027706061433e+1, 7.5043163907103936624165e+2,
      8.1258035174768735759855e+3, 5.2440529172056355429883e+4,
      1.8434070063353677359298e+5, 2.5666493484897117319268e+5 };
  //----------------------------------------------------------------------
  // Coefficients for -4.0 <= x < -1.0
  //----------------------------------------------------------------------
  static const double c[9] =
    { 3.828573121022477169108e-1, 1.107326627786831743809e+1,
      7.246689782858597021199e+1, 1.700632978311516129328e+2,
      1.698106763764238382705e+2, 7.633628843705946890896e+1,
      1.487967702840464066613e+1, 9.999989642347613068437e-1,
      1.737331760720576030932e-8 };
  static const double d[9] = 
    { 8.258160008564488034698e-2, 4.344836335509282083360e+0,
      4.662179610356861756812e+1, 1.775728186717289799677e+2,
      2.953136335677908517423e+2, 2.342573504717625153053e+2,
      9.021658450529372642314e+1, 1.587964570758947927903e+1,
      1.000000000000000000000e+0 };
  //----------------------------------------------------------------------
  // Coefficients for x < -4.0
  //----------------------------------------------------------------------
  static const double e[10] =
    { 1.3276881505637444622987e+2, 3.5846198743996904308695e+4,
      1.7283375773777593926828e+5, 2.6181454937205639647381e+5,
      1.7503273087497081314708e+5, 5.9346841538837119172356e+4,
      1.0816852399095915622498e+4, 1.0611777263550331766871e+3,
      5.2199632588522572481039e+1, 9.9999999999999999087819e-1 };
  static const double f[10] =
    { 3.9147856245556345627078e+4, 2.5989762083608489777411e+5,
      5.5903756210022864003380e+5, 5.4616842050691155735758e+5,
      2.7858134710520842139357e+5, 7.9231787945279043698718e+4,
      1.2842808586627297365998e+4, 1.1635769915320848035459e+3,
      5.4199632588522559414924e+1, 1.0e+0 };
  // a constant :-)
  const double xbig = 701.84e+0;
  // constants from Numerical Recipes
  const int maxit = 100;
  const double euler = 0.577215664901533;
  const double eps = std::numeric_limits<double>::epsilon();
  const double fpmin = std::numeric_limits<double>::min()/eps;

  // ---------------------------------------------------------------------
  // x < 0.
  // ---------------------------------------------------------------------
  if(x < 0.) {
    const double y = fabs(x);
    if(y <= 1.) {
      double sump = a[6]*y + a[0];
      double sumq =      y + b[0];
      for(int i=1; i<6; i++) {
	sump *= y; sump += a[i];
	sumq *= y; sumq += b[i];
      }
      return log(y) - sump/sumq;
    }
    else if(y <= 4.) {
      const double w = 1./y;
      double sump = c[0];
      double sumq = d[0];
      for(int i=1; i<9; i++) {
	sump *= w; sump += c[i];
	sumq *= w; sumq += d[i];
      }
      return -(sump/sumq)*exp(-y);
    }
    else if(x <= xbig) {
      const double w = 1./y;
      double sump = e[0];
      double sumq = f[0];
      for(int i=1; i<10; i++) {
	sump *= w; sump += e[i];
	sumq *= w; sumq += f[i];
      }
      return -w*(1.-w*(sump/sumq))*exp(-y);
    }
    else
      return 0.;
  }
  // ---------------------------------------------------------------------
  // x = 0.
  // ---------------------------------------------------------------------
  else if(x == 0.)
    return -std::numeric_limits<double>::infinity();
  // ---------------------------------------------------------------------
  // x > 0.
  // ---------------------------------------------------------------------
  else {
    int k;
    double fact, prev, sum, term;
    
    if(x < fpmin) return log(x)+euler;
    if(x <= -log(eps)) {
      sum = 0.;
      fact = 1.;
      for(k=1; k<=maxit; k++) {
	fact *= x/k;
	term = fact/k;
	sum += term;
	if (term < eps*sum) break;
      }
      if (k > maxit)
	throw std::runtime_error("Error in soubroutine func_ei: series failed.");
      return sum+log(x)+euler;
    }
    else {
      sum = 0.;
      term = 1.;
      for(k=1; k<=maxit; k++) {
	prev = term;
	term *= k/x;
	if (term < eps) break;
	if (term < prev) sum += term;
	else {
	  sum -= prev;
	  break;
	}
      }
      return exp(x)*(1.+sum)/x;
    }
  }
} // func_ei

// double func_daw(const double x) - calculates value of the Dawson Integral:
//
//   Daw(x) = exp(-x^2) integral (from t=0 to t=x) exp(t^2)
//
// for any real argument x.
//
// Function uses an approximative method from Numerical Recipes.
//
// September 28, 2007, Zdenek Matej
//------------------------------------------------------------------------
double func_daw(const double x)
{
  static const int nmax = 6;
  const double h = 0.4; const double a1 = 2./3.;
  const double a2 = 0.4;  const double a3 = 2./7.;
  int i, n0;
  static bool init = true;
  double d1, d2, e1, e2, sum, x2, xp, xx, ans;
  static double c[nmax];
  
  if (init) {
    init = false;
    for (i=0; i<nmax; i++) c[i] = exp(-(2.*i+1.)*h*(2.*i+1.)*h);
  }
  if (fabs(x) < 0.2) {
    x2 = x*x;
    ans = x*(1.-a1*x2*(1.-a2*x2*(1.-a3*x2)));
  } else {
    xx = fabs(x);
    n0 = 2*int(0.5*xx/h+0.5);
    xp = xx-n0*h;
    e1 = exp(2.*xp*h);
    e2 = e1*e1;
    d1 = n0+1;
    d2 = d1-2.;
    sum = 0.;
    for (i=0; i<nmax; i++, d1 += 2., d2 -= 2., e1 *= e2)
      sum += c[i]*(e1/d1+1./(d2*e1));
    ans = 0.5641895835*((x>=0.) ? 1. : -1.)*exp(-xp*xp)*sum;
  }
  return ans;
} // func_daw

// already defined in the file ReflectionProfile.cpp
#if defined(_MSC_VER) || defined(__BORLANDC__)

double erfc(const double x)// in C99, but not in VC++....
{
   if(x<0.0) return 2.0-erfc(-x);
   if(x<3.8)
   { // Series, Abramowitz & Stegun 7.1.6
      double y=x,y0=x;
      for(int i=1;i<=50;i++)
      {
         y0*=2*x*x/(2*i+1.0);
         y+=y0;
      }
      static const double spi=2/sqrt(M_PI);
      return 1-spi*exp(-x*x)*y;
   }
   double y=1.0,y0=1.0;
   for(int i=1;i<=10;i++)
   {// Asymptotic, Abramowitz & Stegun 7.1.23
      y0*=-(2*i-1)/(2*x*x);
      y+=y0;
   }
   static const double invsqrtpi=1.0/sqrt(M_PI);
   return invsqrtpi*exp(-x*x)/x*y;
}

#endif
