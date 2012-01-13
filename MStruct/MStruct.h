/* 
 * MStruct.h
 * 
 * MStruct++ - Object-Oriented computer program/library for MicroStructure analysis
 * 					   from powder diffraction data.
 * 
 * Copyright (C) 2009-2011  Zdenek Matej
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
 
#ifndef __MSTRUCT__H__
#define __MSTRUCT__H__

extern bool bsavecalc;

#include "RefinableObj/LSQNumObj.h"
#include "RefinableObj/RefinableObj.h"
#include "RefinableObj/GlobalOptimObj.h"
#include "ObjCryst/PowderPattern.h"

#include <vector>
#include <utility>

namespace RotationTB {

/** The auxilliary function calculates the Euler rotation angles
 * phi1, phi and phi2 representing rotation described by
 * suppplied rotation matrix A.
 */
void GetEulerAngles(REAL& phi1, REAL& phi, REAL& phi2,
										const CrystMatrix_REAL& A);

CrystMatrix_REAL GetEulerMatrix(const REAL phi1, const REAL phi, const REAL phi2);

CrystMatrix_REAL MatrixMult(const CrystMatrix_REAL& A, const CrystMatrix_REAL& B);

CrystMatrix_REAL MatrixTranspose(const CrystMatrix_REAL& A);

CrystVector_REAL VectorTransf(const CrystMatrix_REAL& A, const CrystVector_REAL& b);

} // RotationTB

namespace MStruct {

  REAL CalcUnitCellMass(const ObjCryst::Crystal& crystal); // (1e-24 g)

// ReflData
class ReflData {
public:
  long H,K,L;
  REAL x;
  void* data;
public:
  ReflData(long h,long k,long l,REAL val,void* p):H(h),K(k),L(l),
						 x(val),data(p) {};
  virtual ~ReflData() {H=0;K=0;L=0;x=0.;data=0;}
}; // class ReflData

// ReflStore
class ReflStore {
protected:
  std::vector<ReflData> mvData;
public:
  ReflStore() {};
  virtual ~ReflStore() {};
  void add(ReflData& data) {mvData.push_back(data);}
  void add(long h,long k,long l,REAL x,void* data) {
    mvData.push_back(ReflData(h,k,l,x,data));}
  void erase(int i) {mvData.erase(mvData.begin()+i);}
  int find(const ReflData& data)const;
  int find(long h,long k,long l,REAL x)const {
    return find(ReflData(h,k,l,x,0));}
  const ReflData& at(int i)const {
    return mvData.at(i);}
  int size()const {return mvData.size();}
  void clear() {mvData.clear();}
  void print(std::ostream& s)const;
}; // class ReflStore

// LSQNumObj
class LSQNumObj: public ObjCryst::LSQNumObj, public ObjCryst::OptimizationObj {
public:
  LSQNumObj(std::string objName="Unnamed LSQ object"):ObjCryst::LSQNumObj(objName) {};
  virtual ~LSQNumObj() {};
  void Optimize (long &nbSteps, const bool silent=false, const REAL finalcost=0, const REAL maxTime=-1) {};
  void MultiRunOptimize (long &nbCycle, long &nbSteps, const bool silent=false, const REAL finalcost=0, const REAL maxTime=-1) {};
  void Refine (int nbCycle=1,bool useLevenbergMarquardt=false,
                   const bool silent=false, const bool callBeginEndOptimization=true,
                   const float minChi2var=0.01);
  void XMLOutput (ostream &os, int indent=0) const {};
  void XMLInput (istream &is, const ObjCryst::XMLCrystTag &tag) {};
}; // class LSQNumObj

// RandomOptimizationObj
class RandomOptimizationObj: public ObjCryst::OptimizationObj
{
 public:
  RandomOptimizationObj (std::string objName="Unnamed OptimizationObj object");
  void SetLSQNumObj (ObjCryst::LSQNumObj *pLSQNumObj, int nbCycle=1,bool useLevenbergMarquardt=false,
		     const bool silent=false, const bool callBeginEndOptimization=true,
		     const float minChi2var=0.01);
  void SetRefinableObjRandomised (ObjCryst::RefinableObj *pRefinableObj);
  void Optimize (long &nbSteps, const bool silent=false, const REAL finalcost=0, const REAL maxTime=-1);
  void WriteResultsToFile(const string &filename) const;
  void WriteCurrentParamSetToFile(const string &filename, const bool append=true);

  void MultiRunOptimize (long &nbCycle, long &nbSteps, const bool silent=false, const REAL finalcost=0,
			 const REAL maxTime=-1) {};
  void XMLOutput (ostream &os, int indent=0) const {};
  void XMLInput (istream &is, const ObjCryst::XMLCrystTag &tag) {};
 protected:
  ObjCryst::LSQNumObj* mpLSQNumObj;
  int mLSQnbCycle;
  bool mLSQuseLevenbergMarquardt;
  bool mLSQsilent;
  bool mLSQcallBeginEndOptimization;
  REAL mLSQminChi2var;
  ObjCryst::RefinableObj* mpRefinableObjRandomised;
  std::vector< std::pair<unsigned long, REAL> > mSetInfo;

}; // RandomOptimizationObj

/** PowderPatternBackgroundBase : Auxilliary base class for various background models.
 *
 * Purpose of this auxilliary (virtual) class is to reduce duplicit code becouse of many pure
 * virtual methods in the original ObjCryst::PowderPatternComponent.
 * 
 * See a code for deatils of implementation.
 *    i.e. mClcockMaster childs set in the SetParrentPowderPattern(...) method could be
 *    essential in some case when calculation performance should be optimized.
 */
class PowderPatternBackgroundBase: public ObjCryst::PowderPatternComponent
{
public:
	/// Constructor
	PowderPatternBackgroundBase();
	/// Copy constructor
	PowderPatternBackgroundBase(const PowderPatternBackgroundBase &);
	/// Name of this class (MStruct:: PowderPatternBackgroundBase)
	virtual const string& GetClassName()const; 
	/// Set the PowderPattern object which uses this component.
	virtual void SetParentPowderPattern(ObjCryst::PowderPattern &);
	/// Get the calculated powder pattern for this component.
	virtual const CrystVector_REAL& GetPowderPatternCalc()const;
	/// Get the integrated values of the powder pattern.
	virtual pair< const CrystVector_REAL*, const ObjCryst::RefinableObjClock* > GetPowderPatternIntegratedCalc()const;
	/// Get the variance associated to each point of the calculated powder pattern, for this component.
	virtual const CrystVector_REAL& GetPowderPatternCalcVariance()const;
	/// Get the variance associated to each point of the calculated powder pattern, for this component (integrated version).
	virtual pair< const CrystVector_REAL*, const ObjCryst::RefinableObjClock* > GetPowderPatternIntegratedCalcVariance()const;
	/** Does this component have a variance associated with each calculated point ? 
	 *  i.e., do we use maximum likelihood to take into account incomplete models ?
	 */
	virtual bool HasPowderPatternCalcVariance()const;
	/// Use variable slit intensity correction.
	void UseVariableSlitIntensityCorr(const bool b = false);
	/// Get PowderPattern sin(Theta). TODO:: this should be method of the ObjCryst::PowderPattern class
	const CrystVector_REAL& GetPowderPatternSinTheta()const;
	/// Last time when sin(Theta) for each PowderPattern x-point was calcualted. TODO:: this should be method of the ObjCryst::PowderPattern class
	const ObjCryst::RefinableObjClock& GetClockPowderPatternSinTheta()const; 
protected:
	/// Calc the integrated powder pattern.
	virtual void CalcPowderPatternIntegrated()const;
	/// Get the pixel positions separating the integration intervals around reflections.
	virtual const CrystVector_long & GetBraggLimits()const;
	//virtual const CrystVector_long& GetBraggLimits()const;
	/// Set the maximum value for sin(theta)/lambda. (Dummy method - calculation ignores this settings)
	virtual void SetMaxSinThetaOvLambda(const REAL max);
	/** This will be called by the parent PowderPattern object, before calculating the first powder pattern.
	 *  Or maybe it should be called automatically by the object itself...
	 */
	virtual void Prepare();
protected:
	/** Maximum sin(theta)/lambda for all calculations (10 by default).
    *
    * This keeps all data in memory, but only the part which is below
    * the max is calculated.
    *
    * Dummy value - calculation ignores this value.
    *
    */
  REAL mMaxSinThetaOvLambda;
  /** Flag if Variable-Slit intensity correction should be used.
    *
    * The Bragg-Brenaton geometry variable slit intensity correction is applied:
    * 	Calculated intensity is in addition multipled by sin(Theta).  
    */
 	bool  mUseVariableSlitIntensityCorr;
 	/// Vector of calculated Sin(theta) values for each PowderPattern x-point. TODO:: Should be an atribute of the PowderPattern class
 	mutable CrystVector_REAL mPowderPatternSinTheta;
 	/// When sin(Theta) values were calculated last time. TODO:: Should be an atribute of the PowderPattern class
 	mutable ObjCryst::RefinableObjClock mClockPowderPatternSinTheta;
};

/** PowderPatternBackgroundInvX : class to represent an 1/X background.
 *
 * It has no parameters, but it is a scalable component.
 *
 */
class PowderPatternBackgroundInvX: public PowderPatternBackgroundBase
{
public:
	/// Constructor
	PowderPatternBackgroundInvX();
	/// Copy constructor
	PowderPatternBackgroundInvX(const PowderPatternBackgroundInvX &);
	/// Name of this class (MStruct:: PowderPatternBackgroundInvX)
	virtual const string& GetClassName()const;
	/// Symbol for the 1/X type of the background function.
	static const int FUNCTION_OF_X = 0;
	/// Symbol for the 1/sin(Theta) type of the background function.
	static const int FUNCTION_OF_SIN_TH = 1;
	/// Set type of the argument of the 'InvX=1/X' function. (X or sin(Theta))
	void SetXFunctionType(const int type = FUNCTION_OF_X);
protected:
	/// Calc the powder pattern.
	virtual void CalcPowderPattern()const;
	/// Argument 'X' of the 'InvX=1/X' function. (X or sin(Theta))
	int mXFunctionType;
};

/** PowderPatternBackgroundChebyshev : class to represent a Chebyshev polynomial background.
 *
 */
class PowderPatternBackgroundChebyshev: public PowderPatternBackgroundBase
{
public:
	/// Constructor
	PowderPatternBackgroundChebyshev();
	/// Copy constructor
	PowderPatternBackgroundChebyshev(const PowderPatternBackgroundChebyshev &);
	/// Name of this class (MStruct:: PowderPatternBackgroundChebyshev)
	virtual const string& GetClassName()const; 
	/** Set Chebyshev polynomial coefficients.
    *
    * The order of coefficients in the input vector is same as the order of
    * the involved Chebyshev polynomials Tn(x):
    *  b(x) = coef(0)*T0(x) + coef(1)*T1(x) + ... 
    *  
    */
	void SetCoefficients(const CrystVector_REAL &coef);
	/// Get current Chebyshev polynomial coefficients
	const CrystVector_REAL& GetCoefficients()const;
	/// Symbol for the 1/X type of the background function.
	static const int FUNCTION_OF_X = 0;
	/// Symbol for the 1/sin(Theta) type of the background function.
	static const int FUNCTION_OF_SIN_TH = 1;
	/// Set type of the argument of the 'InvX=1/X' function. (X or sin(Theta))
	void SetXFunctionType(const int type = FUNCTION_OF_X);
	/// Get the first derivative values for the LSQ function, for a given parameter.
	virtual const CrystVector_REAL& GetLSQDeriv(const unsigned int, ObjCryst::RefinablePar &);
protected:
	/// Calc values of the Chebyshev polynomials. (Update mChebyshevPolynomials matrix.)
	void CalcChebyshevPolynomials()const;
	/// Calc the powder pattern.
	virtual void CalcPowderPattern()const;
	/// Init parameters and options.
	void Init();
protected:
	/// Chebyshev polynomials coefficients approximating background.
	CrystVector_REAL mChebyshevCoef;
	/// Precalculated values of the used Chebyshev polynomials in the current data points.
	mutable CrystMatrix_REAL mChebyshevPolynomials;
	/// When were values of the used Chebyshev polynomials last computed ?
	mutable ObjCryst::RefinableObjClock mClockChebyshevPolynomialsCalc;
	/// Argument 'X' of the Chebyshev polynomial T_n(X) functions. (X or sin(Theta))
	int mXFunctionType;
};

// PowderPattern
class PowderPattern: public ObjCryst::PowderPattern {
public:
  PowderPattern();
  PowderPattern(const PowderPattern &old):ObjCryst::PowderPattern(old) {};
  const CrystVector_REAL& GetLSQCalc(const unsigned int n) const;
  const CrystVector_REAL& GetLSQObs(const unsigned int n) const;
  const CrystVector_REAL& GetLSQWeight(const unsigned int n) const;
  const CrystVector_REAL& GetLSQDeriv(const unsigned int n,
				      ObjCryst::RefinablePar &par);
  void SetLinearPolarRate(const REAL f);
  void SetIncidenceAngle(const REAL omega);
  REAL GetIncidenceAngle(const REAL omega) const;
  void BeginOptimization(const bool allowApproximations=false, 
			 const bool enableRestraints=false);
	void AddAdditionalLSQObj(ObjCryst::RefinableObj& obj);
	void RemoveAdditionalLSQObj(ObjCryst::RefinableObj& obj);
protected:
	REAL mOmega;
  mutable CrystVector_REAL mLSQCalcNotExcluded;
  mutable CrystVector_REAL mLSQObsNotExcluded;
  mutable CrystVector_REAL mLSQWeightNotExcluded;
  mutable CrystVector_REAL mLSQDerivNotExcluded;
  mutable ObjCryst::RefinableObjClock mClockLSQCalcNotExcluded;
  mutable ObjCryst::RefinableObjClock mClockLSQObsNotExcluded;
  mutable ObjCryst::RefinableObjClock mClockLSQWeightNotExcluded;
  ObjCryst::ObjRegistry< ObjCryst::RefinableObj > mAdditionalLSQObjRegistry;
}; // class PowderPattern

// Texture Correction

// auxiliary computation routines

/* vector cross product (c = a x b) */
void cross(REAL* a, REAL* b, REAL* c);

/* normalized Gauss function */
REAL fgauss(REAL* a, REAL x);

/* calculate rotation matrix for given Euler angles,
   the so-called x-convention is used (see: Euler angles,
   www.mathworld.com */
void rotationMatrix(REAL* A, REAL phi, REAL th, REAL psi);

/* calculate Euler angle from a given rotation matrix,
   the so-called x-convention is used (see: Euler angles,
   www.mathworld.com */
void getAngles(REAL* A, REAL& phi, REAL& th, REAL& psi);

// TextureCalculator
class TextureCalculator {
protected:
  // texture params
  CrystVector_REAL mParams;
  // pointer to a Crystal object
  const ObjCryst::Crystal* mpCrystal;
  // 
  bool mbForceTextureSymmetry;
public:
  TextureCalculator();
  void SetTextureParams(const CrystVector_REAL& params,bool bForceTextureSymmetry=false);
  void SetCrystal(const ObjCryst::Crystal &crystal);
  REAL funcodf(REAL phi,REAL th,REAL psi,REAL* aa=0) const;
  REAL funcodf0(REAL phi,REAL th,REAL psi,REAL* aa=0) const;
  REAL funcodf1(REAL phi,REAL th,REAL psi,REAL* aa=0) const;
  REAL funcodf2(REAL phi,REAL th,REAL psi,REAL* aa=0) const;
  CrystMatrix_REAL GetCrystalOrientations(REAL th0,REAL psi0,CrystVector_REAL vphi0,
					  REAL h,REAL k,REAL l) const;
  REAL CalcCorr(REAL th0,REAL psi0,REAL h,REAL k,REAL l) const;
  CrystMatrix_REAL CalcCorr(const CrystVector_REAL &vth0,
			    const CrystVector_REAL &vpsi0,
			    REAL h,REAL k,REAL l) const;
  void ExportOdf(const string filename) const;
protected:
  // calculate a normalization factor for the ODF function
  REAL CalcOdfNFactor() const;
  // prepre object for calc. of a texture corr.
  void PrepareForCalc(REAL h,REAL k,REAL l) const;
  // auxiliary variables for texture calc.
  mutable int mMultiplicity;
  mutable CrystVector_REAL mvri;
  // tmp memory for odf calc (funcodf)
  mutable CrystVector_REAL mAAMatrix; // aux memory for aa matrix
  mutable CrystVector_REAL mAMatrix; // aux memory for a matrix
  mutable CrystVector_REAL mTMatrix; // aux memory for ta matrix
  //mutable CrystVector_REAL mMAMatrix; // aux memory for ma matrix
  mutable CrystVector_REAL mAVector; // aux memory for a vector
  // normalised vectors of all crystal directions equivalent to set main {HKL} texture axis 
  CrystMatrix_REAL mnmtaHKL;
  // normalised vectors of all crystal directions equivalent to set secondary {HKL} texture axis
  CrystMatrix_REAL mnstaHKL;
  // parameters of the integration grid
  REAL phimin, phimax, phistep;
  int nphi;
  REAL thmin, thmax, thstep;
  int nth;
  REAL psimin, psimax, psistep;
  int npsi;
  // norm. factor of the ODF function (for actual params)
  mutable REAL mOdfNFactor;
}; // class TextureCalculator

class ScatteringCorr: public ObjCryst::ScatteringCorr {
public:
  ScatteringCorr(const ObjCryst::ScatteringData & data);
  
  const CrystVector_REAL& GetCorr(bool needRecalc=true)const;
}; // class ScatteringCorr

// AbsorptionCorr
class AbsorptionCorr:public ScatteringCorr {
public:
  AbsorptionCorr(const ObjCryst::ScatteringData & data);
  virtual ~AbsorptionCorr();
  virtual const string & GetName() const;
  virtual const string & GetClassName() const;
  void SetAbsorptionCorrParams(REAL thickness, REAL depth, REAL absfactor,
			       REAL omega);
protected:
  virtual void CalcCorr() const;
private:
  REAL mOmega;
  REAL mAbsFactor;
  REAL mThickness;
  REAL mDepth;
}; // AbsorptionCorr()

// TextureCorr
class TextureCorr:public ScatteringCorr,
                  public ObjCryst::RefinableObj {
public:
  class TexturePhase {
  public:
    REAL fraction;
    CrystVector_REAL params; // wth,th0,dth,wpsi,psi0,dpsi,hklz(3),hklx(3),tilt(2)
    bool bForceTextureSymmetry;
    TextureCalculator mTextureObj;
    ObjCryst::RefinableObjClock mClockParams;
  public:
    TexturePhase():fraction(0.),params(14),bForceTextureSymmetry(false)
      {params=0.;params(2)=10.*DEG2RAD;params(5)=10.*DEG2RAD;}
    const string & GetName() const
      {const static string name="MStruct::TexturePhase"; return name;}
  };
public:
  TextureCorr(const ObjCryst::ScatteringData & data);
  TextureCorr(const TextureCorr & old);
  virtual ~TextureCorr();
  virtual const string & GetName() const;
  virtual const string & GetClassName() const;
  void AddPhase(const REAL fraction,const REAL thweight,const REAL th0,
		const REAL thwidth,const REAL psiweight,const REAL psi0,
		const REAL psiwidth);
	void AddPhase(const REAL fraction,const CrystVector_REAL params,bool bForceTextureSymmetry=false);
  int GetNbPhase() const;
  virtual void SetCrystal(ObjCryst::Crystal &crystal);
  void SetTextureCorrParams(REAL omega);
protected:
  virtual void CalcCorr() const;
protected:
  REAL mOmega;
  REAL mDivergenceOmega;
  REAL mDivergence2Theta;
  REAL mDivergencePsi_i;
  REAL mDivergencePsi_f;
  REAL mTiltTh;
  REAL mTiltPsi;
  mutable ObjCryst::ObjRegistry<TexturePhase> mPhaseRegistry;
  mutable ObjCryst::RefinableObjClock mClockTiltParams;
private:
  void InitParameters();
}; // TextureCorr()

// HKLIntensityCorr
class HKLIntensityCorr:public ScatteringCorr,
                       public ObjCryst::RefinableObj {
public:
  class IntensityCorrData {
  public:
    int i;
    REAL val; };
public:
  HKLIntensityCorr(const ObjCryst::ScatteringData & data);
  HKLIntensityCorr(const HKLIntensityCorr & old);
  virtual ~HKLIntensityCorr();
  virtual const string & GetName() const;
  virtual const string & GetClassName() const;
  void SetHKLIntensityCorr(int h,int k,int l,REAL val,bool fixed=false);
  void Reset();
  virtual void BeginOptimization(const bool allowApproximations=false,
				 const bool enableRestraints=false);
  const ReflStore& GetReflStore()const;
protected:
  virtual void CalcCorr() const;
protected:
  mutable ObjCryst::RefinableObjClock mClockIndexesCalc;
  ReflStore mReflStore;
}; // HKLIntensityCorr()

// PowderPatternDiffraction
class PowderPatternDiffraction: public ObjCryst::PowderPatternDiffraction {
protected:
  REAL mOmega;
  AbsorptionCorr mCorrAbsorption;
  TextureCorr mCorrTexture;
  HKLIntensityCorr mCorrHKLIntensity;
  mutable ObjCryst::RefinableObjClock mClock2IntensityCorr;
  mutable CrystVector_REAL mIntensityCorrTheta;
public:
  PowderPatternDiffraction();
  PowderPatternDiffraction(const PowderPatternDiffraction &old);
  virtual PowderPatternDiffraction* CreateCopy()const;
  virtual const string& GetClassName()const;
  virtual void SetCrystal(ObjCryst::Crystal &crystal);
  void SetAbsorptionCorrParams(REAL thickness, REAL depth, REAL absfactor,REAL omega);
  void SetTextureCorrParams(REAL omega);
  void AddTextureCorrPhase(REAL fraction,const CrystVector_REAL& params, bool bForceTextureSymmetry=false);
  void SetHKLIntensityCorrParams(int h, int k, int l, REAL val, bool fixed=false);
  void GenerateHKLIntensityCorrForAllReflections(const REAL relIntensity=-1.);
  void PrintHKLIntensityCorr(ostream& s)const;
  void WriteHKLIntensityCorrToFile(const char* filename)const;
  void ReadHKLIntensityCorrFromFile(const char* filename);
  //CrystVector_REAL GetIncidenceIngle(const CrystVector_REAL& theta);
  virtual const CrystVector_REAL& GetLSQDeriv(const unsigned int n, ObjCryst::RefinablePar &par);
  void PrintHKLInfo (ostream &s);
  void PrintHKLInfo2 (ostream &s, const REAL accur=-1.) const;
protected:
  void CalcIntensityCorr () const;
}; // class PowderPatternDiffraction

// pre-definition of class MStruct::ReflectionProfile
class ReflectionProfile;

class ReflectionProfileComponent: virtual public ObjCryst::RefinableObj {
public:
  ReflectionProfileComponent();
  virtual void SetParentReflectionProfile(const ReflectionProfile &);
  virtual const ReflectionProfile& GetParentReflectionProfile()const;
  
  virtual CrystVector_REAL GetProfile(const CrystVector_REAL &x,
				      const REAL xcenter,
				      const REAL h, const REAL k,
				      const REAL l)=0;
  virtual REAL GetApproxFWHM(const REAL xcenter,
			     const REAL h, const REAL k, const REAL l)const=0;
  virtual bool IsRealSpaceType()const;
  virtual bool IsAnisotropic()const;
  virtual REAL GetPositionCorr(const REAL xcenter,
			       const REAL h, const REAL k, const REAL l)const;
protected:
  const ReflectionProfile *mpParentReflectionProfile;
}; // class ReflectionProfileComponent

class SizeBroadeningEffect: public ReflectionProfileComponent {
private:
  REAL mM;
  REAL mSigma;
public:
  SizeBroadeningEffect();
  CrystVector_REAL GetProfile(const CrystVector_REAL &x,
			      const REAL xcenter,
			      const REAL h, const REAL k, const REAL l);
  REAL GetApproxFWHM(const REAL xcenter,
		     const REAL h, const REAL k, const REAL l)const;
  bool IsRealSpaceType()const;
  void SetProfilePar(const REAL m, const REAL sigma);
private:
  void InitParameters();
}; // class SizeBroadeningEffect

#ifdef __DEPRECATED_DoubleComponentBroadeningEffect__

class DoubleComponentBroadeningEffect: public ReflectionProfileComponent {
private:
	REAL mWeight;
public:
	DoubleComponentBroadeningEffect();
	void SetComponents(const ReflectionProfileComponent* effect1, const ReflectionProfileComponent* effect2);
  void SetProfilePar(const REAL weight);
  REAL GetApproxFWHM(const REAL xcenter,
		     const REAL h, const REAL k, const REAL l)const;
  bool IsRealSpaceType()const;
  bool IsAnisotropic()const;
protected:
	ReflectionProfileComponent* mEffect1;
	ReflectionProfileComponent* mEffect2;
}; // DoubleComponentSizeBroadeningEffect

#endif // __DEPRECATED_DoubleComponentBroadeningEffect__

class SizeDistribBroadeningEffect: public ReflectionProfileComponent {
private:
	CrystVector_REAL mD1;
	CrystVector_REAL mD2;
	CrystVector_REAL mDistrib;
	REAL mLSQAlpha;
public:
  SizeDistribBroadeningEffect();
  const string& GetClassName()const;
  CrystVector_REAL GetProfile(const CrystVector_REAL &x,
			      const REAL xcenter,
			      const REAL h, const REAL k, const REAL l);
  REAL GetApproxFWHM(const REAL xcenter,
		     const REAL h, const REAL k, const REAL l)const;
  bool IsRealSpaceType()const;
  void SetDistribution(const CrystVector_REAL d1, const CrystVector_REAL d2,
		       const CrystVector_REAL distrib, const CrystVector_int fixed);
  /// Get vector of values representing the size distribution (unnormalised histogram).
  const CrystVector_REAL & GetDistribution () const;
  /// Get vectors defining lower (D1) and upper (D2) bounds of the distribution bins. 
  void GetDistributionBins (CrystVector_REAL & D1, CrystVector_REAL & D2) const; 
	void ReadDistributionFromFile(const char* filename);
	void WriteDistributionToFile(const char* filename) const;
	void SetLSQConstraintScale(const REAL scale);
	//const CrystVector_REAL& GetLSQCalc(const unsigned int n) const;
	//const CrystVector_REAL& GetLSQObs(const unsigned int n) const;
	//const CrystVector_REAL& GetLSQWeight(const unsigned int n) const;
	//const CrystVector_REAL& GetLSQDeriv(const unsigned int n, ObjCryst::RefinablePar &par);
  unsigned int GetNbLSQConstraints() const;
  void GetLSQConstraint(const unsigned int n,
			std::vector< const ObjCryst::RefinablePar* > &parList, CrystVector_REAL &coef) const;
  virtual unsigned int GetNbLSQRegularizationOperator(const unsigned int LSQfunc) const;
  virtual const ObjCryst::LSQRegularizationOperator & GetLSQRegularizationOperator(const unsigned int nOp,
										   const unsigned int LSQfunc) const;

    void GenerateRandomDistrib();
  /// TODO:: comment here
  virtual void GlobalOptRandomMove(const REAL mutationAmplitude,
				   const ObjCryst::RefParType *type=ObjCryst::gpRefParTypeObjCryst);

  /**   \brief No size regularization used in the LSQ refinement.
   *
   * Type option for a distribution regularization method used during the LSQ refinement.
   * No regularization is done for this option.
   *
   */
  static const int LSQRegOpt_None =               0;

  /**   \brief Simple method of size distribution regularization by minimizing roughly its derivative
   *           is used for the LSQs.
   *
   * Type option for a distribution regularization method used during the LSQ refinement.
   * 
   * The vector repsenting a rough numerical approximation of the distribution derivative in the middle
   * of (D1,D2) intervals were distribution is defined is returned byt the object's GetLSQCalc(...) method.
   * Appropriate GetLSQObs(...) values are set zero hence the derivatives are minimized in advance.
   * The calculated ChiSq value of the derivatives is weighted by mLSQAlpha.
   *
   * Note: Size distributions with "holes" in their definition (if defined by intervasls) are not correctly
   *       handled during the spline interpolation. Define distribution with fixed zero values in the "holes"
   *       to obey this problem.
   *
   */
  static const int LSQRegOpt_DistribDeriv =       1;

  /**   \brief Simple method of size distribution regularization by minimizing roughly the derivative
   *           of the volume weighted size distribution.
   *
   * Type option for a distribution regularization method used during the LSQ refinement.
   * 
   * The vector repsenting a rough numerical approximation of the distribution derivative in the middle
   * of (D1,D2) intervals were distribution is defined is returned byt the object's GetLSQCalc(...) method.
   * Appropriate GetLSQObs(...) values are set zero hence the derivatives are minimized in advance.
   * The calculated ChiSq value of the derivatives is weighted by mLSQAlpha.
   *
   * Note: Size distributions with "holes" in their definition (if defined by intervasls) are not correctly
   *       handled during the spline interpolation. Define distribution with fixed zero values in the "holes"
   *       to obey this problem.
   *
   */
  static const int LSQRegOpt_VolumeDistribDeriv = 2;

  /**   \brief Simple method of size distribution regularization by minimizing roughly the derivative
   *           of the arithmetic and the volume weighted size distribution with the same weight.
   *
   * Type option for a distribution regularization method used during the LSQ refinement.
   * 
   * The vector repsenting a rough numerical approximation of the distribution derivative in the middle
   * of (D1,D2) intervals were distribution is defined is returned byt the object's GetLSQCalc(...) method.
   * Appropriate GetLSQObs(...) values are set zero hence the derivatives are minimized in advance.
   * The calculated ChiSq value of the derivatives is weighted by mLSQAlpha. 
   *
   * The first part of the LSQFunction vector represents the derivative of the arithmetic distribution,
   * the second part is linked with the volume weighted distribution.
   *
   */
  static const int LSQRegOpt_BothDistribDeriv = 3;

  /**   \brief Integral of unsigned size distribution curvature is minimized to regularize the distrubution.
   *
   * Type option for a distribution regularization method used during the LSQ refinement.
   * 
   * The distribution is interpolated by a cubic spline and an integral of the unsigned curvature
   * is calcualted over the range of crystallites sizes (0,Dmax). This integral weighted by mLSQAlpha
   * is then minimized also by LSQs. The object's GetLSQCalc(...) method returns a vector containing only
   * one element: Sqrt( IntegralOfCurvature ).
   *
   * In addition it is assumend for spline interpolation that size distribution f(D) is zero and has zero
   * derivatives at points D=0 and D=Dmax.
   *
   * Note: Size distributions with "holes" in their definition (if defined by intervasls) are not correctly
   *       handled during the spline interpolation. Define distribution with fixed zero values in the "holes"
   *       to obey this problem.
   *
   */
  static const int LSQRegOpt_CurvIntegral =       4;

  /**   \brief  Set a weight and a method type used for the distribution regularization in the LSQ refinement.
   *
   * See a particular option documentation for details of a method used.
   *
   */
  void AddLSQRegularizationMethod(const int option=LSQRegOpt_None, const REAL weight=0.);

  /**   \brief  Build a new distribution.
   *
   * \par NbIntervals intervals of crystallites sizes from \par \Dmin (A) to \par Dmax (A) are generated
   * with "liner", "log" or "sqrt" spacing.
   *
   * 1 - cos( 2*Pi*(D-Dmin)/(Dmax-Dmin) ) like distribution is assigned to the centers of generated
   * intervals. The distribution values in all the intervals are set to be refined.
   *
   */
  void BuildDistribution(const REAL Dmin=10., const REAL Dmax=1.e3, const int NbIntervals=10,
			 const string spacing=string("linear"));

  /// Get integral of the area below the arithmetic distribution.
  REAL GetDistribIntegral() const;
  /// Get integral of the area below the volume weighted distribution.
  REAL GetVolumeDistribIntegral() const;

  /**   \brief  This should be called by any optimization class at the begining of an optimization.
   *
   * It prints information about regularization.
   *
   */
  virtual void BeginOptimization (const bool allowApproximations, const bool enableRestraints);
  /**   \brief  This should be called by any optimization class at the end of an optimization.
   *
   * It prints information about regularization.
   *
   */
  virtual void 	EndOptimization ();
  
  void PrintRegularizationStatistics () const;
  void PrintConstraintsStatistics () const;

protected:

  /// Calculate integrated areas below the arithmetic and the volume weighted distribution
  void CalcDistIntegral () const;
  /// Rebuild the list of Regularization Operators used in the LSQ refinement.
  void RebuildLSQRegOpList ();

  //mutable CrystVector_REAL mLSQCalc;
  //CrystVector_REAL mLSQObs;
  //CrystVector_REAL mLSQWeight;
  //CrystVector_REAL mLSQDeriv;
  REAL mLSQConstraintScale;
  /// List of types of distribution regularization methods used in the LSQ refinement
  vector<int> mLSQRegTypeList;
  /// List of weights of distribution regularization operators used in the LSQ refinement
  vector<REAL> mLSQRegWeightList;
  /** \brief List of LSQ Regularizator Operators - requiring smoothness of the distribution
   * 
   * The LSQRegularizationOperators are normalised "on the fly" (dynamically) in the const method
   * GetLSQRegularizationOperator(...), hence the operators are mutable.
   *
   */
  mutable std::vector< ObjCryst::LSQRegularizationOperator > mLSQRegOpList;
  /** \brief Norms of LSQ Regularizator Operators as their were used last time
   * 
   * The LSQRegularizationOperators are normalised "on the fly" (dynamically) in the const method
   * GetLSQRegularizationOperator(...), the calculated norms as used in that method last time
   * are stored here.
   *
   */
  mutable std::vector< REAL > mLSQRegOpNorm;
  /// Maximum crystallites size for the given distribution
  REAL mDmax;
  
  /// Integrated area of the arithmetic distribution  
  mutable REAL mDistIntegral;
  /// Integrated area of the volume weighted distribution  
  mutable REAL mVolumeDistIntegral;
  /// Last time the integrals of the distribution were calculated
  mutable ObjCryst::RefinableObjClock mClockDistIntegralCalc;

private:
  /// TODO::
  ObjCryst::LSQRegularizationOperator CreateLSQRegOpVolumeDistribDeriv (const REAL weight=0.) const;
  /// TODO::
  ObjCryst::LSQRegularizationOperator CreateLSQRegOpDistribDeriv (const REAL weight=0.) const;
  /**  \brief  Calculate approximative size distribution integrated curvature.
   *
   * Size distribution is interpolated by a cubic spline. Than its unsigned curvature is
   * integrated over the whole range (0,Dmax). For spline interpolation The size distribution
   * is complemented by zero values at points D=0 and D=Dmax and also zero derivatives at these
   * point are assumed.
   *
   * Note: Size distributions with "holes" in their definition (if defined by intervasls) are not correctly
   *       handled during the spline interpolation. Define distribution with fixed zero values in the "holes"
   *       to obey this problem.
   *
   */
  REAL CalcCurvInt() const;

  /// Step for numerical integration of the unsigned curvature of the distribution
  REAL mCurvIntStep;

  /// Number of times the Begin and End Optimization methods were called
  int mBeginEndOptimizationCalled;

}; // class SizeDistribBroadeningEffect

/** RandomSizeDistribBroadeningEffect : This provides reflection profile
 * calculation from spherical crystallites with randomly distrubuted diameter D.
 *
 * Distribution of crystallites size diameter D is defined in the N bins between
 * D=0 and D=Dmax. The bin width is dD=Dmax/N. The distribution is represented by
 * a probability P(i) of finding P crystallites with size in the i-th bin (D lies in
 * the interval i*dD + (0,dD) ). This probability P(i) is distributed uniformly in
 * each bin independently of others bins and it is normalised that
 * sum(i=0..N-1, P(i)*dD)==1 holds.
 *
 */
class RandomSizeDistribBroadeningEffect: public ReflectionProfileComponent {
 public:
  /// Constructor
  RandomSizeDistribBroadeningEffect();
  /// Name for this class ("MStruct::RandomSizeDistribBroadeningEffect")
  virtual const string & GetClassName() const;
  /**   \brief Set new distribution
   *  \param Nbins   : number of bins in which the interval of crystallite size D is divided.
   *  \param Dmax    : crystallite sizes are distributed in the interval (0,Dmax) (in nm)
   *  \param distrib : vector of length Nbins containing initial distribution (optional)
   */
  void SetDistribution(const int Nbins=20, const REAL Dmax=1.e3,
		       const CrystVector_REAL &distrib=CrystVector_REAL(0));
  /**   \brief Generate random distribution
   *
   *  Generate random distribution: a random value P(i) representiing the crystallite size
   *  distrution is choosen from an uniform distribution for each bin independently
   *  on other bins. Aftervards the values P(i) in all bins are normalised so 
   *  sum(i=0..Nbins-1, P(i)*BinWidth)==1 holds.
   */
  void GenerateRandomDistrib();

  /// TODO:: comment here
  virtual void GlobalOptRandomMove(const REAL mutationAmplitude,
				   const ObjCryst::RefParType *type=ObjCryst::gpRefParTypeObjCryst);

  virtual CrystVector_REAL GetProfile(const CrystVector_REAL &x,
				      const REAL xcenter,
				      const REAL h, const REAL k, const REAL l);
  virtual REAL GetApproxFWHM(const REAL xcenter,const REAL h, const REAL k, const REAL l)const;
  virtual bool IsRealSpaceType()const;

 protected:
  /// Initialize RefinableObj parameters
  void Init();

  // data members section

  /// Number of bins between 0 and Dmax representing the crystallites size distribution
  int mNbins;
  /// Width of the bin representing the crystallites size distribution, Width = Dmax/Nbins
  REAL mBinWidth;
  /// Vector representing crystallites size distribution P(i)
  CrystVector_REAL mDistrib;
  /// Clock when reflection profile was calculated last time
  ObjCryst::RefinableObjClock mClockReflProfCalc;

}; // RandomSizeDistribBroadeningEffect

class DislocationBroadeningEffectSvB: public ReflectionProfileComponent {
private:
	REAL mReOrMWilk;
	REAL mRou;
	REAL mCg0;
	mutable REAL mQ1;
	mutable REAL mQ2;
	bool mUseMWilk;
	int mFormula;
	int mArgument;
	REAL mKaganerEta0;
	REAL mKaganerEta1;

	bool mIsIntialised;
public:
	DislocationBroadeningEffectSvB();
	void SetParentReflectionProfile(const ReflectionProfile &);
	CrystVector_REAL GetProfile(const CrystVector_REAL &x,
			      const REAL xcenter,
			      const REAL h, const REAL k, const REAL l);
	REAL GetApproxFWHM(const REAL xcenter,
			     const REAL h, const REAL k, const REAL l)const;
	void PrepareCalcAuxParams(const REAL xcenter, const REAL h, const REAL k, const REAL l,
														double &s0, double &Chkl, double &thkl)const;
	bool IsRealSpaceType()const;
	bool IsAnisotropic()const;
	void SetProfilePar(const REAL reOrMWilk, const REAL rou, const REAL cg0 = 1.,
					   const REAL q1 = 0., const REAL q2 = 0.);
	
	void SetUseMWilk(const bool useMWilk);
	void SetFormula(const int formula, const int arg);

	/*	unsigned int GetNbLSQConstraints() const;
	void GetLSQConstraint(const unsigned int n,
			      std::vector< const ObjCryst::RefinablePar* > &parList,
			      CrystVector_REAL &coef) const;*/
private:
	mutable int mCellType;
	mutable REAL mACr; // a/c ratio
	// length of the dislocations Burgers vector for cubic cells,
	// |b<11-20>|=a for hcp, 1. in other cases 
	mutable REAL mb;
	// clock for auxiliary params (cell type, a/c ratio, Burgers vector length) 
	mutable ObjCryst::RefinableObjClock mClockAuxParams;
private:
	void InitParameters(const bool reinitialize = false);
	void SetAuxParameters()const;
}; // class DislocationBroadeningEffectSvB

/** SizeDistribPowderPatternDiffraction: is a PowderPatternDiffraction class optimised
 * for a refinement of a histogram-like crystallites size representation. It should
 * provide effectively derivatives of the powder diffraction pattern with respect to
 * size distribution parameters (histogram values) and increase performace of the
 * LSQ-refinement especially if solely the size distribution is optimised.
 *
 */
class SizeDistribPowderPatternDiffraction: virtual public MStruct::PowderPatternDiffraction {
 public:
  /// Constructor.
  SizeDistribPowderPatternDiffraction ();
  /// Virtual copy constructor. 
  SizeDistribPowderPatternDiffraction (const SizeDistribPowderPatternDiffraction &);
  /// Creates copy of the object.
  virtual SizeDistribPowderPatternDiffraction * CreateCopy () const;
  /// Name for this class ("MStruct::PowderPattern", ...). 
  virtual const string & GetClassName () const;

  /** \brief Get the first derivative values for the LSQ function, for a given parameter.
   *
   *  The SizeDistrib PowderPatternDiffraction class is optimised to provide effectively
   *  LSQ derivatives with respect to SizeDistrib parameters (histogram values). 
   *  GetLSQDeriv should return these derivatives effectively without unnecessary
   *  repeated numerical recalculations.
   */
  virtual const CrystVector_REAL & GetLSQDeriv (const unsigned int, ObjCryst::RefinablePar &);

 protected:
  /// Prepare everything for an optimization/calculation. 
  virtual void Prepare ();

 protected:
  /// Pointer to the SizeDistrib ReflectionProfile Component Object whose calculation should be optimised.
  SizeDistribBroadeningEffect * mpSizeDistribReflProf;
  /// Clock when all possibly related refinable parameters except the SizeDistrib parameters were mutated.
  ObjCryst::RefinableObjClock mOtherParamsClock;
  /// Clocks when LSQDerivatives with respect to SizeDistrib parameters were calculated.
  std::vector< ObjCryst::RefinableObjClock > mvClockLSQDerivCalculated;
  
  /// Calculated LSQDerivatives with respect to the SizeDistrib parameters
  std::vector< CrystVector_REAL > mvLSQDerivSizeDistrib;
 private:
  /// Values of a Size Fourier coefficients for a zero length A(0) for all SizeDistrib bins 
  CrystVector_REAL mSizeDistribA0;

}; // SizeDistribPowderPatternDiffraction

class FaultsBroadeningEffectFCC: public ReflectionProfileComponent {
protected:
	REAL mAlpha;
	REAL mBeta;
public:
	FaultsBroadeningEffectFCC();
	virtual void SetParentReflectionProfile(const ReflectionProfile &);
	virtual CrystVector_REAL GetProfile(const CrystVector_REAL &x,
			      								  				const REAL xcenter,
			      													const REAL h=0, const REAL k=0, const REAL l=0)=0;
	virtual REAL GetApproxFWHM(const REAL xcenter,
	 									 				 const REAL h=0, const REAL k=0, const REAL l=0) const;
	virtual bool IsRealSpaceType()const;
	virtual bool IsAnisotropic()const;
	virtual void SetProfilePar(const REAL alpha, const REAL beta);
protected:
  mutable const ObjCryst::UnitCell *mpUnitCell;
	mutable CrystVector_int mAbsL0;
	mutable CrystVector_int mCount;
	mutable CrystVector_int mSign;
	mutable ObjCryst::RefinableObjClock mClockAuxParams;
protected:
	void InitParameters();
	void SetAuxParameters()const;
	virtual void GetAllGroupsProfiles(CrystMatrix_REAL &profiles,
																		CrystVector_REAL &shifts,
																		const CrystVector_REAL &x,
			  				  									const REAL xcenter,
			    													const REAL h=0, const REAL k=0, const REAL l=0)=0;
	void PrepareCalcAuxParams(const REAL xcenter, const REAL h, const REAL k, const REAL l,
														double &s0)const;
}; // class FaultsBroadeningEffectFCC

class FaultsBroadeningEffectVelteropFCC: public ReflectionProfileComponent {
private:
	REAL mAlpha;
	REAL mBeta;
public:
	FaultsBroadeningEffectVelteropFCC();
	void SetParentReflectionProfile(const ReflectionProfile &);
	CrystVector_REAL GetProfile(const CrystVector_REAL &x,
			      								  const REAL xcenter,
			      									const REAL h, const REAL k, const REAL l);
	REAL GetApproxFWHM(const REAL xcenter,
			     										const REAL h, const REAL k, const REAL l)const;
	bool IsRealSpaceType()const;
	bool IsAnisotropic()const;
	void SetProfilePar(const REAL alpha, const REAL beta);
private:
  mutable const ObjCryst::UnitCell *mpUnitCell;
	mutable CrystVector_int mAbsL0;
	mutable CrystVector_int mCount;
	mutable CrystVector_int mSign;
	mutable ObjCryst::RefinableObjClock mClockAuxParams;
private:
	void InitParameters();
	void SetAuxParameters()const;
	void GetAllGroupsProfiles(CrystMatrix_REAL &profiles,
														CrystVector_REAL &shifts,
														const CrystVector_REAL &x,
			  				  					const REAL xcenter,
			    									const REAL h, const REAL k, const REAL l);
	void PrepareCalcAuxParams(const REAL xcenter, const REAL h, const REAL k, const REAL l,
														double &s0)const;
}; // class FaultsBroadeningEffectVelteropFCC

class FaultsBroadeningEffectFCCBaloghUngar: public FaultsBroadeningEffectFCC {
protected:
	int mFaultsType;
public:
	FaultsBroadeningEffectFCCBaloghUngar();
	virtual CrystVector_REAL GetProfile(const CrystVector_REAL &x,
			    				  								  const REAL xcenter,
			      													const REAL h=0, const REAL k=0, const REAL l=0);
	void SetProfilePar(const int type, const REAL alpha);
	// alpha(fraction), fwhm and shift in (1/A) units
	void GetSubComponentsPar(const REAL alpha, const int type,
													 const REAL h, const REAL k, const REAL l,
													 REAL weight[4], REAL fwhm[4], REAL shift[4])const;
	// Symbols for stacking faults type
	static const int INTRINSIC = 0;												 
	static const int TWINS = 1;
	static const int EXTRINSIC = 2;
protected:
  // shifts (2theta(rad))
	virtual void GetAllGroupsProfiles(CrystMatrix_REAL &profiles,
																		CrystVector_REAL &shifts,
																		const CrystVector_REAL &x,
			  				  									const REAL xcenter,
			    													const REAL h=0, const REAL k=0, const REAL l=0);
	// Data type for a table of Levente Balogh FCC faulting parameters
	struct ParametersTable {
		static const int Nhkl = 13; // number of hkl reflections stored in the table
		int hkl[Nhkl][3]; // hkl indexes of stored reflections
		REAL weight[Nhkl][4]; // subcomponents weights
		REAL fwhm[4][Nhkl][6]; // fwhm polynomial coefficients for all subcomponents
		REAL shift[4][Nhkl][6]; // shift polynomial coefficients for all subcomponents
	};
	static const ParametersTable TableIntrinsic; // Parameters table for intrinsic stacking faults
	static const ParametersTable TableTwins; // Parameters table for twins
	static const ParametersTable TableExtrinsic; // Parameters table for extrinsic stacking faults
}; // class FaultsBroadeningEffectFCCBaloghUngar

/** FaultsBroadeningEffectWC11m23 : Describes line broadeng from [11-23] displacement faults in WC. 
 *
 * It simulates line broadenig profiles from a WC sample containing [11-23] displacement faults as described
 * by S. Hagege et al.(1980).
 *
 * A simple model assuming stacking of the aAbB tetra-layers is used, hence only a B|a fault is included.
 * Fourier coefficients derived by Zdenek (Jan, 2012) by applying directly the simple approach of Warren.
 *
 * ref: [1] S. Hagege, J. Vicens, G. Nouet, P. Delavignette, phys. stat. sol (a) 61, 675 (1980)
 *      [2] B. E. Warren, X-Ray Diffraction, Addison-Wesley (1969)
 *
 */
class FaultsBroadeningEffectWC11m23: public ReflectionProfileComponent {
private:
  /// Fault probability
  REAL mAlpha;
public:
  /// Constructor
  FaultsBroadeningEffectWC11m23();
  /// Get peak profile (returns Four.Coefs. here)
  virtual CrystVector_REAL GetProfile (const CrystVector_REAL &x, const REAL xcenter, const REAL h, const REAL k, const REAL l);
  /// Get a rough estimate of peak width
  virtual REAL 	GetApproxFWHM (const REAL xcenter, const REAL h, const REAL k, const REAL l);
  /// Is this effect of real space type? (Yes - it returns. Four.Coefs. insteed of intensity)
  virtual bool 	IsRealSpaceType () const;
  /// Is this effect anisotropic? (Yes - it is (hkl) dependent)
  virtual bool 	IsAnisotropic () const;
  /// Get an additional correction for peak position (zero here)
  virtual REAL GetPositionCorr (const REAL xcenter, const REAL h, const REAL k, const REAL l) const;
}; // class FaultsBroadeningEffectWCHagege11m23

class PseudoVoigtBroadeningEffectA: public ReflectionProfileComponent {
private:
  REAL mCagliotiU;
  REAL mCagliotiV;
  REAL mCagliotiW;
  REAL mPseudoVoigtEta0;
  REAL mPseudoVoigtEta1;
  REAL mAsym0;
  REAL mAsym1;
  REAL mAsym2;
  REAL mAsymXMax;
public:
  PseudoVoigtBroadeningEffectA();
  CrystVector_REAL GetProfile(const CrystVector_REAL &x,
			      const REAL xcenter,
			      const REAL h, const REAL k, const REAL l);
  REAL GetApproxFWHM(const REAL xcenter,
		     const REAL h, const REAL k, const REAL l)const;
  bool IsRealSpaceType()const;
  void SetProfilePar (const REAL fwhmCagliotiW, const REAL fwhmCagliotiU=0,
		      const REAL fwhmCagliotiV=0,
		      const REAL eta0=0.5, const REAL eta1=0.);
	void SetAsymXMax(const REAL asymXMax);
private:
  void 	InitParameters ();
}; // class PseudoVoigtBroadeningEffectA

class PseudoVoigtBroadeningEffect: public ReflectionProfileComponent,
                      public virtual ObjCryst::ReflectionProfilePseudoVoigt {
private:
  REAL mAccuracy;
  // REAL mParam;
  REAL mAsymXMax;
public:
  PseudoVoigtBroadeningEffect();
  CrystVector_REAL GetProfile(const CrystVector_REAL &x,
			      const REAL xcenter,
			      const REAL h, const REAL k, const REAL l);
  REAL GetApproxFWHM(const REAL xcenter,
		     const REAL h, const REAL k, const REAL l)const;
  bool IsRealSpaceType()const;
  void SetAsymXMax(const REAL asymXMax);
}; // class PseudoVoigtBroadeningEffect

class HKLPseudoVoigtBroadeningEffectA: public ReflectionProfileComponent {
public:
  class HKLProfilePar {
  public:
    REAL dx;
    REAL fwhm;
    REAL eta;
  };
private:
  ReflStore mReflStore;
public:
  HKLPseudoVoigtBroadeningEffectA();
  virtual ~HKLPseudoVoigtBroadeningEffectA();
  CrystVector_REAL GetProfile(const CrystVector_REAL &x,
			      const REAL xcenter,
			      const REAL h, const REAL k, const REAL l);
  REAL GetApproxFWHM(const REAL xcenter,
		     const REAL h, const REAL k, const REAL l)const;
  bool IsRealSpaceType()const;
  void SetProfilePar(int h,int k,int l,REAL dx,REAL fwhm,REAL eta,bool dx_fixed=true,
		     bool fwhm_fixed=true,bool eta_fixed=true);
  REAL GetPositionCorr(const REAL xcenter,
		       const REAL h, const REAL k, const REAL l)const;
}; // class HKLPseudoVoigtBroadeningEffectA

class ReflectionProfile: public ObjCryst::ReflectionProfile {
public:
  class ReflCalcData {
  public:
    REAL s1, s0, ds;
    int N;
    vector<double> vProfile;
  };
protected:
  const ObjCryst::UnitCell* mpUnitCell;
  const ObjCryst::Radiation* mpRadiation;
  REAL mOmega;
  CrystVector_REAL mStressCoeff;
  REAL mFactor;
  // internal working variables
  int mN;
  vector< complex<double> > mvIn;
  vector< complex<double> > mvOut;
  vector<double> mvProfile;
  double mdL;
  double mds;
  double mQ;
  double ms0;
  double ms1;
  double mLambda;
  int mNbPoints;
  CrystVector_REAL mvL;
  CrystVector_REAL mvs;
  ReflStore mReflStore;
  ObjCryst::RefinableObjClock mClockReflectionProfileCalc;
protected:
  ObjCryst::ObjRegistry<ReflectionProfileComponent> mReflectionProfileComponentRegistry;
  const ObjCryst::PowderPatternDiffraction* mpParentPowderPatternDiffraction;
public:
  ReflectionProfile(const ObjCryst::UnitCell& cell, const ObjCryst::Radiation& radiation);
  ReflectionProfile(const ReflectionProfile &old);
  virtual ~ReflectionProfile();
  ReflectionProfile* CreateCopy()const;
  virtual const string& GetClassName()const;

  void AddReflectionProfileComponent(ReflectionProfileComponent &);
  long GeReflectionProfileComponentNb() const;
  const ReflectionProfileComponent & GetReflectionProfileComponent (long ) const;
  ReflectionProfileComponent & GetReflectionProfileComponent (long );
  const ReflectionProfileComponent & GetReflectionProfileComponent (const string &objName) const;
  ReflectionProfileComponent & GetReflectionProfileComponent (const string &objName);

  virtual void SetParentPowderPatternDiffraction(const ObjCryst::PowderPatternDiffraction &);
  virtual const ObjCryst::PowderPatternDiffraction& GetParentPowderPatternDiffraction()const;
  //void SetProfilePar(const REAL m, const REAL sigma);

  bool IsAnisotropic()const;
  CrystVector_REAL GetProfile (const CrystVector_REAL &x, const REAL xcenter,
			       const REAL h, const REAL k, const REAL l);
  REAL GetApproxFWHM (const REAL xcenter,
		      const REAL h, const REAL k, const REAL l)const;
  REAL GetFullProfileWidth (const REAL relativeIntensity, const REAL xcenter,
			    const REAL h, const REAL k, const REAL l);
  REAL GetIntegralWidth (const REAL xcenter,const REAL h,const REAL k,const REAL l); 
  void XMLOutput (ostream &os, int indent=0) const {};
  void XMLInput (istream &is, const ObjCryst::XMLCrystTag &tag) {};
  void SetIncidenceAngle(const REAL omega);
  REAL GetIncidenceAngle(const REAL xcenter) const;
  virtual REAL GetPositionCorr (const REAL xcenter,
		                const REAL h, const REAL k, const REAL l) const;
protected:
  
private:
  void InitParameters();
  bool PrepareForCalc(const CrystVector_REAL &x, const REAL xcenter,
		      REAL h, REAL k, REAL l);
 private:
  /// Clock when phenomenological stress correction parameters (q) were chnged last time.
  ObjCryst::RefinableObjClock mClockStressCorrQ;

}; // class ReflectionProfile

class DoubleComponentReflectionProfile: public ObjCryst::ReflectionProfile {
public:
	DoubleComponentReflectionProfile();
	~DoubleComponentReflectionProfile();
	DoubleComponentReflectionProfile(const DoubleComponentReflectionProfile &old);
	const string& GetClassName() const;
	DoubleComponentReflectionProfile* CreateCopy()const;
  CrystVector_REAL GetProfile(const CrystVector_REAL &x, const REAL xcenter,
												      const REAL h, const REAL k, const REAL l);
	REAL GetFullProfileWidth(const REAL relativeIntensity, const REAL xcenter,
			    								 const REAL h, const REAL k, const REAL l);
	bool IsAnisotropic()const;
	void XMLOutput (ostream &os, int indent=0) const {};
  void XMLInput (istream &is, const ObjCryst::XMLCrystTag &tag) {};
  
  void SetComponents(ObjCryst::ReflectionProfile &comp1, ObjCryst::ReflectionProfile &comp2);
  void SetProfileParams(const REAL weight);
protected:
	REAL mWeight;
	ObjCryst::ReflectionProfile* mComponent1;
	ObjCryst::ReflectionProfile* mComponent2;
	void InitParameters();
}; // class DoubleComponentReflectionProfile

/** ReflectionPositionCorrBase : Auxilliary base class for various reflection shift corrections.
 *
 * This is a base (virtual) class to describe diffraction peak shifts due to various
 * effects.
 * 
 * The essential GetPositionCorr(center,h,k,l) method should return the appropriate
 * reflection position correction in 2Theta(rad).
 * 
 * TODO:: Currently, it is based on the MStruct::ReflectionProfileComponent class.
 *   A general class describing reflection shift effects shoud be designed and implemented.
 * 
 **/
class ReflectionPositionCorrBase: public ReflectionProfileComponent
{
public:
	/// Constructor
	ReflectionPositionCorrBase();
	/// Copy constructor
	ReflectionPositionCorrBase(const  ReflectionPositionCorrBase & );
	/// Name of this class (MStruct::ReflectionPositionCorrBase)
	virtual const string& GetClassName()const;
	/// Get reflection profile (delta function in this case)
  CrystVector_REAL GetProfile(const CrystVector_REAL &x,
			      				const REAL xcenter,
			      				const REAL h, const REAL k, const REAL l);
	/// Get an approxiamative value (guess) of the FWHM (0.rad in this case)
  REAL GetApproxFWHM(const REAL xcenter,
		     		   const REAL h, const REAL k, const REAL l)const;
	/// Is this effect of real space type ? (NO in this case)
  	bool IsRealSpaceType()const;
  /// Get 2Theta(rad) correction for a given reflection
  	virtual REAL GetPositionCorr(const REAL xcenter,
		       					 const REAL h, const REAL k, const REAL l)const;
}; // ReflectionPositionCorrBase

/** RefractionPositionCorr : This calculates reflection position correction
 *   due to the refraction effect.
 * 
 * In the reflection geometry, if the incidence or exit angle of x-rays
 * are small (close to the critical angle of total external reflection) 
 * rerefraction effect causes a nonnegligible 2Theta shift of diffraction
 * lines [1-3].
 * 
 * This calculates the diffractions lines shift due to the refraction effect
 * according to the Snell's law:
 *   
 *         n_re(film) * cos(alpha_film) = cos(alpha_vac) ,
 * 
 * where n_re(film) is a real part of a material refraction coefficient and
 * alpha are the incidence angles in the film and vacuum respectively.
 * 
 * It can be written:
 * 
 *         n(material) ~ 1 + chi0(material)/2 ,
 * 
 * where chi0 is a material dielectric permittivity coefficient. Permittivity
 * is connected with electron density end hence also with the structure
 * factor:
 * 
 * 				chi0(matrial) = -4*pi*rel^2/Vel*K^2 * F(hkl=0) ,
 * 
 * where 'rel' is the classical electron radius, Vel is the UnitCell volume,
 * K = 2*pi/Lambda and F(hkl) is the structure facor.
 * 
 * Moreover it is often a case that due to some porosity or other effects
 * the film density is a little bit different from (theoretical/structural)
 * density. If it can be written for densities:
 * 
 *        rou(film) = nr * rou(material),
 * 
 * where 'nr' is a relative film density, than also material permittivity
 * can be multiplied with the same number:
 * 
 * 			  chi0(film) = nr * chi0(material) ~ nr * F(hkl=0).
 * 
 * This relative density 'nr' is a refinable parameter of this effect.
 * 
 * Material permittivity chi0(material) can be calculated from the crystal
 * structure or set directly.
 * 
 * [1] Hart ...
 * [2] G.Lim,W.Parrish,C.Ortiz,M.Bellotto,M.Hart, J. Mater.Res. (1987) 2, 471
 * [3] P.Colombi,P.Zanola,E.Bontempi,L.E.Depero, Spectrochimica Acta B (2007) 62, 554
 *   
 * blab, bla, bla
 */
class RefractionPositionCorr: public ReflectionPositionCorrBase
{
public:
	/// Constructor
	RefractionPositionCorr();
	/// Copy constructor
	RefractionPositionCorr(const RefractionPositionCorr & );
	/// Name of this class (MStruct::RefractionPositionCorr)
	virtual const string& GetClassName()const;
	/** \brief Get 2Theta(rad) correction (same for all reflections)
	 * 
	 * For dielectric permittivity and the index of refraction it is written:
	 * 
	 *     n = 1 - delata - ii * beta ~= 1 + chi0/2 , hence
	 * 
	 *     delta ~= - Re(chi0/2)   and    beta ~= - Im(chi0/2) ,
	 *   
	 * where chi0 is the model permittivity - it means permittivity
	 * multiplied by the relative density model parameter.
	 * 
	 * The 2Theta correction is calcualted according to [1-2]:
	 * 
	 *  delta(2Theta) ~= alpha - 1/sqrt(2) * { (alpha-2*delta) + 
	 *          + [ (alpha-2*delta)^2 + 4*beta^2 ]^(1/2) }^(1/2)
	 * 
	 *  
	 * [1] M.F.Toney, S.Brennan, Phys.Rev.B 39 (1989), 7963
	 * [2] T.Noma,K.Takada,A.Iida,X-Ray Spectrom. 28 (1999), 433-439
	 * 
	 */ 
	virtual REAL GetPositionCorr(const REAL xcenter,
							       					 const REAL h, const REAL k, const REAL l)const;
	/** \brief Set Crystal object used to calculate index of refraction
	 * 
	 * The supplied crystal structure is used to calculate chi0 value
	 * of the material:
	 * 
	 *     chi0 = -4*pi*rel/(Vcell*K^2) * F(hkl=0) ,
	 * 
	 * where rel is the classical electron radius rel = 2.8179e-5 A,
	 * Vcell is the unit cell volume, K = 2*pi/Lambda and F(hkl=0) is
	 * the structure factor for the reflection (h,k,l)=(0,0,0).
	 *  
	 * The chi0 value is computed after the Crystal is set and it is considered
	 * to not change (during further calculations) if the \param \fixed
	 * is \val true. In opossite case after any modification of the crystal structure
	 * or the unit cell chi0 value is recalculated. 
	 * 
	 * This method is complementary to the \name SetChi0(...) methods, which
	 * can be used to directly set the chi0 value or \name SetChemicalFormula(...) method
	 * calculating the chi0 vlaue from a chemical formula and material absolute density.
	 * 
	 */ 
	void SetCrystal(const ObjCryst::Crystal & , const bool fixed=true);
	/** \brief Set chi0 value used to calculate index of refraction
	 * 
	 * This method is complementary to the \name SetCrystal(...) or
	 * \name SetChemicalFormula(..) methods, which set a crystal structure and a unit cell
	 * or chemical formula and absolute material density used to calculate
	 * material refraction index.
	 * 
	 */
	void SetChi0(const complex<REAL> & chi0);
	/** \brief Set chemical formula amd absolute density of material used to calculate index of refraction
	 * 
	 * This method is complementary to the \name SetChi0(...) or
	 * \name SetCrystal(..) methods, which set a chi0 value or crystal structure
	 * and a unit cell used to calculate material refraction index.
	 * 
	 */
	void SetChemicalFormula(const string & formula, const REAL density);
	/** Returns chi0 value used to calculate index of refraction
	 * 
	 *       n ~= 1 + chi0/2
	 * 
	 * The value was either set directly or calculated from a Crystal structure.
	 * 
	 * If \param forceReCalc is \val \true than a recalculation of the chi0 value
	 * is forced.
	 *
	 */
	const complex< REAL > & GetChi0(const bool forceReCalc=false)const;
	/** \brief Set relative density
	 * 
	 * Chi0 value (supplied directly to the Object or calculated from the Crystal structure
	 * or chemical formuala and absolute density) is multiplied by this relative density
	 * factor to obtain a 'model' chi0 value used to calculate refraction index
	 * of the material. This should correct for various porosity effects etc.
	 * 
	 */
	void SetParams(const REAL density);
	
	/// Chi0 value set directly
	static const int CHI0_VALUE        = 0;
		/// Chi0 calculated from crystal structure
	static const int CHI0_CRYSTAL      = 1;
	/// Chi0 calculated from chemical formula and density
	static const int CHI0_CHEM_FORMULA = 2;
	/** \brief Class representing chemical formula
	 * 
	 * This class encapsulates on object representing chemical formula (string, ex.: "PbSO4")
	 * by a list of atomic scattering powers (ex.: Pb,S,O) and their counts (ex.: 1,1,4). This
	 * can be used eg. for forward and resonant scatting factors calculation.  
	 * 
	 */
	class ChemicalFormula {
	public:
		/// Constructor
		ChemicalFormula();
		/// Copy constructor
		ChemicalFormula(const ChemicalFormula & );
		/// Name for this class ("MStruct::RefractionPositionCorr::ChemicalFormula") 
		virtual const string & GetClassName() const;
		/// Set object name
		void SetName(const string & name);
		/// Set chemical formula (Some information is printed into the stream \param s.)
		void SetFormula(const string & formula, ostream & s = cout);
		/// Get chemical formula represented by this object
		const string & GetFormula() const;
		/// Get list of scattering atom powers of this object
		const vector< ObjCryst::ScatteringPowerAtom > & GetScatteringPowerAtomList() const;
		/// Get list of counts of scattering atom powers of this object
		const vector< int > & GetScatteringPowerAtomCountList() const;
		/// Get clock when the chemical formula was changed for the last time.
		const ObjCryst::RefinableObjClock & GetClock() const;
	protected:
		/// Object name
		string mName;
		/// Chemical formula (ex.: "PbSO4")
		string mFormula;
		/// The vector of scattering atom powers (ex.: Pb,S,O)
		vector< ObjCryst::ScatteringPowerAtom > mvScatt;
		/// The vector of scattering atom powers counts in the chemical formula (ex.: 1,1,4)
		vector< int > mvScattCount;
		/// When the chemical formula was set/changed for the last time?
		ObjCryst::RefinableObjClock mClock;
	}; // ChemicalFormula
protected:
	/// Relative density of material
	REAL mDensity;
	/// Pointer to the Crystal object (can be NULL)
	const ObjCryst::Crystal * mpCrystal;
	/// Absolute density (g/cm^3) of material (can be undefined)
	mutable REAL mAbsDensity;
	/// Chemical formula of material (can be undefined)
	ChemicalFormula mFormula;
	/// Choice flag how chi0 value is specified (if a value is specified directly
	/// or calculated from Crystal structure or chemical formula and absolute density)
	int mChi0ValueChoice;
	/// Flag if the Crystal (structure + unit cell) should be considered as nonchanging
	bool mConsiderCrystalFixed;
	/// Chi0 value (directly specified or computed from Crystal or chemical formula and absolute density)
	mutable complex< REAL > mChi0;
	/// When the material chi0 value was set/calculated for the last time? 
	mutable ObjCryst::RefinableObjClock mClockChi0;
	/// Initialization flag - currently used to signalize if chi0 value has been computed
	mutable int mInitializationFlag;
	/// Initialise object parameters (internal auxilliary function)
	void InitParameters();
}; // RefractionPositionCorr

/** XECsObj : This calculates x-ray elastic constants (XECs).
 * 
 * This is a base (purly virtual) class for XECs encapsulation.
 * 
 * X-ray elastic constants (XECs) S1 and S2 appear in many models
 * as scaling quantities between elastic lattice deformation measured
 * by diffraction and stress values in the sample [1,2] (See Class
 * ResidualStressPositionCorrection). XECs reflects elastic properties
 * of the material under study and an elastic grain interaction model
 * desribing stress/deformation state of crystallites in the sample.
 * XECs are, in general, hkl dependent.
 * 
 * This class shoud be a base class for various of such XECs computation
 * models (isotropic/anisotropic, various grain interaction type).
 * 
 * [1] I.C.Noyan, J.B.Cohen, Residual Stress (1987), Berlin: Springer Verlang
 * [2] M.Dopita, D.Rafaja, Z.Kristallogr. Suppl. 23 (2006) 67-72
 * 
 */
class XECsObj : virtual public ObjCryst::RefinableObj
{
public:
	/// Constructor
	XECsObj();
	/// Copy constructor
	XECsObj(const XECsObj & );
	/// Destructor
	virtual ~XECsObj();
	/// Name of this class (MStruct::XECsObj)
	virtual const string& GetClassName()const;
	/// Get XECs (s1 and s2) for the given reflection (in 1/GPa)
	virtual void GetXECs(REAL & s1, REAL & s2, const REAL h=0., const REAL k=0., const REAL l=0.)const = 0;
};

/** ResidualStressPositionCorrection : This calculates reflection position correction
 *   due to the residual stress in the sample.
 * 
 * A simple model of a rotationally symmetrical biaxial surface parallel state of the residual 
 * stress is assumed (quasi-isotropic specimen, no shear stress).
 * 
 * Relationship between a deformation 'e' measured by diffraction and the residual stress
 * in the specimen 'sigma' can be than written (e.g. [1-4]) as:
 * 
 *          e = 1/2 s2 * sigma * sin(psi)^2 + 2 s1 * sigma ,
 * 
 * where 'sigma' is the stress value (tensile(+), copressive(-)), 'psi' is the inclination
 * angle of the diffracting lattice planes from the sample surface and 's1' and 's2'
 * are the x-ray elastic constants (XECs).
 * 
 * A distance d(hkl) of diffracting lattice planes is affected by the stress state in the sample
 * and hence it is different from the stress free (relaxed) unit cell defined lattice planes distance
 * d0(hkl). Deformation measured by the diffracting family of lattice planes (hkl) is:
 * 
 *         e(hkl) = ( d(hkl) - d0(hkl) ) / d0(hkl).
 * 
 * Calculated 2Theta reflection position correction is conected whit this deformation:
 * 
 *         Corr ( 2Theta(rad) ) = -2 * e(hkl) * tan(Theta) .   
 * 
 * Inclination angle 'psi' can be in the case of a coplanar (a)symmetric diffraction
 * given as: 
 *
 *         psi(hkl) = Theta(hkl) - Omega ,
 * 
 * where Omega is the incidence angle. 
 *
 * XEconstants depend on the elastic properties of the studied material and on the grain
 * interaction model introduced in averaging procedure of the stress state influence on
 * grains in the specimen. XECs values are supplied to this object by a XECsObj class.
 * 
 * [1] I.C.Noyan, J.B.Cohen, Residual Stress (1987), Berlin: Springer Verlang
 * [2] M.Dopita, D.Rafaja, Z.Kristallogr. Suppl. 23 (2006) 67-72
 * [3] H.P.Klug, L.E.Alexander, X-ray Diffract. Procedures, chap.11, Stress Meas. in Metals
 * [4] I.Kraus, N.Ganev, chap.19, X-ray analysis of inhomog. stress state, in Defect and
 *     Microstructure analysis by Diffraction, edited by R.L.Snyder, J.Fiala, H.J.Bunge  
 * 
 */
class ResidualStressPositionCorrection: public ReflectionPositionCorrBase
{
public:
	/// Constructor
	ResidualStressPositionCorrection();
	/// Copy constructor
	ResidualStressPositionCorrection(const ResidualStressPositionCorrection & );
	/// Destructor
	virtual ~ResidualStressPositionCorrection();
	/// Name of this class (MStruct::ResidualStressPositionCorrection)
	virtual const string& GetClassName()const;
	/// Get 2Theta(rad) correction for a given reflection
	virtual REAL GetPositionCorr(const REAL xcenter,
							       					 const REAL h, const REAL k, const REAL l)const;
	/// Set XECs object which will be used to get XECs.
	void SetXECsObj(XECsObj & );
	/// Set value of the Stress parameter (GPa)
	void SetParams(const REAL stress);
protected:
	/// Stress value
	REAL mStress;
	/// Pointer to the XECs object used to get XECs
	XECsObj * pXECsObj;
	/// Initialise object parameters (internal auxilliary function)
	void InitParameters();
}; // ResidualStressPositionCorrection

/** XECsIsotropic : This calculates x-ray elastic constants (XECs)
 * 	  in the case of an isotropic material characterised by its
 *    Young's modulus and Poisson's ratio.
 * 
 * XECs are in this model independent of the crystallographic direction
 * and can be simply calculated as [1-3]:
 * 
 *   1/2*s2 = (1+ni)/E      ,   s1 = -ni/E  ,
 * 
 * where 'E' is Young's modulus and 'ni' is Poisson's ratio.
 * 
 * [1] I.C.Noyan, J.B.Cohen, Residual Stress (1987), Berlin: Springer Verlang
 * [2] H.P.Klug, L.E.Alexander, X-ray Diffract. Procedures, chap.11, Stress Meas. in Metals
 * [3] I.Kraus, N.Ganev, chap.19, X-ray analysis of inhomog. stress state, in Defect and
 *     Microstructure analysis by Diffraction, edited by R.L.Snyder, J.Fiala, H.J.Bunge  
 * 
 */
class XECsIsotropic : public XECsObj
{
public:
	/// Constructor
	XECsIsotropic();
	/// Copy constructor
	XECsIsotropic(const XECsIsotropic & );
	/// Name of this class (MStruct::XECsIsotropic)
	virtual const string& GetClassName()const;
	/// Get XECs (s1 and s2) for the given reflection (in 1/GPa)
	virtual void GetXECs(REAL & s1, REAL & s2, const REAL h=0., const REAL k=0., const REAL l=0.)const;
	/// Set material elastic parameters (Young's modulus (GPa), Poisson's ratio)
	void SetElasticParameters(const REAL E, const REAL ni);
protected:
	/// Young's modulus of the material (GPa)
	REAL mE;
	/// Poisson's ratio of the material
	REAL mNi;
	/// Initialise object parameters
	void InitParameters();
};

/** XECsReussVoigt : This calculates x-ray elastic constants (XECs)
 * 	  in the case of an anisotropic material characterised by its
 *    single crystal elastic constants.
 * 
 * The XECs are calculated as weighted average of XECs int the case
 * of the Reuss and the Voigt grain interaction model [1-3]:
 * 
 *   Sn = (1-w) * Sn_Reuss(hkl) + w * Sn_Voigt , 
 * 
 * where Sn are XECs (n=1,2). As it is indicated above, in the case
 * of the Reuss model XECs are dependent on the crystallographic
 * direction (hkl).
 * 
 * XECs are computed according to simplified equations from [2]
 * for any Laue group symmetry of the material. An Unitcell object
 * should be supplied and a list of required single-crystal
 * stifness constants tensor components is automatically created.
 * 
 * The Reuss/Voigt model weigth and the single crystal stiffness
 * constants are the model parameters.
 * 
 * [1] H.Behnken, V.Hauk, Z.Metall.(IJMR) 77 (1986)H.9 620−626
 * [2] J.C.Popa, J.Appl.Cryst. 33 (2000) 103-107
 * [3] M.Dopita, D.Rafaja, Z.Kristallogr. Suppl. 23 (2006) 67-72 
 * 
 */
class XECsReussVoigt : public XECsObj
{
public:
	/// Constructor
	XECsReussVoigt();
	/// Copy constructor
	XECsReussVoigt(const XECsReussVoigt & );
	/// Destructor
	virtual ~XECsReussVoigt();
	/// Name of this class (MStruct::XECsReussVoigt)
	virtual const string& GetClassName()const;
	/// Get XECs (s1 and s2) for the given reflection (in 1/GPa)
	virtual void GetXECs(REAL & s1, REAL & s2, const REAL h=0., const REAL k=0., const REAL l=0.)const;
	/** \brief Set UnitCell object
	 * 
	 * The UnitCell Object is required for computation of directional cosines of
	 * crystallographic directions. XECs are only slightly dependent on
	 * the UnitCell metrics and hence if \param fixed is \val true, than
	 * an approxiamtion of a fixed UnitCell is used: for a given hkl XECs
	 * are calculated with actual UnitCell parameters some time after the UnitCell
	 * was set and than XECs are usually not recalculated even thought UnitCell
	 * parameters can be changing/modified (e.g. becouse of refinement).
	 * 
	 * \return a flag indicating recognised (Laue group) stiffness\compliance
	 * tensor symmetry.
	 *  
	 */
	int SetUnitCell(const ObjCryst::UnitCell & , const bool fixed=true);
	/** \brief Set single crystall Stiffness constants
	 * 
	 * An UnitCell Object should be set before elastic
	 * constants are supplied. If UnitCell is set, a material symmetry
	 * can be recognised and a list of required stiffnes constant
	 * names (e.g. {"C11","C12","C44"} for a cubic symmetry) created. The list
	 * is returned by \name GetStiffnessConstantsNames() method.
	 * 
	 * A map of constats names (e.g. "C11") and values (in GPa) of single crystal
	 * stiffness constants should be supplied. This map should contain all
	 * stiffnessconstants listed in the list returned by
	 * \name GetStiffnessConstantsNames() method. Otherwise they are considered to
	 * be zero. Constants names should be same as those in the list. Redundant
	 * values are ignored (without warning).
	 * 
	 * The Stiffness constant matrix is initialised after the constants
	 * were set and can be returned by \name GetCijMatrix() method.
	 *
	 */
	void SetStiffnessConstants(const map< string, REAL > & );
	/// Return a vector of stiffness constants names (e.g. {"C11","C12","C44"} for a cubic symmetry)
	/// See description of \name SetStiffnessConstants(...) method.
	const vector< string > & GetStiffnessConstantsNames()const;
	/// Get matrix (6x6) of Stiffness constants in the Voigt notation (in GPa)
	const CrystMatrix_REAL & GetCijMatrix()const;
	/// Flags for Stiffness tensor settings for different Laue groups
	static const int LAUE_UNDEFINED         =  0; // undefined
	static const int LAUE_TRICLINIC         =  1; // triclinic -1
	static const int LAUE_MONOCLINIC_AXIS_C =  2; // monoclinic 2/m unique axis c
	static const int LAUE_MONOCLINIC_AXIS_B =  3; // monoclinic 2/m unique axis b
	static const int LAUE_ORTHORHOMBIC      =  4; // orthorhombic 2/mmm
	static const int LAUE_TETRAGONAL_LOW    =  5; // tetragonal 4/m
	static const int LAUE_TETRAGONAL_HIGH   =  6; // tetragonal 4/mmm
	static const int LAUE_TRIGONAL_LOW      =  7; // trigonal -3
	static const int LAUE_TRIGONAL_HIGH     =  8; // trigonal -3m
	static const int LAUE_HEXAGONAL         =  9; // hexagonal 6/m, 6/mmm
	static const int LAUE_CUBIC             = 10; // cubic m-3, m-3m
	/// Set model parameters (Reuss-Voigt model weigth)
	void SetParams(const REAL weight);
protected:
	/// Reuss/Voigt model weight (0..Reuss, 1..Voigt)
	REAL mWeight;
	/// Pointer to the material UnitCell Object
	const ObjCryst::UnitCell * mpUnitCell;
	/// Flag if the UnitCell lattice parameters should be considered as nonchanging
	bool mConsiderUnitCellFixed;
	/// Flag denoting Stiffness tensor settings for a given UnitCell lattice symmetry
	int mGroupType;
	/// Vector of single crystal Stiffness constants names (e.g. {"C11","C12","C44"} for a cubic symmetry)
	vector< string > mvCijConstantsNames;
	/// Vector of single crystal Stiffness constants values
	CrystVector_REAL mCijValues;
	/// When single crystal Stiffness constants values were changed for the last time?
	ObjCryst::RefinableObjClock mClockCijValues;
	/// Single crystal Stiffness constants matrix
	mutable CrystMatrix_REAL mCijMatrix;
	/** \brief Single crystal Compliance constants matrix
	 * 
	 * \note Internally Popa notation is used. In this notation it is written:
	 * 
	 *  (e11,e22,e33,e23,e31,e12)^T = Sij * (s11,s22,s33,2*s23,2*s31,2*s12) ,
	 * 
	 * hence definition of the Complience tensor is a little bit different from
	 * the common Voigt notation:
	 *  
	 *  (e11,e22,e33,2*e23,2*e31,2*e12)^T = Sij * (s11,s22,s33,s23,s31,s12)
	 *   
	 * \name mSijMatrix is written in the internal (Popa) notation.
	 */
	mutable CrystMatrix_REAL mSijMatrix;
	/// Voigt S1 XEC (Calculated when Cij matrix is set) 
	mutable REAL mS1Voigt;
	/// Voigt S2 XEC (Calculated when Cij matrix is set) 
	mutable REAL mS2Voigt;
	/// When Cij matrix was set, Sij matrix and Voigt XECs recalculated for the last time?
	mutable ObjCryst::RefinableObjClock mClockMatricesCalc;
	/// Calculate matrices (Cij, Sij, S1Voigt, S2Voigt) (from cij values)
	void CalcMatrices()const;
	/// Initialise object parameters (not Cij constants parameters)
	void InitParameters();
	/// Auxilliary class for Stiffness/Complience tensor indices manipulation
	class MC_indices {
		public:
			/// Constructor
			MC_indices(int nb, const int indices[][2])
			{
				for (int i=0; i<nb; i++) {
					indices_.push_back( make_pair(indices[i][0],indices[i][1]) );
					ostringstream s;
					s << indices[i][0] << indices[i][1];
					indices_str_.push_back(s.str());
				}
			}
			/// Returns vector of tensor constants indices (e.g. ((1,1),(1,2),(4,4)) for cubic symmetry)
			const vector< pair<int,int> > & indices() const { return indices_; }
			/// Returns vector of tensor constants indices as strings (e.g. ("11","12","44") for cubic symmetry)
			const vector< string > & indices_str() const { return indices_str_; }
		private:
			/// Tensor indices pairs (e.g. (1,1))
			vector< pair<int,int> > indices_;
			/// Tensor indices pairs as strings (e.g. "11")
			vector< string > indices_str_;
	};
	/// Arrays of tensor indices for different Laue group symmetries
	static const int MCij_a_triclinic[][2];
	static const int MCij_a_monoclinic_axis_c[][2];
	static const int MCij_a_monoclinic_axis_b[][2];
	static const int MCij_a_orthorhombic[][2];
	static const int MCij_a_tetragonal_low[][2];
	static const int MCij_a_tetragonal_high[][2];
	static const int MCij_a_trigonal_low[][2];
	static const int MCij_a_trigonal_high[][2];
	static const int MCij_a_hexagonal[][2];
	static const int MCij_a_cubic[][2];
	/// Auxilliary structures for tensor indices manipulation different Laue group symmetries
	static const MC_indices MCij_undefined;
	static const MC_indices MCij_triclinic;
	static const MC_indices MCij_monoclinic_axis_c;
	static const MC_indices MCij_monoclinic_axis_b;
	static const MC_indices MCij_orthorhombic;
	static const MC_indices MCij_tetragonal_low;
	static const MC_indices MCij_tetragonal_high;
	static const MC_indices MCij_trigonal_low;
	static const MC_indices MCij_trigonal_high;
	static const MC_indices MCij_hexagonal;
	static const MC_indices MCij_cubic;
	/// Current material constants tensor indices (should correspod to material symmetry flag \name mGroupType)
	MC_indices mMCij;
	/// Auxilliary structure, where already calculated Reuss XECs constants are stored
	mutable ReflStore mReflStore;
	/// When the first data were saved in the Reflection data store
	mutable ObjCryst::RefinableObjClock mClockReflStore;
};

// an interpolation routine
double interp1(const vector<double> &vx, const vector<double> &vy,
	       const double x);

// a convolution routine
double conv2(const vector<double> &fx,const vector<double> &fy,
	     const vector<double> &gx,const vector<double> &gy,
	     const double x);


/** TextureOdfBase : The base class for representation of the ODF. This
 * class can be supplied for example to the OdfNumCalculator. This can be
 * used for example to calculation of pole figers or inverse pole figures.
 * You should derive your own class to represent your peculiar ODF model.
 * For a proper implementation you shold create appropriate GetOdfValue,
 * IsOdfSymmetrised and IsOdfNormalized methods.
 */
class TextureOdfBase: virtual public ObjCryst::RefinableObj
{
public:
	/** This method should return the value of the model ODF
	 * for the selected orientation of crystallites defined by
	 * the first parameter which is the Eulerian rotation matrix (3x3)
	 * representing rotations in the (phi1,phi,phi2) notation.
	 * Values of the appropriate angles (phi1,phi,phi2) have to
	 * be supplied only if the AnglesInitialised parameter is true.
	 * If only angles are known function
	 * A = RotationTB::GetEulerMatrix(phi1,phi,phi2) can be used to supply
	 * proper parameters to the method. Function
	 * RotationTB::GetEulerAngles(phi1,phi,phi2,A) can be used
	 * to get the these angles in the opposite case.
	 */
	virtual REAL GetOdfValue(const CrystMatrix_REAL& A,
								  	  		 REAL phi1 = 0., REAL phi = 0., REAL phi2 = 0.,
												   const bool AnglesInitialised = false) const;
	/// Is the ODF properly symmetrised (has the correct crystal symmetry)
	virtual bool IsOdfSymmetric() const;
	/// Is the ODF properly normalised.
	virtual bool IsOdfNormalised() const;
};

/** TextureModelFiber : The representation of the ODF fiber texture model
 * with few types of the fiber texture shape functions.
 */
 class TextureModelFiber: public TextureOdfBase
 {
 protected:
 	/// hkl indexes of the fiber texture
 	CrystVector_REAL mHKL;
 	/// Fiber texture shape function type
 	int mShapeFuncType;
 	/// Texture shape function width parameter
 	REAL mWidthPar;
 	/// Texture shape function shape parameter
 	REAL mShapePar;
 	/// Texture hkl tilt parameters phi1, phi and phi2
 	CrystVector_REAL mTiltParams;
 	/// Pointer to the crystal unit cell
 	const ObjCryst::UnitCell* mpUnitCell;
public:
	/// Constructor
	TextureModelFiber();
	/// Destructor
	virtual ~TextureModelFiber();
	/// Set unit cell of the crystal structure
	void SetUnitCell(const ObjCryst::UnitCell* unitcell);
	/** Set the fiber texture hkl indexes
	 * Use forceFriedelLaw = true to force considering Friedel mates as equivalent
	 * (add a center of symetry).
	 */
	void SetFiberHKL(const int h, const int k, const int l, const bool forceFriedelLaw=false);
	/// Get the ODF value
	REAL GetOdfValue(const CrystMatrix_REAL& A,
				  	  		 REAL phi1 = 0., REAL phi = 0., REAL phi2 = 0.,
								   const bool AnglesInitialised = false) const;
	/// Is the ODF properly symmetrised (has the correct crystal symmetry)
	bool IsOdfSymmetric() const;
	/// Is the ODF properly normalised.
	bool IsOdfNormalised() const;
protected:
	/// Normalised orthonormal coordinates of all equivalent fiber HKL directions 
 	CrystMatrix_REAL mHKLDirections;
 	/// Normalised coordinates of the tlted normal sample direction
 	mutable CrystVector_REAL mTiltDirection;
 	/// Normalisation factor for the used shape function
 	mutable REAL mFuncNorm;
	/// Internal method for the object parameters and options initialisation
	void Init();
	/// Calculate the normalisation factor for the used shape function
	REAL CalcShapeFuncNorm() const;
	/// Clock for the texture shape parameters (width, shape)
	ObjCryst::RefinableObjClock mClockFuncParams;
	/// Clock to catch the event when the fiber HKL indexes was set last time
	mutable ObjCryst::RefinableObjClock mClockHKL;
	/// Clock for Tilt parameters
	ObjCryst::RefinableObjClock mClockTiltParams;
	/// Clock storing time when the coordinates of the reference direction were calculated
	mutable ObjCryst::RefinableObjClock mClockTiltDirection;
	/// Clock storing time when the shape function normalisation factor was calculated
	mutable ObjCryst::RefinableObjClock mClockNormCalculated;
 };
 
/** TextureOdfNumCalculator: Class implementing various methods mainly
 * for texture numerical calculations with the supplied ODF model.
 * Beside such tasks (as a pole figure or inverse pole figure simulation)
 * the class can be also used e.g. for the odf export in a text format.
 */
class TextureOdfNumCalculator: virtual public ObjCryst::RefinableObj
{
protected:
	/// Pointer to the current ODF model
	const TextureOdfBase* mpOdfModel;
	/// Pointer to the current unit cell used for crystallographic calculations
	const ObjCryst::UnitCell* mpUnitCell;
	/** Parameters of the ODF grid representation (for export, normalisation) - 
	 * (3x3) matrix with the structure:
	 *      { { phi1_min, phi1_max, phi1_step },
	 * 			  {  phi_min,  phi_max,  phi_step },
	 *        { phi2_min, phi2_max, phi2_step } }
	 */
	CrystMatrix_REAL mOdfGridParams;
public:
	/// Constructor
	TextureOdfNumCalculator();
	/// Destructor
	virtual ~TextureOdfNumCalculator();
	/// Set the current unit cell used for crystallographic calculations
	void SetUnitCell(const ObjCryst::UnitCell* unitcell); 
	/// Set the current ODF model. 
	void SetOdfModel(const TextureOdfBase* OdfModel);
	/// Get normalisation factor of the current ODF
	REAL GetOdfNorm() const;
	/** Exports the current ODF model in the Panalytical XPert like
	 * text format - not exactly the same used by XPert-Texture,
	 * however with almost similar data structure.
	 */
	void ExportOdfXPert(ostream& os) const;
	// pf calculation (for one point, multiple points)
	/** Calculate the integral (projection) of the ODF function over the crystallites
	 * having their crystal direction (hkl) parallel to the sample direction y =
	 * ( cos(phi)*sin(psi), sin(phi)*sin(phi), cos(psi) ).
	 *
	 * If a special crystal direction (hkl indexes are integer values) is supplied
	 * the calculation speed is optimised using the rotation symmetry of the (hkl)
	 * crystal axis. To the contrary, the sample symmetry has not been yet implemented
	 * and used here for the calcualtion speed optimization.
	 *
	 * Futhermore, if a special crystal direction (hkl indexes are integer values) is supplied
	 * and the current used ODF model is not properly symmetrised (IsOdfSymmetric returns false)
	 * the projection is averaged over all symmetry equivalent (hkl) directions and the Friedel
	 * Symmetry is forced. If the non-integer (hkl) indexes are supplied, the equivalent directions
	 * are not generated and if the ODF model is not correctly symmetric the calculated
	 * projection doesn't respect crystal symmetry in any way.
	 */
	REAL CalcOdfProjection(const REAL psi, const REAL phi,
								 				 const REAL h, const REAL k, const REAL l) const;
	/** Calculate the Pole figure for the (hkl) crystal direction in the sample
	 * directions specified by angles in input vectors psi (tilt) and phi (rotation).
	 * The output matrix contains calculated Pole Figure projection values in the grid
	 * defined by these input vectors. The appropriate calculated values of the projection
	 * for a constant psi tilt and for all phi rotations are stored in an apprpriate row
	 * of the output matrix.
	 */
	CrystMatrix_REAL CalcPFProjection(const CrystVector_REAL& psi, const CrystVector_REAL& phi,
												 						const REAL h, const REAL k, const REAL l) const;
	// ipf calcualtion (for one point and multiple points)
	/** Exports the supplied Pole Figure data in the Panalytical XPert like
	 * text format - not exactly the same used by XPert-Texture,
	 * however with almost similar data structure. (It is suitable
	 * to use this method immediately after using this object to
	 * calculate the Pole Figure data by the CalcPFProjection or
	 * CalcOdfProjection method. In this way more correct information
	 * is written in the output stream but still in most cases information
	 * in the header part of the exported stream is only ballast.)
	 * The optional parameter "scale" is the multiplication factor for
	 * the exported the Pole figure projection values.
	 */
	void ExportPFProjectionXPert(ostream& os,	const CrystVector_REAL& psi, const CrystVector_REAL& phi,
															 const CrystMatrix_REAL& pf_data, const REAL scale = 1.) const;
protected:
	/// Prepare the object for the odf projection calculation
	void PrepareForOdfProjectionCalc(const REAL h, const REAL k, const REAL l,
																	 const bool forceFriedelLaw=false) const;
	/** Auxilliary inverse rotation matrix(ces -all equivalent)
	 * used during the odf projection calculation to transform the given HKL
	 * into the normal direction (001).
	 */
	mutable std::vector<CrystMatrix_REAL> mAuxInvRotMatHKL;
	/// Indexes of the last HKL direction for which the ODF projection calculation was prepared
	mutable CrystVector_REAL mLastOdfProjectionHKL;
	/// Clocks to control the state of the object for the ODF projection calculation
	mutable ObjCryst::RefinableObjClock mClockOdfProjectionCalcPrepared;
	/** Default step for the ODF projection integration procedure - it is not ensured
	 * that exactly this step is used by the algorithm. The atualy used integration
	 * step has the approxiamtively same value but it is adjusted to fix the periodcity
	 * of the angular interval of rotation around the current integration axis
	 */ 
	REAL mSetOdfProjectionIntegStep;
	/// The actual ODF projection integration step adjusted by the prepare method
	mutable REAL mOdfProjectionIntegStep;
	/// The actual number of intervals for ODF projection integration
	mutable int mOdfProjectionIntegPointsNb;
};

/** TextureModelHKL : structure to represent a texture model
 * of a multiple component {HKL} texture. It should be able
 * to effectively calculate projection of the encapsulated
 * ODF function into the any supplied direction of a hkl
 * diffraction vector. 
 *
 * This keeps an ODF function describing the model for multiple components.
 * It also keeps fractional values of all {HKL} components and other
 * parameters for each component. The odf function and model parameters 
 * are supplied to an enclosed TextureCalculator object during
 * the calculation of the ODF projection.
 *
 */
class TextureModelHKL: virtual public ObjCryst::RefinableObj
{
public:
  /** Type of the texture model function - number of crystalities inclined
   * from a defined texture axis by an angle Omega is decreasing with this
   * angle like the well known Normal Distribution (Gauss) function.
   * The only one parameter of this model function is its FWHM.
   * 
   * ref: http://mathworld.wolfram.com/NormalDistribution.html
   */
  static const int TEXTURE_MODEL_HKL_GAUSS = 0; 
  /** Type of the texture model function - number of crystalities inclined
   * from a defined texture axis by an angle Omega is decreasing with this
   * angle like the function
   *                         ~ exp( - ((1-cos(Omega/N))/X)^n )
   *  
   * described e.g. in ref. Simek (JApplCryst.). This model function has
   * two parameters: the parameter X - conected with the function FWHM
   * and the shape parameter n. The constant N is equal to 1 if this function
   * is used for description of the tilt of the crystallites z-axis from
   * the prescribed z-{HKL} texture axis (Omega is running from 0 to PI)
   * and N = 2 in the case that the function is used for description
   * of the rotation around texture {HKL} axis (deviation angle Omega of
   * the crystallites x-axis from the selected x-{HKL} axis is running
   * from 0 to 2*PI).   
   * 
   * ref: D.Simek, R.Kuzel, D.Rafaja, J.Appl.Cryst. (2006) 39, p. 487-501
   */
	static const int TEXTURE_MODEL_HKL_SIMEK = 1;
	/** Type of the texture model function - number of crystalities inclined
   * from a defined texture axis by an angle Omega is decreasing with this
   * angle like the function
   *                         ~ | cos(Omega/N) |^n
   *  
   * described in ref. Birkholz (JApplCryst.). The only one parameter
   * of this model is the texture degree n. Constant N is equal to 2
   * if it is used for desribtion of tilt part of the {HKL} texture,
   * N = 4 if it is used for the rotational part (see description of
   * the SIMEK's-like model) and N = 1 in the case it is used for description
   * of the tilt part of the {HKL} texture symmetry and the symmetry of
   * the rotational part is set to AXIAL (this is done so, to be consistent
   * with the original Birkholz's model.
   * 
   * ref: M.Birkholz, J.Appl.Cryst. (2007) 40, p. 735-742
   */
  static const int TEXTURE_MODEL_HKL_BIRKHOLZ = 2;
	/** Type of the texture model function - this is just a simple constant
	 * value function indicating no dependence of the texture on the rotation
	 * around the {HKL} texture axis. Can be used only for the description of
	 * the phi1 or phi2 rotational part. (Note that the model is different
	 * from the classical definition by the idea phi2(Euler) = phi2 - phi1)
   */
	static const int TEXTURE_MODEL_HKL_AXIAL = 3;
protected:
	/// Number of {HKL} texture components
	int mNbComponents;
	/// Texture model type (phi1) for each {HKL} texture component
	std::vector< int* > mComponentTextureModels_phi1;
	/// Texture model type (phi) for each {HKL} texture component
	std::vector< int* > mComponentTextureModels_phi;
	/// Texture model type (phi2) for each {HKL} texture component
	std::vector< int* > mComponentTextureModels_phi2;
	/// Vector of fractions of all {HKL} texture components
	std::vector< REAL* > mComponentFractions;
	/// HKLz vectors for all components
	std::vector< CrystVector_REAL* > mComponentHKLz;
	/// HKLx vectors for all components (dummy values for axial textures)
	std::vector< CrystVector_REAL* > mComponentHKLx;
	/// All texture parameters for all {HKL} components 
	std::vector< CrystVector_REAL* > mComponentTextureParams;
	/// The odf function for texture calculation
	REAL funcodf(REAL phi,REAL th,REAL psi,REAL* aa=0) const;
	/// The TextureCalaculator object
	TextureCalculator mTextureCalculator;
	/// Euler phi1 angle of the whole texture function
	REAL mGlobalPhi1;
	/// Euler phi1 angle of the whole texture function
	REAL mGlobalPhi;
	/// Euler phi2 angle of the whole texture function
	REAL mGlobalPhi2;
	/// Texture parameters clock
	ObjCryst::RefinableObjClock mClockParams;
	/** Pointer to a Crystal object. It is used to get informatiion about the crystal
	 * symmetry and to help with some crystallogrphic transformations.
	 */
  const ObjCryst::Crystal* mpCrystal;
	/** This is used for considering the crystal symmetry during odf function generation. If true,
	 * a center of symmetry will be added to crystal (to force considering Friedel mates HKL
	 * crystal orientations as equivalent). This attribute affects only the symmetry of the odf
	 * function on its own.
	 */
	bool mbforceFriedelSymmetry;
	/// All symmetry equivalent HKLz vectors for all components - in the orthogonal coord.
	std::vector< std::vector< CrystVector_REAL > > mComponentEquivHKLz;
	/// HKLx vectors for all components (dummy values for axial textures)
	std::vector< std::vector< CrystVector_REAL > > mComponentEquivHKLx;
	/** Reference systems (properly rotated) for all texture components. E1 coordinates of
	 * the vectors of the reference basis are stored in the rows (e.g. the i-th row: coordinates
	 * of the i-th vector of the basis)
	 * */ 
	std::vector< CrystMatrix_REAL > mComponentReferenceSystems; 
public:
	/// Constructor
	TextureModelHKL();
	/// Copy constructor
	TextureModelHKL(const TextureModelHKL& old);
	/// Destructor
	virtual ~TextureModelHKL();
	/// Get class name (MStruct::TextureModelHKL)
	const string & GetClassName() const;
	/// Initialize basic parameters - used in constructor
	void InitParameters();
	/** Initialize parameters for the texture component no. comp_nb.
	 * It is assumed that parameters were set first. If you are
	 * changing component model etc., you should deregister all component
	 * parameters first (including the fraction parameter).
	 */
	void InitComponentParameters(const int comp_nb);
	/** Pointer to a Crystal object, which is used to get informatiion about the crystal
	 * symmetry and to help with some crystallogrphic transformations. If value of the another
	 * parameter - forceFriedelSymmetry - is set to true, center of symmetry is added to crystal
	 * symmetry and this will affect the symmetry of the odf function.  
	 */
	void SetCrystal(const ObjCryst::Crystal& crystal, const bool forceFriedelSymmetry=false);
	/** Prepare the object for the odf function calculation - generete all equivalent
	 * crystal orientations for all HKL components, transform their coordinates to
	 * the normal orthogonal system E1; create texture reference systems (properly rotated
	 * base) for all components; calculate the odf normalization factor.
	 */
	void PrepareForOdfCalc();
	void AddComponent(const int texture_model_phi=0,
										const int texture_model_phi2=3,
										const int texture_model_phi1=3,
										const REAL fraction = 1.,
										const CrystVector_REAL& HKLz=CrystVector_REAL(0),
										const CrystVector_REAL& HKLx=CrystVector_REAL(0),
										const CrystVector_REAL& params=CrystVector_REAL(0));
	/// Remove the selected component
	void RemoveComponent(const int comp_nb);
protected:
	/// Auxiliary structure to store refences to a part of a texture component parameters
	struct ComponentPartParams{
			/// pointers to parameters of parts (subsequently phi1,phi,phi2) of HKL texture component
			std::vector< REAL* > params;
			/// type of the model function
			int model;
	};
	/** Auxiliary function - parses a vector of a HKL texture component parameters
	 * and return a vector of ComponentPartParams encapsulating poiters to parameters
	 * for each part (subsequently phi1, phi, phi2) of the HKL texture component.
	 */
	std::vector< ComponentPartParams > SplitHKLComponentParams(
																												 CrystVector_REAL& params,
																												 const int texture_model_phi,
																												 const int texture_model_phi2,
																							 					 const int texture_model_phi1);
};

/// Auxiliary structure to store very simple information about a peak.
struct PeakParams {
	REAL xmax;
	REAL ymax;
	REAL intensity;
	REAL fwhm;
	REAL asym;
};

// auxiliary routine to calculate simple peak parmeters of numerical data representing peak-like function 
PeakParams CalcPeakParams(const CrystVector_REAL &x, const CrystVector_REAL &y, const int nbApproxPoints = 5);
/** Calculate simple characteristic parameters of a peak-like function, which is represented
 * numericaly in the x and y input data arrays. (Please consider that the data spacing and also
 * the range significantelly affects precision of the calculated values.) nbApproxPoints is
 * the number of points used for a function approximation by a quadratic polynomial around
 * the maximum function value to get a better guass of the maximum.
 */

extern const ObjCryst::RefParType *gpRefParTypeScattDataCorrHKLIntensity;
extern const ObjCryst::RefParType *gpRefParTypeMaterialData;
extern const ObjCryst::RefParType *gpRefParTypeScattDataProfileSizeDistrib;

// Global register for all XECs Objects
extern ObjCryst::ObjRegistry<XECsObj> gXECsObjRegistry;

class NiftyStaticGlobalObjectsInitializer_MStruct {
public:
  NiftyStaticGlobalObjectsInitializer_MStruct() {
    if (mCount++ == 0) {
      gpRefParTypeScattDataCorrHKLIntensity
	=new ObjCryst::RefParType(ObjCryst::gpRefParTypeObjCryst,
				  "HKL Intensity Corr.");
      gpRefParTypeMaterialData
	=new ObjCryst::RefParType(ObjCryst::gpRefParTypeObjCryst,
				  "Material Data");
      gpRefParTypeScattDataProfileSizeDistrib
	=new ObjCryst::RefParType(ObjCryst::gpRefParTypeScattDataProfileWidth,
				  "Crystallites Size Distrib.");	  
    }
  }
  ~NiftyStaticGlobalObjectsInitializer_MStruct() {
    if (--mCount == 0) {
      delete gpRefParTypeScattDataCorrHKLIntensity;
      gpRefParTypeScattDataCorrHKLIntensity=0;
      delete gpRefParTypeMaterialData;
      gpRefParTypeMaterialData=0;
      delete gpRefParTypeScattDataProfileSizeDistrib;
      gpRefParTypeScattDataProfileSizeDistrib=0;
    }
  }
private:
  static long mCount;
};
static NiftyStaticGlobalObjectsInitializer_MStruct
               NiftyStaticGlobalObjectsInitializer_MStruct_counter;

} // namespace MStruct

#endif /* __MSTRUCT__H__ */
