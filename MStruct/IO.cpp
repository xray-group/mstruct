/* 
 * IO.cpp
 * 
 * MStruct++ - Object-Oriented computer program/library for MicroStructure analysis
 * 					   from powder diffraction data.
 * 
 * Copyright (C) 2009-2014  Zdenek Matej, Charles University in Prague
 * Copyright (C) 2014-2018  Zdenek Matej, MAX IV Laboratory, Lund University
 * Copyright (C) 2016-2018  Milan Dopita, Jan Endres, Charles University in Prague
 * Copyright (C) 2017-2018  Jiri Wollman, Charles University in Prague
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License`
 * along with MStruct++. If not, see <http://www.gnu.org/licenses/>.
 * 
 */
 
#ifdef MSCV
#pragma warning(disable : 4996)
#endif // MSVC

#include "MStruct.h"
#include "IO.h"
#include "boost/format.hpp"

using namespace ObjCryst;

namespace MStruct {

////////////////////////////////////////////////////////////////////////
//
//    I/O PowderPatternDiffraction
//
////////////////////////////////////////////////////////////////////////

void PowderPatternDiffraction::XMLOutput(ostream &os,int indent)const
{
   VFN_DEBUG_ENTRY("MStruct::PowderPatternDiffraction::XMLOutput():"<<this->GetName(),5)
   for(int i=0;i<indent;i++) os << "  " ;
   XMLCrystTag tag("PowderPatternCrystal");
   tag.AddAttribute("Name",this->GetName());
   tag.AddAttribute("Crystal",this->GetCrystal().GetName());
   {
      stringstream ss;
      ss<<this->IsIgnoringImagScattFact();
      tag.AddAttribute("IgnoreImagScattFact",ss.str());
   }
   os <<tag<<endl;
   indent++;

   if(mFreezeLatticePar)
   {
      XMLCrystTag t("FrozenLatticePar");
      t.AddAttribute("a", (boost::format("%f")%mFrozenLatticePar(0)).str() );
      t.AddAttribute("b", (boost::format("%f")%mFrozenLatticePar(1)).str() );
      t.AddAttribute("c", (boost::format("%f")%mFrozenLatticePar(2)).str() );
      t.AddAttribute("alpha", (boost::format("%f")%(mFrozenLatticePar(3)*180/M_PI)).str() );
      t.AddAttribute("beta" , (boost::format("%f")%(mFrozenLatticePar(4)*180/M_PI)).str() );
      t.AddAttribute("gamma", (boost::format("%f")%(mFrozenLatticePar(5)*180/M_PI)).str() );
      t.SetIsEmptyTag(true);
      for(int i=0;i<indent;i++) os << "  " ;
      os<<t<<endl;
   }

   if(mpReflectionProfile!=0) mpReflectionProfile->XMLOutput(os,indent);

   this->GetPar(&mGlobalBiso).XMLOutput(os,"globalBiso",indent);
   os <<endl;

   if(mCorrTextureMarchDollase.GetNbPhase()>0)
   {
      mCorrTextureMarchDollase.XMLOutput(os,indent);
   }

   mCorrTextureEllipsoid.XMLOutput(os,indent);

   // Geometry - Omega
   {
	   XMLCrystTag t("Geometry");
	   t.AddAttribute("Omega", (boost::format("%f")%(mOmega*180./M_PI)).str() );
	   t.SetIsEmptyTag(true);
	   for(int i=0;i<indent;i++) os << "  " ;
	   os<<t<<endl;
   }
   
   mCorrAbsorption.XMLOutput(os,indent);
   
   //mCorrTexture.XMLOutput(os,indent);
   
   //mCorrHKLIntensity.XMLOutput(os,indent);
	   
   #if 0
   if(mFhklObsSq.numElements()>0)
   {
      XMLCrystTag tag2("FhklObsSq");
      for(int i=0;i<indent;i++) os << "  " ;
      os <<tag2<<endl;

      for(long j=0;j<this->GetNbRefl();j++)
      {
         for(int i=0;i<=indent;i++) os << "  " ;
         os << mIntH(j) <<" "
            << mIntK(j) <<" "
            << mIntL(j) <<" "
            << mFhklObsSq(j) <<endl;
      }

      tag2.SetIsEndTag(true);
      for(int i=0;i<indent;i++) os << "  " ;
      os <<tag2<<endl;
   }
   #else
   if(mpLeBailData!=0) mpLeBailData->XMLOutput(os,indent);
   #endif

   indent--;
   tag.SetIsEndTag(true);
   for(int i=0;i<indent;i++) os << "  " ;
   os <<tag<<endl;
   VFN_DEBUG_EXIT("MStruct::PowderPatternDiffraction::XMLOutput():"<<this->GetName(),5)
}

////////////////////////////////////////////////////////////////////////
//
//    I/O PowderPattern
//
////////////////////////////////////////////////////////////////////////

void PowderPattern::XMLInput(istream &is,const XMLCrystTag &tagg)
{
   VFN_DEBUG_ENTRY("MStruct::PowderPattern::XMLInput():"<<this->GetName(),5)
   for(unsigned int i=0;i<tagg.GetNbAttribute();i++)
   {
      if("Name"==tagg.GetAttributeName(i)) this->SetName(tagg.GetAttributeValue(i));
   }
   while(true)
   {
      XMLCrystTag tag(is);
      if(("PowderPattern"==tag.GetName())&&tag.IsEndTag())
      {
         this->UpdateDisplay();
         VFN_DEBUG_EXIT("MStruct::PowderPattern::Exit():"<<this->GetName(),5)
         return;
      }
      if("Radiation"==tag.GetName()) mRadiation.XMLInput(is,tag);
      if("MaxSinThetaOvLambda"==tag.GetName())
      {
         is>>mMaxSinThetaOvLambda;
         XMLCrystTag junk(is);
      }
      if("Par"==tag.GetName())
      {
         for(unsigned int i=0;i<tag.GetNbAttribute();i++)
         {
            if("Name"==tag.GetAttributeName(i))
            {
               if(("2ThetaZero"==tag.GetAttributeValue(i)) ||("Zero"==tag.GetAttributeValue(i)))
               {
                  this->GetPar(&mXZero).XMLInput(is,tag);
                  break;
               }
               if("2ThetaDisplacement"==tag.GetAttributeValue(i))
               {
                  this->GetPar(&m2ThetaDisplacement).XMLInput(is,tag);
                  break;
               }
               if("2ThetaTransparency"==tag.GetAttributeValue(i))
               {
                  this->GetPar(&m2ThetaTransparency).XMLInput(is,tag);
                  break;
               }
               if("TOF-DIFC"==tag.GetAttributeValue(i))
               {
                  this->GetPar(&mDIFC).XMLInput(is,tag);
                  break;
               }
               if("TOF-DIFA"==tag.GetAttributeValue(i))
               {
                  this->GetPar(&mDIFA).XMLInput(is,tag);
                  break;
               }
            }
         }
         continue;
      }
      if("Option"==tag.GetName())
      {
         for(unsigned int i=0;i<tag.GetNbAttribute();i++)
            if("Name"==tag.GetAttributeName(i))
               mOptionRegistry.GetObj(tag.GetAttributeValue(i)).XMLInput(is,tag);
         continue;
      }
      if("PowderPatternBackground"==tag.GetName())
      {
         PowderPatternBackground *comp=new PowderPatternBackground;
         comp->SetParentPowderPattern(*this);
         comp->XMLInput(is,tag);
         continue;
      }
      if("PowderPatternBackgroundInvX"==tag.GetName())
      {
         PowderPatternBackgroundInvX *comp=new PowderPatternBackgroundInvX;
         comp->SetParentPowderPattern(*this);
         comp->XMLInput(is,tag);
         continue;
      }
      if("PowderPatternBackgroundChebyshev"==tag.GetName())
      {
         PowderPatternBackgroundChebyshev *comp=new PowderPatternBackgroundChebyshev;
         comp->SetParentPowderPattern(*this);
         comp->XMLInput(is,tag);
         continue;
      }
      if("PowderPatternCrystal"==tag.GetName())
      {
         PowderPatternDiffraction *comp=new PowderPatternDiffraction;
         comp->SetParentPowderPattern(*this);
         comp->XMLInput(is,tag);
         continue;
      }
      if("PowderPatternComponent"==tag.GetName())
      {
         REAL scale=1.0;
         string name;
         for(unsigned int i=0;i<tag.GetNbAttribute();i++)
         {
            if("Scale"==tag.GetAttributeName(i))
            {
               stringstream ss(tag.GetAttributeValue(i));
               ss>>scale;
               continue;
            }
            if("Name"==tag.GetAttributeName(i)) name=tag.GetAttributeValue(i);
         }
         this->AddPowderPatternComponent(gPowderPatternComponentRegistry.GetObj(name));
         mScaleFactor(mPowderPatternComponentRegistry.GetNb()-1)=scale;
         VFN_DEBUG_MESSAGE("->Adding Component :"<<name<<"with scale="<<scale,8);
         continue;
      }
      if("ExcludeX"==tag.GetName())
      {
         float min,max;
         is>>min>>max;
         if(this->GetRadiation().GetWavelengthType()==WAVELENGTH_TOF)
            this->AddExcludedRegion(min,max);
         else this->AddExcludedRegion(min*DEG2RAD,max*DEG2RAD);
         XMLCrystTag end(is);
         continue;
      }
      if("IobsSigmaWeightList"==tag.GetName())
      {
         // Old version, just for 2theta-intensity pattern (no TOF)
         VFN_DEBUG_ENTRY("Loading Iobs-Sigma-Weight List...",8);
         REAL min,step;
         for(unsigned int i=0;i<tag.GetNbAttribute();i++)
         {
            if("TThetaMin"==tag.GetAttributeName(i))
            {
               stringstream ss(tag.GetAttributeValue(i));
               ss>>min;
               VFN_DEBUG_MESSAGE("2Theta min="<<min,8);
               min*=DEG2RAD;
            }
            if("TThetaStep"==tag.GetAttributeName(i))
            {
               stringstream ss(tag.GetAttributeValue(i));
               ss>>step;
               VFN_DEBUG_MESSAGE("2Theta step="<<step<<tag.GetAttributeValue(i),8);
               step*=DEG2RAD;
            }
         }
         while(0==isgraph(is.peek())) is.get();
         if(is.peek()=='<')
         {
            cout <<"MStruct::PowderPattern::XMLInput(): no data point in the powder pattern !"<<endl;
            XMLCrystTag junk(is);
            VFN_DEBUG_EXIT("Loading Iobs-Sigma-Weight List...",8);
            continue;
         }
         mNbPoint=0;
         mPowderPatternObs.resize(500);
         mPowderPatternObsSigma.resize(500);
         mPowderPatternWeight.resize(500);
         do
         {
            is >>mPowderPatternObs(mNbPoint)
               >>mPowderPatternObsSigma(mNbPoint)
               >>mPowderPatternWeight(mNbPoint);
            mNbPoint++;
            VFN_DEBUG_MESSAGE("Point #"<<mNbPoint,5);
            if(mNbPoint==(unsigned long)mPowderPatternObs.numElements())
            {
               mPowderPatternObs.resizeAndPreserve(mNbPoint+500);
               mPowderPatternObsSigma.resizeAndPreserve(mNbPoint+500);
               mPowderPatternWeight.resizeAndPreserve(mNbPoint+500);
            }
            while(0==isgraph(is.peek())) is.get();
         }
         while(is.peek()!='<');//until next tag
         this->SetPowderPatternPar(min,step,mNbPoint);
         mClockPowderPatternPar.Click();

         XMLCrystTag junk(is);
         VFN_DEBUG_EXIT("Loading Iobs-Sigma-Weight List...",8);
         continue;
      }
      if("XIobsSigmaWeightList"==tag.GetName())
      {
         VFN_DEBUG_ENTRY("Loading X-Iobs-Sigma-Weight List...",8);
         while(0==isgraph(is.peek())) is.get();
         if(is.peek()=='<')
         {
            cout <<"MStruct::PowderPattern::XMLInput(): no data point in the powder pattern !"<<endl;
            XMLCrystTag junk(is);
            VFN_DEBUG_EXIT("Loading Iobs-Sigma-Weight List...",8);
            continue;
         }
         mNbPoint=0;
         mX.resize(500);
         mPowderPatternObs.resize(500);
         mPowderPatternObsSigma.resize(500);
         mPowderPatternWeight.resize(500);
         do
         {
            is >>mX(mNbPoint)
               >>mPowderPatternObs(mNbPoint)
               >>mPowderPatternObsSigma(mNbPoint)
               >>mPowderPatternWeight(mNbPoint);
            mNbPoint++;
            VFN_DEBUG_MESSAGE("Point #"<<mNbPoint,5);
            if(mNbPoint==(unsigned long)mPowderPatternObs.numElements())
            {
               mX.resizeAndPreserve(mNbPoint+500);
               mPowderPatternObs.resizeAndPreserve(mNbPoint+500);
               mPowderPatternObsSigma.resizeAndPreserve(mNbPoint+500);
               mPowderPatternWeight.resizeAndPreserve(mNbPoint+500);
            }
            while(0==isgraph(is.peek())) is.get();
         }
         while(is.peek()!='<');//until next tag
         mX.resizeAndPreserve(mNbPoint);
         if(this->GetRadiation().GetWavelengthType()!=WAVELENGTH_TOF)
            mX*=DEG2RAD;
         this->SetPowderPatternX(mX);

         XMLCrystTag junk(is);
         VFN_DEBUG_EXIT("Loading X-Iobs-Sigma-Weight List...",8);
         continue;
      }
   }
}

////////////////////////////////////////////////////////////////////////
//
//    I/O AbsorptionCorr
//
////////////////////////////////////////////////////////////////////////

void AbsorptionCorr::XMLOutput(ostream &os,int indent)const
{
    VFN_DEBUG_ENTRY("MStruct::AbsorptionCorr::XMLOutput():"<<this->GetName(),11)
    for(int i=0; i<indent; i++) os << "  ";
    XMLCrystTag tag("AbsorptionCorr");
    tag.AddAttribute("Name", this->GetName());
	tag.AddAttribute("Thickness", (boost::format("%f")%(mThickness/10.)).str());
	tag.AddAttribute("Depth", (boost::format("%f")%(mDepth/10.)).str());
	tag.AddAttribute("AbsorptionFactor", (boost::format("%f")%(mAbsFactor*1e8)).str());
	tag.SetIsEmptyTag(true);
    os << tag << endl;
    
    VFN_DEBUG_EXIT("MStruct::AbsorptionCorr::XMLOutput():"<<this->GetName(),11)
}

////////////////////////////////////////////////////////////////////////
//
//    I/O PowderPatternBackgroundInvX
//
////////////////////////////////////////////////////////////////////////

void PowderPatternBackgroundInvX::XMLOutput(ostream &os,int indent)const
{
    VFN_DEBUG_ENTRY("PowderPatternBackgroundInvX::XMLOutput():"<<this->GetName(),11)
    for(int i=0; i<indent; i++) os << "  ";
    XMLCrystTag tag("PowderPatternBackgroundInvX");
    tag.AddAttribute("Name", this->GetName());
    os << tag << endl;
    indent++;

    mXFunctionType.XMLOutput(os,indent);
    os << endl;

    indent--;
    tag.SetIsEndTag(true);
    for(int i=0; i<indent; i++) os << "  ";
    os << tag << endl;
    VFN_DEBUG_EXIT("PowderPatternBackgroundInvX::XMLOutput():"<<this->GetName(),11)
}

void PowderPatternBackgroundInvX::XMLInput(istream &is,const XMLCrystTag &tagg)
{
	VFN_DEBUG_ENTRY("PowderPatternBackgroundInvX::XMLInput():"<<this->GetName(),11)
	for(unsigned int i=0;i<tagg.GetNbAttribute();i++)
	{
		if("Name"==tagg.GetAttributeName(i)) this->SetName(tagg.GetAttributeValue(i));
	}
	while(true)
	{
		XMLCrystTag tag(is);
		if(("PowderPatternBackgroundInvX"==tag.GetName())&&tag.IsEndTag())
		{
			this->UpdateDisplay();
			VFN_DEBUG_EXIT("PowderPatternBackgroundInvX::Exit():"<<this->GetName(),11)
			return;
		}
		if("Option"==tag.GetName())
		{
			for(unsigned int i=0;i<tag.GetNbAttribute();i++)
				if("Name"==tag.GetAttributeName(i))
					mOptionRegistry.GetObj(tag.GetAttributeValue(i)).XMLInput(is,tag);
			continue;
		}
	}
}

////////////////////////////////////////////////////////////////////////
//
//    I/O PowderPatternBackgroundChebyshev
//
////////////////////////////////////////////////////////////////////////

void PowderPatternBackgroundChebyshev::XMLOutput(ostream &os,int indent)const
{
    VFN_DEBUG_ENTRY("PowderPatternBackgroundChebyshev::XMLOutput():"<<this->GetName(),11)
    for(int i=0; i<indent; i++) os << "  ";
    XMLCrystTag tag("PowderPatternBackgroundChebyshev");
    tag.AddAttribute("Name", this->GetName());
    tag.AddAttribute("Polynomial_degree", std::to_string(mChebyshevCoef.numElements()-1));
    os << tag << endl;
    indent++;

    mXFunctionType.XMLOutput(os,indent);
    os << endl;
    
    // save coefficients of the Chebyshev polynomials
    const REAL* coefs = (const REAL*) mChebyshevCoef.data();
    for(int i=0; i<mChebyshevCoef.numElements(); i++) {
        ostringstream ss;
        ss << "Background_Coef_" << i;
        this->GetPar(&coefs[i]).XMLOutput(os,ss.str(),indent);
        os<<endl;
    }
    
    indent--;
    tag.SetIsEndTag(true);
    for(int i=0; i<indent; i++) os << "  ";
    os << tag << endl;
    VFN_DEBUG_EXIT("PowderPatternBackgroundChebyshev::XMLOutput():"<<this->GetName(),11)
}

void PowderPatternBackgroundChebyshev::XMLInput(istream &is,const XMLCrystTag &tagg)
{
	VFN_DEBUG_ENTRY("PowderPatternBackgroundChebyshev::XMLInput():"<<this->GetName(),11)
	for(unsigned int i=0;i<tagg.GetNbAttribute();i++)
	{
		if("Name"==tagg.GetAttributeName(i)) this->SetName(tagg.GetAttributeValue(i));
		if("Polynomial_degree"==tagg.GetAttributeName(i)) {
			CrystVector_REAL tcoefs(atoi(tagg.GetAttributeValue(i).c_str())+1);
			tcoefs = 0.;
			this->SetCoefficients(tcoefs);
		}
	}
	while(true)
	{
		XMLCrystTag tag(is);
		if(("PowderPatternBackgroundChebyshev"==tag.GetName())&&tag.IsEndTag())
		{
			this->UpdateDisplay();
			VFN_DEBUG_EXIT("PowderPatternBackgroundChebyshev::Exit():"<<this->GetName(),11)
			return;
		}
		if("Option"==tag.GetName())
		{
			for(unsigned int i=0;i<tag.GetNbAttribute();i++)
				if("Name"==tag.GetAttributeName(i))
					mOptionRegistry.GetObj(tag.GetAttributeValue(i)).XMLInput(is,tag);
			continue;
		}
		if("Par"==tag.GetName())
		{
			for(unsigned int i=0;i<tag.GetNbAttribute();i++)
			{
				if("Name"==tag.GetAttributeName(i))
				{
					this->GetPar(tag.GetAttributeValue(i)).XMLInput(is,tag);
					break;
				}
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////
//
//    I/O SizeBroadeningEffect
//
////////////////////////////////////////////////////////////////////////

void SizeBroadeningEffect::XMLOutput(ostream &os,int indent)const
{
    VFN_DEBUG_ENTRY("SizeBroadeningEffect::XMLOutput():"<<this->GetName(),11)
    for(int i=0; i<indent; i++) os << "  ";
    XMLCrystTag tag("SizeLn");
    tag.AddAttribute("Name", this->GetName());
    os << tag << endl;
    indent++;

	this->GetPar("M").XMLOutput(os,"M",indent);
	os << endl;
	this->GetPar("Sigma").XMLOutput(os,"Sigma",indent);
    os << endl;
	
    indent--;
    tag.SetIsEndTag(true);
    for(int i=0; i<indent; i++) os << "  ";
    os << tag << endl;
    VFN_DEBUG_EXIT("SizeBroadeningEffect::XMLOutput():"<<this->GetName(),11)
}

////////////////////////////////////////////////////////////////////////
//
//    I/O PseudoVoigtBroadeningEffectA
//
////////////////////////////////////////////////////////////////////////
void PseudoVoigtBroadeningEffectA::XMLOutput(ostream &os,int indent)const
{
    VFN_DEBUG_ENTRY("PseudoVoigtBroadeningEffectA::XMLOutput():"<<this->GetName(),11)
    for(int i=0; i<indent; i++) os << "  ";
    XMLCrystTag tag("pVoigtA");
    tag.AddAttribute("Name", this->GetName());
    os << tag << endl;
    indent++;

	this->GetPar("U").XMLOutput(os,"U",indent);
	os << endl;
	this->GetPar("V").XMLOutput(os,"V",indent);
    os << endl;
	this->GetPar("W").XMLOutput(os,"W",indent);
	os << endl;
	this->GetPar("Eta0").XMLOutput(os,"Eta0",indent);
    os << endl;
	this->GetPar("Eta1").XMLOutput(os,"Eta1",indent);
	os << endl;
	this->GetPar("Asym0").XMLOutput(os,"Asym0",indent);
    os << endl;
	this->GetPar("Asym1").XMLOutput(os,"Asym1",indent);
	os << endl;
	this->GetPar("Asym2").XMLOutput(os,"Asym2",indent);
    os << endl;
	
    indent--;
    tag.SetIsEndTag(true);
    for(int i=0; i<indent; i++) os << "  ";
    os << tag << endl;
    VFN_DEBUG_EXIT("PseudoVoigtBroadeningEffectA::XMLOutput():"<<this->GetName(),11)
}

////////////////////////////////////////////////////////////////////////
//
//    I/O RefractionPositionCorr
//
////////////////////////////////////////////////////////////////////////

void RefractionPositionCorr::XMLOutput(ostream &os,int indent)const
{
    VFN_DEBUG_ENTRY("RefractionPositionCorr::XMLOutput():"<<this->GetName(),11)
    for(int i=0; i<indent; i++) os << "  ";
    XMLCrystTag tag("RefractionCorr");
    tag.AddAttribute("Name", this->GetName());
    os << tag << endl;
    indent++;

	// --- options ---
    mChi0ValueChoice.XMLOutput(os,indent);
    os << endl;
	
	// --- refinable parameters ---
	this->GetPar("Density").XMLOutput(os,"relDensity",indent);
	os << endl;
	
    indent--;
    tag.SetIsEndTag(true);
    for(int i=0; i<indent; i++) os << "  ";
    os << tag << endl;
    VFN_DEBUG_EXIT("RefractionPositionCorr::XMLOutput():"<<this->GetName(),11)
}

////////////////////////////////////////////////////////////////////////
//
//    I/O XECsReussVoigt
//
////////////////////////////////////////////////////////////////////////
void XECsReussVoigt::XMLOutput(ostream &os,int indent)const
{
    VFN_DEBUG_ENTRY("XECsReussVoigt::XMLOutput():"<<this->GetName(),11)
    for(int i=0; i<indent; i++) os << "  ";
    XMLCrystTag tag("XECsReussVoigt");
    tag.AddAttribute("Name", this->GetName());
    os << tag << endl;
    indent++;
	
	// Stiffness constants
	const vector< string > & cij_names = this->GetStiffnessConstantsNames();
	for(int i=0; i<cij_names.size(); i++) {
		for(int i=0; i<indent; i++) os << "  ";
		XMLCrystTag tag("StiffnessConstant");
		tag.AddAttribute("Name", cij_names[i]);
		os << tag << mCijValues(i);
		tag.SetIsEndTag(true);
		os << tag << endl;
	}
	
	// model weight (0..Reuss,1..Voigt)
	this->GetPar("RV_weight").XMLOutput(os,"RV_weight",indent);
	os << endl;
	
    indent--;
    tag.SetIsEndTag(true);
    for(int i=0; i<indent; i++) os << "  ";
    os << tag << endl;
    VFN_DEBUG_EXIT("XECsReussVoigt::XMLOutput():"<<this->GetName(),11)
}

////////////////////////////////////////////////////////////////////////
//
//    I/O ResidualStressPositionCorrection
//
////////////////////////////////////////////////////////////////////////

void ResidualStressPositionCorrection::XMLOutput(ostream &os,int indent)const
{
    VFN_DEBUG_ENTRY("ResidualStressPositionCorrection::XMLOutput():"<<this->GetName(),11)
    for(int i=0; i<indent; i++) os << "  ";
    XMLCrystTag tag("StressSimple");
    tag.AddAttribute("Name", this->GetName());
    os << tag << endl;
    indent++;

	// XECs object
	if(pXECsObj != NULL)
		pXECsObj->XMLOutput(os,indent);
	
	// stress value
	this->GetPar(&mStress).XMLOutput(os,"Stress",indent);
	os << endl;
	
    indent--;
    tag.SetIsEndTag(true);
    for(int i=0; i<indent; i++) os << "  ";
    os << tag << endl;
    VFN_DEBUG_EXIT("ResidualStressPositionCorrection::XMLOutput():"<<this->GetName(),11)
}

////////////////////////////////////////////////////////////////////////
//
//    I/O ReflectionProfile
//
////////////////////////////////////////////////////////////////////////

void ReflectionProfile::XMLOutput(ostream &os,int indent)const
{
    VFN_DEBUG_ENTRY("MStruct::ReflectionProfile::XMLOutput():"<<this->GetName(),11)
    for(int i=0; i<indent; i++) os << "  ";
    XMLCrystTag tag("ReflectionProfile");
    tag.AddAttribute("Name", this->GetName());
    os << tag << endl;
    indent++;

	// XMLOutput for all components
	for(int icomp=0; icomp<this->GetReflectionProfileComponentNb(); icomp++) {
		this->GetReflectionProfileComponent(icomp).XMLOutput(os,indent);
		os << endl;
	}
    
    indent--;
    tag.SetIsEndTag(true);
    for(int i=0; i<indent; i++) os << "  ";
    os << tag << endl;
    VFN_DEBUG_EXIT("MStruct::ReflectionProfile::XMLOutput():"<<this->GetName(),11)
}

////////////////////////////////////////////////////////////////////////
//
//    Global functions
//
////////////////////////////////////////////////////////////////////////

void XMLCrystFileLoadAllObject(const string & filename)
{
	VFN_DEBUG_ENTRY("MStruct::XMLCrystFileLoadAllObject(filename,)",5)
	ifstream is(filename.c_str());
	if(is.fail()) throw ObjCrystException("XMLCrystFileLoadAllObject()   failed input");
	XMLCrystFileLoadAllObject(is);
	(*fpObjCrystInformUser)("Finished loading XML file:"+filename);
	VFN_DEBUG_EXIT("MStruct::XMLCrystFileLoadAllObject(filename,)",5)
}

void XMLCrystFileLoadAllObject(istream &is)
{
	VFN_DEBUG_ENTRY("MStruct::XMLCrystFileLoadAllObject(istream)",5)
	is.imbue(std::locale::classic());
	XMLCrystTag tag;
	do {is>>tag;} while(("ObjCryst"!=tag.GetName()) && (false==is.eof()));

	while(true)
	{
		XMLCrystTag tag(is);
		if(true==is.eof()) break;
		if(tag.GetName()=="Crystal")
		{
			Crystal* obj = new Crystal;
			obj->XMLInput(is,tag);
			(*fpObjCrystInformUser)("XML: finished reading Crystal object:"+obj->GetName());
		}
		if(tag.GetName()=="PowderPattern")
		{
			MStruct::PowderPattern* obj = new MStruct::PowderPattern;
			obj->XMLInput(is,tag);
			(*fpObjCrystInformUser)("XML: finished reading Powder Pattern object:"+obj->GetName());
		}
		if(tag.GetName()=="DiffractionDataSingleCrystal")
		{
			DiffractionDataSingleCrystal* obj = new DiffractionDataSingleCrystal;
			obj->XMLInput(is,tag);
			(*fpObjCrystInformUser)("XML: finished reading Single Crystal Diffraction object:"+obj->GetName());
		}
		if(tag.GetName()=="GlobalOptimObj")
		{
			MonteCarloObj* obj = new MonteCarloObj;
			obj->XMLInput(is,tag);
			(*fpObjCrystInformUser)("XML: finished reading Global Optimization object:"+obj->GetName());
		}
	}
	(*fpObjCrystInformUser)("Finished loading XML");
	VFN_DEBUG_EXIT("MStruct::XMLCrystFileLoadAllObject(istream)",5)
}

} // namespace MStruct
/* ------------------------------------------------------------------------------------------------ */

//void XMLOutput(ostream &os, int indent=0) const;
//void XMLInput(istream &is,const ObjCryst::XMLCrystTag &tag);