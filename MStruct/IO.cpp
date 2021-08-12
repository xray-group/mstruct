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

using namespace ObjCryst;

namespace MStruct {


////////////////////////////////////////////////////////////////////////
//
//    PowderPatternBackgroundInvX
//
////////////////////////////////////////////////////////////////////////

void PowderPatternBackgroundInvX::XMLOutput(ostream &os,int indent)const
{
    VFN_DEBUG_ENTRY("PowderPatternBackgroundInvX::XMLOutput():Begin"<<this->GetName(),11)
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
    VFN_DEBUG_ENTRY("PowderPatternBackgroundInvX::XMLOutput():Begin"<<this->GetName(),11)
}

////////////////////////////////////////////////////////////////////////
//
//    PowderPatternBackgroundChebyshev
//
////////////////////////////////////////////////////////////////////////

void PowderPatternBackgroundChebyshev::XMLOutput(ostream &os,int indent)const
{
    VFN_DEBUG_ENTRY("PowderPatternBackgroundChebyshev::XMLOutput():Begin"<<this->GetName(),11)
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
    VFN_DEBUG_ENTRY("PowderPatternBackgroundChebyshev::XMLOutput():Begin"<<this->GetName(),11)
}

////////////////////////////////////////////////////////////////////////
//
//    SizeBroadeningEffect
//
////////////////////////////////////////////////////////////////////////

void SizeBroadeningEffect::XMLOutput(ostream &os,int indent)const
{
    VFN_DEBUG_ENTRY("SizeBroadeningEffect::XMLOutput():Begin"<<this->GetName(),11)
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
    VFN_DEBUG_ENTRY("SizeBroadeningEffect::XMLOutput():Begin"<<this->GetName(),11)
}

////////////////////////////////////////////////////////////////////////
//
//    PseudoVoigtBroadeningEffectA
//
////////////////////////////////////////////////////////////////////////
void PseudoVoigtBroadeningEffectA::XMLOutput(ostream &os,int indent)const
{
    VFN_DEBUG_ENTRY("PseudoVoigtBroadeningEffectA::XMLOutput():Begin"<<this->GetName(),11)
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
    VFN_DEBUG_ENTRY("PseudoVoigtBroadeningEffectA::XMLOutput():Begin"<<this->GetName(),11)
}

////////////////////////////////////////////////////////////////////////
//
//    RefractionPositionCorr
//
////////////////////////////////////////////////////////////////////////

void RefractionPositionCorr::XMLOutput(ostream &os,int indent)const
{
    VFN_DEBUG_ENTRY("RefractionPositionCorr::XMLOutput():Begin"<<this->GetName(),11)
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
    VFN_DEBUG_ENTRY("RefractionPositionCorr::XMLOutput():Begin"<<this->GetName(),11)
}

////////////////////////////////////////////////////////////////////////
//
//    XECsReussVoigt
//
////////////////////////////////////////////////////////////////////////
void XECsReussVoigt::XMLOutput(ostream &os,int indent)const
{
    VFN_DEBUG_ENTRY("XECsReussVoigt::XMLOutput():Begin"<<this->GetName(),11)
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
    VFN_DEBUG_ENTRY("XECsReussVoigt::XMLOutput():Begin"<<this->GetName(),11)
}

////////////////////////////////////////////////////////////////////////
//
//    ResidualStressPositionCorrection
//
////////////////////////////////////////////////////////////////////////

void ResidualStressPositionCorrection::XMLOutput(ostream &os,int indent)const
{
    VFN_DEBUG_ENTRY("ResidualStressPositionCorrection::XMLOutput():Begin"<<this->GetName(),11)
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
    VFN_DEBUG_ENTRY("ResidualStressPositionCorrection::XMLOutput():Begin"<<this->GetName(),11)
}

////////////////////////////////////////////////////////////////////////
//
//    ReflectionProfile
//
////////////////////////////////////////////////////////////////////////

void ReflectionProfile::XMLOutput(ostream &os,int indent)const
{
    VFN_DEBUG_ENTRY("MStruct::ReflectionProfile::XMLOutput():Begin"<<this->GetName(),11)
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
    VFN_DEBUG_ENTRY("MStruct::ReflectionProfile::XMLOutput():Begin"<<this->GetName(),11)
}

} // namespace MStruct
/* ------------------------------------------------------------------------------------------------ */

//void XMLOutput(ostream &os, int indent=0) const;