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

} // namespace MStruct
/* ------------------------------------------------------------------------------------------------ */
