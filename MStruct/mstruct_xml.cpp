/* 
 * mstruct_xml.cpp
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
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MStruct++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MStruct++. If not, see <http://www.gnu.org/licenses/>.
 * 
 */
 
/*
 *  MStruct refinement program - XML API.
 *
 */

#include "config.h"
#include "MStruct.h"
#include "IO.h"
#include "ObjCryst/IO.h"

//#include <fstream>
//#include <iomanip>
//#include <cstring>

using std::cout;
using std::flush;

int main (int argc, char *argv[])
{
  if(argc>=2 && (string(argv[1])==string("-v") || string(argv[1])==string("--version"))) // print version
  {
    // print version and license information
    cout << "version: " << mstruct_version_str << "\n";
    cout << "mstruct Copyright\n";
    cout << "(C) 2009-2018 Charles University in Prague\n";
    cout << "(C) 2014-2022 Zdenek Matej, MAX IV Laboratory, Lund University\n";
    cout << "e-mail: <zdenek.matej@maxiv.lu.se>\n";
    cout << "web: <https://xray.cz/mstruct/>\n";
    cout << "License GNU GPL: <https://gnu.org/licenses/gpl.html>.\n";
    cout << "This program comes with ABSOLUTELY NO WARRANTY;\n";
    cout << "This is free software, and you are welcome to redistribute it.\n";
    cout << flush;
    return 0;
  }

  if(argc<2) { // not enough arguments: print usage
    cout << "Not enough arguments.\n";
    cout << "Usage: mstruct_xml sample.xml\n";
    cout << flush;
    return 0;
  }   

  int level=12;
  if(argc>=3) { // debug level hase been supplied
    level=atoi(argv[2]);
  }
  cout << "Setting debug level: " << level << "\n";
  VFN_DEBUG_GLOBAL_LEVEL(level);
	   
  MStruct::XMLCrystFileLoadAllObject(argv[1]);

  // Get PowderPattern object
  ObjCryst::RefinableObj &obj = ObjCryst::gRefinableObjRegistry.GetObj("pattern0","MStruct::PowderPattern");
  MStruct::PowderPattern &data = dynamic_cast<MStruct::PowderPattern&>(obj);

  // Prepare data
  data.Prepare();
  //data.FitScaleFactorForRw();

  data.SavePowderPattern("pattern0_xml.dat");

  ObjCryst::XMLCrystFileSaveGlobal("xray_out.xml");

  return 0;
}
