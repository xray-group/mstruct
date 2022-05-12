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
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * MStruct++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MStruct++. If not, see <https://www.gnu.org/licenses/>.
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

#include <boost/program_options.hpp>
namespace po = boost::program_options;

//#include <fstream>
#include <iostream>
#include <iomanip>
//#include <cstring>

using std::cout;
using std::cerr;
using std::flush;
using std::string;

int main (int argc, char *argv[])
{
  string input_file;
  string output_file("xray_out.xml");
  string output_dat("pattern0_xml.dat");
  int niter = -1;
  bool fit_scale_factor = false;
  
  try {
    
    po::options_description desc("Allowed options");
    desc.add_options()
      ("help,h", "print help message")
      ("version,v", "print version message")
      ("input,i", po::value<string>(), "input file")
      ("output,o", po::value<string>(&output_file), "output file (xray_out.xml)")
      ("output-data,O", po::value<string>(&output_dat), "output data file (pattern0_xml.dat)")
      ("debug-level", po::value<int>(), "debug level")
      ("niteraction,n", po::value<int>(&niter), "number of refinement iteractions")
      ("fit-scale-factor,n", po::value<bool>(&fit_scale_factor), "fit scale factor minimising Rw (before refinement)")
      ;

    po::positional_options_description p;
    p.add("input", 1);
    p.add("output", 1);
    
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
	      options(desc).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help")) {
      cout << "Usage: mstruct_xml input.xml\n";
      cout << desc << "\n";
      return 0;
    }

    if (vm.count("debug-level")) {
      int level = vm["debug-level"].as<int>();
      cout << "Setting debug level: " << level << "\n";
      VFN_DEBUG_GLOBAL_LEVEL(level);
    }

    if (vm.count("version")) {
      // print version and license information
      cout << "version: " << mstruct_version_str << "\n";
      cout << "mstruct Copyright\n";
      cout << "(C) 2009-2018 Charles University in Prague\n";
      cout << "(C) 2014-2022 Zdenek Matej, MAX IV Laboratory, Lund University\n";
      cout << "e-mail: <zdenek.matej@maxiv.lu.se>\n";
      cout << "web: <https://xray.cz/mstruct/>\n";
      cout << "License GNU GPL: <https://www.gnu.org/licenses/gpl-2.0.html>.\n";
      cout << "This program comes with ABSOLUTELY NO WARRANTY;\n";
      cout << "This is free software, and you are welcome to redistribute it.\n";
      cout << flush;
      return 0;
    }

    if (vm.count("input")) {
      input_file = vm["input"].as<string>();
    } else {
      cout << "Not enough arguments.\n";
      cout << "Usage: mstruct_xml sample.xml\n";
      cout << flush;                                                                                                                                                                                                                                                         
      return 1;
    }

    if (vm.count("niterations") && (niter>0) && (vm.count("fit-scale-factor")==0)) {
      fit_scale_factor = true; // switch default if refinement enabled
    }
      
  }
  catch(std::exception& e) {
    cerr << "error: " << e.what() << "\n";
    return 1;
  }
  catch(...) {
    cerr << "Exception of unknown type!\n";
  }
	   
  MStruct::XMLCrystFileLoadAllObject(input_file.c_str());

  // Get PowderPattern object
  ObjCryst::RefinableObj &obj = ObjCryst::gRefinableObjRegistry.GetObj("pattern0","MStruct::PowderPattern");
  MStruct::PowderPattern &data = dynamic_cast<MStruct::PowderPattern&>(obj);

  // Prepare data
  data.Prepare();
  if (fit_scale_factor)
    data.FitScaleFactorForRw();

  data.SavePowderPattern(output_dat.c_str());

  cout << "Output file: " << output_file << "\n";
  ObjCryst::XMLCrystFileSaveGlobal(output_file.c_str());

  return 0;
}
