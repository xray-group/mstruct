/* 
 * IO.h
 * 
 * MStruct++ - Object-Oriented computer program/library for MicroStructure analysis
 * 					   from powder diffraction data.
 * 
 * Copyright (C) 2009-2014  Zdenek Matej, Charles University in Prague
 * Copyright (C) 2014-2024  Zdenek Matej, MAX IV Laboratory, Lund University
 * Copyright (C) 2016-2019  Milan Dopita, Jan Endres, Charles University in Prague
 * Copyright (C) 2017-2018  Jiri Wollman, Charles University in Prague
 * 
 * Special credit to Vincent Favre Nicolin & ObjCryst++ for a part of I/O code.
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
 * along with MStruct++. If not, see <https://www.gnu.org/licenses/>.
 * 
 */
 
#ifndef __IO_MSTRUCT__H__
#define __IO_MSTRUCT__H__

#include <iostream>
#include <string>

namespace MStruct {

/** \brief Load all 'top' objects from a file (Crystal, PowderPattern, DiffDataSingleCrystal
*  and GlobalOptimObj objects). All objects are directly allocated, and can be accessed through
* their respective global registry (eg gCrystalRegistry fro a Crysta, etc...)
*
* \param file: the filename from which the objects will be loaded.
*/
void XMLCrystFileLoadAllObject(const std::string & file);
/** \brief Load all 'top' objects from a file (Crystal, PowderPattern, DiffDataSingleCrystal
*  and GlobalOptimObj objects). All objects are directly allocated, and can be accessed through
* their respective global registry (eg gCrystalRegistry fro a Crysta, etc...)
*
* \param file: the filename from which the objects will be loaded.
*/
void XMLCrystFileLoadAllObject(std::istream &is);

} // namespace MStruct

#endif /* __IO_MSTRUCT__H__ */
