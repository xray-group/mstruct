/* 
 * mstruct_tests.h
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
 
#ifndef __MSTRUCT_TESTS__
#define __MSTRUCT_TESTS__

#include <iostream>
#include <sstream>
#include <vector>
#include <functional>

namespace MStruct {

int mstruct_test1(int argc, char* argv[], std::istream &ccin);
int mstruct_test2(int argc, char* argv[], std::istream &ccin);
int mstruct_test3(int argc, char* argv[], std::istream &ccin);
int mstruct_test4(int argc, char* argv[], std::istream &ccin);
int mstruct_test5(int argc, char* argv[], std::istream &is);
int mstruct_test6(int argc, char* argv[], std::istream &is);
int mstruct_test7(int argc, char* argv[], std::istream &is);
int mstruct_test8(int argc, char* argv[], std::istream &is);
int mstruct_test9(int argc, char* argv[], std::istream &is);
int mstruct_test10(int argc, char* argv[], std::istream &is);
int mstruct_test11(int argc, char* argv[], std::istream &is);
int mstruct_test12(int argc, char* argv[], std::istream &is);

} // namespace MStruct

// unary function to convert char to lowercase
struct char2lower : public std::unary_function<char, void> {
	void operator() (char& c) {
       c = tolower(c); // C - function
	}
};

// unary function to convert char to lowercase
struct char2upper : public std::unary_function<char, void> {
	void operator() (char& c) {
       c = toupper(c); // C - function
	}
};

// The function reads all empty and commented (//) lines from the input
// stream 'is' until an uncommented and non-empty line is found. This last read line is returned
// as content of the 'iss' input-string-stream which can be used for a standard input.
// If the 'load_objects' parameter is true, than all found 'objects' are loaded and poiters to
// them are returned in the output vector.
std::vector< void* > read_line ( std::istringstream &iss, std::istream &is = std::cin, const bool load_objects=true);


#endif /* __MSTRUCT_TESTS__ */
