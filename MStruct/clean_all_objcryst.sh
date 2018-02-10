#!/bin/bash

# absolute path to the MStruct directory
mstruct_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

answer='';
echo -n "Clean all ObjCryst + MStruct project? (y/n): ";
read answer;
if [ $answer = 'y' ] || [ $answer = 'Y' ];
then

	echo "Deleting ObjCryst + MStruct project";
	
	cd ${mstruct_dir}; pwd; make -f gnu.mak clean;
	cd ${mstruct_dir}/doc; pwd; make -f gnu.mak clean;
	cd ${mstruct_dir}/../cctbx; pwd; make -f gnu.mak clean;
	cd ${mstruct_dir}/../newmat; pwd; rm -f *.o *.a;
	cd ${mstruct_dir}/../static-libs/lib; rm -f *.a;
	cd ${mstruct_dir}/../ObjCryst; pwd; make -f gnu.mak clean;
	cd ${mstruct_dir}/../ObjCryst/doc; pwd; make -f gnu.mak clean;
	
else
	echo "Nothing done.";
fi

answer='';
echo -n "Build MStruct project (mstruct_am)? (y/n): ";
read answer;
if [ $answer = 'y' ] || [ $answer = 'Y' ];
then

	echo "Building ObjCryst + MStruct project";
	
	cd ${mstruct_dir}; pwd;
	make -f gnu.mak wxcryst=0 opengl=0 fftw=0 shared-wxgtk=1 shared-glut=1 shared-newmat=0 shared-fftw=0 unicode=0 shared=1 clapack=1 debug=0 zdenek=1 mstruct_am
	
fi

answer='';
echo -n "Generate Documentation for ObjCryst + MStruct project? (y/n): ";
read answer;
if [ $answer = 'y' ] || [ $answer = 'Y' ];
then

	echo "Generating Documentation for ObjCryst + MStruct project";
	
	cd ${mstruct_dir}; pwd;
	
	make -f gnu.mak wxcryst=0 opengl=0 fftw=0 shared-wxgtk=1 shared-glut=1 shared-newmat=0 shared-fftw=0 unicode=0 shared=1 clapack=1 debug=0 zdenek=1 doc

fi

exit 0;
