OBJCRYSTLIB_DIR=./../../libobjcryst/build/fast-x86_64
MSTRUCTLIB_DIR=./../../libmstruct/build/fast-x86_64
CXXFLAGS=-D__ZDENEK__ -DREAL=double -I$(OBJCRYSTLIB_DIR) -I$(OBJCRYSTLIB_DIR)/ObjCryst -I$(MSTRUCTLIB_DIR) -I$(MSTRUCTLIB_DIR)/MStruct
LDFLAGS=-L$(OBJCRYSTLIB_DIR) -L$(MSTRUCTLIB_DIR)
LIBS=-lObjCryst -lMStruct -llapack -lfftw3 -lgsl

test7: test7.cpp
	g++ $(CXXFLAGS) $(LDFLAGS) -o test7 test7.cpp $(LIBS)

testlib: testlib.cpp
	g++ $(CXXFLAGS) $(LDFLAGS) -o testlib testlib.cpp $(LIBS)

clean:
	rm -f testlib test7 *.o

