OBJCRYSTLIB_DIR=./../../libobjcryst/build/fast-x86_64
MSTRUCTLIB_DIR=./../../libmstruct/build/fast-x86_64
CXXFLAGS=-D__ZDENEK__ -I$(OBJCRYSTLIB_DIR) -I$(OBJCRYSTLIB_DIR)/ObjCryst -I$(OBJCRYSTLIB_DIR)/MStruct
LDFLAGS=-L$(OBJCRYSTLIB_DIR) -L$(MSTRUCTLIB_DIR)
LIBS=-lObjCryst -lMStruct -llapack -lfftw3 -lgsl

testlib: testlib.cpp
	g++ $(CXXFLAGS) $(LDFLAGS) -o testlib testlib.cpp $(LIBS)

test7: test7.cpp
	g++ $(CXXFLAGS) $(LDFLAGS) -o test7 test7.cpp $(LIBS)

clean:
	rm -f testlib test7 *.o

