OBJCRYSTLIB_DIR=./../build/fast-x86_64
CXXFLAGS=-D__ZDENEK__ -I$(OBJCRYSTLIB_DIR) -I$(OBJCRYSTLIB_DIR)/ObjCryst
LDFLAGS=-L$(OBJCRYSTLIB_DIR)
LIBS=-lObjCryst -llapack

testlib: testlib.cpp
	g++ $(CXXFLAGS) $(LDFLAGS) -o testlib testlib.cpp $(LIBS)

clean:
	rm -f testlib *.o

