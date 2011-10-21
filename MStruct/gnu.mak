BUILD_DIR=$(CURDIR)/..
include $(BUILD_DIR)/ObjCryst/rules.mak

OBJ= $@.o

%.o : %.cpp
	@$(MAKEDEPEND)
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -c $< -o $@

-include $(OBJ:.o=.dep)

mstruct_am: mstruct_am.o MStruct.o mstruct_test1.o libCrystVector libQuirks libRefinableObj libcctbx libCryst libnewmat
	${LINKER} ${CRYST_LDFLAGS} -o $@ ${filter-out %.a %.so lib%, $^} ${LOADLIBES} -lfftw3 -lblas -llapack

# target for making everything
.PHONY : all
all: mstruct_am

default: all

#target to make documentation (requires doxygen)
#also makes tags file, although it is not related to doxygen
.PHONY : doc
doc::
	cd $(BUILD_DIR)/MStruct/doc; $(MAKE) -f gnu.mak doc

.PHONY : tidy
tidy::
	@${RM} core *.o *.obj *.dep

.PHONY : clean
clean:: tidy
	@${RM} core *.a *.exe mstruct_am

.PHONY : clean_all
clean_all::
	$(MAKE) -f gnu.mak -C ${DIR_CRYST} clean
	$(MAKE) -f gnu.mak -C ${DIR_CRYST}/doc clean
	$(MAKE) -f gnu.mak -C $(BUILD_DIR)/MStruct clean
	$(MAKE) -f gnu.mak -C $(BUILD_DIR)/MStruct/doc clean
	$(MAKE) -f gnu.mak -C ${DIR_CRYST}/../cctbx clean
	cd $(BUILD_DIR)/newmat; ${RM} *.o *.a
	cd $(BUILD_DIR)/static-libs/lib; ${RM} *.o *.a

cvsignore:
	cp -f ${DIR_CRYST}/.cvsignore ./

