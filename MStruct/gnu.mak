BUILD_DIR=$(CURDIR)/..
include $(BUILD_DIR)/ObjCryst/rules.mak

OBJ= $@.o

%.o : %.cpp
	@$(MAKEDEPEND)
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -c $< -o $@

-include $(OBJ:.o=.dep)

mstruct_am: mstruct_am.o MStruct.o mstruct_test1.o libCrystVector libQuirks libRefinableObj libcctbx libCryst libnewmat
	${LINKER} ${CRYST_LDFLAGS} -o $@ ${filter-out %.a %.so lib%, $^} ${LOADLIBES}

# target for making everything
.PHONY : all
all: mstruct_am

default: all

.PHONY : tidy
tidy::
	@${RM} core *.o *.obj *.dep

.PHONY : clean
clean:: tidy
	@${RM} core *.a *.exe mstruct_am

cvsignore:
	cp -f ${DIR_CRYST}/.cvsignore ./

