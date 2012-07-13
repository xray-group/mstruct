BUILD_DIR=$(CURDIR)/../..
include $(BUILD_DIR)/ObjCryst/rules.mak

doc: ../*.h ${DIR_CRYST}/*/*.h
	cd ${DIR_DOC}; $(MAKE) -f gnu.mak doc
	doxygen Doxyfile

clean:
	@${RM} -fr html/* latex/* *.tag
	
