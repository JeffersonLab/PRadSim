# Makefile for the ESEPP event generator

include Makefile.arch
-include ../MyConfig.mk


ESEPPO       = esepp.$(ObjSuf)
ESEPPS       = esepp.$(SrcSuf)
ESEPP        = esepp$(ExeSuf)


OBJS          = $(ESEPP)

PROGRAMS      = $(ESEPP)


.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)


$(ESEPP): $(ESEPPO)
	$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@
	$(MT_EXE)
	@echo "$@ done"


.$(SrcSuf).$(ObjSuf):
	$(CXX)  $(CXXFLAGS) -c $<

clean:
	rm $(ESEPPO) $(ESEPP)
