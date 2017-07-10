# Makefile for the EDS program
# (c) A. V. Gramolin, June 2017


include Makefile.arch
-include ../MyConfig.mk


EDSO       = eds.$(ObjSuf)
EDSS       = eds.$(SrcSuf)
EDS        = eds$(ExeSuf)


OBJS          = $(EDS)

PROGRAMS      = $(EDS)


.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)


$(EDS):         $(EDSO)
		$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@
		$(MT_EXE)
		@echo "$@ done"


.$(SrcSuf).$(ObjSuf):
	$(CXX)  $(CXXFLAGS) -c $<
