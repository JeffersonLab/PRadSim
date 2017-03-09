# Makefile for the MOLLER program
# (c) A. V. Gramolin, April 2014


include Makefile.arch
-include ../MyConfig.mk


MOLLERO       = moller.$(ObjSuf)
MOLLERS       = moller.$(SrcSuf)
MOLLER        = moller$(ExeSuf)


OBJS          = $(MOLLER)

PROGRAMS      = $(MOLLER)


.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)


$(MOLLER):      $(MOLLERO)
		$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@
		$(MT_EXE)
		@echo "$@ done"


.$(SrcSuf).$(ObjSuf):
	$(CXX)  $(CXXFLAGS) -c $<

clean:
	rm $(MOLLERO) $(MOLLER)
