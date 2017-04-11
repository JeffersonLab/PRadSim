# Makefile for the ESEPP event generator

include Makefile.arch

EPNORCO     = epnorc.$(ObjSuf)
EPNORCS     = norc.cxx ep.h
EPNORC      = epnorc
EDNORCO     = ednorc.$(ObjSuf)
EDNORCS     = norc.cxx ed.h
EDNORC      = ednorc

OBJS        = $(EPNORCO) $(EDNORCO)

PROGRAMS    = $(EPNORC) $(EDNORC)

.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)

all: $(PROGRAMS)

$(EPNORC): $(EPNORCO)
	$(LD) $(LDFLAGS) $(LIBS) -o $@ $<
	@echo "$@ done"

$(EPNORCO): $(EPNORCS)
	$(CXX) -c $(CXXFLAGS) -DELASTIC_EP -o $@ $<

$(EDNORC): $(EDNORCO)
	$(LD) $(LDFLAGS) $(LIBS) -o $@ $<
	@echo "$@ done"

$(EDNORCO): $(EDNORCS)
	$(CXX) -c $(CXXFLAGS) -DELASTIC_ED -o $@ $<

clean:
	rm -rf $(EPNORCO) $(EPNORC) $(EDNORCO) $(EDNORC)
