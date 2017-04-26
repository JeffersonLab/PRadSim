# Makefile for the ESEPP event generator

include Makefile.arch

EENORCO     = eenorc.$(ObjSuf)
EENORCS     = norc.cc ee.hh
EENORC      = eenorc
EPNORCO     = epnorc.$(ObjSuf)
EPNORCS     = norc.cc ep.hh
EPNORC      = epnorc
EDNORCO     = ednorc.$(ObjSuf)
EDNORCS     = norc.cc ed.hh
EDNORC      = ednorc

OBJS        = $(EENORCO) $(EPNORCO) $(EDNORCO)

PROGRAMS    = $(EENORC) $(EPNORC) $(EDNORC)

.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)

all: $(PROGRAMS)

$(EENORC): $(EENORCO)
	$(LD) $(LDFLAGS) $(LIBS) -o $@ $<
	@echo "$@ done"

$(EENORCO): $(EENORCS)
	$(CXX) -c $(CXXFLAGS) -DELASTIC_EE -o $@ $<

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
	rm -rf $(EENORCO) $(EENORC) $(EPNORCO) $(EPNORC) $(EDNORCO) $(EDNORC)
