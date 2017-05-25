########################################################################
#   This Makefile shows how to compile all C++, C and Fortran
#   files found in $(SRCDIR) directory.
#   Linking is done with g++. Need to have $ROOTSYS defined
########################################################################

########################################################################
MYOS        := $(shell uname)
ARCH        := $(shell uname -m)
USER        := $(shell whoami)
MYHOST      := $(shell hostname -s)

########################################################################
EXECFILE    := moller

########################################################################
SRCDIR      := .
INCDIR      := .
OBJDIR      := obj.$(ARCH)

########################################################################
# Compiler
AR          := ar
CXX         := g++
FF          := gfortran
LD          := g++

########################################################################
# Flags
ifeq ($(ARCH),i686)
    MODE    := -m32
else
    MODE    := -m64
endif
INCDIRS     := $(patsubst %,-I%,$(subst :, ,$(INCDIR)))
CFLAGS      := -Wall -fPIC -O3 -g $(MODE) $(INCDIRS)
CXXFLAGS    := -Wall -fPIC -O3 -g $(MODE) $(INCDIRS)
FFLAGS      := -Wall -fPIC -O3 -g $(MODE) $(INCDIRS)
CXXFLAGS    += -std=c++11
FFLAGS      += -cpp -ffixed-line-length-none
GPPFLAGS    := -MM
LDFLAGS     := -O3 -g $(MODE)

########################################################################
# Generate obj file lis
FSOURCES    := $(wildcard $(SRCDIR)/*.[Ff])
CSOURCES    += $(wildcard $(SRCDIR)/*.[Cc][Cc])
CSOURCES    += $(wildcard $(SRCDIR)/*.[Cc][XxPp][XxPp])
SOURCES     := $(FSOURCES) $(CSOURCES)
# header files
HEADERS     := $(foreach n,$(subst :, ,$(INCDIR)),$(wildcard $(n)/*.hh))
# add .o to all the source files
OBJS        := $(addsuffix .o, $(basename $(SOURCES)))
OBJS        := $(patsubst  $(SRCDIR)/%.o,$(OBJDIR)/%.o,$(OBJS))
DEPS        := $(subst .o,.d,$(OBJS))

########################################################################
# Libs
SYSLIBS     := -lstdc++ -lgfortran

########################################################################
# ROOT configure
ROOTCFLAGS  := $(shell root-config --cflags)
ROOTLIBS    := $(shell root-config --libs) -lMathMore -lFoam
ROOTGLIBS   := $(shell root-config --glibs) -lMathMore -lFoam

CXXFLAGS    += $(ROOTCFLAGS)
LIBS        := $(SYSLIBS) $(ROOTLIBS)
GLIBS       := $(SYSLIBS) $(ROOTGLIBS)

########################################################################
# You can specify the .SUFFIXES
.SUFFIXES: .c .C .cc .CC .cpp .cxx .f .F
.PHONY: all clean tes
VPATH       := $(SRCDIR)

########################################################################
all: exe

########################################################################
# Make the $(TARGET).d file and include it.
$(OBJDIR)/%.d: %.c
	@echo Making dependency for file $< ......
	@set -e;
	$(CXX) $(GPPFLAGS) $(CXXFLAGS) $< | \
	sed 's!$*\.o!$(OBJDIR)/& $@!' > $@; \
	[ -s $@ ] || rm -f $@

$(OBJDIR)/%.d: %.C
	@echo Making dependency for file $< ......
	@set -e;
	$(CXX) $(GPPFLAGS) $(CXXFLAGS) $< | \
	sed 's!$*\.o!$(OBJDIR)/& $@!' > $@; \
	[ -s $@ ] || rm -f $@

$(OBJDIR)/%.d: %.cc
	@echo Making dependency for file $< ......
	@set -e;
	$(CXX) $(GPPFLAGS) $(CXXFLAGS) $< | \
	sed 's!$*\.o!$(OBJDIR)/& $@!' > $@; \
	[ -s $@ ] || rm -f $@

$(OBJDIR)/%.d: %.CC
	@echo Making dependency for file $< ......
	@set -e;
	$(CXX) $(GPPFLAGS) $(CXXFLAGS) $< | \
	sed 's!$*\.o!$(OBJDIR)/& $@!' > $@; \
	[ -s $@ ] || rm -f $@

$(OBJDIR)/%.d: %.cpp
	@echo Making dependency for file $< ......
	@set -e;
	$(CXX) $(GPPFLAGS) $(CXXFLAGS) $< |
	sed 's!$*\.o!$(OBJDIR)/& $@!' > $@;
	[ -s $@ ] || rm -f $@

$(OBJDIR)/%.d: %.cxx
	@echo Making dependency for file $< ......
	@set -e;
	$(CXX) $(GPPFLAGS) $(CXXFLAGS) $< | \
	sed 's!$*\.o!$(OBJDIR)/& $@!' > $@; \
	[ -s $@ ] || rm -f $@

#$(OBJDIR)/%.d: %.f
#	@echo Making dependency for file $< ......
#	@set -e;
#	$(FF) -cpp $(GPPFLAGS) $(FFLAGS) $< | \
#	sed 's!$*\.o!$(OBJDIR)/& $@!' > $@; \
#	[ -s $@ ] || rm -f $@

#$(OBJDIR)/%.d: %.F
#	@echo Making dependency for file $< ......
#	@set -e;
#	$(FF) -cpp $(GPPFLAGS) $(FFLAGS) $< | \
#	sed 's!$*\.o!$(OBJDIR)/& $@!' > $@; \
#	[ -s $@ ] || rm -f $@

ifneq ($(DEPS),)
-include $(DEPS)
endif

########################################################################
exe: dir $(OBJS)
	@$(LD) $(LDFLAGS) -o $(EXECFILE) $(OBJS) $(LIBS) $(OTHERLIBS)
	@echo "Linking $(EXECFILE) ...... done!"

########################################################################
$(OBJDIR)/%.o: %.c
	@echo Compiling $< ......
	@$(CXX) -c $< -o $@  $(CXXFLAGS)

$(OBJDIR)/%.o: %.C
	@echo Compiling $< ......
	@$(CXX) -c $< -o $@  $(CXXFLAGS)

$(OBJDIR)/%.o: %.cc
	@echo Compiling $< ......
	@$(CXX) -c $< -o $@  $(CXXFLAGS)

$(OBJDIR)/%.o: %.CC
	@echo Compiling $< ......
	@$(CXX) -c $< -o $@  $(CXXFLAGS)

$(OBJDIR)/%.o: %.cpp
	@echo Compiling $< ......
	@$(CXX) -c $< -o $@  $(CXXFLAGS)

$(OBJDIR)/%.o: %.cxx
	@echo Compiling $< ......
	@$(CXX) -c $< -o $@  $(CXXFLAGS)

$(OBJDIR)/%.o: %.f
	@echo Compiling $< ......
	@$(FF) -c $< -o $@  $(FFLAGS)

$(OBJDIR)/%.o: %.F
	@echo Compiling $< ......
	@$(FF) -c $< -o $@  $(FFLAGS)

dir:
	@if [ ! -d $(OBJDIR) ] ; then mkdir -p $(OBJDIR) ;fi

########################################################################
clean: dir
	@rm -f $(OBJDIR)/*
	@rm -f $(EXECFILE)
	@rm -f *~ *# */*~ */*#

distclean: clean

test:
	@echo \\MYOS\:$(MYOS) \\ARCH\:$(ARCH)
	@echo \\CFLAGS\:$(CFLAGS)
	@echo \\CXXFLAGS\:$(CXXFLAGS)
	@echo \\FFLAGS\:$(FFLAGS)
	@echo \\LDFLAGS\:$(LDFLAGS)
	@echo \\SYSLIBS\:$(SYSLIBS)
	@echo \\fsources\: $(FSOURCES)
	@echo \\sources\: $(SOURCES)
	@echo \\headers\: $(HEADERS)
	@echo \\objs\: $(OBJS)
	@echo \\dependencies: \$(DEPS)

help: tes

env: tes
