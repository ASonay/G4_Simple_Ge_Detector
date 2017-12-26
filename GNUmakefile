# $Id: GNUmakefile 66357 2012-12-18 09:06:10Z gcosmo $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := GeDet
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../..
endif

ifdef DMXENV_GPS_USE
  CPPFLAGS += -DDMXENV_GPS_USE
endif

include $(G4INSTALL)/config/architecture.gmk

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk

CPPFLAGS += $(shell $(ROOTSYS)/bin/root-config --cflags)
EXTRALIBS += $(shell $(ROOTSYS)/bin/root-config --libs)
LDLIBS += $(shell $(ROOTSYS)/bin/root-config --libs) 

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*
