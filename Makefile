#---------------------------------------------------
#QTINCDIR =$(wildcard /usr/include/qt4/ ${QTDIR}/include/)
#QTINCLUDEDIRS = -I. -I/usr/include/qt4/ -I/usr/include/qt4/QtCore -I/usr/include/qt4/QtGui -I${QTDIR}/include/QtCore -I${QTDIR}/include/QtGui
#QTINCLUDEDIRS = -I$(QTINCDIR)/ -I$(QTINCDIR)/QtCore -I$(QTINCDIR)/QtGui 
#QTLIBS =   -L/usr/lib -lQtGui -lQtCore -lpthread 

# Detecting OS type
OS := $(shell uname)

# get version tag using git describe

VERSION := $(shell git describe --tags)

DEBUG         = 1
DEBUGON = -g 
COMFLAGS = -Wall -fopenmp -flto=auto
COMPILER = c++
#CXXFLAGS      = `${ROOTSYS}/bin/root-config --cflags` -fPIC -flto=auto -O2 -I src
ifeq ($(OS),Darwin)
  # Flags for OSX
  CXXFLAGS      = `${ROOTSYS}/bin/root-config --cflags` -fPIC -O2 $(DEBUGON) $(COMFLAGS) -I src -Wall -Wno-unused-variable -Wno-sign-compare -Wno-unused-function -Wno-unused-but-set-variable -Wno-reorder $(QTINCLUDEDIRS) -DVERSION=\"$(VERSION)\"
else
  # Flags for others
  CXXFLAGS      = `${ROOTSYS}/bin/root-config --cflags` -g  -march=native -std=c++17  -fPIC -O2 $(DEBUGON) $(COMFLAGS) -I src -Wl,--no-as-needed -Wall -Wno-unused-variable -Wno-sign-compare -Wno-unused-function -Wno-unused-but-set-variable -Wno-reorder $(QTINCLUDEDIRS) -DVERSION=\"$(VERSION)\"
endif
#LDFLAGS       = `${ROOTSYS}/bin/root-config --libs` -lPythia6 -lEG -lEGPythia6 -lCore  -lCint -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lGeom -lpthread -lm -ldl -rdynamic -lHist $(QTLIBS)
LDFLAGS       = `${ROOTSYS}/bin/root-config --libs` -lPythia6 -lEG -lEGPythia6 -lGeom -lMinuit -lgfortran $(QTLIBS) $(COMFLAGS)
LD	      = $(COMPILER)
CXX	      = $(COMPILER)
CC 	      = $(COMPILER)
FC            = gfortran

%.o: %.cc
		$(COMPILER) ${CXXFLAGS} -c $< -o $@

%.o: %.f
		gfortran  -c $< -o $@

TRGTS =         $(addprefix $(BIN)/,nuwro kaskada myroot glue event1.so nuwro2neut nuwro2nuance nuwro2rootracker\
                dumpParams test_beam_rf test_makehist test_nucleus test_beam \
                fsi niwg ladek_topologies test mb_nce_run ganalysis reweight_to reweight_along whist\
								mktabular mktabular2d nuwro_metropolis\
                )

DIS=    charge.o LeptonMass.o parameters.o grv94_bodek.o dis_cr_sec.o  dis_nc.o dis_cc_neutron.o delta.o dis2res.o \
	dis_cc_proton.o fragmentation.o fragmentation_nc.o fragmentation_cc.o singlepion.o \
	disevent.o resevent2.o singlepionhadr.o alfa.o res_kinematics.o res_xsec.o

ESPP_OBJS=$(patsubst %.cc,%.o,$(wildcard src/espp/*.cc)) src/e_spp_event.o 
SF_OBJS = $(patsubst %.cc,%.o,$(wildcard src/sf/*.cc))
HYBRID_OBJS = $(patsubst %.cc,%.o,$(wildcard src/hybrid/*.cc)) src/resevent_hybrid.o
GUI_OBJS = $(patsubst %.cc,%.o,$(wildcard src/gui/*.cc))
GUI_OBJS += $(patsubst src/gui/C%.cc,src/gui/moc_C%.o,$(wildcard src/gui/C*.cc))

DIS_OBJS= $(addprefix src/dis/,$(DIS))
BIN=bin
#BIN=.

EVENT_OBJS =  $(addprefix src/, event1.o event1dict.o pdg.o particle.o generatormt.o dirs.o rew/rewparams.o)


all:            $(TRGTS)

$(BIN)/whist: src/whist.o $(EVENT_OBJS)
		$(LINK.cc) $^ -o $@

$(BIN)/nuwro:   $(addprefix src/, event1.o event1dict.o generatormt.o particle.o pauli.o cohevent2.o cohdynamics2.o qelevent1.o hipevent.o\
	    mecdynamics.o mecevent.o mecevent_tem.o mecevent_Nieves.o mecevent_SuSA.o mecevent_common.o e_el_event.o e_el_sigma.o\
	    mecdynamics2.o mecevent2.o rew/rewparams.o\
        qel_sigma.o kinsolver.o kinematics.o pdg.o target_mixer.o nucleus.o  sfevent.o ff.o dirs.o rpa_2013.o\
        nucleus_data.o isotopes.o elements.o rew/PythiaQuiet.o\
        nuwro.o beam.o nd280stats.o beamHist.o coh.o fsi.o pitab.o scatter.o kaskada7.o Interaction.o input_data.o data_container.o  main.o) \
        $(SF_OBJS) $(DIS_OBJS) $(ESPP_OBJS) $(HYBRID_OBJS)
		$(LINK.cc) $^ -fwhole-program -o $@

$(BIN)/nuwro_metropolis:   $(addprefix src/, event1.o event1dict.o generatormt.o particle.o pauli.o cohevent2.o cohdynamics2.o qelevent1.o hipevent.o\
	    mecdynamics.o mecevent.o mecevent_tem.o mecevent_Nieves.o mecevent_SuSA.o mecevent_common.o e_el_event.o e_el_sigma.o\
	    mecdynamics2.o mecevent2.o rew/rewparams.o\
        qel_sigma.o kinsolver.o kinematics.o pdg.o target_mixer.o nucleus.o  sfevent.o ff.o dirs.o rpa_2013.o\
        nucleus_data.o isotopes.o elements.o rew/PythiaQuiet.o\
       beam.o nd280stats.o beamHist.o coh.o fsi.o pitab.o scatter.o kaskada7.o Interaction.o input_data.o data_container.o  nuwro_metropolis_main.o nuwro_metropolis.o) \
        $(SF_OBJS) $(DIS_OBJS) $(ESPP_OBJS) $(HYBRID_OBJS)
		$(LINK.cc) $^ -fwhole-program -o $@

$(BIN)/mktabular:   $(addprefix src/, event1.o event1dict.o generatormt.o particle.o pauli.o cohevent2.o cohdynamics2.o qelevent1.o hipevent.o\
	    mecdynamics.o mecevent.o mecevent_tem.o mecevent_Nieves.o mecevent_SuSA.o mecevent_common.o e_el_event.o e_el_sigma.o\
	    mecdynamics2.o mecevent2.o rew/rewparams.o\
        qel_sigma.o kinsolver.o kinematics.o pdg.o target_mixer.o nucleus.o  sfevent.o ff.o dirs.o rpa_2013.o\
        nucleus_data.o isotopes.o elements.o rew/PythiaQuiet.o\
        nuwro.o beam.o nd280stats.o beamHist.o coh.o fsi.o pitab.o scatter.o kaskada7.o Interaction.o input_data.o data_container.o  mktabular.o) \
        $(SF_OBJS) $(DIS_OBJS) $(ESPP_OBJS) $(HYBRID_OBJS)
		$(LINK.cc) $^ -fwhole-program -o $@

$(BIN)/mktabular2d:   $(addprefix src/, event1.o event1dict.o generatormt.o particle.o pauli.o cohevent2.o cohdynamics2.o qelevent1.o hipevent.o\
	    mecdynamics.o mecevent.o mecevent_tem.o mecevent_Nieves.o mecevent_SuSA.o mecevent_common.o e_el_event.o e_el_sigma.o\
	    mecdynamics2.o mecevent2.o rew/rewparams.o\
        qel_sigma.o kinsolver.o kinematics.o pdg.o target_mixer.o nucleus.o  sfevent.o ff.o dirs.o rpa_2013.o\
        nucleus_data.o isotopes.o elements.o rew/PythiaQuiet.o\
        nuwro.o beam.o nd280stats.o beamHist.o coh.o fsi.o pitab.o scatter.o kaskada7.o Interaction.o input_data.o data_container.o  mktabular2d.o) \
        $(SF_OBJS) $(DIS_OBJS) $(ESPP_OBJS) $(HYBRID_OBJS)
		$(LINK.cc) $^ -fwhole-program -o $@

$(BIN)/kaskada:   $(addprefix src/, scatter.o generatormt.o particle.o event1.o event1dict.o kaskada7.o Interaction.o input_data.o data_container.o dirs.o\
				  pdg.o nucleus.o kaskada.o fsi.o pitab.o nucleus_data.o isotopes.o elements.o rew/rewparams.o)
		$(LINK.cc) $^ -o $@

$(BIN)/myroot:  $(EVENT_OBJS) src/myroot.o
		$(LINK.cc) $^ -o $@

$(BIN)/event1.so: $(EVENT_OBJS)
	        $(LD) -shared  $(LDFLAGS) $^ $(LIBS) -o $@

#$(BIN)/event1.a: $(EVENT_OBJS)
#	         ar r $@ $^  

$(BIN)/glue:    $(EVENT_OBJS) src/glue.o
		$(LINK.cc) $^ -o $@

$(BIN)/nuwro2neut:  $(EVENT_OBJS) src/nuwro2neut.o 
		$(LINK.cc) $^ -o $@

$(BIN)/nuwro2nuance: $(EVENT_OBJS) src/nuwro2nuance.o
		 $(LINK.cc) $^ -o $@

$(BIN)/nuwro2rootracker: $(EVENT_OBJS) src/nuwro2rootracker.o
		 $(LINK.cc) $^ -o $@

$(BIN)/fsi:   src/scatter.o src/generatormt.o src/particle.o src/event1.o src/event1dict.o src/kaskada7.o src/Interaction.o src/input_data.o src/data_container.o src/pdg.o src/dirs.o  src/nucleus.o  src/nucleus_data.o src/isotopes.o src/elements.o\
       src/fsi.o src/pitab.o src/calculations.o src/simulations.o src/vivisection.o src/plots.o  src/mplots.o  src/dirs.o src/fsi_main.o  src/rew/rewparams.o
		$(LINK.cc) $^ -o $@

$(BIN)/mb_nce_run:   src/mb_nce.o src/mb_nce_run.o src/event1.o src/event1dict.o src/mb_nce_fit.o src/pdg.o src/scatter.o src/generatormt.o src/dirs.o src/particle.o
		$(LINK.cc) $^ -o $@

$(BIN)/niwg:   src/scatter.o src/generatormt.o src/particle.o src/event1.o src/event1dict.o src/kaskada7.o src/Interaction.o src/input_data.o src/data_container.o src/pdg.o src/dirs.o  src/nucleus.o  src/nucleus_data.o src/isotopes.o src/elements.o\
        src/fsi.o src/pitab.o src/calculations.o src/niwg_ccqe.o src/niwg_tech.o src/niwg_ccpi.o src/niwg.o  
		$(LINK.cc) $^ -o $@

$(BIN)/ladek_topologies: src/event1.o src/event1dict.o src/pdg.o src/particle.o  src/generatormt.o src/ladek_topologies.o src/dirs.o \
          src/fsi.o src/pitab.o
		$(LINK.cc) $^ -o $@

$(BIN)/test: src/event1.o src/event1dict.o src/pdg.o src/particle.o  src/generatormt.o src/test.o src/dirs.o src/nucleus.o src/nucleus_data.o src/isotopes.o src/elements.o
		$(LINK.cc) $^ -o $@
		
$(BIN)/ganalysis: $(addprefix src/, \
		event1.o event1dict.o generatormt.o particle.o pauli.o cohevent2.o cohdynamics2.o qelevent1.o mecdynamics.o mecevent.o hipevent.o\
	    mecdynamics2.o mecevent2.o mecevent_tem.o mecevent_Nieves.o mecevent_SuSA.o mecevent_common.o e_el_event.o e_el_sigma.o rew/PythiaQuiet.o\
        qel_sigma.o kinsolver.o kinematics.o pdg.o target_mixer.o nucleus.o  sfevent.o ff.o dirs.o rpa_2013.o nucleus_data.o isotopes.o elements.o \
        nuwro.o beam.o nd280stats.o beamHist.o coh.o fsi.o pitab.o scatter.o kaskada7.o Interaction.o input_data.o data_container.o ganalysis.o rew/rewparams.o) \
        $(SF_OBJS) $(DIS_OBJS) $(ESPP_OBJS) $(HYBRID_OBJS)
		$(LINK.cc) $^ -o $@

$(BIN)/reweight_to: $(addprefix src/, \
		event1.o event1dict.o generatormt.o particle.o pauli.o cohevent2.o cohdynamics2.o qelevent1.o mecdynamics.o mecevent.o hipevent.o\
	    mecdynamics2.o mecevent2.o mecevent_tem.o mecevent_Nieves.o mecevent_SuSA.o mecevent_common.o e_el_event.o e_el_sigma.o\
        qel_sigma.o kinsolver.o kinematics.o pdg.o target_mixer.o nucleus.o  sfevent.o ff.o dirs.o rpa_2013.o nucleus_data.o isotopes.o elements.o \
        nuwro.o beam.o nd280stats.o beamHist.o coh.o fsi.o pitab.o scatter.o kaskada7.o Interaction.o input_data.o data_container.o\
        rew/rewparams.o rew/Reweighters.o rew/rewQEL.o rew/rewRES.o rew/rewNorm.o rew/reweight_to.o rew/PythiaQuiet.o) \
        $(SF_OBJS) $(DIS_OBJS) $(ESPP_OBJS) $(HYBRID_OBJS)
		$(LINK.cc)  $^ -o $@ 

$(BIN)/reweight_along: $(addprefix src/, \
		event1.o event1dict.o generatormt.o particle.o pauli.o cohevent2.o cohdynamics2.o qelevent1.o mecdynamics.o mecevent.o hipevent.o\
	    mecdynamics2.o mecevent2.o mecevent_tem.o mecevent_Nieves.o mecevent_SuSA.o mecevent_common.o e_el_event.o e_el_sigma.o\
        qel_sigma.o kinsolver.o kinematics.o pdg.o target_mixer.o nucleus.o  sfevent.o ff.o dirs.o rpa_2013.o nucleus_data.o isotopes.o elements.o \
        nuwro.o beam.o nd280stats.o beamHist.o coh.o fsi.o pitab.o scatter.o kaskada7.o Interaction.o input_data.o data_container.o\
        rew/rewparams.o rew/Reweighters.o rew/rewQEL.o rew/rewRES.o rew/rewNorm.o rew/reweight_along.o rew/PythiaQuiet.o) \
        $(SF_OBJS) $(DIS_OBJS) $(ESPP_OBJS) $(HYBRID_OBJS)
		$(LINK.cc)  $^ -o $@ 


$(BIN)/dumpParams:      src/dumpParams.o src/dirs.o
		$(LINK.cc) $^ -o $@

$(BIN)/test_nucleus:   src/generatormt.o src/nucleus.o src/test_nucleus.o src/pdg.o src/dirs.o  src/nucleus_data.o src/isotopes.o src/elements.o
		$(LINK.cc) $^ -o $@

$(BIN)/test_beam:	src/generatormt.o src/pdg.o src/test_beam.o 
		$(LINK.cc) $^ -o $@

$(BIN)/test_beam_rf:      src/test_beam_rf.o src/particle.o  src/generatormt.o 
		$(LINK.cc) $^ -o $@

$(BIN)/test_makehist:    src/test_makehist.o src/nd280stats.o
		$(LINK.cc) $^ -o $@

$(BIN)/test_balancer:       src/test_balancer.cc  src/generatormt.o
		$(LINK.cc) $^ -o $@

clean:;         @rm -f          *.o *.d src/event1dict.* src/event1dict_rdict.pcm       core src/dis/*.o src/dis/*.d src/sf/*.o src/sf/*.d src/*.o src/*.d\
		src/gui/*.o src/gui/*.d src/gui/moc_* src/rew/*.o src/espp/*.o src/hybrid/*.o


distclean:;     @rm -f $(TRGTS) *.o *.d src/event1dict.* {src,bin}/event1dict_rdict.pcm core src/dis/*.o src/dis/*.d src/sf/*.o src/sf/*.d src/*.o src/*.d\
		src/gui/*.o src/gui/*.d src/gui/moc_* src/rew/*.o src/espp/*.o src/hybrid/*.o *.root *.root.txt


src/event1dict.h src/event1dict.cc:  src/params_all.h src/params.h src/event1.h src/event1LinkDef.h src/event1.o
		@echo "Generating dictionary ..."
		cd src;${ROOTSYS}/bin/rootcint -f event1dict.cc -c event1.h event1LinkDef.h;cd ..
		cp src/event1dict_rdict.pcm bin

src/params_all.h:  src/params.xml src/params.h src/params.sed Makefile
		@echo "Building params_all.h"
		@xmllint --dtdvalid src/params.dtd src/params.xml --noout
		@echo "#define PARAMS_ALL()\\">src/params_all.h
		@sed -f src/params.sed src/params.xml >> src/params_all.h 
		@echo "" >> src/params_all.h


%.d: %.cc
	@echo Making dependencies for $<
	@$(SHELL) -ec '$(CC) -MM -MT "$@ $<" $(CXXFLAGS) $< \
	| sed s/.cc\:/.o:/ > $@;\
	[ -s $@ ] || rm -f $@'


ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
 -include $(patsubst %.cc,%.d,$(wildcard  src/dis/*.cc src/sf/*.cc src/*.cc src/gui/*.cc)) src/event1dict.d
endif
endif

doxygen:
	@doxygen doc/doxygen/Doxyfile

# DO NOT DELETE
