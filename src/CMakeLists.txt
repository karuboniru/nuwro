# ---------------------------------------------------------------------------
# params_all.h  — auto-generated from params.xml via params.sed
# Generated into the build directory so the source tree stays clean.
# ---------------------------------------------------------------------------
add_custom_command(
  OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/params_all.h
  COMMAND ${CMAKE_COMMAND}
    -DOUTPUT=${CMAKE_CURRENT_BINARY_DIR}/params_all.h
    -DSED_FILE=${CMAKE_CURRENT_SOURCE_DIR}/params.sed
    -DXML_FILE=${CMAKE_CURRENT_SOURCE_DIR}/params.xml
    -P ${CMAKE_SOURCE_DIR}/cmake/generate_params.cmake
  DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/params.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/params.sed
  COMMENT "Generating params_all.h"
)
add_custom_target(params_header
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/params_all.h)

# ---------------------------------------------------------------------------
# ROOT dictionary  — event1Dict.cxx + event1Dict_rdict.pcm
# ROOT_GENERATE_DICTIONARY calls rootcling (ROOT 6) internally.
# The OPTIONS -I flag makes rootcling see the local headers.
# ---------------------------------------------------------------------------
ROOT_GENERATE_DICTIONARY(event1Dict
  event1.h
  LINKDEF event1LinkDef.h
  OPTIONS
    -I${CMAKE_CURRENT_SOURCE_DIR}
    -I${CMAKE_CURRENT_BINARY_DIR}   # for generated params_all.h
)

# ---------------------------------------------------------------------------
# Source lists  (mirror the Makefile groupings)
# ---------------------------------------------------------------------------

# --- DIS ---
set(DIS_SOURCES
  dis/charge.cc         dis/LeptonMass.cc       dis/parameters.cc
  dis/grv94_bodek.cc    dis/dis_cr_sec.cc        dis/dis_nc.cc
  dis/dis_cc_neutron.cc dis/delta.cc             dis/dis2res.cc
  dis/dis_cc_proton.cc  dis/fragmentation.cc     dis/fragmentation_nc.cc
  dis/fragmentation_cc.cc dis/singlepion.cc      dis/disevent.cc
  dis/resevent2.cc      dis/singlepionhadr.cc    dis/alfa.cc
  dis/res_kinematics.cc dis/res_xsec.cc
)

# --- MEC ---
set(MEC_SOURCES
  mecdynamics.cc         mecevent.cc
  mecevent_tem.cc        mecevent_Nieves.cc
  mecevent_SuSA.cc       mecevent_2020Valencia.cc
  mecevent_common.cc     mecdynamics2.cc
  mecevent2.cc
)

# --- Spectral functions (all .cc in sf/) ---
file(GLOB SF_SOURCES CONFIGURE_DEPENDS sf/*.cc)

# --- Hybrid RES model (all .cc in hybrid/ + resevent_hybrid.cc) ---
file(GLOB HYBRID_CC CONFIGURE_DEPENDS hybrid/*.cc)
set(HYBRID_SOURCES ${HYBRID_CC} resevent_hybrid.cc)

# --- ESPP ---
set(ESPP_SOURCES
  espp/crossection.cc
  espp/DM.cc
)

# --- Reweighting ---
set(REW_SOURCES
  rew/Reweighters.cc
  rew/rewQEL.cc
  rew/rewRES.cc
  rew/rewNorm.cc
  rew/rewparams.cc
  rew/PythiaQuiet.f
)

# --- FSI (subset used by the standalone kaskada tool) ---
set(FSI_SOURCES
  scatter.cc          kaskada7.cc         Interaction.cc
  hyperon_cascade.cc  hyperon_interaction.cc
  input_data.cc       data_container.cc
  nucleus.cc          nucleus_data.cc
  isotopes.cc         elements.cc
  fsi.cc              pitab.cc
)

# --- Core event structure (shared / static library) ---
set(EVENT_SOURCES
  event1.cc
  generatormt.cc
  particle.cc
  pdg.cc
  dirs.cc
  ${CMAKE_CURRENT_BINARY_DIR}/event1Dict.cxx
)

# --- Physics common to nuwro, reweight_to, reweight_along ---
set(COMMON_SOURCES
  pauli.cc
  cohevent2.cc        cohdynamics2.cc
  qelevent1.cc        qel_sigma.cc
  kinsolver.cc        kinematics.cc
  target_mixer.cc     nucleus.cc
  sfevent.cc          ff.cc
  rpa_2013.cc         nucleus_data.cc
  isotopes.cc         elements.cc
  beam.cc             nd280stats.cc
  beamHist.cc         coh.cc
  fsi.cc              pitab.cc
  scatter.cc          kaskada7.cc
  Interaction.cc      input_data.cc
  data_container.cc   hyperon_cascade.cc
  hypevent.cc         hyperon_interaction.cc
  e_el_sigma.cc
  lepevent.cc         nu_e_el_sigma.cc
  generatormt.cc      rew/rewparams.cc
  rew/PythiaQuiet.f
)

# ---------------------------------------------------------------------------
# Common include dirs and link libraries used by all targets in this dir
# ---------------------------------------------------------------------------
set(NUWRO_INCLUDE_DIRS
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}   # for params_all.h
  ${ROOT_INCLUDE_DIRS}          # needed by rootcling at dict-generation time
)

# Use ROOT imported targets so CMake knows library locations and can compute
# RPATH correctly.  ${ROOT_LIBRARIES} is a list of raw -lFoo flags that
# produce empty RPATH entries; named targets do not have that problem.
set(NUWRO_ROOT_LIBS
  ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist
  ROOT::Physics ROOT::MathCore
  ROOT::EG ROOT::Geom ROOT::Minuit
)

set(NUWRO_LINK_LIBS
  ${NUWRO_ROOT_LIBS}
  Pythia6::Pythia6
  Pythia6::EGPythia6
  ${NUWRO_FORTRAN_RUNTIME}
)

# ---------------------------------------------------------------------------
# Libraries
# ---------------------------------------------------------------------------

# Shared event library — loaded by ROOT macros (myroot, analysis scripts)
add_library(event SHARED ${EVENT_SOURCES})
add_dependencies(event params_header)
# Use generator expressions so the exported target points to the installed
# include/ directory rather than the build tree.
target_include_directories(event
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    $<INSTALL_INTERFACE:include>
)
# ROOT headers/libs are PUBLIC because event1.h exposes ROOT types.
# Pythia6 and the Fortran runtime are not used by the event library sources
# (event1.cc, particle.cc, pdg.cc, dirs.cc, generatormt.cc, event1Dict.cxx)
# so they are not linked here at all.
target_link_libraries(event PUBLIC ${NUWRO_ROOT_LIBS})
target_link_options(event PRIVATE $<$<PLATFORM_ID:Linux>:-Wl,--no-undefined>)

# Static event library — linked into executables
add_library(event1 STATIC ${EVENT_SOURCES})
add_dependencies(event1 params_header)
target_include_directories(event1 PUBLIC ${NUWRO_INCLUDE_DIRS})
target_link_libraries(event1 PUBLIC ${NUWRO_LINK_LIBS})

# DIS
add_library(dis1 STATIC ${DIS_SOURCES})
add_dependencies(dis1 params_header)
target_include_directories(dis1 PRIVATE ${NUWRO_INCLUDE_DIRS})
target_link_libraries(dis1 PRIVATE event1)

# Spectral functions
add_library(sf1 STATIC ${SF_SOURCES})
add_dependencies(sf1 params_header)
target_include_directories(sf1 PRIVATE ${NUWRO_INCLUDE_DIRS})
target_link_libraries(sf1 PRIVATE event1)

# Hybrid RES model
add_library(hybrid1 STATIC ${HYBRID_SOURCES})
add_dependencies(hybrid1 params_header)
target_include_directories(hybrid1 PRIVATE ${NUWRO_INCLUDE_DIRS})
target_link_libraries(hybrid1 PRIVATE event1)

# MEC
add_library(mec1 STATIC ${MEC_SOURCES})
add_dependencies(mec1 params_header)
target_include_directories(mec1 PRIVATE ${NUWRO_INCLUDE_DIRS})
target_link_libraries(mec1 PRIVATE event1)

# ESPP
add_library(espp1 STATIC ${ESPP_SOURCES})
add_dependencies(espp1 params_header)
target_include_directories(espp1 PRIVATE ${NUWRO_INCLUDE_DIRS})
target_link_libraries(espp1 PRIVATE event1)

# FSI (standalone subset)
add_library(fsi1 STATIC ${FSI_SOURCES} pdg.cc)
add_dependencies(fsi1 params_header)
target_include_directories(fsi1 PRIVATE ${NUWRO_INCLUDE_DIRS})
target_link_libraries(fsi1 PRIVATE event1)

# Common physics
add_library(common1 STATIC ${COMMON_SOURCES})
add_dependencies(common1 params_header)
target_include_directories(common1 PRIVATE ${NUWRO_INCLUDE_DIRS})
target_link_libraries(common1 PRIVATE event1 mec1 dis1 sf1 hybrid1 espp1)

# Reweighting
add_library(rew1 STATIC ${REW_SOURCES})
add_dependencies(rew1 params_header)
target_include_directories(rew1 PRIVATE ${NUWRO_INCLUDE_DIRS})
target_link_libraries(rew1
  PUBLIC  event1 common1 fsi1 dis1 sf1
  PRIVATE ${NUWRO_FORTRAN_RUNTIME}
)

# ---------------------------------------------------------------------------
# Helper macro: set common target properties in one call
# ---------------------------------------------------------------------------
macro(nuwro_target_setup _tgt)
  add_dependencies(${_tgt} params_header)
  target_include_directories(${_tgt} PRIVATE ${NUWRO_INCLUDE_DIRS})
endmacro()

# ---------------------------------------------------------------------------
# Executables
# ---------------------------------------------------------------------------

add_executable(nuwro main.cc nuwro.cc)
nuwro_target_setup(nuwro)
target_link_libraries(nuwro PRIVATE
  event1 common1 mec1 dis1 espp1 sf1 hybrid1 ${NUWRO_LINK_LIBS})

add_executable(kaskada
  kaskada.cc        scatter.cc         kaskada7.cc
  Interaction.cc    hyperon_interaction.cc hyperon_cascade.cc
  input_data.cc     data_container.cc  nucleus.cc
  nucleus_data.cc   isotopes.cc        elements.cc
  rew/rewparams.cc  ff.cc              fsi.cc
  pitab.cc
)
nuwro_target_setup(kaskada)
target_link_libraries(kaskada PRIVATE event1 ${NUWRO_LINK_LIBS})

add_executable(myroot myroot.cc ${CMAKE_CURRENT_BINARY_DIR}/event1Dict.cxx)
nuwro_target_setup(myroot)
target_link_libraries(myroot PRIVATE event1 ${NUWRO_LINK_LIBS} ROOT::Rint)

add_executable(nuwro2neut       nuwro2neut.cc)
nuwro_target_setup(nuwro2neut)
target_link_libraries(nuwro2neut PRIVATE event1 ${NUWRO_LINK_LIBS})

add_executable(nuwro2nuance     nuwro2nuance.cc)
nuwro_target_setup(nuwro2nuance)
target_link_libraries(nuwro2nuance PRIVATE event1 ${NUWRO_LINK_LIBS})

add_executable(nuwro2rootracker nuwro2rootracker.cc)
nuwro_target_setup(nuwro2rootracker)
target_link_libraries(nuwro2rootracker PRIVATE event1 ${NUWRO_LINK_LIBS})

add_executable(glue glue.cc)
nuwro_target_setup(glue)
target_link_libraries(glue PRIVATE event1 ${NUWRO_LINK_LIBS})

add_executable(dumpParams dumpParams.cc dirs.cc)
nuwro_target_setup(dumpParams)
# dumpParams only needs ROOT core (no event lib needed)
target_link_libraries(dumpParams PRIVATE ${ROOT_LIBRARIES})

add_executable(reweight_to rew/reweight_to.cc)
nuwro_target_setup(reweight_to)
target_link_libraries(reweight_to PRIVATE
  event1 common1 fsi1 dis1 sf1 rew1 ${NUWRO_LINK_LIBS})

add_executable(reweight_along rew/reweight_along.cc)
nuwro_target_setup(reweight_along)
target_link_libraries(reweight_along PRIVATE
  event1 common1 fsi1 dis1 sf1 rew1 ${NUWRO_LINK_LIBS})

add_executable(whist whist.cc)
nuwro_target_setup(whist)
target_link_libraries(whist PRIVATE event1 ${NUWRO_LINK_LIBS})

# Optional / disabled test executables (uncomment to build):
# add_executable(test_nucleus generatormt.cc nucleus.cc test_nucleus.cc
#                pdg.cc dirs.cc nucleus_data.cc isotopes.cc elements.cc)
# target_link_libraries(test_nucleus PRIVATE event1 ${NUWRO_LINK_LIBS})

# add_executable(test_beam generatormt.cc pdg.cc test_beam.cc)
# target_link_libraries(test_beam PRIVATE event1 ${NUWRO_LINK_LIBS})

# add_executable(test_beam_rf test_beam_rf.cc particle.cc generatormt.cc)
# target_link_libraries(test_beam_rf PRIVATE event1 ${NUWRO_LINK_LIBS})

# add_executable(test_makehist test_makehist.cc nd280stats.cc)
# target_link_libraries(test_makehist PRIVATE event1 ${NUWRO_LINK_LIBS})

# ---------------------------------------------------------------------------
# Installation
# ---------------------------------------------------------------------------

install(TARGETS
    nuwro kaskada myroot
    nuwro2neut nuwro2nuance nuwro2rootracker
    glue dumpParams reweight_to reweight_along whist
  RUNTIME DESTINATION bin
)

install(TARGETS event
  EXPORT  NuWroTargets
  LIBRARY DESTINATION lib
)

# ---------------------------------------------------------------------------
# Generator library  (libnuwrogen.so)
# Provides libnuwro::Generator — fluent MH-based event generation API.
# Self-contained: all physics statics are linked privately so consumers only
# need to link this one target (plus NuWro::event for ROOT I/O if needed).
# ---------------------------------------------------------------------------
add_library(nuwrogen SHARED libnuwro.cc)
add_dependencies(nuwrogen params_header)
target_include_directories(nuwrogen
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    $<INSTALL_INTERFACE:include>
)
# ROOT types appear in the public headers (event1.h, params.h) → PUBLIC.
# All physics and Pythia6/Fortran are implementation details → PRIVATE.
target_link_libraries(nuwrogen
  PUBLIC  ${NUWRO_ROOT_LIBS}
  PRIVATE event1 common1 mec1 dis1 espp1 sf1 hybrid1
          Pythia6::Pythia6 Pythia6::EGPythia6 ${NUWRO_FORTRAN_RUNTIME}
)
target_link_options(nuwrogen PRIVATE $<$<PLATFORM_ID:Linux>:-Wl,--no-undefined>)

install(TARGETS nuwrogen
  EXPORT  NuWroTargets
  LIBRARY DESTINATION lib
)

install(EXPORT NuWroTargets
  NAMESPACE   NuWro::
  DESTINATION lib/cmake/NuWro
)

# ROOT 6 PCM + rootmap files — needed for dictionary reflection and auto-load
foreach(_dest lib bin)
  install(FILES
      ${CMAKE_CURRENT_BINARY_DIR}/libevent1Dict_rdict.pcm
      ${CMAKE_CURRENT_BINARY_DIR}/libevent1Dict.rootmap
    DESTINATION ${_dest}
    OPTIONAL      # OPTIONAL: ROOT 5 does not produce these files
  )
endforeach()

install(DIRECTORY ${CMAKE_SOURCE_DIR}/data/
  DESTINATION data
)

install(FILES
    params.h       vec.h           dirs.h
    particle.h     generatormt.h
    jednostki.h    pdg.h           vect.h
    vecrand.h      event1.h
    libnuwro.h
  DESTINATION include
)
install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/params_all.h
  DESTINATION include
)
