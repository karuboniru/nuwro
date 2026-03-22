cmake_minimum_required(VERSION 3.16 FATAL_ERROR)

project(NuWro
  DESCRIPTION "Monte Carlo neutrino event generator"
  LANGUAGES CXX Fortran
)

# ---------------------------------------------------------------------------
# Version from git
# ---------------------------------------------------------------------------
find_package(Git QUIET)
if(Git_FOUND)
  execute_process(
    COMMAND ${GIT_EXECUTABLE} describe --tags
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE NUWRO_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
  )
endif()
if(NOT NUWRO_VERSION OR NUWRO_VERSION STREQUAL "")
  set(NUWRO_VERSION "unknown")
endif()
message(STATUS "NuWro version: ${NUWRO_VERSION}")

# ---------------------------------------------------------------------------
# Build type
# ---------------------------------------------------------------------------
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY
    STRINGS Debug Release RelWithDebInfo MinSizeRel)
endif()

# ---------------------------------------------------------------------------
# CERN ROOT
# Use ROOT's own CMake config (find ROOT_DIR from ROOTSYS if not set).
# ---------------------------------------------------------------------------
if(NOT ROOT_DIR AND DEFINED ENV{ROOTSYS})
  set(ROOT_DIR "$ENV{ROOTSYS}/share/root/cmake" CACHE PATH "ROOT CMake config dir" FORCE)
endif()

find_package(ROOT 6 REQUIRED
  CONFIG
  COMPONENTS Core RIO Tree Hist Physics MathCore EG Geom Minuit
)
# include(${ROOT_USE_FILE})

message(STATUS "ROOT version    : ${ROOT_VERSION}")
message(STATUS "ROOT include    : ${ROOT_INCLUDE_DIRS}")
message(STATUS "ROOT lib dir    : ${ROOT_LIBRARY_DIR}")

# ---------------------------------------------------------------------------
# Pythia6 / EGPythia6
# Searched in ROOTEGPythia6_ROOT (env or cmake var), then ROOT's lib dir.
# ---------------------------------------------------------------------------
if(NOT DEFINED ROOTEGPythia6_ROOT AND DEFINED ENV{ROOTEGPythia6_ROOT})
  set(ROOTEGPythia6_ROOT "$ENV{ROOTEGPythia6_ROOT}"
    CACHE PATH "Path to the EGPythia6 / Pythia6 installation prefix")
endif()

if(ROOTEGPythia6_ROOT)
  find_library(Pythia6_LIBRARY   NAMES Pythia6
    PATHS "${ROOTEGPythia6_ROOT}/lib" NO_DEFAULT_PATH REQUIRED)
  find_library(EGPythia6_LIBRARY NAMES EGPythia6
    PATHS "${ROOTEGPythia6_ROOT}/lib" NO_DEFAULT_PATH REQUIRED)
  find_path(Pythia6_INCLUDE_DIR  NAMES TPythia6.h
    PATHS "${ROOTEGPythia6_ROOT}/include" NO_DEFAULT_PATH)
else()
  # Fall back to ROOT's own library directory (some ROOT builds bundle Pythia6)
  find_library(Pythia6_LIBRARY   NAMES Pythia6
    HINTS ${ROOT_LIBRARY_DIR} REQUIRED)
  find_library(EGPythia6_LIBRARY NAMES EGPythia6
    HINTS ${ROOT_LIBRARY_DIR} REQUIRED)
  find_path(Pythia6_INCLUDE_DIR  NAMES TPythia6.h
    HINTS ${ROOT_INCLUDE_DIRS})
endif()

message(STATUS "Pythia6 library  : ${Pythia6_LIBRARY}")
message(STATUS "EGPythia6 library: ${EGPythia6_LIBRARY}")

# Create imported targets so subdirectories can link them by name
foreach(_lib Pythia6 EGPythia6)
  if(NOT TARGET Pythia6::${_lib})
    add_library(Pythia6::${_lib} UNKNOWN IMPORTED)
    set_target_properties(Pythia6::${_lib} PROPERTIES
      IMPORTED_LOCATION "${${_lib}_LIBRARY}")
  endif()
endforeach()

if(Pythia6_INCLUDE_DIR)
  foreach(_lib Pythia6 EGPythia6)
    set_target_properties(Pythia6::${_lib} PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${Pythia6_INCLUDE_DIR}")
  endforeach()
endif()

# ---------------------------------------------------------------------------
# Fortran runtime
# When C++ executables link static libs containing Fortran objects, CMake
# does not automatically pull in the Fortran runtime — list them explicitly.
# ---------------------------------------------------------------------------
set(NUWRO_FORTRAN_RUNTIME ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})

# ---------------------------------------------------------------------------
# Global compiler options
# ---------------------------------------------------------------------------
add_compile_options(
  -fPIC
  -Wall
  -Wno-unused-variable
  -Wno-sign-compare
  -Wno-unused-function
  -Wno-unused-but-set-variable
  -Wmissing-braces
)

if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
  add_link_options(-Wl,--no-as-needed)
endif()

add_compile_definitions(VERSION="${NUWRO_VERSION}")

# ---------------------------------------------------------------------------
# RPATH: keep full rpath in the build tree; trim to install prefix at install
# ---------------------------------------------------------------------------
set(CMAKE_SKIP_BUILD_RPATH        FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH           "${CMAKE_INSTALL_PREFIX}/lib")

# ---------------------------------------------------------------------------
# Sources
# ---------------------------------------------------------------------------
add_subdirectory(src)

# ---------------------------------------------------------------------------
# Generate setup.sh using configure_file so variables are expanded at
# configure time (not hardcoded paths at build time).
# ---------------------------------------------------------------------------
execute_process(
  COMMAND root-config --bindir
  OUTPUT_VARIABLE ROOT_BINDIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

if(EXISTS "${ROOT_BINDIR}/thisroot.sh")
  set(_thisroot_line "source ${ROOT_BINDIR}/thisroot.sh")
else()
  set(_thisroot_line "")
endif()

if(ROOTEGPythia6_ROOT)
  set(_pythia6_env
    "export LD_LIBRARY_PATH=${ROOTEGPythia6_ROOT}/lib:\${LD_LIBRARY_PATH}")
else()
  set(_pythia6_env "")
endif()

configure_file(
  "${CMAKE_SOURCE_DIR}/cmake/setup.sh.in"
  "${CMAKE_BINARY_DIR}/setup.sh"
  @ONLY
)
install(PROGRAMS "${CMAKE_BINARY_DIR}/setup.sh" DESTINATION .)
