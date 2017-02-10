########################################################
# cmake file for building Marlin example Package
# @author Jan Engels, Desy IT
CMAKE_MINIMUM_REQUIRED(VERSION 2.6 FATAL_ERROR)
########################################################

# project name
PROJECT( PolyaFit )

### DEPENDENCIES ############################################################

FIND_PACKAGE( ILCUTIL REQUIRED COMPONENTS ILCSOFT_CMAKE_MODULES )

FIND_PACKAGE( Boost REQUIRED COMPONENTS  )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )

# load default settings from ILCSOFT_CMAKE_MODULES
INCLUDE( ilcsoft_default_settings )

FIND_PACKAGE( ROOT REQUIRED COMPONENTS TMVA XMLIO Minuit MLP MathMore MathCore Minuit2 RooFit RooFitCore) # minimum required Marlin version
INCLUDE_DIRECTORIES( ${ROOT_INCLUDE_DIRS} )
LINK_LIBRARIES( ${ROOT_LIBRARIES} )
ADD_DEFINITIONS( ${ROOT_DEFINITIONS} )

### LIBRARY #################################################################

# definitions to pass to the compiler
ADD_DEFINITIONS( "-Wall -ansi -pedantic" )
ADD_DEFINITIONS( "-Wno-long-long" )

#include directories
INCLUDE_DIRECTORIES( ./include )

#add library
AUX_SOURCE_DIRECTORY( ./src library_sources )

# display some variables and write them to cache
DISPLAY_STD_VARIABLES()

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cpp)

file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.h)


SET (Fit ./src/fit.cxx)
SET (FitSim ./src/fitSim.cxx)
SET (Draw ./src/draw.cxx)
SET (Compar ./src/comparArnaud.cxx)

#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../bin/")

ADD_EXECUTABLE( Fit ${Fit} ${headers} ${sources} )
ADD_EXECUTABLE( FitSim ${FitSim} ${headers} ${sources} )
ADD_EXECUTABLE( Draw ${Draw} ${headers} ${sources} )
ADD_EXECUTABLE( Compar ${Compar} ${headers} ${sources} )

TARGET_LINK_LIBRARIES( Fit ${ROOT_COMPONENT_LIBRARIES} ${Boost_LIBRARIES} )
TARGET_LINK_LIBRARIES( FitSim ${ROOT_COMPONENT_LIBRARIES} ${Boost_LIBRARIES} )
TARGET_LINK_LIBRARIES( Draw ${ROOT_COMPONENT_LIBRARIES} ${Boost_LIBRARIES} )
TARGET_LINK_LIBRARIES( Compar ${ROOT_COMPONENT_LIBRARIES} ${Boost_LIBRARIES} )

INSTALL(TARGETS Fit RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/bin)
INSTALL(TARGETS FitSim RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/bin)
INSTALL(TARGETS Draw RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/bin)
INSTALL(TARGETS Compar RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/bin)
