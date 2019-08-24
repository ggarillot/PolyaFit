########################################################
# cmake file for building Marlin example Package
# @author Jan Engels, Desy IT
CMAKE_MINIMUM_REQUIRED(VERSION 2.6 FATAL_ERROR)
########################################################

# project name
PROJECT( PolyaFit )

FIND_PACKAGE( ROOT REQUIRED COMPONENTS XMLIO Minuit MathMore MathCore Minuit2 )
INCLUDE_DIRECTORIES( SYSTEM ${ROOT_INCLUDE_DIRS} )
LINK_LIBRARIES( ${ROOT_LIBRARIES} )
ADD_DEFINITIONS( ${ROOT_DEFINITIONS} )

### LIBRARY #################################################################

# definitions to pass to the compiler
ADD_DEFINITIONS( "-Wall -ansi -pedantic -std=c++11 -O2 -Weffc++" )
ADD_DEFINITIONS( "-Wno-long-long" )

# include directories
INCLUDE_DIRECTORIES( ./include )

# add library
AUX_SOURCE_DIRECTORY( ./src library_sources )


file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cpp)

file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.h ${PROJECT_SOURCE_DIR}/include/*.hpp)

add_library( MyLib ${sources} ${headers} )


file( GLOB APP_SOURCES ${PROJECT_SOURCE_DIR}/exec/*.cxx )
foreach( testsourcefile ${APP_SOURCES} )
	get_filename_component(testname ${testsourcefile} NAME_WE)
	add_executable( ${testname} ${testsourcefile} )
	target_link_libraries( ${testname} MyLib ${ROOT_COMPONENT_LIBRARIES} )
	install(TARGETS ${testname} RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/bin)
endforeach( testsourcefile ${APP_SOURCES} )
