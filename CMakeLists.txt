# CMakeLists.txt for event package. It creates a library with dictionary and a main program
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(event)

# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS  Core Cint RIO Net Hist Graf Graf3d Gpad Tree Rint Postscript Matrix Physics MathCore Thread Gui pthread)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

include_directories(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
add_definitions(${ROOT_CXX_FLAGS})

include_directories(include)
include_directories(InputTree)

SET(CMAKE_CXX_FLAGS "-std=c++0x")

set(INCLUDE
    include/CentralityDetectorEvent.h     
    include/CentralityEventContainer.h
    include/CentralitySlicesFinder.h
    include/CentralityGetter.h
    include/NA61DataEvent.h
    include/NA61DataEventNew.h
    include/CentralityManager.h
    include/CentralitySlice.h
    include/TreeInterface.h
    include/CentralityContainerNormalizer.h
    Glauber/GlauberFitter.h 
    
    Glauber/AliGlauberMC.h
    Glauber/AliGlauberNucleon.h     
    Glauber/AliGlauberNucleus.h  

    InputTree/DataTreeEvent.h
    InputTree/DataTreeMCTrack.h
    InputTree/DataTreePSDModule.h
    InputTree/DataTreePSDSection.h
    InputTree/DataTreeTrack.h
    InputTree/DataTreeTOFSegment.h
    InputTree/DataTreeV0Candidate.h    
    InputTree/DataTreeTrigger.h
    InputTree/DataTreeBPD.h
    )

ROOT_GENERATE_DICTIONARY(G__Centrality ${INCLUDE} LINKDEF EventLinkDef.h)

#---Create a shared library with geneated dictionary

set(SOURCES 
    src/CentralityDetectorEvent.cxx 
    src/CentralityEventContainer.cxx 
    src/CentralitySlicesFinder.cxx
    src/CentralityGetter.cxx
    src/NA61DataEvent.cxx
    src/NA61DataEventNew.cxx
    src/CentralityManager.cxx
    src/CentralitySlice.cxx
    src/TreeInterface.cxx
    src/CentralityContainerNormalizer.cxx
 
    Glauber/GlauberFitter.cxx 

    Glauber/AliGlauberMC.cxx
    Glauber/AliGlauberNucleon.cxx     
    Glauber/AliGlauberNucleus.cxx  
    
    
    InputTree/DataTreeEvent.cxx
    InputTree/DataTreeMCTrack.cxx
    InputTree/DataTreePSDModule.cxx
    InputTree/DataTreePSDSection.cxx
    InputTree/DataTreeTrack.cxx
    InputTree/DataTreeTOFSegment.cxx
    InputTree/DataTreeV0Candidate.cxx 
    InputTree/DataTreeTrigger.cxx
    InputTree/DataTreeBPD.cxx
    )




# add_library(Glauber SHARED ${GLAUBER_SOURCES} )
# target_link_libraries(Glauber ${ROOT_LIBRARIES})    
    
add_library(Centrality SHARED ${SOURCES} G__Centrality.cxx)
target_link_libraries(Centrality ${ROOT_LIBRARIES})


#---Create  a main program using the library
add_executable(Main main.C)
target_link_libraries(Main Centrality)
