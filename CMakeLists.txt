cmake_minimum_required( VERSION 3.20.3 )                                    # Using the current cmake version as a reqiremnet. Not really matter

project( tadpole )															# Project name, can be used as target name too

file( GLOB SOURCES ${PROJECT_SOURCE_DIR}/src/*.cpp )						# Specify all of the source files in the src subdirectory in SOURCE variable

include_directories( "inc" )												# All header files of this project are placed in the inc subdirectory
include_directories( "../include/boost_1_77_0" )							# Include the path of BGL libraries and all the ohter boost libraries.
include_directories( "/opt/homebrew/include" )							

set( CMAKE_CXX_STANDARD 17 )												# Using c++17 standard, it's like -std=c++17 option
set( CMAKE_CXX_STANDARD_REQUIRED ON )

set( CMAKE_BUILD_TYPE Debug )												# Specify debug flag, it's like -g option

set( CMAKE_CXX_FLAGS "-L/opt/homebrew/lib -lxlsxwriter" )				# Specify some other flags

set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin )				# Setting the output directory of the generated executable to the "bin"
																			# directory which is a subdirectory of the current directory the
																			# CMakeLists.txt file is placed in

find_package(fmt)															

# Build a executable target (1st param: Target name; others: source files)

add_executable( 
    test_dirtyarchieve 
    ${PROJECT_SOURCE_DIR}/src/test_dirtyarchieve.cpp 
)						
target_link_libraries( test_dirtyarchieve PRIVATE fmt::fmt-header-only )       

add_executable( 
    examine_test_case
    ${PROJECT_SOURCE_DIR}/src/examine_test_case.cpp
)						
target_link_libraries( examine_test_case PRIVATE fmt::fmt-header-only )       

add_executable( 
    test_outter_routing
    ${PROJECT_SOURCE_DIR}/src/test_outter_routing.cpp 
)						
target_link_libraries( test_outter_routing PRIVATE fmt::fmt-header-only )       

add_executable( 
    test_sample1
    ${PROJECT_SOURCE_DIR}/src/test_sample1.cpp 
)						
target_link_libraries( test_sample1 PRIVATE fmt::fmt-header-only )       

add_executable( 
    test_heuristic1_59
    ${PROJECT_SOURCE_DIR}/src/test_heuristic1_59.cpp 
)						
target_link_libraries( test_heuristic1_59 PRIVATE fmt::fmt-header-only )       

add_executable( 
    test_heuristic1_99
    ${PROJECT_SOURCE_DIR}/src/test_heuristic1_99.cpp 
)						
target_link_libraries( test_heuristic1_99 PRIVATE fmt::fmt-header-only )       

add_executable( 
    test_heuristic1_256
    ${PROJECT_SOURCE_DIR}/src/test_heuristic1_256.cpp 
)						
target_link_libraries( test_heuristic1_256 PRIVATE fmt::fmt-header-only )       

add_executable( 
    test_heuristic1_500
    ${PROJECT_SOURCE_DIR}/src/test_heuristic1_500.cpp 
)						
target_link_libraries( test_heuristic1_500 PRIVATE fmt::fmt-header-only )       

add_executable( 
    experiment 
    ${PROJECT_SOURCE_DIR}/src/experiment.cpp 
)						
target_link_libraries( experiment PRIVATE fmt::fmt-header-only )       

add_executable( 
    demo_1 
    ${PROJECT_SOURCE_DIR}/src/demo_1.cpp 
)						
add_executable( 
    demo_2 
    ${PROJECT_SOURCE_DIR}/src/demo_2.cpp 
)						
add_executable( 
    demo_3 
    ${PROJECT_SOURCE_DIR}/src/demo_3.cpp 
)						
add_executable( 
    demo_4 
    ${PROJECT_SOURCE_DIR}/src/demo_4.cpp 
)						
add_executable( 
    demo_5 
    ${PROJECT_SOURCE_DIR}/src/demo_5.cpp 
)						

add_executable( 
    experiment2 
    ${PROJECT_SOURCE_DIR}/src/experiment2.cpp 
)						
target_link_libraries( experiment2 PRIVATE fmt::fmt-header-only )       

add_executable( 
    demo_6 
    ${PROJECT_SOURCE_DIR}/src/demo_6.cpp 
)						
target_link_libraries( demo_6 PRIVATE fmt::fmt-header-only )       

