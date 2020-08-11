cmake_minimum_required(VERSION 3.10)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Set compiler for unix system
if(UNIX)
	set(CMAKE_C_COMPILER "gcc")
	set(CMAKE_CXX_COMPILER "g++")
endif()

# set output folder for *.exe or *.dll
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# set the project name
project(SYNCOM)

# Compile dynamic linked library or executables
set(SYNCOM_API 1)
	
# Set include library path
include_directories("${PROJECT_SOURCE_DIR}/include")

# Add the excutable from the src folder
if (SYNCOM_API)
	add_library(SYNCOM_API SHARED
		src/error.cpp
		src/printOut_api.cpp
		src/readIn_api.cpp
		src/setting.cpp
		src/strainSolver_api.cpp
		src/stressSolver_api.cpp
		src/SYNCOM_API.cpp
		include/rapidxml-1.13/rapidxml_print.hpp
		include/rapidxml-1.13/rapidxml_utils.hpp
		include/rapidxml-1.13/rapidxml.hpp
		)
# Add the excutable from the src folder		
else()
	add_executable(SYNCOM 
		src/error.cpp
		src/printOut.cpp
		src/readIn.cpp
		src/setting.cpp
		src/strainSolver.cpp
		src/stressSolver.cpp
		src/SYNCOM.cpp
		include/rapidxml-1.13/rapidxml_print.hpp
		include/rapidxml-1.13/rapidxml_utils.hpp
		include/rapidxml-1.13/rapidxml.hpp
		)
endif()