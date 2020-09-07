# SynCOM makefile with output SynCOM.exe

# Set dependency directories
VPATH = ../src/
INC = ../include/

LFLAGS = -shared -static -static-libstdc++ 

CFLAGS = -c -O3 -g -w -Wall -static -std=gnu++0x -static-libstdc++ -I$(INC)

CC = g++

all: SYNCOM_API.dll
 
SYNCOM_API.dll: SynCOM_API.o setting.o readIn_api.o error.o stressSolver_api.o strainSolver_api.o printOut_api.o
	$(CC) $(LFLAGS) -o SynCOM_API.dll SynCOM_API.o setting.o readIn_api.o error.o \
				stressSolver_api.o strainSolver_api.o printOut_api.o 

SynCOM_API.o: SynCOM_API.h SynCOM_API.cpp
	$(CC) $(CFLAGS) $(VPATH)SynCOM_API.cpp

error.o: error.h error.cpp
	$(CC) $(CFLAGS) $(VPATH)error.cpp

printOut_api.o: printOut_api.h printOut_api.cpp strainSolver_api.h  strainSolver_api.cpp \
		stressSolver_api.h stressSolver_api.cpp error.h error.cpp
	$(CC) $(CFLAGS) $(VPATH)printOut_api.cpp

readIn_api.o: readIn_api.h readIn_api.cpp setting.h setting.cpp error.h error.cpp 
	$(CC) $(CFLAGS) $(VPATH)readIn_api.cpp

setting.o: setting.h setting.cpp error.h error.cpp
	$(CC) $(CFLAGS) $(VPATH)setting.cpp

strainSolver_api.o: strainSolver_api.h strainSolver_api.cpp setting.h setting.cpp error.h error.cpp
	$(CC) $(CFLAGS) $(VPATH)strainSolver_api.cpp

stressSolver_api.o: stressSolver_api.h stressSolver_api.cpp setting.h setting.cpp error.h error.cpp
	$(CC) $(CFLAGS) $(VPATH)stressSolver_api.cpp

clean:
	-del *.o

