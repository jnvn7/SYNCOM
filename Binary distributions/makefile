# SynCOM makefile with output SynCOM.exe

# Set dependency directories
VPATH = src/
INC = include/

LFLAGS = -static -static-libstdc++ 

CFLAGS = -c -O3 -g -w -Wall -static -std=gnu++0x -static-libstdc++ -I$(INC)

CC = g++

all: SYNCOM
 
SYNCOM: SynCOM.o setting.o readIn.o error.o stressSolver.o strainSolver.o printOut.o
	$(CC) $(LFLAGS) -o SynCOM SynCOM.o setting.o readIn.o error.o stressSolver.o strainSolver.o printOut.o 

SynCOM.o: SynCOM.cpp
	$(CC) $(CFLAGS) $(VPATH)SynCOM.cpp

error.o: error.h error.cpp
	$(CC) $(CFLAGS) $(VPATH)error.cpp

printOut.o: printOut.h printOut.cpp strainSolver.h  strainSolver.cpp stressSolver.h \
		stressSolver.cpp error.h error.cpp
	$(CC) $(CFLAGS) $(VPATH)printOut.cpp

readIn.o: readIn.h readIn.cpp setting.h setting.cpp error.h error.cpp resource.h
	$(CC) $(CFLAGS) $(VPATH)readIn.cpp

setting.o: setting.h setting.cpp error.h error.cpp
	$(CC) $(CFLAGS) $(VPATH)setting.cpp

strainSolver.o: strainSolver.h strainSolver.cpp setting.h setting.cpp error.h error.cpp
	$(CC) $(CFLAGS) $(VPATH)strainSolver.cpp

stressSolver.o: stressSolver.h stressSolver.cpp setting.h setting.cpp error.h error.cpp
	$(CC) $(CFLAGS) $(VPATH)stressSolver.cpp

clean:
	-del *.o

