# makefile for MoorDyn-Syncom compiled as MoorDynSC.dll

# search for dependencies in directories
VPATH = ../src/
INC = ../include/

LFLAGS = -shared -static -static-libgcc -static-libstdc++ -lws2_32 -DMoorDyn_EXPORTS

CFLAGS = -c -O3 -g -w -Wall -static -static-libgcc -static-libstdc++ -std=gnu++0x -DMoorDyn_EXPORTS -DUSEGL \
		-I$(INC)


all: MoorDynSC.dll
 
MoorDynSC.dll: MoorDyn.o Line.o Connection.o Misc.o kiss_fft.o \
		SynCOM.o SC_readIn_api.o SC_error.o SC_stressSolver_api.o
	g++ $(LFLAGS) -o MoorDynSC.dll MoorDyn.o Line.o Connection.o Misc.o kiss_fft.o \
		SynCOM.o SC_readIn_api.o SC_error.o SC_stressSolver_api.o -lopengl32

MoorDyn.o: MoorDyn.cpp MoorDyn.h Line.h Line.cpp Connection.h Connection.cpp QSlines.h Misc.h Misc.cpp \
		SynCOM.h SynCOM.cpp SC_readIn_api.h SC_readIn_api.cpp SC_stressSolver_api.h SC_stressSolver_api.cpp
	g++ $(CFLAGS) $(VPATH)MoorDyn.cpp
	
kiss_fft.o: kiss_fft.h kiss_fft.c
	g++ $(CFLAGS) $(VPATH)kiss_fft.c
	
Line.o: Line.h Line.cpp Connection.h Connection.cpp QSlines.h Misc.h
	g++ $(CFLAGS) $(VPATH)Line.cpp

Connection.o: Line.h Line.cpp Connection.h Connection.cpp QSlines.h Misc.h Misc.cpp
	g++ $(CFLAGS) $(VPATH)Connection.cpp
	
Misc.o: Misc.h Misc.cpp
	g++ $(CFLAGS) $(VPATH)Misc.cpp

SynCOM.o: SynCOM.h SynCOM.cpp SC_readIn_api.h SC_readIn_api.cpp SC_error.h SC_error.cpp \
		SC_stressSolver_api.h SC_stressSolver_api.cpp 
	g++ $(CFLAGS) $(VPATH)SynCOM.cpp

SC_readIn_api.o: SC_readIn_api.h SC_readIn_api.cpp SC_error.h SC_error.cpp 
	g++ $(CFLAGS) $(VPATH)SC_readIn_api.cpp

SC_error.o: SC_error.h SC_error.cpp
	g++ $(CFLAGS) $(VPATH)SC_error.cpp

SC_stressSolver_api.o: SC_stressSolver_api.h SC_stressSolver_api.cpp SC_error.h SC_error.cpp
	g++ $(CFLAGS) $(VPATH)SC_stressSolver_api.cpp

clean:
	-del *.o

