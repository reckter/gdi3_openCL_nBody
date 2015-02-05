#eated with genmake.pl v1.1 on Wed Aug  8 16:08:07 2012
# genmake.pl home: http://muquit.com/muquit/software/
# Copryright: GNU GPL (http://www.gnu.org/copyleft/gpl.html)
##----------------------------------------------------------------------------
rm=/bin/rm -f
CC=g++
DEFS=  
PROGNAME= nBody
INCLUDES=-I. -I../common 
LIBS=-lOpenCL


DEFINES=$(INCLUDES) $(DEFS)
CFLAGS=-g -O3 -fopenmp $(DEFINES)

SRCS = main.cpp OpenCLManager.cpp

OBJS = main.o OpenCLManager.o


.c.o:
	$(rm) $@
	$(CC) $(CFLAGS) -c $*.cpp

all: $(PROGNAME)

$(PROGNAME) : $(SRCS) 
	$(CC) $(CFLAGS) $(SRCS) -o $(PROGNAME) $(LIBS)

clean:
	$(rm) $(OBJS) $(PROGNAME) core *~
