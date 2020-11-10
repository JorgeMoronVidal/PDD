# Date Created 2013-01-02
# Adapted from http://www.cs.swarthmore.edu/~newhall/unixhelp/howto_makefiles.html
# License MIT
#
# 'make depend' uses makedepend to automatically generate dependencies 
#               (dependencies are added to end of Makefile)
# 'make'        build executable file 'mycc'
# 'make clean'  removes all .o and executable files
#

# define the C/C++ compiler to use
CC = mpicxx

# define compile-time flags NOTE: last two flags required for C++ 11 spec on osx
CFLAGS = -g -Wall -fopenmp

# define any directories containing header files other than /usr/include
INCLUDES = -I/usr/local/include

# define library paths in addition to /usr/lib
#   if I wanted to include libraries not in /usr/lib specify
#   their path using -Lpath, something like:
LFLAGS = -L/usr/local/lib

# define any libraries to link into executable:
#   To ling libusb-1.0 :
#   LIBS = -lusb-1.0
LIBS = -lgsl -lgslcblas -lm 

# define the source files, default here is main.cpp
SRCS = main.cpp PDDSparseGM.cpp interface.cpp subdomain.cpp BVP.cpp stencil.cpp scalarfunction.cpp scalarfunctionN.cpp vectorfunction.cpp matrixfunction.cpp rbfunction.cpp boundary.cpp lookuptable.cpp poissontest.cpp rectangle.cpp FKACSolver.cpp GMSolver.cpp

# define the object files 
#
# This uses Suffix Replacement within a macro:
#   $(name:string1=string2)
#         For each word in 'name' replace 'string1' with 'string2'
# Below we are replacing the suffix .c of all words in the macro SRCS
# with the .o suffix
#
OBJS = $(SRCS:.c=.o)

# define the executable file  name
MAIN =  main



#
# The following part of the makefile is generic; it can be used to 
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
#

.PHONY: depend clean

all:    $(MAIN)
	@echo  Compiled!

$(MAIN): $(OBJS) 
	$(CC) $(CFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) $(LFLAGS) $(LIBS)

# this is a suffix replacement rule for building .o's from .c's
# it uses automatic variables $<: the name of the prerequisite of
# the rule(a .cpp file) and $@: the name of the target of the rule (a .o file) 
# (see the gnu make manual section about automatic variables)
.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@

clean:
	$(RM) *.o *~ $(MAIN)

depend: $(SRCS)
	makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it

