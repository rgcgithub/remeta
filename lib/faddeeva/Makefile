#
# Makefile to compile PGEN library
#
AR         ?= ar
CXX        ?= g++
CXXFLAGS    = -O3 -Wall -std=c++11
CFLAGS      = 

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
 CXXFLAGS += -arch x86_64
endif

faddeeva.a: Faddeeva.o
	${AR} rcs $@ $^ 

%.o: %.cc
	${CXX} ${CXXFLAGS} -o $@ -c $<

clean:
	rm -f *.o *.a
