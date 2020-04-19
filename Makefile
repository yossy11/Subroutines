CXX=gfortran
CXXFLAGS=-Iinclude -Wall
DEBUG_CXXFLAGS=$(CXXFLAGS) -g
RELEASE_CXXFLAGS=-Wall -DENDEBUG -O2
LDFLAGS=-lm

all: debug

TARGETS=Subroutine
DEBUG_TARGETS=$(TARGETS:%=Debug/%)

debug:Debug $(DEBUG_TARGETS)

Debug:
	mkdir Debug

Debug/%.o : src/%.for
	$(CXX) -c $< -o $@

OBJS1=Subroutine.o
Debug/test : $(OBJS1:%=Debug/%)
	$(CXX) -o $@ $^ $(LDFLAGS)

clean: 
	rm -rf Debug
