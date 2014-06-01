# name of the library
LIBNAME = EnvelopeCode

#Necessary to use shell built-in commands
SHELL=bash

# figure out os
UNAME := $(shell uname)
BASEDIR = $(shell pwd)

$(shell mkdir -p lib)
$(shell mkdir -p obj)

USERINCLUDES += -I$(ROOTSYS)/include
USERINCLUDES += -I$(ROOFITSYS)/include
USERINCLUDES += -I$(BOOSTPATH)/include
USERLIBS += -L$(ROOFITSYS)/lib -lRooFit -lRooFitCore -lRooStats
USERLIBS += $(shell root-config --glibs) -lTMVA -lMinuit -lFoam -lRooFit -lRooFitCore -lMathMore -lRooStats
USERLIBS += -L$(BOOSTPATH)/lib -lboost_regex -lboost_program_options
CXXFLAGS = -Wall -g -Wl,-rpath=$(BASEDIR)/lib 
LDFLAGS = -shared -W 

CXX=g++
LD=g++

INCLUDE +=  -I $(BASEDIR)/interface
INCLUDE += $(USERINCLUDES)
CXXFLAGS += $(INCLUDE)

LIBS += $(USERLIBS)

# this is where you  would normaly have the lib, bin and source directories
LIBDIR = $(BASEDIR)/lib
EXEDIR = $(BASEDIR)/bin
MACRODIR = $(BASEDIR)/src
SRCDIR = $(BASEDIR)/src
OBJDIR = $(BASEDIR)/obj
INCDIR = $(BASEDIR)/interface
TESTDIR = $(BASEDIR)/test
DOCDIR= $(BASEDIR)/docs
OBJ_EXT=o
TEST_EXT=cpp

ROOT_DICT=obj/RootDict.cxx
ROOT_OBJ=$(subst cxx,$(OBJ_EXT),$(ROOT_DICT))
SRCS=$(wildcard $(BASEDIR)/src/*.cc)
EXES=$(wildcard $(BASEDIR)/test/*.cpp)
HEADERS=$(wildcard $(BASEDIR)/interface/*.h)
OBJS=$(subst $(SRCDIR), $(OBJDIR),$(subst cc,$(OBJ_EXT),$(SRCS)))
BINS=$(subst $(TESTDIR), $(EXEDIR),$(subst .$(TEST_EXT),,$(EXES)))

all: $(LIBDIR)/lib$(LIBNAME).so

$(OBJDIR)/%.$(OBJ_EXT): $(SRCDIR)/%.cc $(ROOT_DICT) 
	@echo Making object $@
	@$(CXX) $(CXXFLAGS) -fPIC -c $<  -o $@

$(LIBDIR)/lib$(LIBNAME).so: $(OBJS) $(ROOT_OBJ)
	@echo Building shared library $@
	@$(LD) $(LDFLAGS) -o $(LIBDIR)/lib$(LIBNAME).so $^ $(LIBS) 
	
$(ROOT_DICT): $(SRCS) $(HEADERS)
	@echo Making dictionary $@
	@rootcint -f $@ -c -L$(ROOFITSYS)/lib -I$(ROOFITSYS)/include $(HEADERS)
	@$(CXX) $(CXXFLAGS) -fPIC -c $(ROOT_DICT) -o $(ROOT_OBJ)

vars:
	@echo "LIBS: " $(LIBS)
	@echo "CXXFLAGS: " $(CXXFLAGS)
	@echo "Header files: " $(HEADERS)
	@echo "Source files: " $(SRCS) 
	@echo "Root H files: " $(ROOTHEADERS)
	@echo "Root S files: " $(ROOTCLASSES) 
	@echo "Object files: " $(OBJS)
	@echo "Executables:  " $(TARGETS)

clean:
	@echo "Cleaning: " 
	@echo "  " $(OBJS)
	@echo "  " $(LIBDIR)/lib$(LIBNAME).so
	@echo "  " $(ROOT_DICT)
	@rm -rf $(OBJS) $(LIBDIR)/lib$(LIBNAME).so $(BINS) $(ROOT_DICT)


