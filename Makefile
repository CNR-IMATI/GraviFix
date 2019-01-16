egg_TARGET := gravifix
egg_SOURCES := Gravifix.cpp 

ImatiSTLdir := ./ImatiSTL-4.2-2



DEFINES := -DIS64BITPLATFORM 
INCLUDES := -I$(ImatiSTLdir)/include/ImatiSTL  -I$(ImatiSTLdir)/include/Kernel  -I$(ImatiSTLdir)/include/TMesh

egg_SOURCE_PATH := ./ 

egg_FLAGS := 
egg_RELEASE_FLAGS := -O3
egg_DEBUG_FLAGS := -g

egg_CFLAGS := 
egg_CXXFLAGS := -fpermissive -std=c++11 -Wno-write-strings -Wno-int-to-pointer-cast $(DEFINES) $(INCLUDES)


egg_LDFLAGS :=  -L$(ImatiSTLdir)/make/
egg_LDLIBS :=  -lImatiSTL

include ./eggmakelib/engine.mk
