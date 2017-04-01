#export CFLAGS=`root-config --cflags`
#export LDFLAGS=`root-config --ldflags --libs --lRooFit --lRooFitCore -lMinuit`

ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)
ROOTGLIBS = $(shell root-config --glibs)
CXXFLAGS += $(ROOTCFLAGS)
MODE=GPU
LIBS = $(ROOTLIBS)
GLIBS = $(ROOTGLIBS)
CODEDIR = ./
CURRDIR = $(shell pwd)
INCLUDES = $(CXXFLAGS) -D USEROOT -D ROO_COMPLEX -I$(CODEDIR) -D FIXPHI
OBJS = pwa.o DPFPWAPoint.o DPFPWAPdf.o DPFAngular.o DPFPropogator.o common_tools.o phikk_structure.o phipipi_structure.o fitproxy.o PWA_CTRL.o

ifeq ($(MODE),GPU)
	OBJS+=gpuCode.o cu_DPFPropogator.o kernel_calEva.o
endif

EXEC = pwa

$(EXEC):$(OBJS)
	g++ -o $@ $^ `root-config --libs` -L`root-config --libdir` -lRooFit -D ROO_COMPLEX -fopenmp  -lgomp -lpthread -lRooFitCore -lcudart -L/usr/local/cuda/lib64

ifeq ($(MODE),GPU)
gpuCode.o:kernel_calEva.o cu_DPFPropogator.o
	nvcc $^ -o gpuCode.o -arch=sm_20 -dlink
endif

pwa.o:pwa.C
	g++ -c $< $(INCLUDES) -D PWAINI="$(CURRDIR)/pwa_ini.C"

ifeq ($(MODE),GPU)
cu_DPFPropogator.o:cu_DPFPropogator.cu cu_DPFPropogator.h 
	nvcc -dc $<

kernel_calEva.o:kernel_calEva.cu cu_DPFPropogator.h kernel_calEva.h
	nvcc -dc $<
endif

DPFPWAPdf.o:$(CODEDIR)/DPFPWAPdf.cc $(CODEDIR)/DPFPWAPdf.h $(CODEDIR)/cu_PWA_PARAS.h $(CODEDIR)/kernel_calEva.h
	g++ -c $< $(INCLUDES) -D $(MODE) -fopenmp `root-config --libs` -L`root-config --libdir` -lRooFit -lRooFitCore

DPFPWAPoint.o:$(CODEDIR)/DPFPWAPoint.cc $(CODEDIR)/DPFPWAPoint.h
	g++ -c $< $(INCLUDES)

DPFAngular.o:$(CODEDIR)/DPFAngular.cc $(CODEDIR)/DPFAngular.h
	g++ -c $< $(INCLUDES)

DPFPropogator.o:$(CODEDIR)/DPFPropogator.cc $(CODEDIR)/DPFPropogator.h
	g++ -c $< $(INCLUDES)

common_tools.o:$(CODEDIR)/common_tools.cc $(CODEDIR)/common_tools.h
	g++ -c $< $(INCLUDES)

phikk_structure.o:$(CODEDIR)/phikk_structure.cc $(CODEDIR)/phikk_structure.h
	g++ -c $< $(INCLUDES)

phipipi_structure.o:$(CODEDIR)/phipipi_structure.cc $(CODEDIR)/phipipi_structure.h
	g++ -c $< $(INCLUDES)

fitproxy.o:$(CODEDIR)/fitproxy.cc $(CODEDIR)/fitproxy.h
	g++ -c $< $(INCLUDES)

PWA_CTRL.o:$(CODEDIR)/PWA_CTRL.C $(CODEDIR)/PWA_CTRL.H
	g++ -c $< $(INCLUDES)

clean:
	rm -f $(EXEC) $(OBJS)
