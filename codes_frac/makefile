#export CFLAGS=`root-config --cflags`
#export LDFLAGS=`root-config --ldflags --libs --lRooFit --lRooFitCore -lMinuit`

ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)
ROOTGLIBS = $(shell root-config --glibs)
CXXFLAGS += $(ROOTCFLAGS)
LIBS = $(ROOTLIBS)
GLIBS = $(ROOTGLIBS)
CODEDIR = ./
CURRDIR = $(shell pwd)
INCLUDES = $(CXXFLAGS) -D USEROOT -D ROO_COMPLEX -I$(CODEDIR) -D FIXPHI
OBJS = pwa.o DPFPWAPoint.o DPFPWAPdf.o DPFAngular.o DPFPropogator.o common_tools.o phikk_structure.o phipipi_structure.o fitproxy.o PWA_CTRL.o gpuCode.o
EXEC = pwa

$(EXEC):$(OBJS)
	g++ -o $@ $^ `root-config --libs` -L`root-config --libdir` -lRooFit -D ROO_COMPLEX -fopenmp  -lgomp -lpthread -lRooFitCore -lcudart -L/usr/local/cuda/lib64

gpuCode.o:kernel_calEva.o cu_DPFPropogator.o
	nvcc -arch=sm_20 -dlink $^ -o gpuCode.o

pwa.o:pwa.C
	g++ -c $< $(INCLUDES) -D PWAINI="$(CURRDIR)/pwa_ini.C"

cu_DPFPropogator.o:cu_DPFPropogator.cu
	nvcc -dc $<

kernel_calEva.o:kernel_calEva.cu
	nvcc -dc $<

DPFPWAPdf.o:$(CODEDIR)/DPFPWAPdf.cc $(CODEDIR)/DPFPWAPdf.h
	g++ -c $< $(INCLUDES) -fopenmp `root-config --libs` -L`root-config --libdir` -lRooFit -lRooFitCore

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
