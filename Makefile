
CFLAGS	= -g 
EXE	= Laplacian3D
CC	= g++
NVCC	= g++
INC	= 
LIBS	= 
NVARCH	= 
NVFLAGS	= $(CFLAGS) 
TARGS	= main.o setAnalysis.o setModelData.o

all:	$(EXE)

%.o:	%.cpp
	$(NVCC) $(NVFLAGS) -c $? 

$(EXE): $(TARGS) 
	$(NVCC) $(NVFLAGS) $(TARGS) -o $@

clean:
	rm -rf *.o
