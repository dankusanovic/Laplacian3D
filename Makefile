
CFLAGS	= -g 
EXE	= Laplacian3D
CC	= g++
NVCC	= g++
INC	= -I /usr/include/mumps_seq
LPATH	= -L/usr/lib/gcc/x86_64-linux-gnu/4.8 
LIBS	= -lgfortran -lblas -lpthread -lmpiseq_seq -lmumps_common_seq -ldmumps_seq 
NVARCH	= 
NVFLAGS	= $(CFLAGS) 
TARGS	= Main.o setAnalysis.o setModelData.o allocateModelData.o setModelDofs.o setForceVector.o setStiffnessMatrix.o setBoundaries.o MUMPSSolver.o freeModelData.o 

all:	$(EXE)

%.o:	%.cpp
	$(NVCC) $(NVFLAGS) $(INC) -c $? 

$(EXE): $(TARGS) 
	$(NVCC) $(NVFLAGS) $(LPATH) $(LIBS) $(TARGS) -o $@

clean:
	rm -rf *.o
