NVCC= nvcc
NVCCFLAGS= -O3 -std=c++14 -arch=sm_70
CXXFLAGS= -Xcompiler="-fopenmp -march=native"

all: main

main: main.cu
	$(NVCC) $(NVCCFLAGS) $(CXXFLAGS) main.cu -o main

clean:
	rm -f main
