
#include <iostream>
#include "simple.h"
// #include <cuda.h>

// __global__ void add(int a, int b, int *c)
// {
// 	*c = a + b; 
// }

int test_cuda(){
	// cudaSetDevice (2);
	int a,b,c;
	int *dev_c;
	a=3;
	b=4;
	std::cout << "\n\n\n\nhello world\n\n\n" << std::endl;
	// cudaSetDevice (2);
	// cudaMalloc((void**)&dev_c, sizeof(int)); 
	// add<<<1,1>>>(a,b,dev_c);
	// cudaMemcpy(&c, dev_c, sizeof(int), cudaMemcpyDeviceToHost); 
	// printf("%d + %d :is: %d\n", a, b, c);
	// cudaFree(dev_c);
	return 0;
}