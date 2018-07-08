#include "mex.h"
#include "gpu/mxGPUArray.h"

/*
 * Device code
 */


__global__ void MatMulKernel(double* A, double* B, double* C, dim3 dimsA, dim3 dimsB)//[1]  
{
	// Each thread computes one element of C  
	// by accumulating results into Cvalue  

	double Cvalue = 0;
	int row = blockIdx.y * blockDim.y + threadIdx.y;
	int col = blockIdx.x * blockDim.x + threadIdx.x;
	for (int e = 0; e < dimsA.y; ++e)
		Cvalue += A[row * dimsA.y + e] * B[e * dimsB.y + col];
	C[row * dimsB.y + col] = Cvalue;
}

void MatrixMultiplication_CUDA(const double* A, const double* B, double* C)
{
	dim3 dimsA(135, 135);// the size of matrix A which you have to modify  (4,3)
	dim3 dimsB(135, 3);// the size of matirx B which you have to modify   

					 //copy memory from host to devices  
	unsigned int size_A = dimsA.x * dimsA.y;
	unsigned int mem_size_A = sizeof(double) * size_A;
	double *d_A;
	cudaMalloc(&d_A, mem_size_A);
	cudaMemcpy(d_A, A, mem_size_A, cudaMemcpyHostToDevice);
	double *d_B;
	unsigned int size_B = dimsB.x * dimsB.y;
	unsigned int mem_size_B = sizeof(double) * size_B;
	cudaMalloc(&d_B, mem_size_B);
	cudaMemcpy(d_B, B, mem_size_B, cudaMemcpyHostToDevice);
	unsigned int mem_size_C = sizeof(double)* dimsA.x*dimsB.y;
	double *d_C;
	cudaMalloc(&d_C, mem_size_C);

	//dimBlock represents the threads'size within block which you have to modify[2]  
	dim3 dimBlock(3, 2);
	dim3 dimGrid(dimsB.y / dimBlock.x, dimsA.x / dimBlock.y);//[3]  

	MatMulKernel << <dimGrid, dimBlock >> >(d_A, d_B, d_C, dimsA, dimsB);
	// Read C from device memory  
	cudaMemcpy(C, d_C, mem_size_C,
		cudaMemcpyDeviceToHost);
	// Free device memory  
	cudaFree(d_A);
	cudaFree(d_B);
	cudaFree(d_C);

}





/*
 * Host code
 */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, mxArray const *prhs[])
{

    /*
     * Call the kernel using the CUDA runtime API. We are using a 1-d grid here,
     * and it would be possible for the number of elements to be too large for
     * the grid. For this example we are not guarding against this possibility.
     */
    //N = (int)(mxGPUGetNumberOfElements(A));
    //blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;
    //TimesTwo<<<blocksPerGrid, threadsPerBlock>>>(d_A, d_B, N);

	int M0 = mxGetM(prhs[0]);    //得到arr0的行数
    int N0 = mxGetN(prhs[0]);    //得到arr0的列数
 //   double* pArr0 = (double*)mxGetPr(prhs[0]);//得到arr0的指针

    int M1 = mxGetM(prhs[1]);
    int N1 = mxGetN(prhs[1]);
//    double* pArr1 = (double*)mxGetPr(prhs[1]);


	
	
	double AA[] = { 6,2,3,
		8,3,5,
		7,2,4,
		8.3,2,5 };
	double BB[] = { 1,2,3,
		4,5,6,
		7,8,9 };
		//AA=[6,2,3;8,3,5;7,2,4;8.3,2,5 ];BB= [ 1,2,3;4,5,6;7,8,9 ];

	double* pArr0= (double*)mxGetPr(prhs[0]);
	double* pArr1= (double*)mxGetPr(prhs[1]);
	double*CC = new double[135*3];
//	MatrixMultiplication_CUDA(AA, BB, CC);
	MatrixMultiplication_CUDA(pArr0, pArr1, CC);
	
	plhs[0] = mxCreateDoubleMatrix(3, 135, mxREAL);
	double* pRe =(double*)mxGetPr(plhs[0]);
	for(int i=0;i<3*135;i++)
	{
	pRe[i]=CC[i];
	}
	
    /* Wrap the result up as a MATLAB gpuArray for return. */
    //plhs[0] = mxGPUCreateMxArrayOnGPU(B);

    /*
     * The mxGPUArray pointers are host-side structures that refer to device
     * data. These must be destroyed before leaving the MEX function.
     */

}
