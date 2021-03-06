
#include "mex.h"    //必须有这个
#include <Eigen/Sparse>

//y===mcount(x,L)

//调用形式 re=SUM(arr0,arr1),将两个矩阵相加赋值给结果矩阵。
//nlhs：输出参数个数
//plhs：输出参数列表
//nrhs：输入参数个数
//prhs：输入参数列表
void mexFunction(int nlhs,mxArray *plhs[], int nrhs,const mxArray *prhs[])    //相当于一般的main()了
{
    int M0 = mxGetM(prhs[0]);    //得到arr0的行数
    int N0 = mxGetN(prhs[0]);    //得到arr0的列数
    double* pArr0 = (double*)mxGetPr(prhs[0]);    //得到arr0的指针

    int M1 = mxGetM(prhs[1]);
    int N1 = mxGetN(prhs[1]);
    double* pArr1 = (double*)mxGetPr(prhs[1]);

	mwIndex *ir=mxGetIr(prhs[0]);
	mwIndex *jc=mxGetJc(prhs[0]);
	double *pr=	mxGetPr(prhs[0]);
	
   // if(M0!=N0||M1!=N1)
    //    mexErrMsgTxt("两个矩阵行列应该相等");

    plhs[0] = mxCreateDoubleMatrix(M0, N1, mxREAL);    //创建一个M0行，N0列的矩阵
    double* pRe =(double*)mxGetPr(plhs[0]);
	

	for(int i=0;i<100;i++)
	{
		pRe[i]=ir[i];
		pRe[M0+i]=jc[i];
		pRe[2*M0+i]=pr[i];
	}

/*    for(int i=0;i<M0;i++)
    {
        for (int j=0;j<N0;j++)
        {
            //pRe[i*N0+j]=pArr0[i*N0+j]+pArr1[i*N0+j];    //两个矩阵逐个相加给结果矩阵
			pRe[i*N0+j]=pArr0[i*N0+j];
        }
    }*/	
	

/*	
	for(int i=0;i<M0;i++)
    {
        for (int j=0;j<N1;j++)
        {
            //pRe[i*N0+j]=pArr0[i*N0+j]+pArr1[i*N0+j];    //两个矩阵逐个相加给结果矩阵
			//pRe[i*N0+j]=pArr0[i*N0+j];
			pRe[j*M0+i]=0;
			for (int k=0;k<N0;k++)
			{
				pRe[j*M0+i]+=pArr0[k*M0+i]*pArr1[j*M1+k];
				//pRe[i*N1+j]+=1;
			}
        }
    }
*/
	
}
