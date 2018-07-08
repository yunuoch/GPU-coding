#include "mex.h"

#include <Eigen/Sparse>
#include <vector>

typedef Eigen::SparseMatrix<float> SparseMat;
typedef Eigen::Triplet<float> T;
std::vector<T> coefTri;


void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray*prhs[])

{
//x

int M1 = mxGetM(prhs[1]);
int N1 = mxGetN(prhs[1]);
double* pArr1 = (double*)mxGetPr(prhs[1]);

//Eigen::MatrixXd xx = MatrixXd::zeros(M1,3);  


	Eigen::VectorXf rx(M1);										// right side
	Eigen::VectorXf ry(M1);										
	Eigen::VectorXf rz(M1);	
	Eigen::VectorXf lrx(M1);										// right side
	Eigen::VectorXf lry(M1);										
	Eigen::VectorXf lrz(M1);	

	for (int i = 0; i < M1; i++)
	{
		rx(i)=pArr1[i];
		ry(i)=pArr1[M1+i];
		rz(i)=pArr1[2*M1+i];
	}



	
	
double  *pr;

mwIndex     *ir, *jc;

int      row, col, total=0, number_of_columns;

int      starting_row_index, stopping_row_index,

current_row_index;

/* Get the starting positions of the data in the sparse array. */

pr = mxGetPr(prhs[0]);

ir = mxGetIr(prhs[0]);

jc = mxGetJc(prhs[0]);

/* Display the nonzero elements of the sparse array. */

number_of_columns = mxGetN(prhs[0]);

for (col=0; col<number_of_columns; col++)

{

starting_row_index = jc[col];

stopping_row_index = jc[col+1];

if (starting_row_index == stopping_row_index)

continue;

else
	
{

for (current_row_index = starting_row_index;

current_row_index < stopping_row_index;

current_row_index++)
{
mexPrintf("(%d,%d) = %g\n", ir[current_row_index]+1,col+1, pr[total]);

//coefTri.emplace_back(T(ir[current_row_index]+1,col+1, pr[total]));//
coefTri.emplace_back(T(ir[current_row_index],col, pr[total]));

total++;

}

}

}

SparseMat L(mxGetM(prhs[0]), mxGetN(prhs[0]));
L.setFromTriplets(coefTri.begin(), coefTri.end());

    for(int itr=0; itr<1024; itr++)
	{
	lrx=-L*rx;lry=-L*ry;lrz=-L*rz;
	rx=lrx;ry=lry;rz=lrz;
	}

/*	for (int i = 0; i < M1; i++)
	{
		rx(i)=lrx(i);
		ry(i)=lry(i);
		rz(i)=lrz(i);
	}
	
	rx=lrx;ry=lry;rz=lrz;*/
	plhs[0] = mxCreateDoubleMatrix(M1, N1, mxREAL);    //创建一个M1行，N1列的矩阵
    double* pRe =(double*)mxGetPr(plhs[0]);
	for (mwIndex i=0; i < M1;i++)//                        int is not okay
	{
		pRe[i]=rx(i);
		pRe[M1+i]=ry(i);
		pRe[2*M1+i]=rz(i);
	}//display x
	//pRe[0]=rx(0);pRe[1]=rx(1);pRe[2]=rx(2);pRe[3]=ry(0);pRe[4]=ry(1);pRe[5]=rz(0);
	//pRe[3]=rx(3);pRe[M1]=rx(1);pRe[2*M1]=rx(1);





//Eigen::MatrixXd xI=L*rx;



//mexPrintf("(%d,%d) = %g\n", ir[current_row_index]+1,col+1, pr[total]);

}
