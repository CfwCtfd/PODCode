#include "matrix_def.h"
void matmul(double** A, double* B, int dim1A,int dim2A,int dim2B, double** AxB)
{
    double sum;
    for(int i=0;i<dim1A;i++){
        for( int k=0; k<dim2B; k++){
           sum = 0;
           for(int j=0; j<dim2A; j++){
                 sum += A[i][j]*B[j*dim2B+k];		   
             }
              AxB[i][k] = sum;
	     }
         }	
}

void transpose(double* A, double* At, int rowA, int colA)
{
	double temp;
	for(int i=0;i<rowA;++i){
		for(int j=0;j<colA;++j){
			At[j*colA+i]=A[i*colA+j];
						
		}
	}
}

void printmat(double** A,int row,int col)
{
	
  for (int  i=0; i<row; ++i)
 {  
    for (int j=0; j<col; ++j)
   {
    cout<<A[i][j]<<"\t";
   }
  cout<<"\n";
 }
 }
 
 void printvecInMat(double *A,int row,int col,int nfiles)
 {
  for(int i=0;i<row; ++i){
	 for(int j=0; j<col;++j){
		 cout<<A[i*nfiles+j]<<"\t";
	 }
	 cout<<"\n";
  }
 }