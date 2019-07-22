# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <cmath>
# include <string>
# include <sstream>
# include<vector>
# include <complex>
using namespace std;
# include "linpack_d.hpp"
# include "blas1_d.hpp"
# include "blas0.hpp"
# include "matrix_imp.h"
class POD{

fstream readdata,griddata;
ofstream eigval,tempmodes,logeig;
std::string line,filename,first_file,metric_file; 
std::stringstream metrics;
int nfiles,filepos;
int col,row,ncomp,Ntot,nmodes,nx;
int max_files, no_modes_for_calib;
double temp;
double **snap, **phivectp;
double *correlt,*sval,*u,*v,*lambda,*area,*uTranspose,*average;
std::vector<double> data;
std::vector<double> metric;

public:
POD( std::string firstfile, std::string metricfile,  int no_modes_for_calib=6, int max_files=1000 );
~POD();                     // destructor
   void ReadSnapData( );
   void AverageOutSnap();
   void readmetric();
   double correlation(int , int );
   void r8mat_svd_linpack ( );
   void write_eigen_ric();
   void write_temp_modes();
   void find_spatial_modes();
   void write_spatial_modes();
   void write_modes_calib();
   void file_name_inc_nowrap ( string *filename );
   void file_info( );
   void computePOD();
};
POD::POD(std::string firstfile, std::string metricfile,  int r, int y)
{
	first_file=firstfile;
	metrics.str("");
        metrics<<metricfile;
	cout<<"First File is"<<"\t"<<first_file<<endl;
	cout<<"Metric File name is"<<"\t"<<metrics.str()<<endl;
	no_modes_for_calib =r;
	cout<<"no of modes for calib=\t"<<r<<endl;
        max_files=y;	
}
POD::~POD(void)
{
	data.clear();
}

void POD::computePOD()
{
    file_info();	
    cout<<"The total number of files="<<nfiles<<endl;
    snap=new double*[Ntot];
  for(int i = 0; i < Ntot; ++i)  snap[i] = new double[nfiles];
    ReadSnapData( ); 
    cout<<"Data has been read to snapshot matrix"<<endl;
    average=new double[Ntot];
    cout<<"Data is being averaged"<<endl;
    AverageOutSnap();
    nx = Ntot/ncomp;
    cout<<"The number of grid points="<<nx<<"\n";
    area=new double[nx];
    readmetric();
    correlt=new double[nfiles*nfiles];
    sval=new double[nfiles*nfiles];
    u=new double[nfiles*nfiles];
    v=new double[nfiles*nfiles];
    uTranspose=new double[nfiles*nfiles];
    cout<<"Determining the correlation Matrix...."<<endl;
 //to find the correlation matrix
  for (int i=0;i<nfiles;i++){
    for(int j=0;j<nfiles;j++ ){
     correlt[i*nfiles+j] = correlation( i, j );
     }
  }
    cout<<"Finished determining the correlation matrix"<<endl;
    cout<<"The first 5x5 components of the correlation Matrix is:"<<endl;
    printvecInMat(correlt,5,5,nfiles);	
  //To find the Singular Value Decomposition
    r8mat_svd_linpack (  );
    lambda = new double [nfiles];
 for (int i = 0;i<nfiles;++i){
	 for(int j=0;j<nfiles;++j){
		 if(i==j){
        	 lambda[i] = sval[i*nfiles+j];
		}
	 }
	// cout<<"lambda["<<i<<"]=\t"<<lambda[i]<<"\n";
}
    cout<<"Writing the Eigen information"<<endl;
    write_eigen_ric();
 // to Normalise the temporal modes
    transpose(u,uTranspose,nfiles,nfiles);
    cout<<"Writing the temporal modes"<<endl;
    write_temp_modes();
    phivectp = new double*[data.size()];
 for(int i = 0; i < data.size(); ++i) 
  phivectp[i] = new double[nfiles];
  find_spatial_modes();
  cout<<"Writing the spatial modes"<<endl;
  write_spatial_modes();
  write_modes_calib();	
 //delete the arrays
    for(int i = 0; i < data.size(); ++i){ 
     delete [] snap[i];
     delete [] phivectp[i];
  }
     delete [] snap;
     delete [] phivectp;
     delete sval;
     delete correlt;
     delete u;
     delete v;
     delete uTranspose;
     delete lambda;
     delete area;
    delete average;
  cout << "POD exectuted on " << __DATE__;
  cout << " at " << __TIME__ << ".\n";
 }

 void POD::file_info()
{
std::string filename;
filename = first_file;
nfiles=0;
readdata.open(first_file.c_str());
if(!readdata.is_open()){
  cout << first_file << "not found \n";
  exit (0); 
}
 data.clear();
 getline ( readdata, line );
 getline ( readdata, line );
 filepos=readdata.tellg();
 getline ( readdata, line );
 std::stringstream ss(line);
ncomp=0;
while(ss>>temp){
ncomp++;
}
readdata.seekg (filepos, readdata.beg);
while(true)
{
readdata>>temp;
if(readdata.eof()){break;}
data.push_back(temp);
}
   readdata.close();
   Ntot=data.size();
while(nfiles<max_files){
readdata.open(filename.c_str()); 
 if(!readdata.is_open())
 break;
 readdata.close();
 ++nfiles;
   file_name_inc_nowrap ( &filename );
  }
}
 
 
void POD::ReadSnapData( )
{
std::string filename;
filename = first_file;
for(int fiter=1;fiter<=nfiles;fiter++){
filepos=0;
temp=0;	
readdata.open(filename.c_str());
if(!readdata.is_open()){
  cout << filename << "not found \n";
  exit (0); 
}
 data.clear();
 getline ( readdata, line );
 getline ( readdata, line );
 filepos=readdata.tellg();
 getline ( readdata, line );
 std::stringstream ss(line);
ncomp=0;
while(ss>>temp){
ncomp++;
}
readdata.seekg (filepos, readdata.beg);
while(true)
{
readdata>>temp;
if(readdata.eof()){break;}
data.push_back(temp);
}
   readdata.close();
   cout<<filename<<"\n";
   Ntot=data.size();
   cout<<"The total number of components="<<ncomp<<endl;
   cout<<"The total number of solution componenets="<<data.size()<<endl;
file_name_inc_nowrap(&filename);
for (int i=0;i<Ntot;i++)snap[i][fiter-1] = data[i];
}
}

void POD::AverageOutSnap()
{
	for(int i=0;i<Ntot;i++){
		double sum=0;
		for(int j=0;j<nfiles;j++){
		  sum+=snap[i][j];
		}
		sum /= nfiles;
		average[i]=sum;
	}
	/*To subtract the average vector from the snapshots*/
	
	for(int j=0;j<nfiles;j++){
		for(int i=0;i<Ntot;i++){
                snap[i][j] -=average[i];
		}
	}
	
}

 void POD::readmetric()
 {
	ifstream griddata; 
	double temp;
	griddata.open(metrics.str().c_str());
        while(true)
        {
	  griddata>>temp;
	  if(griddata.eof()){break;}
  	  metric.push_back(temp);
        }
          if(metric.size()!=nx){
          cout<<"Some problem with the grid file---->exiting the program"<<"\n";
	   exit(0);
	}
for(int i=0;i<metric.size();++i)
	area[i]=metric[i];
      	
}

double POD::correlation(int i, int j)
{
   int upj,dpj;
   double u1,u2,v1,v2;
   double sum=0;
   for(int k=0;k<nx;++k){
   upj = k*ncomp;
   dpj =1+k*ncomp;
   u1=snap[upj][i];
   u2=snap[upj][j];
   v1=snap[dpj][i];
   v2=snap[dpj][j];
   sum+= (u1*u2+v1*v2)*area[k];
   }	   
   sum /=	nfiles;
    return sum;
}

 void POD::write_eigen_ric()
{
ofstream eigval,logeig;
eigval.open("../results_pod/lambda.dat", ios::trunc );
logeig.open("../results_pod/RIC.dat",ios::trunc);
double sumeig=0;
 for (int i = 0;i<nfiles;++i){
		eigval<<std::setprecision(20)<<i<<"\t"<<log10(lambda[i]/lambda[0])<<"\n";
		sumeig+=lambda[i]; 
	 
}
double ricsum=0;
for (int i = 0;i<nfiles;++i){
	 	  ricsum+=lambda[i];
		  logeig<<std::setprecision(20)<<ricsum/sumeig*100<<"\n"; 
}
logeig.close();
eigval.close(); 	
return;	
}

void POD::write_temp_modes()
 {
    ofstream tempmodes,templstm;
    std::stringstream ss1,ss2;
    
    

    for(int fiter=1;fiter<=nfiles;fiter++){
    ss1.str("");
   if(fiter<=9)
    ss1<<"../results_pod/avectp_000"<<fiter<<".dat";
   if(fiter>=10  && fiter<=99)
   ss1<<"../results_pod/avectp_00"<<fiter<<".dat";
   if(fiter>=100 && fiter<=nfiles)
   ss1<<"../results_pod/avectp_0"<<fiter<<".dat";
 

   tempmodes.open(ss1.str().c_str());
    if(!tempmodes.is_open()){
       cout << ss1.str() << "not found \n";
       exit (0); 
    } 

     for(int i=0;i<nfiles;++i){
	   uTranspose[i*nfiles+(fiter-1)] *=sqrt(lambda[fiter-1]*nfiles);
	   tempmodes<<std::setprecision(20)<<uTranspose[i*nfiles+(fiter-1)] <<"\n";
         }
	tempmodes.close();
        	   
 }
      ss2<<"../results_pod/temp_lstm"<<".dat";
      templstm.open(ss2.str().c_str());
      if(!templstm.is_open()){
       cout << ss2.str() << "not found \n";
       exit (0); 
      } 

    for(int i=0;i<nfiles;++i){
          for(int j=0;j<no_modes_for_calib;++j){
            
           templstm<<std::setprecision(15)<<uTranspose[i*nfiles+j] <<"\t";

        }
       templstm<<"\n";
    }




 templstm.close();

 return;
 }
 
 void POD::find_spatial_modes()
{
matmul(snap,uTranspose,Ntot,nfiles,nfiles, phivectp); 
	for(int j=0;j<nfiles;++j){
	  for(int i=0;i<Ntot;++i){
	            phivectp[i][j] = phivectp[i][j]/(nfiles*lambda[j]);
        	}
	}
      return;
}

 void POD::write_spatial_modes()
 {
   std::stringstream ss1;
   ofstream tempmodes;
   for(int fiter=1;fiter<=nfiles;fiter++){
    ss1.str("");
   if(fiter<=9)
    ss1<<"../results_pod/phivectp_000"<<fiter<<".dat";
   if(fiter>=10  && fiter<=99)
   ss1<<"../results_pod/phivectp_00"<<fiter<<".dat";
   if(fiter>=100 && fiter<=nfiles)
   ss1<<"../results_pod/phivectp_0"<<fiter<<".dat";
   tempmodes.open(ss1.str().c_str());
   if(!tempmodes.is_open()){
    cout << ss1.str() << "not found \n";
    exit (0); 
   } 
   for(int i=0;i<data.size();++i)
   {tempmodes<<std::setprecision(20)<<"\t"<<phivectp[i][fiter-1] ;
	   if((i+1)%2==0) tempmodes<<"\n";}
   tempmodes.close();	   
    }
    cout<<"writing the average vector"<<endl;
    ss1.str("");
    ss1<<"../results_pod/phivectp_000.dat";
    tempmodes.open(ss1.str().c_str());
    for(int i=0;i<data.size();++i){
	tempmodes<<std::setprecision(20)<<"\t"<<average[i] ;
	   if((i+1)%2==0) tempmodes<<"\n";    
    }
       tempmodes.close();
}

void POD::write_modes_calib(){
   std::stringstream ss1,ss2;
   ofstream eigen,at;
   ss1.str("");
   ss2.str("");
   ss1<<"../results_pod/lambda_saved.txt";
   ss2<<"../results_pod/avectp_saved.txt";
   eigen.open(ss1.str().c_str());
   at.open(ss2.str().c_str());
   for(int j=0;j<no_modes_for_calib;j++){
   eigen<<lambda[j]<<"\n";
     for(int i=0;i<nfiles;i++)
      {
      at<<uTranspose[i*nfiles+j]<<std::setprecision(20)<<"\n";
      }
     }
       eigen.close();
      at.close();
}
  
void POD::file_name_inc_nowrap ( string *filename )
//****************************************************************************80
//
//  Purpose:
//
//    FILE_NAME_INC_NOWRAP increments a partially numeric file name.
//
//  Discussion:
//
//    It is assumed that the digits in the name, whether scattered or
//    connected, represent a number that is to be increased by 1 on
//    each call.  If this number is all 9's on input, the output number
//    is all 0's.  Non-numeric letters of the name are unaffected.
//
//    If the (nonempty) name contains no digits, or all the digits are
//    9, then the empty string is returned.
//
//    If the empty string is input, the routine stops.
//
//  Example:
//
//      Input            Output
//      -----            ------
//      "a7to11.txt"     "a7to12.txt"  (typical case.  Last digit incremented)
//      "a7to99.txt"     "a8to00.txt"  (last digit incremented, with carry.)
//      "a8to99.txt"     "a9to00.txt"
//      "a9to99.txt"     " "
//      "cat.txt"        " "
//      " "              STOP!
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 November 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, string *FILENAME, the filename to be incremented.
//
{
  char c;
  int carry;
  int change;
  int i;
  int lens;

  lens = (*filename).length ( );

  if ( lens <= 0 )
  {
    cerr << "\n";
    cerr << "FILE_NAME_INC_NOWRAP - Fatal error!\n";
    cerr << "  The input string is empty.\n";
    exit ( 1 );
  }

  change = 0;
  carry = 0;

  for ( i = lens - 1; 0 <= i; i-- )
  {
    c = (*filename)[i];

    if ( '0' <= c && c <= '9' )
    {
      change = change + 1;
      carry = 0;

      if ( c == '9' )
      {
        carry = 1;
        c = '0';
        (*filename)[i] = c;
      }
      else
      {
        c = c + 1;
        (*filename)[i] = c;
        return;
      }
    }
  }
//
//  Unsatisfied carry.  The input digits were all 9.  Return blank.
//
  if ( carry == 1 )
  {
    for ( i = lens - 1; 0 <= i; i-- )
    {
      (*filename)[i] = ' ';
    }
  }

//
//  No digits were found.  Return blank.
//
  if ( change == 0 )
  {
    for ( i = lens - 1; 0 <= i; i-- )
    {
      (*filename)[i] = ' ';
    }
  }

  return;
}

void POD::r8mat_svd_linpack ( )
//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_SVD_LINPACK gets the SVD of a matrix using a call to LINPACK.
//
//  Discussion:
//
//    The singular value decomposition of a real MxN matrix A has the form:
//
//      A = U * S * V'
//
//    where 
//
//      U is MxM orthogonal,
//      S is MxN, and entirely zero except for the diagonal;
//      V is NxN orthogonal.
//
//    Moreover, the nonzero entries of S are positive, and appear
//    in order, from largest magnitude to smallest.
//
//    This routine calls the LINPACK routine DSVDC to compute the
//    factorization.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 September 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns in the matrix A.
//
//    Input, double A[M*N], the matrix whose singular value
//    decomposition we are investigating.
//
//    Output, double U[M*M], S[M*N], V[N*N], the factors
//    that form the singular value decomposition of A.
//
{
  double *a_copy;
  double *e;
  int i;
  int info;
  int j;
  int lda;
  int ldu;
  int ldv;
  int job;
  int lwork;
  double *sdiag;
  double *work;
	
  int m = nfiles;
  int n=nfiles;
	
//
//  The correct size of E and SDIAG is min ( m+1, n).
//
  a_copy = new double[m*n];
  e = new double[m+n];
  sdiag = new double[m+n];
  work = new double[m];
//
//  Compute the eigenvalues and eigenvectors.
//
  job = 11;
  lda = m;
  ldu = m;
  ldv = n;
//
//  The input matrix is destroyed by the routine.  Since we need to keep
//  it around, we only pass a copy to the routine.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    { 
      a_copy[i+j*m] = correlt[i+j*m];
    }
  }
  info = dsvdc ( a_copy, lda, m, n, sdiag, e, u, ldu, v, ldv, work, job );
 
  if ( info != 0 )
  {
    cout << "\n";
    cout << "R8MAT_SVD_LINPACK - Failure!\n";
    cout << "  The SVD could not be calculated.\n";
    cout << "  LINPACK routine DSVDC returned a nonzero\n";
    cout << "  value of the error flag, INFO = " << info << "\n";
    return;
  }
//
//  Make the MxN matrix S from the diagonal values in SDIAG.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      if ( i == j )
      {
        sval[i+j*m] = sdiag[i];
      }
      else
      {
        sval[i+j*m] = 0.0;
      }
    }
  }
//
//  Note that we do NOT need to transpose the V that comes out of LINPACK!
//
  delete [] a_copy;
  delete [] e;
  delete [] sdiag;
  delete [] work;
  return;
}
