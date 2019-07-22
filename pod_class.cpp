#include <iostream>
#include <stdlib.h>
#include<stdio.h>
#include<vector>
using namespace std;
class POD{
private:
int length;
int *data;
std::string filename; 
std::vector<double> data1;
public:
	POD( int len );             // simple constructor
	POD(const POD& p1);	    // copy constructor

       POD( std::string first_file, int max_files=1000, int no_modes_for_calib=6 );

       // POD( std::string first_file);
        POD& operator=(const POD& p1);
	//output(int *data_main){};
      	~POD();                     // destructor

};
POD::POD(std::string first_file, int y, int r)
{
	filename=first_file;
	cout<<"First File is"<<"\t"<<filename<<endl;
	//std::cout << "hi\n" ;
}

 POD::POD(int len)
{
	cout<<"normal constructor called\n";
	length=len;
	//data=new int [len];
	data1.resize(length);
}	
 POD::POD(const POD& p1)   //copy constructor
 {
  //      To copy data when, for example, used as an argument of another function
          /*n=v1.n;
          v=new T [n];
          for(int i=0;i<n;i++)
          v[i]=v1.v[i];*/

	cout<<"copy constructor called\n";

	length=p1.length;
	data=new int [length];

	for(int i=0;i<length;i++)
	data[i]=p1.data[i];
  
  }

POD& POD::operator=(const POD& p1)	//  p=p1 
{


	if(&p1==this)
	return *this;

	if(data)
	delete [] data;

	length=p1.length;
	
	data=new int [length];

	for(int i=0;i<length;i++)
	data[i]=p1.data[i];

	//cout<<"n="<<n<<endl;

	/*if(v)		// If v exists
	delete [] v;

	n=p1.n;
	v=new T [p1.n];
	
	
	for(int i=0;i<p1.n;i++)
	{
		v[i]=p1.v[i];
	}*/

	return *this;
}
POD::~POD(void)
{
//	if(data)
//	delete [] data;
	data1.clear();
}
void function(POD pod)
{
	cout<<"In function\n";
}
int main(int argc, char *argv[])
{
	POD pod(10);
	//cout<<pod.length<<endl;
	std::string test;
	test="uv_0001.dat";
	POD pod1(test);
	cout<<"sucessfully called the constructur"<<endl;
	//int *data_main;
	//pod.output(data_main);
	//function(pod);
	return 0;
	
}

