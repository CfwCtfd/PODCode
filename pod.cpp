#include "pod_def.h"
int main(int argc, char *argv[])
{
//constructor pod("first file name", "metric file", No of modes for calibration)
 POD pod("uv_0001.dat", "area.dat",6);
 //pod.input();
 pod.computePOD();
 return 0;	
}
