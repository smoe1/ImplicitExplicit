#include <fstream>
#include "tensors.h"
using namespace std;

// This is a user-supplied routine that sets the
// auxiliary arrays at all the points "xpts"
//
void AuxFunc(const dTensor2& xpts, 
	     dTensor2& auxvals)
{
    int i;
    int numpts=xpts.getsize();
    double x,gamma,mass_ratio,debye,cs_light,larmor_radius;
    char buffer[256];
    ifstream read_file("param.data", ios::in);

    // Parameters
    read_file >> gamma;
    read_file.getline(buffer,256);
    read_file >> mass_ratio;
    read_file.getline(buffer,256);
    read_file >> debye;
    read_file.getline(buffer,256);
    read_file >> cs_light;
    read_file.getline(buffer,256);
    read_file >> larmor_radius;
    read_file.getline(buffer,256);
    read_file.close();
 
    for (i=1; i<=numpts; i++)
    {
	x = xpts.get(i);
		
	auxvals.set(i,1, gamma );
	auxvals.set(i,2, mass_ratio );
	auxvals.set(i,3, debye );
	auxvals.set(i,4, cs_light );
	auxvals.set(i,5, larmor_radius );
    }

}
