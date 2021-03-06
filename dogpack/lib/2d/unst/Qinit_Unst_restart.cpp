#include "dogdefs.h"
#include "DogParams.h"

// Initial data from the file qXXXX_restart.dat (where XXXX=nstart)
double QinitRestartASCII_Unst(char* fname, dTensor3& q)
{

    const int NumElems = q.getsize(1);
    const int     meqn = q.getsize(2);
    const int     kmax = q.getsize(3);

    char fname1[1024];
    snprintf(fname1,1024,"%s.dat",fname);
    FILE* restart_file = fopen(fname1,"r");

    if(restart_file==NULL)
    {
        printf("Qinit_Unst_restart: failed to open file %s\n",fname1);
        exit(1);
    }

    // read in old data (written in Output.cpp)
    double t;
    int garbage = fscanf(restart_file,"%lf",&t);
    for (int i=1; i<=NumElems; i++)
    for (int m=1; m<=meqn; m++)
    for (int k=1; k<=kmax; k++)
    {
        double tmp;
        garbage = fscanf(restart_file,"%lf",&tmp);
        q.set(i,m,k, tmp );
    }
    fclose(restart_file);

    return t;

}

double QinitRestart_Unst(int nstart, char* varname, 
             dTensor3& q, const char* outputdir)
{  

    char tmp[5];
    tmp[0]='0';
    tmp[1]='0';
    tmp[2]='0';
    tmp[3]='0';
    tmp[4]='\0';
    if( nstart < 10 )
    {

        char* chartmp = new char[1];
        sprintf(chartmp, "%d", nstart);

        tmp[3] = chartmp[0];      
        delete chartmp;

    }
    else if (nstart<100)
    {

        char* chartmp = new char[2];
        sprintf(chartmp, "%d", nstart);

        tmp[2]=chartmp[0]; 
        tmp[3]=chartmp[1];
        delete chartmp;

    }
    else if (nstart<1000)
    {
        char* chartmp = new char[3];
        sprintf(chartmp, "%d", nstart);

        tmp[1]=chartmp[0]; 
        tmp[2]=chartmp[1];
        tmp[3]=chartmp[2];
        delete chartmp;

    }
    else
    {
        char* chartmp = new char[4];
        sprintf(chartmp, "%d", nstart);

        tmp[0]=chartmp[0]; 
        tmp[1]=chartmp[1];
        tmp[2]=chartmp[2];
        tmp[3]=chartmp[3]; 
        delete chartmp;

    }
    char fname[1024];
    snprintf(fname,1024,"%s/%s%s_restart",outputdir,varname,tmp);

    double retval;
    switch(dogParams.get_datafmt())
    {

        case ASCII:
            double QinitRestartASCII_Unst(char* fname, dTensor3& q);
            retval = QinitRestartASCII_Unst(fname, q);
            break;
        default:
            printf("Qinit_Unst_restart.cpp: your chosen datafmt is not supported\n");
            exit(1);
            break;

    }
    return retval;
}

double Qinit_Unst_restart(int nstart, 
              dTensor3& q, 
              const char* outputdir)
{
    char* qname = new char[1];
    qname[0]='q';
    return QinitRestart_Unst(nstart, qname, q, outputdir);
    delete qname;
}
