#include "dogdefs.h"
#include "IniDocument.h"
#include "DogParams.h"
#include "DogParamsUnst2.h"
#include "edge_data_Unst.h"
#include "mesh.h"
#include "DogStateUnst2.h"
#include "SparseCholesky.h"
#include "ext_time.h" /* for get_utime and timeval_diff */

int RunDogpack_Unst(const char* outputdir)
{

  // ------------------------------------------------------------
  // Function definitions
  void RunMeshCopyScript(const char* outputdir);
  void RunMatlabCopyScript(const char* outputdir);
  void ParseArguments(int argc,char**argv,string& outputdir);
  void GetParams_Unst(string inputdir,string& time_stepping_method);
  void L2Project_Unst(const int istart, 
		      const int iend, 
		      const int QuadOrder,		    
		      const int BasisOrder_fout,
		      const mesh& Mesh, 
		      dTensor3* fout, 
		      void (*Func)(const dTensor2&,dTensor2&));
  void Output_Unst(const mesh& Mesh, 
		   const dTensor3& aux, 
		   const dTensor3& q, 
		   double t, 
		   int nframe, 
		   const char* outputdir);
  void InitApp_Unst(IniDocument& ini_doc,const mesh& Mesh);
  // ------------------------------------------------------------
  
  // Get current time
  timeval start_time = get_utime();

  // Output title information
  printf("\n"
	 "   ------------------------------------------------   \n"
	 "   | DoGPack: The Discontinuous Galerkin Package  |   \n"
	 "   | Developed by the research group of           |   \n"
	 "   |            James A. Rossmanith               |   \n"
	 "   |            Department of Mathematics         |   \n"
	 "   |            Iowa State University             |   \n"
	 "   ------------------------------------------------   \n\n");

  // Get parameters from parameters.ini
  dogParams.init();  

  // allocate state variables
  dogStateUnst2.init();
  dogStateUnst2.set_initial_dt(dogParams.get_initial_dt());
  
  // Check to see if a mesh has been generated,
  //   if YES, then read in basic mesh parameters
  FILE* mesh_file = fopen("Unstructured_Mesh/mesh_output/mesh_params.dat","r");
  if(mesh_file==NULL)
    {
      printf(" ERROR: file not found:"
	     " 'Unstructured_Mesh/mesh_output/mesh_params.dat' \n"
	     "   In order to run DoGPack in unstructured grid mode \n"
	     "   you must first generate a mesh. You can do this by \n"
	     "   following these steps: \n \n"
	     "      (1) Type: $DOGPACK/scripts/create_unst2_dir \n"
	     "      (2) Type: cd Unstructured_Mesh \n"
	     "      (3) Modify the following files as desired:  \n"
	     "            - input2D.data  \n"
	     "            - SignedDistance.cpp  \n"
	     "            - GridSpacing.cpp  \n"
	     "            - MeshPreProcess.cpp  \n"
	     "            - MeshPostProcess1.cpp  \n"
	     "            - MeshPostProcess2.cpp  \n"
	     "      (4) Run mesh generator by typing: mesh2d.exe  \n"
	     "      (5) Visualize mesh using the 'plotmesh2.m' MATLAB script  \n \n");
      exit(1);
    }

  // Read-in in basic mesh parameters
  int NumElems,NumPhysElems,NumGhostElems,NumNodes;
  int NumPhysNodes,NumBndNodes,NumEdges,NumBndEdges;
  char buffer[1024];

  fscanf(mesh_file,"%i",&NumElems);  
  fgets(buffer, sizeof buffer, mesh_file);
  fscanf(mesh_file,"%i",&NumPhysElems);  
  fgets(buffer, sizeof buffer, mesh_file);
  fscanf(mesh_file,"%i",&NumGhostElems);  
  fgets(buffer, sizeof buffer, mesh_file);
  fscanf(mesh_file,"%i",&NumNodes);  
  fgets(buffer, sizeof buffer, mesh_file);
  fscanf(mesh_file,"%i",&NumPhysNodes);  
  fgets(buffer, sizeof buffer, mesh_file);
  fscanf(mesh_file,"%i",&NumBndNodes);  
  fgets(buffer, sizeof buffer, mesh_file);
  fscanf(mesh_file,"%i",&NumEdges);  
  fgets(buffer, sizeof buffer, mesh_file);
  fscanf(mesh_file,"%i",&NumBndEdges);  
  fgets(buffer, sizeof buffer, mesh_file);
  fclose(mesh_file);
    
  // Initialize 2d unstructured parameters
  dogParamsUnst2.init(NumElems,
		      NumPhysElems,
		      NumGhostElems,
		      NumNodes,
		      NumPhysNodes,
		      NumBndNodes,
		      NumEdges,
		      NumBndEdges);  

  // Create and read-in entire unstructured mesh
  mesh Mesh(NumElems,NumPhysElems,NumNodes,NumPhysNodes,
            NumBndNodes,NumEdges,NumBndEdges);
  char mesh_dir[] = "Unstructured_Mesh/mesh_output";
  Mesh.InputMesh(mesh_dir);

  // create help file for plotting purposes
  char qhelp[1024];
  snprintf(qhelp,1024,"%s/qhelp.dat",outputdir);
  dogParams.write_qhelp(qhelp);  
  dogParamsUnst2.write_qhelp(qhelp);

  // Copy mesh into output directory
  RunMeshCopyScript(outputdir);

  // Copy matlab files for Cholesky factorization into output directory
  RunMatlabCopyScript(outputdir);

  // Get application parameters 
  InitApp_Unst(ini_doc,Mesh);

  // Get sub-factor value
  char subfname[1024];
  snprintf(subfname,1024,"%s/matlab_output/subfactor.dat",outputdir);
  FILE* subread = fopen(subfname,"r");
  void filenotfound(char*);
  if (subread==NULL)
    { filenotfound(subfname); }
  int subtype;
  fscanf(subread,"%i",&subtype);
  fclose(subread);
  assert_eq(subtype,dogParams.get_space_order());

  // Dimension arrays
  const int kmax = dogParams.get_kmax();
  const int space_order = dogParams.get_space_order();  
  const int meqn = dogParams.get_meqn();
  const int maux = dogParams.get_maux();
  dTensor3   q(NumElems,meqn,kmax);
  dTensor3 aux(NumElems,maux,kmax);
  int SubNumPhysNodes;
  int SubNumBndNodes;
  int SubFactor;
  
  if (space_order==1)
    { SubFactor = 1; }
  else
    { SubFactor = Mesh.get_SubFactor(); }

  assert_eq(SubFactor,space_order);
  
  if (space_order==1)
    {
      SubNumPhysNodes = NumPhysNodes;
      SubNumBndNodes  = NumBndNodes;
    }
  else
    {
      SubNumPhysNodes = Mesh.get_SubNumPhysNodes();
      SubNumBndNodes  = Mesh.get_SubNumBndNodes();
    }

  // Find out boundary condition type
  void filenotfound(char*);
  char bcfname[1024];
  snprintf(bcfname,1024,"%s/matlab_output/boundary_condition.dat",outputdir);
  FILE* bcread = fopen(bcfname,"r");
  if (bcread==NULL)
    { filenotfound(bcfname); }
  int bctype;
  fscanf(bcread,"%i",&bctype);
  fclose(bcread);

  // Read in redundant node information
  char fname[1024];
  if (SubFactor==1)
    {
      snprintf(fname,1024,"%s/mesh_output/mesh_rdnt_node.dat",outputdir);
    }
  else
    {
      snprintf(fname,1024,"%s/mesh_output/submesh_rdnt_node.dat",outputdir);
    }
  FILE* read1 = fopen(fname,"r");
  if (read1==NULL)
    { filenotfound(fname); }
  int NumRdntBndNodes;
  fscanf(read1,"%i",&NumRdntBndNodes);
  iTensor2 rdnt_bnd_node(NumRdntBndNodes,2);
  for (int i=1; i<=NumRdntBndNodes; i++)
    {
      int tmp1,tmp2;
      fscanf(read1,"%i %i",&tmp1,&tmp2);
      rdnt_bnd_node.set(i,1, tmp1 );
      rdnt_bnd_node.set(i,2, tmp2 );
    }
  fclose(read1);
  
  // -------------------------- 
  // Start new computation
  // -------------------------- 
  
  // Variables
  dTensor1 rhs(SubNumPhysNodes);
  dTensor1 phi(SubNumPhysNodes);
  dTensor2  E1(NumElems,kmax);
  dTensor2  E2(NumElems,kmax);

  // Get Cholesky factorization matrix R
  SparseCholesky R(SubNumPhysNodes);
  R.read(outputdir);

  // Create right-hand side vector
  void Rhs2D_unst(const int bctype,
		  const int space_order,
		  const int NumRdntBndNodes,
		  const iTensor2& rdnt_bnd_node,
		  const mesh& Mesh,
		  dTensor1& rhs);
  Rhs2D_unst(bctype,
	     space_order,
	     NumRdntBndNodes,
	     rdnt_bnd_node,
	     Mesh,
	     rhs);
  
  // Call Poisson solver  
  void PoissonSolver2D_unst(const int bctype,
			    const int space_order,
			    const int NumRdntBndNodes,
			    const iTensor2& rdnt_bnd_node,
			    const mesh& Mesh,			    
			    const SparseCholesky& R,
			    const dTensor1& rhs,
			    dTensor1& phi,
			    dTensor2& E1,
			    dTensor2& E2);
  PoissonSolver2D_unst(bctype,
		       space_order,
		       NumRdntBndNodes,
		       rdnt_bnd_node,
		       Mesh,
		       R,
		       rhs,
		       phi,
		       E1,
		       E2);
  
  // Output data to file
  void DumpSolution(const int space_order,
		    const mesh& Mesh,
		    const dTensor1& phi,
		    const dTensor2& E1,
		    const dTensor2& E2,
		    dTensor3& q);
  DumpSolution(space_order,Mesh,phi,E1,E2,q);
  Output_Unst(Mesh,aux,q,0.0,0,outputdir); 
  
  // Get current time
  timeval end_time = get_utime();

  // Output elapsed time
  double diff_utime = timeval_diff(end_time, start_time);
  printf(" Total elapsed time in seconds = %11.5f\n\n", diff_utime);
  
  return 0;
}

// Dump electric potential and fields into variable "q"
void DumpSolution(const int space_order,
		  const mesh& Mesh,
		  const dTensor1& phi,
		  const dTensor2& E1,
		  const dTensor2& E2,
		  dTensor3& q)
{
  const int NumPhysElems = Mesh.get_NumPhysElems();
  const int meqn = q.getsize(2);
  const int kmax = q.getsize(3);
  
  switch(space_order)
    {
    case 1:
      assert(kmax==1);
      for (int i=1; i<=NumPhysElems; i++)
	{
	  int n1 = Mesh.get_tnode(i,1);
	  int n2 = Mesh.get_tnode(i,2);
	  int n3 = Mesh.get_tnode(i,3);

	  double phi1 = phi.get(n1);
	  double phi2 = phi.get(n2);
	  double phi3 = phi.get(n3);

	  q.set(i,1,1, onethird*(phi1+phi2+phi3)              );
	}
      break;

    case 2:
      assert(kmax==3);
      for (int i=1; i<=NumPhysElems; i++)
	{
	  int n1 = Mesh.get_node_subs(i,1);
	  int n2 = Mesh.get_node_subs(i,2);
	  int n3 = Mesh.get_node_subs(i,3);
	  int n4 = Mesh.get_node_subs(i,4);
	  int n5 = Mesh.get_node_subs(i,5);
	  int n6 = Mesh.get_node_subs(i,6);

	  double phi1 = phi.get(n1);
	  double phi2 = phi.get(n2);
	  double phi3 = phi.get(n3);
	  double phi4 = phi.get(n4);
	  double phi5 = phi.get(n5);
	  double phi6 = phi.get(n6);

	  q.set(i,1,1, onethird*(phi2+phi4+phi5) );
	  q.set(i,1,2, (sq2/60.0)*(-3.0*phi1+4.0*phi2+6.0*phi3-8.0*phi4+4.0*phi5-3.0*phi6) );
	  q.set(i,1,3, (sq2*sq3/60.0)*(-3.0*phi1-4.0*phi2+4.0*phi5+3.0*phi6) );
	}
      break;

    case 3:
      assert(kmax==6);
      for (int i=1; i<=NumPhysElems; i++)
	{
	  int n1  = Mesh.get_node_subs(i,1);
	  int n2  = Mesh.get_node_subs(i,2);
	  int n3  = Mesh.get_node_subs(i,3);
	  int n4  = Mesh.get_node_subs(i,4);
	  int n5  = Mesh.get_node_subs(i,5);
	  int n6  = Mesh.get_node_subs(i,6);
	  int n7  = Mesh.get_node_subs(i,7);
	  int n8  = Mesh.get_node_subs(i,8);
	  int n9  = Mesh.get_node_subs(i,9);
	  int n10 = Mesh.get_node_subs(i,10);

	  double phi1  = phi.get(n1);
	  double phi2  = phi.get(n2);
	  double phi3  = phi.get(n3);
	  double phi4  = phi.get(n4);
	  double phi5  = phi.get(n5);
	  double phi6  = phi.get(n6);
	  double phi7  = phi.get(n7);
	  double phi8  = phi.get(n8);
	  double phi9  = phi.get(n9);
	  double phi10 = phi.get(n10);	  

	  double A1  =  onethird*0.0025*( 40.0*(phi1+phi4+phi10) 
					  + 90.0*(phi2+phi3+phi5+phi7+phi8+phi9) + 540.0*phi6 );
	  double A2  = -onethird*0.05*osq2*(phi1+phi10-2.0*phi4+9.0*(phi2+phi5+phi8+phi9)-18.0*(phi3+phi7) );
	  double A3  = -onethird*0.05*osq2*sq3*(phi1-phi10 + 9.0*(phi2+phi5-phi8-phi9));
	  double A4  =  oneseventh*0.075*osq7*(phi2+phi5+4.0*phi1+17.0*(phi7+phi9)
					       -12.0*(phi4+phi10)+30.0*phi6-23.0*(phi3+phi8));
	  double A5  =  oneseventh*0.05*sq3*osq7*(2.0*(phi8-phi10)+10.0*phi1-15.0*phi2+12.0*(phi4-phi3)
						  +18.0*phi7-3.0*phi9+20.0*phi5-30.0*phi6);
	  double A6  =  oneseventh*0.05*sq3*sq5*(2.0*(phi1-phi5-phi8+phi10)+3.0*(phi2+phi9)-6.0*phi6);
	  
	  q.set(i,1,1,  A1  );
	  q.set(i,1,2,  A2  );
	  q.set(i,1,3,  A3  );
	  q.set(i,1,4,  A4  );
	  q.set(i,1,5,  A5  );
	  q.set(i,1,6,  A6  );
	}
      break;

    default:
      printf("\n");
      printf(" ERROR in DumpSolution.cpp: space_order value not supported.\n");
      printf("       space_order = %i\n",space_order);
      printf("\n");
      exit(1);
      break;	  
    }

  for (int i=1; i<=NumPhysElems; i++)
    for (int k=1; k<=kmax; k++)
      {
	q.set(i,2,k, E1.get(i,k)  );
	q.set(i,3,k, E2.get(i,k)  );
      }

}

void filenotfound(char* fname)
{
  printf("\n");
  printf("  Error in RunDogpack_Unst.cpp.\n");
  printf("  File not found: %s.\n",fname);
  printf("\n");
  exit(1);
}
