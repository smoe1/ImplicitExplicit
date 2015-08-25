#include "tensors.h"
#include "mesh.h"

// Create right-hand side vector
void Rhs2D_unst(const int bctype,
		const int space_order,
		const int NumRdntBndNodes,
		const iTensor2& rdnt_bnd_node,
		const mesh& Mesh,
		dTensor1& rhs)
{
  // Constants
  const int kmax = ((space_order*(space_order+1))/2);
  const int NumElems = Mesh.get_NumElems();
  const int NumPhysElems = Mesh.get_NumPhysElems();

  // Project rhs function onto DG basis functions
  void RhsFunc(const dTensor2& xpts, 
	       const dTensor2& NOT_USED_1,
	       const dTensor2& NOT_USED_2, 
	       dTensor2& rhs);
  void L2Project_Unst(const int istart, 
		      const int iend, 
		      const int QuadOrder, 
		      const int BasisOrder_qin,
		      const int BasisOrder_auxin,
		      const int BasisOrder_fout,
		      const mesh& Mesh, 
		      const dTensor3* qin, 
		      const dTensor3* auxin, 
		      dTensor3* fout, 
		      void (*Func)(const dTensor2&,const dTensor2&,
				   const dTensor2&,dTensor2&));

  dTensor3 rhs_dg(NumElems,
		  1,
		  kmax);
  L2Project_Unst(1,
		 NumElems,
		 space_order,
		 space_order,
		 space_order,
		 space_order,
		 Mesh,
		 &rhs_dg,
		 &rhs_dg,
		 &rhs_dg,
		 &RhsFunc);

  // Project DG solution onto CG basis functions
  void L2Project_DG2CG_Unst(const int istart, 
			    const int iend, 
			    const int QuadOrder,		    
			    const int cg_order,
			    const int dg_comp,
			    const mesh& Mesh,
			    const dTensor3& fdg,
			    dTensor1& fcg);
  L2Project_DG2CG_Unst(1,
		       NumPhysElems,
		       space_order+1,
		       space_order,
		       1,
		       Mesh,
		       rhs_dg,
		       rhs);

  // If bctype==1: Use right-hand side to also enforce Dirichlet boundary conditions
  switch (bctype)
    {
    case 1:
      // If bctype==1: Use right-hand side to also enforce Dirichlet boundary conditions
      if(space_order==1)
	{
	  for (int i=1; i<=Mesh.get_NumBndNodes(); i++)
	    {
	      int j = Mesh.get_bnd_node(i);
	      // ------------------------------
	      // -- TO DO --
	      // double x = Mesh.get_node(j,1);
	      // double y = Mesh.get_node(j,2);
	      // double val = BCFunc(x,y);
	      // rhs.set(j, val );
	      // ------------------------------
	      rhs.set(j, 0.0 );
	    }
	}
      else
	{
	  for (int i=1; i<=Mesh.get_SubNumBndNodes(); i++)
	    {
	      int j = Mesh.get_sub_bnd_node(i);
	      // ------------------------------
	      // -- TO DO --
	      // double x = Mesh.get_sub_node(j,1);
	      // double y = Mesh.get_sub_node(j,2);
	      // double val = BCFunc(x,y);
	      // rhs.set(j, val );
	      // ------------------------------
	      rhs.set(j, 0.0 );
	    }
	}
      break;

    case 2:
      // If bctype==2: Add redundant node rhs values to primary node rhs values
      for (int i=1; i<=NumRdntBndNodes; i++)
	{        
	  int j = rdnt_bnd_node.get(i,1);
	  int k = rdnt_bnd_node.get(i,2);
	  
	  rhs.set(k, rhs.get(k) + rhs.get(j) );
	}
      break;

    default:
      printf("\n");
      printf(" ERROR in Rhs2D_unst.cpp: bctype = %i is currently unsupported\n",bctype);
      printf("\n");
      exit(1);
      break;
    }
}
