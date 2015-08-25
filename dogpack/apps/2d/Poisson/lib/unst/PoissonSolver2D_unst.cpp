#include "dogdefs.h"
#include "edge_data_Unst.h"
#include "mesh.h"
#include "SparseCholesky.h"

// Poisson's equation:
//
//     -\Delta phi = rhs.
//
// Electric field is defined as (E1, E2) = - \grad phi.
//
void PoissonSolver2D_unst(const int bctype,
			  const int space_order,
			  const int NumRdntBndNodes,
			  const iTensor2& rdnt_bnd_node,
			  const mesh& Mesh,			  
			  const SparseCholesky& R,
			  const dTensor1& rhs,
			  dTensor1& phi,
			  dTensor2& E1,
			  dTensor2& E2)
{
  // Solve linear system using cholesky factorization + forward/backward substitution
  const int N = rhs.getsize();
  dTensor1 rhs_tmp(N);
  R.ForwardSubs(rhs,rhs_tmp);
  R.BackwardSubs(rhs_tmp,phi);

  // Apply periodic boundary conditions if necessary
  if (bctype==2)
    {
      for (int i=1; i<=NumRdntBndNodes; i++)
	{        
	  int j = rdnt_bnd_node.get(i,1);
	  int k = rdnt_bnd_node.get(i,2);
	  phi.set(j, phi.get(k) );
	}
    }

  // Compute electric field as gradient of the electric potential
  void ComputeEfield(const int space_order,
		     const mesh& Mesh,
		     const dTensor1& phi,
		     dTensor2& E1,
		     dTensor2& E2);
  ComputeEfield(space_order,Mesh,phi,E1,E2);

}
