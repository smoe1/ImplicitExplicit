#include "debug.h"
#include <cmath>
#include "DogParams.h"
#include "GEMparams.h"
#include "tensors.h"
#include "Boundaries.h"
#include "assert.h"

// for left/right boundaries
inline int getXreverse(int component)
{
  switch(component)
  {
    case 4:
    case 2:
      return -1;
    case 6:
    case 5:
    case 3:
    case 1:
      return 1;
    default:
      unsupported_value_error(component);
  }
}
// for top/bottom boundaries
inline int getYreverse(int component)
{
  switch(component)
  {
    case 4:
    case 3:
      return -1;
    case 6:
    case 5:
    case 2:
    case 1:
      return 1;
    default:
      unsupported_value_error(component);
  }
}

int get_idx_of_first_zero(const int* idx_arr)
{
  int i=0;
  while(idx_arr[i]) i++;
  return i;
}

// In the following methods meq_ind is an optional parameter
// which supplies a subset of the equation indices to which to
// apply the desired boundary condition.

// equivalent to saying that all symmetries are even
// effectively enforces normal partial derivatives to be zero
//
void SetLeftMirrorBCs(dTensorBC4& T, const int* meq_ind)
{
    int meqi_hgh;
    if(meq_ind) meqi_hgh = get_idx_of_first_zero(meq_ind);
    else meqi_hgh = T.getsize(3);
    //
    int mbc = T.getmbc();
    for (int i1=0; i1>=(1-mbc); i1--)
    for (int i2=(1-mbc); i2<=(T.getsize(2)+mbc); i2++)
    for (int meqi=1; meqi<=meqi_hgh; meqi++)
    {
      int meq; if(meq_ind) meq = meq_ind[meqi-1]; else meq = meqi;
      if(meq>dogParams.get_meqn()) invalid_value_error(meq); // continue;
      for (int i4=1; i4<=T.getsize(4); i4++)
          T.set(i1,i2,meq,i4,  getXreverse(i4)*T.get(1-i1,i2,meq,i4) );
    }
}
void SetRightMirrorBCs(dTensorBC4& T, const int* meq_ind)
{
    int meqi_hgh;
    if(meq_ind) meqi_hgh = get_idx_of_first_zero(meq_ind);
    else meqi_hgh = T.getsize(3);
    //
    int mbc = T.getmbc();
    int m1 = T.getsize(1);
    //
    for (int i1=(m1+1); i1<=(m1+mbc); i1++)
    for (int i2=(1-mbc); i2<=(T.getsize(2)+mbc); i2++)
    for (int meqi=1; meqi<=meqi_hgh; meqi++)
    {
      int meq; if(meq_ind) meq = meq_ind[meqi-1]; else meq = meqi;
      if(meq>dogParams.get_meqn()) invalid_value_error(meq); // continue;
      for (int i4=1; i4<=T.getsize(4); i4++)
          T.set(i1,i2,meq,i4,  getXreverse(i4)*T.get(2*m1+1-i1,i2,meq,i4) );
    }
}
void SetBottomMirrorBCs(dTensorBC4& T, const int* meq_ind)
{
    int meqi_hgh;
    if(meq_ind) meqi_hgh = get_idx_of_first_zero(meq_ind);
    else meqi_hgh = T.getsize(3);
    //
    int mbc = T.getmbc();
    for (int i1=(1-mbc); i1<=(T.getsize(1)+mbc); i1++)
    for (int i2=0; i2>=(1-mbc); i2--)
    for (int meqi=1; meqi<=meqi_hgh; meqi++)
    {
      int meq; if(meq_ind) meq = meq_ind[meqi-1]; else meq = meqi;
      if(meq>dogParams.get_meqn()) invalid_value_error(meq); // continue;
      for (int i4=1; i4<=T.getsize(4); i4++)
          T.set(i1,i2,meq,i4,  getYreverse(i4)*T.get(i1,1-i2,meq,i4) );
    }
}
void SetTopMirrorBCs(dTensorBC4& T, const int* meq_ind)
{
    int meqi_hgh;
    if(meq_ind) meqi_hgh = get_idx_of_first_zero(meq_ind);
    else meqi_hgh = T.getsize(3);
    //
    int mbc = T.getmbc();
    int m2 = T.getsize(2);
    //
    for (int i1=(1-mbc); i1<=(T.getsize(1)+mbc); i1++)
    for (int i2=(m2+1); i2<=(m2+mbc); i2++)
    for (int meqi=1; meqi<=meqi_hgh; meqi++)
    {
      int meq; if(meq_ind) meq = meq_ind[meqi-1]; else meq = meqi;
      if(meq>dogParams.get_meqn()) invalid_value_error(meq); // continue;
      for (int i4=1; i4<=T.getsize(4); i4++)
        T.set(i1,i2,meq,i4,  getYreverse(i4)*T.get(i1,2*m2+1-i2,meq,i4) );
    }
}

// periodic boundary conditions
//
void SetLeftPeriodicBCs(dTensorBC4& T)
{
    int mbc = T.getmbc();
    int m1 = T.getsize(1);
    int m2 = T.getsize(2);
    //
    for (int i1=0; i1>=(1-mbc); i1--)
    for (int i2=(1-mbc); i2<=(m2+mbc); i2++)
    for (int i3=1; i3<=T.getsize(3); i3++)
    for (int i4=1; i4<=T.getsize(4); i4++)
        T.set(i1,i2,i3,i4,  T.get(i1+m1,i2,i3,i4) );
}
void SetRightPeriodicBCs(dTensorBC4& T)
{
    int mbc = T.getmbc();
    int m1 = T.getsize(1);
    int m2 = T.getsize(2);
    //
    for (int i1=(m1+1); i1<=(m1+mbc); i1++)
    for (int i2=(1-mbc); i2<=(m2+mbc); i2++)
    for (int i3=1; i3<=T.getsize(3); i3++)
    for (int i4=1; i4<=T.getsize(4); i4++)
        T.set(i1,i2,i3,i4,  T.get(i1-m1,i2,i3,i4) );
}
void SetBottomPeriodicBCs(dTensorBC4& T)
{
    int mbc = T.getmbc();
    int m2 = T.getsize(2);
    //
    for (int i1=(1-mbc); i1<=(T.getsize(1)+mbc); i1++)
    for (int i2=0; i2>=(1-mbc); i2--)
    for (int i3=1; i3<=T.getsize(3); i3++)
    for (int i4=1; i4<=T.getsize(4); i4++)
        T.set(i1,i2,i3,i4,  T.get(i1,i2+m2,i3,i4) );
}
void SetTopPeriodicBCs(dTensorBC4& T)
{
    int mbc = T.getmbc();
    int m2 = T.getsize(2);
    //
    for (int i1=(1-mbc); i1<=(T.getsize(1)+mbc); i1++)
    for (int i2=(m2+1); i2<=(m2+mbc); i2++)
    for (int i3=1; i3<=T.getsize(3); i3++)
    for (int i4=1; i4<=T.getsize(4); i4++)
        T.set(i1,i2,i3,i4,  T.get(i1,i2-m2,i3,i4) );
}

// used for zero-order extrapolation
//
void LeftCopyBCs(dTensorBC4& T, const int* meq_ind)
{
    int meqi_hgh;
    if(meq_ind) meqi_hgh = get_idx_of_first_zero(meq_ind);
    else meqi_hgh = T.getsize(3);
    //
    int mbc = T.getmbc();
    for (int i1=0; i1>=(1-mbc); i1--)
    for (int i2=(1-mbc); i2<=(T.getsize(2)+mbc); i2++)
    for (int meqi=1; meqi<=meqi_hgh; meqi++)
    {
      int meq; if(meq_ind) meq = meq_ind[meqi-1]; else meq = meqi;
      if(meq>dogParams.get_meqn()) invalid_value_error(meq); // continue;
      // copy data from neighbor
      for (int i4=1; i4<=T.getsize(4); i4++)
          T.set(i1,i2,meq,i4,  T.get(1,i2,meq,i4) );
    }
}
void RightCopyBCs(dTensorBC4& T, const int* meq_ind)
{
    int meqi_hgh;
    if(meq_ind) meqi_hgh = get_idx_of_first_zero(meq_ind);
    else meqi_hgh = T.getsize(3);
    //
    int mbc = T.getmbc();
    int m1 = T.getsize(1);
    //
    for (int s1=1; s1<=mbc; s1++)
    for (int i2=(1-mbc); i2<=(T.getsize(2)+mbc); i2++)
    for (int meqi=1; meqi<=meqi_hgh; meqi++)
    {
      int meq; if(meq_ind) meq = meq_ind[meqi-1]; else meq = meqi;
      if(meq>dogParams.get_meqn()) invalid_value_error(meq); // continue;
      // copy data from neighbor
      for (int i4=2; i4<=T.getsize(4); i4++)
          T.set(m1+s1,i2,meq,i4,  T.get(m1,i2,meq,i4) );
    }
}
void BottomCopyBCs(dTensorBC4& T, const int* meq_ind)
{
    int meqi_hgh;
    if(meq_ind) meqi_hgh = get_idx_of_first_zero(meq_ind);
    else meqi_hgh = T.getsize(3);
    //
    int mbc = T.getmbc();
    //
    // BOTTOM BOUNDARY
    for (int i1=(1-mbc); i1<=(T.getsize(1)+mbc); i1++)
    for (int i2=0; i2>=(1-mbc); i2--)
    for (int meqi=1; meqi<=meqi_hgh; meqi++)
    {
      int meq; if(meq_ind) meq = meq_ind[meqi-1]; else meq = meqi;
      if(meq>dogParams.get_meqn()) invalid_value_error(meq); // continue;
      // copy data from neighbor
      for (int i4=1; i4<=T.getsize(4); i4++)
        T.set(i1,i2,meq,i4,  T.get(i1,1,meq,i4) );
    }
}
void TopCopyBCs(dTensorBC4& T, const int* meq_ind)
{
    int meqi_hgh;
    if(meq_ind) meqi_hgh = get_idx_of_first_zero(meq_ind);
    else meqi_hgh = T.getsize(3);
    //
    int mbc = T.getmbc();
    int m2 = T.getsize(2);
    //
    for (int i1=(1-mbc); i1<=(T.getsize(1)+mbc); i1++)
    for (int s2=1; s2<=mbc; s2++)
    for (int meqi=1; meqi<=meqi_hgh; meqi++)
    {
      int meq; if(meq_ind) meq = meq_ind[meqi-1]; else meq = meqi;
      if(meq>dogParams.get_meqn()) invalid_value_error(meq); // continue;
      // copy data from neighbor
      for (int i4=1; i4<=T.getsize(4); i4++)
        T.set(i1,m2+s2,meq,i4,  T.get(i1,m2,meq,i4) );
    }
}

// BC Helper Function
void setLeftRotationalBCs(dTensorBC4& q, const iTensor1* symmetryType=0)
{
  int mbc = q.getmbc();
  int mx  = q.getsize(1);
  int my  = q.getsize(2);
  if(symmetryType) assert_eq(symmetryType->getsize(),q.getsize(3));

  // LEFT BOUNDARY
  for (int ell=1; ell<=q.getsize(3); ell++)
  {
    int flip = (symmetryType && symmetryType->get(ell)==FLIP) ? -1 : 1;
    for (int i=0; i>=(1-mbc); i--)
    {
      int k = my+mbc;
      for (int j=(1-mbc); j<=(my+mbc); j++,k--)
      {
          switch(dogParams.get_space_order())
          {
          case 3:
          q.set(i,j,ell,6,  flip*q.get(1-i,k,ell,6) );
          q.set(i,j,ell,5,  flip*q.get(1-i,k,ell,5) );
          q.set(i,j,ell,4,  flip*q.get(1-i,k,ell,4) );
          case 2:
          q.set(i,j,ell,3, -flip*q.get(1-i,k,ell,3) );
          q.set(i,j,ell,2, -flip*q.get(1-i,k,ell,2) );
          case 1:
          q.set(i,j,ell,1,  flip*q.get(1-i,k,ell,1) );
          }
      }
    }
  }
}

// BC Helper Function
void setRightRotationalBCs(dTensorBC4& q, const iTensor1* symmetryType=0)
{
  int mbc = q.getmbc();
  int mx  = q.getsize(1);
  int my  = q.getsize(2);
  if(symmetryType) assert_eq(symmetryType->getsize(),q.getsize(3));

  // RIGHT BOUNDARY
  for (int ell=1; ell<=q.getsize(3); ell++)
  {
    for (int i=1; i<=mbc; i++)
    {
      int flip = (symmetryType && symmetryType->get(ell)==FLIP) ? -1 : 1;
      int k = my+mbc;
      for (int j=(1-mbc); j<=(my+mbc); j++,k--)
      {
        switch(dogParams.get_space_order())
        {
          case 3:
          q.set(mx+i,j,ell,6,  flip*q.get(mx+1-i,k,ell,6) );
          q.set(mx+i,j,ell,5,  flip*q.get(mx+1-i,k,ell,5) );
          q.set(mx+i,j,ell,4,  flip*q.get(mx+1-i,k,ell,4) );
          case 2:
          q.set(mx+i,j,ell,3, -flip*q.get(mx+1-i,k,ell,3) );
          q.set(mx+i,j,ell,2, -flip*q.get(mx+1-i,k,ell,2) );
          case 1:
          q.set(mx+i,j,ell,1,  flip*q.get(mx+1-i,k,ell,1) );
        }
      }
    }
  }
}

// BC Helper Function
void BCHelper(BoundaryID bdID, const iTensor1& symmetryType, dTensorBC4& q)
{
    int i,j,k,ell,m;
    //int msy = symm.getsize();
    //int mfl = flip.getsize();
    int mbc = q.getmbc();
    int mx  = q.getsize(1);
    int my  = q.getsize(2);    

    switch( bdID )
    {
      case LEFT: // LEFT BOUNDARY
        for (ell=1; ell<=symmetryType.getsize(); ell++)
        {
          int flip = (symmetryType.get(ell)==FLIP) ? -1 : 1;
          for (i=0; i>=(1-mbc); i--)
            for (j=(1-mbc); j<=(my+mbc); j++)
              {
                switch(dogParams.get_space_order())
                {
                case 3:
                q.set(i,j,ell,6,  flip*q.get(1-i,j,ell,6) );
                q.set(i,j,ell,5,  flip*q.get(1-i,j,ell,5) );
                q.set(i,j,ell,4, -flip*q.get(1-i,j,ell,4) );
                case 2:
                q.set(i,j,ell,3,  flip*q.get(1-i,j,ell,3) );
                q.set(i,j,ell,2, -flip*q.get(1-i,j,ell,2) );
                case 1:
                q.set(i,j,ell,1,  flip*q.get(1-i,j,ell,1) );
                }
              }
        }
        break;

      case RIGHT: // RIGHT BOUNDARY
        for (ell=1; ell<=symmetryType.getsize(); ell++)
          for (i=1; i<=mbc; i++)
            for (j=(1-mbc); j<=(my+mbc); j++)
              {
                int flip = (symmetryType.get(ell)==FLIP) ? -1 : 1;
                switch(dogParams.get_space_order())
                {
                case 3:
                q.set(mx+i,j,ell,6,  flip*q.get(mx+1-i,j,ell,6) );
                q.set(mx+i,j,ell,5,  flip*q.get(mx+1-i,j,ell,5) );
                q.set(mx+i,j,ell,4, -flip*q.get(mx+1-i,j,ell,4) );
                case 2:
                q.set(mx+i,j,ell,3,  flip*q.get(mx+1-i,j,ell,3) );
                q.set(mx+i,j,ell,2, -flip*q.get(mx+1-i,j,ell,2) );
                case 1:
                q.set(mx+i,j,ell,1,  flip*q.get(mx+1-i,j,ell,1) );
                }
              }
        break;

      case BOTTOM: // BOTTOM BOUNDARY
        for (ell=1; ell<=symmetryType.getsize(); ell++)
          for (j=0; j>=(1-mbc); j--)
            for (i=(1-mbc); i<=(mx+mbc); i++)
              {
                int flip = (symmetryType.get(ell)==FLIP) ? -1 : 1;
                switch(dogParams.get_space_order())
                {
                case 3:
                q.set(i,j,ell,6,  flip*q.get(i,1-j,ell,6) );
                q.set(i,j,ell,5,  flip*q.get(i,1-j,ell,5) );
                q.set(i,j,ell,4, -flip*q.get(i,1-j,ell,4) );
                case 2:
                q.set(i,j,ell,3, -flip*q.get(i,1-j,ell,3) );
                q.set(i,j,ell,2,  flip*q.get(i,1-j,ell,2) );
                case 1:
                q.set(i,j,ell,1,  flip*q.get(i,1-j,ell,1) );
                }
              }
        break;

      case TOP: // TOP BOUNDARY
        for (ell=1; ell<=symmetryType.getsize(); ell++)
          for (j=1; j<=mbc; j++)
            for (i=(1-mbc); i<=(mx+mbc); i++)
              {
                int flip = (symmetryType.get(ell)==FLIP) ? -1 : 1;
                switch(dogParams.get_space_order())
                {
                case 3:
                q.set(i,my+j,ell,6,  flip*q.get(i,my+1-j,ell,6) );
                q.set(i,my+j,ell,5,  flip*q.get(i,my+1-j,ell,5) );
                q.set(i,my+j,ell,4, -flip*q.get(i,my+1-j,ell,4) );
                case 2:
                q.set(i,my+j,ell,3, -flip*q.get(i,my+1-j,ell,3) );
                q.set(i,my+j,ell,2,  flip*q.get(i,my+1-j,ell,2) );
                case 1:
                q.set(i,my+j,ell,1,  flip*q.get(i,my+1-j,ell,1) );
                }
              }
        break;
    }
}

void SetYbndValues(
    dTensorBC4& q, dTensorBC4& aux,
    const int enforced_symmetry,
    const BoundaryConditionType BCs,
    const iTensor1& y_cond_wall_symmetryType,
    const iTensor1& y_mirror_symmetryType,
    const iTensor1& y_forced_symmetryType)
{
  // top and bottom walls
  if (enforced_symmetry&Y) // Y-symmetry enforced
  {
    assert_eq(gemParams.get_B_guide(), 0.);
    assert(gemParams.get_enforcing_flip_symmetry());

    // mirror bottom boundary
    BCHelper(BOTTOM,y_mirror_symmetryType,q); //bottom boundary
    SetBottomMirrorBCs(aux);

    // top boundary
    //
    switch(BCs)
    {
      case CONDUCTING_WALL:
        BCHelper(TOP,y_cond_wall_symmetryType,q);
        SetTopMirrorBCs(aux);
        break;
      case OPEN_BOUNDARY:
        SetTopMirrorBCs(q);
        SetTopMirrorBCs(aux);
        break;
      case COPY_BOUNDARY:
        TopCopyBCs(q);
        TopCopyBCs(aux);
        break;
      case FORCED:
        unsupported_value_error(BCs);
        // subtract off desired values/slopes
        // apply mirroring
        BCHelper(TOP,y_forced_symmetryType,q);
        SetTopMirrorBCs(aux);
        // add desired values/slopes back on
        unsupported_value_error(BCs);
        break;
      default:
        unsupported_value_error(BCs);
    }
  }
  else // Y-symmetry not enforced
  {
    // top and bottom boundaries
    //
    switch(BCs)
    {
      case CONDUCTING_WALL:
        BCHelper(BOTTOM,y_cond_wall_symmetryType,q);
        BCHelper(TOP,   y_cond_wall_symmetryType,q);
        SetBottomMirrorBCs(aux);
        SetTopMirrorBCs(aux);
        break;
      case OPEN_BOUNDARY:
        SetTopMirrorBCs(q);
        SetTopMirrorBCs(aux);
        SetBottomMirrorBCs(q);
        SetBottomMirrorBCs(aux);
        break;
      case COPY_BOUNDARY:
        TopCopyBCs(q);
        TopCopyBCs(aux);
        BottomCopyBCs(q);
        BottomCopyBCs(aux);
        break;
      case PERIODIC_ON_VERTICALLY_DOUBLED_DOMAIN:
        //BCHelper(TOP,y_periodic_symmetryType,q);
        unsupported_value_error(BCs);
        SetTopMirrorBCs(aux); // wrong BCs
      case FORCED:
        unsupported_value_error(BCs);
        break;
      default:
        unsupported_value_error(BCs);
        break;
    }
  }
}

void SetXbndValues(
    dTensorBC4& q, dTensorBC4& aux,
    const int enforced_symmetry,
    const BoundaryConditionType BCs,
    const iTensor1& x_mirror_symmetryType,
    const iTensor1& rotational_symmetries)
{
    // left and right boundaries
    //
    if (enforced_symmetry&X)
    {
      if(gemParams.get_enforcing_rotational_symmetry())
      {
        setLeftRotationalBCs(q,&rotational_symmetries);
        setLeftRotationalBCs(aux);
        if(true){ // correct
        setRightRotationalBCs(q,&rotational_symmetries);
        setRightRotationalBCs(aux);
        } else { // incorrect (creates shock, tests positivity limiters)
        BCHelper(RIGHT,x_mirror_symmetryType,q);
        SetRightMirrorBCs(aux);
        }
      }
      else
      {
        // mirror left boundary
        BCHelper(LEFT, x_mirror_symmetryType, q);
        SetLeftMirrorBCs(aux);

        switch(BCs)
        {
          case PERIODIC:
            // mirror right boundary
            BCHelper(RIGHT,x_mirror_symmetryType,q);
            SetRightMirrorBCs(aux);
            break;
          case OPEN_BOUNDARY:
            SetRightMirrorBCs(q);
            SetRightMirrorBCs(aux);
            break;
          default:
            unsupported_value_error(BCs);
        }
      }
    }
    else // X-symmetry not enforced
    {
      switch(BCs)
      {
       case PERIODIC:
        assert(!gemParams.get_enforcing_rotational_symmetry());
        // periodic left boundary
        SetLeftPeriodicBCs(q);
        SetLeftPeriodicBCs(aux);

        // periodic right boundary
        SetRightPeriodicBCs(q);
        SetRightPeriodicBCs(aux);
        break;
       case OPEN_BOUNDARY:
        SetLeftMirrorBCs(q);
        SetLeftMirrorBCs(aux);
        SetRightMirrorBCs(q);
        SetRightMirrorBCs(aux);
        break;
       default:
        unsupported_value_error(BCs);
      }
    }
}

void wrap_periodic_ghost_cell_values_at_botm_boundary(
    dTensorBC4& q, int* rho_indices, int* y_indices, double y_period_2)
{
  if (gemParams.get_enforced_symmetry()&Y) return;
  //
  // Walk along the bottom boundary and subtract a multiple of the
  // period from the average value to make it close to the
  // average value of the neighboring element
  //
  const int mx  = q.getsize(1);
  const int my  = q.getsize(2);
  const double y_period = 2*y_period_2;
  const int mbc = q.getmbc();
  //
  int ell=0;
  for (int idx=0; ell=y_indices[idx]; idx++)
  {
    int rho_ = rho_indices[idx];
    for (int i=(1-mbc); i<=(mx+mbc); i++)
    {
      double qval_neighbor = q.get(i,my,ell,1);
      for (int j=1; j<=mbc; j++)
      {
        //double qval = q.get(i,1-j,ell,1);
        //double rho_val = q.get(i,1-j,rho_,1);
        //qval -= rho_val*y_period;
        //assert_le(abs(qval-qval_neighbor),rho_val*y_period_2);
        //q.set(i,1-j,ell,1,  qval);
        q.set(i,1-j,ell,1,  qval_neighbor);
        for(int k=2;k<=q.getsize(4);k++) q.set(i,1-j,ell,k,  0.);
      }
    }
  }
}

void wrap_periodic_ghost_cell_values_at_top__boundary(
    dTensorBC4& q, int* rho_indices, int* y_indices, double y_period_2)
{
  //
  // Walk along the top boundary and add a multiple of the
  // period to the average value to make it close to the
  // average value of the neighboring element
  //
  const int mx  = q.getsize(1);
  const int my  = q.getsize(2);
  const double y_period = 2*y_period_2;
  const int mbc = q.getmbc();
  //
  int ell=0;
  for (int idx=0; ell=y_indices[idx]; idx++)
  {
    int rho_ = rho_indices[idx];
    // const double bi = q.get(1,my,  ell,1);
    // const double bo = q.get(1,my+1,ell,1);
    // // find the value of n for which |n*y_period+bo-bi| <= y_period/2
    // double diff = bo-bi;
    // // the same shift should work for all points
    // // and should not be bigger than the y_period
    // double shift = 0.;
    // if(diff>y_period_2)
    // {
    //   shift = -y_period;
    // }
    // else if(diff<-y_period_2)
    // {
    //   shift = y_period;
    // }
    // diff+=shift;
    // assert_le(diff, y_period_2);
    // assert_ge(diff,-y_period_2);

    for (int i=(1-mbc); i<=(mx+mbc); i++)
    {
      double qval_inside = q.get(i,my,ell,1);
      for (int j=1; j<=mbc; j++)
      {
        double qval = q.get(i,my+j,ell,1);
        double rho_val = q.get(i,my+j,rho_,1);
        qval += rho_val*y_period;
        assert_le(fabs(qval-qval_inside),rho_val*y_period_2);
        q.set(i,my+j,ell,1,  qval);
      }
    }
  }
}

void wrap_periodic_ghost_cell_values_at_left_boundary(
    dTensorBC4& q, int* rho_indices, int* x_indices, double x_period_2)
{
  if (gemParams.get_enforced_symmetry()&X) return;
  //
  // Walk along the left boundary and subtract a multiple of the
  // period from the average value to make it close to the
  // average value of the neighboring element
  //
  const int mx  = q.getsize(1);
  const int my  = q.getsize(2);
  const double x_period = 2*x_period_2;
  const int mbc = q.getmbc();
  //
  int ell=0;
  for (int idx=0; ell=x_indices[idx]; idx++)
  {
    int rho_ = rho_indices[idx];
    for (int j=(1-mbc); j<=(my+mbc); j++)
    {
      double qval_neighbor = q.get(mx,j,ell,1);
      for (int i=1; i<=mbc; i++)
      {
        //double qval = q.get(1-i,j,ell,1);
        //double rho_val = q.get(1-i,j,rho_,1);
        //qval -= rho_val*x_period;
        //assert_le(abs(qval-qval_neighbor),rho_val*x_period_2);
        //q.set(1-i,j,ell,1,  qval);
        q.set(1-i,j,ell,1,  qval_neighbor);
        for(int k=2;k<=q.getsize(4);k++)
        {
          q.set(1-i,j,ell,k,  0.);
        }
      }
    }
  }
}
void wrap_periodic_ghost_cell_values_at_rght_boundary(
    dTensorBC4& q, int* rho_indices, int* x_indices, double x_period_2)
{
  //
  // Walk along the right boundary and add a multiple of the
  // period to the average value to make it close to the
  // average value of the neighboring element
  //
  const int mx  = q.getsize(1);
  const int my  = q.getsize(2);
  const double x_period = 2*x_period_2;
  const int mbc = q.getmbc();
  //
  int ell=0;
  for (int idx=0; ell=x_indices[idx]; idx++)
  {
    int rho_ = rho_indices[idx];
    // const double bi = q.get(mx,  1,ell,1);
    // const double bo = q.get(mx+1,1,ell,1);
    // // find the value of n for which |n*x_period+bo-bi| <= x_period/2
    // double diff = bo-bi;
    // // the same shift should work for all points
    // // and should not be bigger than the x_period
    // double shift = 0.;
    // if(diff>x_period_2)
    // {
    //   shift = -x_period;
    // }
    // else if(diff<-x_period_2)
    // {
    //   shift = x_period;
    // }
    // diff+=shift;
    // assert_le(diff, x_period_2);
    // assert_ge(diff,-x_period_2);

    // I do not think that this is high-order-accurate.
    // Should we shift using average or cell-center values?
    for (int j=(1-mbc); j<=(my+mbc); j++)
    {
      double qval_neighbor = q.get(mx,j,ell,1);
      for (int i=1; i<=mbc; i++)
      {
        // double qval = q.get(mx+i,j,ell,1);
        // double rho_val = q.get(mx+i,j,rho_,1);
        // qval += rho_val*x_period;
        // //if(j>=3 && j<=my && (qval-qval_neighbor)>.20*rho_val*x_period_2)
        // //{
        // //  dprint(qval);
        // //  dprint(qval_neighbor);
        // //  dprint(i);
        // //  dprint(j);
        // //}
        // q.set(mx+i,j,ell,1,  qval);
        q.set(mx+i,j,ell,1,  qval_neighbor);
        for(int k=2;k<=q.getsize(4);k++) q.set(mx+i,j,ell,k,  0.);
        // qval_neighbor = qval;
      }
    }
  }
}
