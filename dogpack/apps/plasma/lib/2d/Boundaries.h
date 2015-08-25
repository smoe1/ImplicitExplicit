#ifndef _BOUNDARIES_H_
#define _BOUNDARIES_H_
#include "tensors.h"

// 1st bit says low/high
// 2nd bit says hor/vert
enum BoundaryID{LEFT=0, RIGHT=1, BOTTOM=2, TOP=3};

enum SymmetryType{SYMM, FLIP};

enum BoundaryConditionType
{
  CONDUCTING_WALL=1,
  PERIODIC=2,
  PERIODIC_ON_VERTICALLY_DOUBLED_DOMAIN=3,
  OPEN_BOUNDARY=4,
  COPY_BOUNDARY=5,
  FORCED=6
};

// 1st bit says whether symmetry is enforced in X axis,
// 2nd bit says whether symmetry is enforced in Y axis.
// with such short names I ought to wrap this in a struct
// or a namespace, e.g.
//
// struct EnforcedSymmetryType { enum e {NONE=0,X=1,Y=2,XY=3}};
//   or
// namespace EnforcedSymmetryType { enum e {NONE=0,X=1,Y=2,XY=3}};
//
enum EnforcedSymmetryType
{
  NONE=0,
  X=1,
  Y=2,
  XY=3
};

// void SetBndValues(
//     dTensorBC4& q, dTensorBC4& aux,
//     const int enforced_symmetry,
//     const BoundaryConditionType BCs,
//     const iTensor1& x_mirror_symmetryType,
//     const iTensor1& y_cond_wall_symmetryType,
//     const iTensor1& y_mirror_symmetryType,
//     const iTensor1& y_periodic_symmetryType,
//     const iTensor1& rotational_symmetries);

void SetYbndValues(
    dTensorBC4& q, dTensorBC4& aux,
    const int enforced_symmetry,
    const BoundaryConditionType BCs,
    const iTensor1& y_cond_wall_symmetryType,
    const iTensor1& y_mirror_symmetryType,
    const iTensor1& y_forced_symmetryType);

void SetXbndValues(
    dTensorBC4& q, dTensorBC4& aux,
    const int enforced_symmetry,
    const BoundaryConditionType BCs,
    const iTensor1& x_mirror_symmetryType,
    const iTensor1& rotational_symmetries);

void SetLeftMirrorBCs  (dTensorBC4& T, const int* meq_ind=0);
void SetBottomMirrorBCs(dTensorBC4& T, const int* meq_ind=0);
void SetRightMirrorBCs (dTensorBC4& T, const int* meq_ind=0);
void SetTopMirrorBCs   (dTensorBC4& T, const int* meq_ind=0);

void LeftCopyBCs  (dTensorBC4& T, const int* meq_ind=0);
void RightCopyBCs (dTensorBC4& T, const int* meq_ind=0);
void BottomCopyBCs(dTensorBC4& T, const int* meq_ind=0);
void TopCopyBCs   (dTensorBC4& T, const int* meq_ind=0);

void wrap_periodic_ghost_cell_values_at_rght_boundary(
    dTensorBC4& q, int* rho_indices, int* x_indices, double x_period_2);
void wrap_periodic_ghost_cell_values_at_left_boundary(
    dTensorBC4& q, int* rho_indices, int* x_indices, double x_period_2);
void wrap_periodic_ghost_cell_values_at_top__boundary(
    dTensorBC4& q, int* rho_indices, int* y_indices, double y_period_2);
void wrap_periodic_ghost_cell_values_at_botm_boundary(
    dTensorBC4& q, int* rho_indices, int* y_indices, double y_period_2);

#endif
