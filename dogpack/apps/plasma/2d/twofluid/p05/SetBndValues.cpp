// lib headers
#include "tensors.h"
#include "DogParams.h"
#include "DogParamsCart2.h"
// local headers
#include "GEMparams.h"
#include "Components.h"
// plasma/lib headers
#include "Boundaries.h"
#include "debug.h"

void setRotationalSymmetries(iTensor1& symmetryType)
{
    symmetryType.setall(SYMM);
    int flip_indices[] = {
        _M1_i , _M2_i,
        _B1   , _B2  ,
        _x_i  , _y_i , 0};
    for(int i=0;flip_indices[i];i++)
      if(flip_indices[i] <= symmetryType.getsize())
        symmetryType.set(flip_indices[i],FLIP);
}

// LEFT AND RIGHT BOUNDARIES (symmetry boundaries)
void setLeftRightSymmetryBoundaries(iTensor1& symmetryType)
{
    int i;
    // initialize to symmetric
    for(i=1;i<=dogParams.get_meqn();i++) symmetryType.set(i,SYMM);
    int flip_indices[] = {_M1_i, _B2, _psi, _x_i, 0};
    for(i=0;flip_indices[i];i++)
      if(flip_indices[i]<=symmetryType.getsize())
        symmetryType.set(flip_indices[i],FLIP);
}

// BOTTOM AND TOP BOUNDARIES (symmetry boundaries)

// mirror boundary for enforced symmetry,
// top also for periodic and enforced symmetry
void setBottomTopSymmetryBoundaries(iTensor1& symmetryType)
{
    int i;
    for(i=1;i<=dogParams.get_meqn();i++) symmetryType.set(i,SYMM);
    int flip_indices[] = {_M2_i, _B1, _psi, _y_i, 0};
    for(i=0;flip_indices[i];i++)
      if(flip_indices[i]<=symmetryType.getsize())
        symmetryType.set(flip_indices[i],FLIP);
}

// CONDUCTING_WALL
void setBottomTopConductingWallBoundaries(iTensor1& symmetryType)
{
    int i;
    // initialize to symmetric
    for(i=1;i<=dogParams.get_meqn();i++) symmetryType.set(i,SYMM);
    int flip_indices[] = {_M2_i, _B2, _E3, _y_i, 0};
    for(i=0;flip_indices[i];i++)
      if(flip_indices[i]<=symmetryType.getsize())
        symmetryType.set(flip_indices[i],FLIP);
}

// PERIODIC_ON_VERTICALLY_DOUBLED_DOMAIN
void setBottomTopPeriodicBoundaries(iTensor1& symmetryType)
{
    int i;
    for(i=1;i<=dogParams.get_meqn();i++) symmetryType.set(i,SYMM);
    int flip_indices[] = {_M2_i, _M3_i, _B2, _E3, _y_i, 0};
    for(i=0;flip_indices[i];i++)
      if(flip_indices[i]<=symmetryType.getsize())
        symmetryType.set(flip_indices[i],FLIP);
}

void SetBndValues(dTensorBC4& q, dTensorBC4& aux)
{
    const int& meqn = dogParams.get_meqn();
    // LEFT AND RIGHT BOUNDARIES
    //
    // (x=const) mirror boundary (enforced_symmetry: X or XY)
    iTensor1 x_mirror_symmetryType(meqn);
    setLeftRightSymmetryBoundaries(x_mirror_symmetryType);

    // TOP AND BOTTOM BOUNDARIES
    //
    // (y=const) mirror boundary (enforced_symmetry: Y or XY)
    iTensor1 y_mirror_symmetryType(meqn);
    setBottomTopSymmetryBoundaries(y_mirror_symmetryType);
    //
    // (y=const) (CONDUCTING_WALL)
    iTensor1 y_cond_wall_symmetryType(meqn);
    setBottomTopConductingWallBoundaries(y_cond_wall_symmetryType);
    //
    // (y=const) top boundary (PERIODIC_ON_VERTICALLY_DOUBLED_DOMAIN)
    iTensor1 y_periodic_symmetryType(meqn);
    setBottomTopPeriodicBoundaries(y_periodic_symmetryType);

    // used when imposing rotational symmetry (for nonzero guide field)
    //
    iTensor1 rotational_symmetries(meqn);
    setRotationalSymmetries(rotational_symmetries);

    SetYbndValues(
        q, aux,
        gemParams.get_enforced_symmetry(),
        gemParams.get_yBCs(),
        y_cond_wall_symmetryType,
        y_mirror_symmetryType,
        y_periodic_symmetryType);

    SetXbndValues(
        q, aux,
        gemParams.get_enforced_symmetry(),
        gemParams.get_xBCs(),
        x_mirror_symmetryType,
        rotational_symmetries);

    if(dogParams.get_meqn() < _y_i) return;

    int xy_indices[] = {_x_i, _y_i, 0};
    LeftCopyBCs  (q, xy_indices);
    BottomCopyBCs(q, xy_indices);
    RightCopyBCs (q, xy_indices);
    TopCopyBCs   (q, xy_indices);
    #if 0
    // coordinates are periodic, i.e. multivalued,
    // so when using periodicity to set coordinate
    // boundary values we have to add or subtract one period.
    int x_indices[] = {_x_i, 0};
    int y_indices[] = {_y_i, 0};
    int rho_indices[] = {_rho_i, 0};
    const double x_period_2 = dogParamsCart2.get_xhigh();
    const double y_period_2 = dogParamsCart2.get_yhigh();
    wrap_periodic_ghost_cell_values_at_left_boundary(
      q, rho_indices, x_indices, x_period_2);
    wrap_periodic_ghost_cell_values_at_botm_boundary(
      q, rho_indices, y_indices, y_period_2);
    wrap_periodic_ghost_cell_values_at_rght_boundary(
      q, rho_indices, x_indices, x_period_2);
    wrap_periodic_ghost_cell_values_at_top__boundary(
      q, rho_indices, y_indices, y_period_2);
    #endif
}
