#include<iostream>
#include "tensors.h"
#include "DogParams.h"
#include "GEMparams.h"
#include "Components.h"

using namespace std;

void setRotationalSymmetries(iTensor1& symmetryType)
{
    symmetryType.setall(SYMM);
    int flip_indices[] = {
        _M1_i , _M2_i, _N13_i, _N23_i,
        _M1_e , _M2_e,
        _B1   , _B2   , _E1   , _E2  , 0};
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
    int flip_indices[] = {
        _M1_i , _N12_i, _N13_i,
        _M1_e ,
        _B2   , _B3   , _E1   , _psi  , 0};
    for(i=0;flip_indices[i];i++) symmetryType.set(flip_indices[i],FLIP);
}

// BOTTOM AND TOP BOUNDARIES (symmetry boundaries)

// mirror boundary for enforced symmetry,
// top also for periodic and enforced symmetry
void setBottomTopSymmetryBoundaries(iTensor1& symmetryType)
{
    int i;
    for(i=1;i<=dogParams.get_meqn();i++) symmetryType.set(i,SYMM);
    int flip_indices[] = {
        _M2_i , _N12_i, _N23_i,
        _M2_e ,
        _B1   , _B3   , _E2   , _psi  , 0};
    for(i=0;flip_indices[i];i++)
      if(flip_indices[i] <= symmetryType.getsize())
        symmetryType.set(flip_indices[i],FLIP);
}

// CONDUCTING_WALL
void setBottomTopConductingWallBoundaries(iTensor1& symmetryType)
{
    int i;
    // initialize to symmetric
    for(i=1;i<=dogParams.get_meqn();i++) symmetryType.set(i,SYMM);
    int flip_indices[] = {
        _M2_i , _N12_i, _N23_i,
        _M2_e ,
        _B2   , _E1   , _E3   , _phi  , 0};
    for(i=0;flip_indices[i];i++)
      if(flip_indices[i] <= symmetryType.getsize())
        symmetryType.set(flip_indices[i],FLIP);
}

// PERIODIC_ON_VERTICALLY_DOUBLED_DOMAIN
void setBottomTopPeriodicBoundaries(iTensor1& symmetryType)
{
    int i;
    for(i=1;i<=dogParams.get_meqn();i++) symmetryType.set(i,SYMM);
    int flip_indices[] = {
        _M2_i , _M3_i , _N12_i, _N13_i,
        _M2_e , _M3_e ,
        _B2   , _B3   , _E2   , _E3   , 0};
    for(i=0;flip_indices[i];i++)
      if(flip_indices[i] <= symmetryType.getsize())
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
}
