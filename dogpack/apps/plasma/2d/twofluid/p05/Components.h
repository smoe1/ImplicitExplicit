#define _COMPONENTS_H_

enum ComponentID{
    _rho_i=1,  // rho (ion)
    _M1_i =2,  // 1-momentum (ion)
    _M2_i =3,  // 2-momentum (ion)
    _M3_i =4,  // 3-momentum (ion)
    _N_i  =5,  // energy (ion)
    _B1   =6, // 1-magnetic field
    _B2   =7, // 2-magnetic field
    _E3   =8, // 3-electric field
    _psi  =9, // B-field cleaning
    _x_i  =10, // x-coordinate tracking
    _y_i  =11, // y-coordinate tracking
    _entropy_i = 12, // entropy tracking
};

