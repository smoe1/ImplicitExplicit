#define COMPONENTS_H

enum ComponentID{
    _rho_i=1,  // rho (ion)
    _M1_i =2,  // 1-momentum (ion)
    _M2_i =3,  // 2-momentum (ion)
    _M3_i =4,  // 3-momentum (ion)
    _N11_i=5,  // energy_11 (ion)
    _N12_i=6,  // energy_12 (ion)
    _N13_i=7,  // energy_13 (ion)
    _N22_i=8,  // energy_22 (ion)
    _N23_i=9,  // energy_23 (ion)
    _N33_i=10, // energy_33 (ion)
    _B1   =11, // 1-magnetic field
    _B2   =12, // 2-magnetic field
    _E3   =13, // 3-electric field
    _psi  =14, // B-field cleaning
    _entropy_i = 15, // entropy tracking
};

