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
    _Q111_i=11,
    _Q112_i=12,
    _Q113_i=13,
    _Q122_i=14,
    _Q123_i=15,
    _Q133_i=16,
    _Q222_i=17,
    _Q223_i=18,
    _Q233_i=19,
    _Q333_i=20,
    _B1   =21, // 1-magnetic field
    _B2   =22, // 2-magnetic field
    _E3   =23, // 3-electric field
    _psi  =24, // B-field cleaning
    _entropy_i = 25, // entropy tracking
};

