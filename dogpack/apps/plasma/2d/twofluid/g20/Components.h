#define _COMPONENTS_H_

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
    _rho_e=21, // rho (elc)
    _M1_e =22, // 1-momentum (elc)
    _M2_e =23, // 2-momentum (elc)
    _M3_e =24, // 3-momentum (elc)
    _N11_e=25, // energy_11 (elc)
    _N12_e=26, // energy_12 (elc)
    _N13_e=27, // energy_13 (elc)
    _N22_e=28, // energy_22 (elc)
    _N23_e=29, // energy_23 (elc)
    _N33_e=30, // energy_33 (elc)
    _Q111_e=31,
    _Q112_e=32,
    _Q113_e=33,
    _Q122_e=34,
    _Q123_e=35,
    _Q133_e=36,
    _Q222_e=37,
    _Q223_e=38,
    _Q233_e=39,
    _Q333_e=40,
    _B1   =41, // 1-magnetic field
    _B2   =42, // 2-magnetic field
    _B3   =43, // 3-magnetic field
    _E1   =44, // 1-electric field
    _E2   =45, // 2-electric field
    _E3   =46, // 3-electric field
    _psi  =47, // B-field cleaning
    _phi  =48, // E-field cleaning
    _entropy_i = 49, // entropy tracking
    _entropy_e = 50, // entropy tracking

    // declare primitive variable equivalents

    _u1_i  = _M1_i,
    _u2_i  = _M2_i,
    _u3_i  = _M3_i,
    _P11_i = _N11_i,
    _P12_i = _N12_i,
    _P13_i = _N13_i,
    _P22_i = _N22_i,
    _P23_i = _N23_i,
    _P33_i = _N33_i,

    _u1_e  = _M1_e,
    _u2_e  = _M2_e,
    _u3_e  = _M3_e,
    _P11_e = _N11_e,
    _P12_e = _N12_e,
    _P13_e = _N13_e,
    _P22_e = _N22_e,
    _P23_e = _N23_e,
    _P33_e = _N33_e,
};


