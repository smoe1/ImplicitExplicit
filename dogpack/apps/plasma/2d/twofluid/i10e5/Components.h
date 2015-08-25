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
    _rho_e=11, // rho (elc)
    _M1_e =12, // 1-momentum (elc)
    _M2_e =13, // 2-momentum (elc)
    _M3_e =14, // 3-momentum (elc)
    _N_e  =15, // energy (elc)
    _B1   =16, // 1-magnetic field
    _B2   =17, // 2-magnetic field
    _B3   =18, // 3-magnetic field
    _E1   =19, // 1-electric field
    _E2   =20, // 2-electric field
    _E3   =21, // 3-electric field
    _psi  =22, // B-field cleaning
    _phi  =23, // E-field cleaning
    _entropy_i = 24, // entropy tracking
    _entropy_e = 25, // entropy tracking
};

