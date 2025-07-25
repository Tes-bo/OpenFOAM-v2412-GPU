inline scalar Cp
(
    const scalar p,
    const scalar T
) const
{
    #ifdef __clang__
    volatile const scalar cv = Cv(p, T);
    volatile const scalar cpmcv = EquationOfState::CpMCv(p, T);
    return cv + cpmcv;
    #else
    return Cv(p, T) + EquationOfState::CpMCv(p, T);
    #endif
}

inline scalar Hs
(
    const scalar p,
    const scalar T
) const
{
    #ifdef __clang__
    volatile const scalar es = Es(p, T);
    volatile const scalar rho = EquationOfState::rho(p, T);
    return es + p/rho;
    #else
    return Es(p, T) + p/EquationOfState::rho(p, T);
    #endif
}

inline scalar Ha
(
    const scalar p,
    const scalar T
) const
{
    #ifdef __clang__
    volatile const scalar ea = Ea(p, T);
    volatile const scalar rho = EquationOfState::rho(p, T);
    return ea + p/rho;
    #else
    return Ea(p, T) + p/EquationOfState::rho(p, T);
    #endif
}
