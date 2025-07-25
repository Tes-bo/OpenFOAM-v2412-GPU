/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::cubicEqn

Description
    Container to encapsulate various operations for
    cubic equation of the forms with real coefficients:

    \f[
        a*x^3 + b*x^2 + c*x + d = 0
          x^3 + B*x^2 + C*x + D = 0
    \f]

    The following two substitutions into the above forms are used:

    \f[
        x = t - B/3
        t = w - P/3/w
    \f]

    This reduces the problem to a quadratic in \c w^3:

    \f[
        w^6 + Q*w^3 - P = 0
    \f]

    where \c Q and \c P are given in the code.

    The properties of the cubic can be related to the properties of this
    quadratic in \c w^3.

    If the quadratic eqn has two identical real roots at zero, three identical
    real roots exist in the cubic eqn.

    If the quadratic eqn has two identical real roots at non-zero, two identical
    and one distinct real roots exist in the cubic eqn.

    If the quadratic eqn has two complex roots (a complex conjugate pair),
    three distinct real roots exist in the cubic eqn.

    If the quadratic eqn has two distinct real roots, one real root and two
    complex roots (a complex conjugate pair) exist in the cubic eqn.

    The quadratic eqn is solved for the most numerically accurate value
    of \c w^3. See the \link quadraticEqn.H \endlink for details on how to
    pick a value. This single value of \c w^3 can yield up to three cube
    roots for \c w, which relate to the three solutions for \c x.

    Only a single root, or pair of conjugate roots, is directly evaluated; the
    one, or ones with the lowest relative numerical error. Root identities are
    then used to recover the remaining roots, possibly utilising a quadratic
    and/or linear solution. This seems to be a good way of maintaining the
    accuracy of roots at very different magnitudes.

    Reference:
    \verbatim
        Kahan's algo. to compute 'p' using fused multiply-adds (tag:JML):
            Jeannerod, C. P., Louvet, N., & Muller, J. M. (2013).
            Further analysis of Kahan's algorithm for the accurate
            computation of 2× 2 determinants.
            Mathematics of Computation, 82(284), 2245-2264.
            DOI:10.1090/S0025-5718-2013-02679-8
    \endverbatim

See also
    Test-cubicEqn.C

SourceFiles
    cubicEqnI.H
    cubicEqn.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_cubicEqn_H
#define Foam_cubicEqn_H

#include "Roots.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class cubicEqn Declaration
\*---------------------------------------------------------------------------*/

class cubicEqn
:
    public VectorSpace<cubicEqn, scalar, 4>
{
public:

    //- Component labeling enumeration
    enum components { A, B, C, D };


    // Constructors

        //- Default construct
        cubicEqn() = default;

        //- Construct initialized to zero
        inline cubicEqn(const Foam::zero);

        //- Construct from components
        inline cubicEqn
        (
            const scalar a,
            const scalar b,
            const scalar c,
            const scalar d
        );


    // Member Functions

    // Access

        scalar a() const noexcept { return this->v_[A]; }
        scalar b() const noexcept { return this->v_[B]; }
        scalar c() const noexcept { return this->v_[C]; }
        scalar d() const noexcept { return this->v_[D]; }

        scalar& a() noexcept { return this->v_[A]; }
        scalar& b() noexcept { return this->v_[B]; }
        scalar& c() noexcept { return this->v_[C]; }
        scalar& d() noexcept { return this->v_[D]; }


    // Evaluate

        //- Evaluate the cubic equation at x
        inline scalar value(const scalar x) const;

        //- Evaluate the derivative of the cubic equation at x
        inline scalar derivative(const scalar x) const;

        //- Estimate the error of evaluation of the cubic equation at x
        inline scalar error(const scalar x) const;

        //- Return the roots of the cubic equation with no particular order
        //  if discriminant > 0: return three distinct real roots
        //  if discriminant < 0: return one real root and one complex root
        //                       (one member of the complex conjugate pair)
        //  if discriminant = 0: return two identical roots,
        //                       and one distinct real root
        //  if identical zero Hessian: return three identical real roots
        //  where the discriminant = - 4p^3 - 27q^2.
        Roots<3> roots() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cubicEqnI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
