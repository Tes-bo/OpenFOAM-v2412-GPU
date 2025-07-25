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
    Foam::quadraticEqn

Description
    Container to encapsulate various operations for
    quadratic equation of the forms with real coefficients:

    \f[
        a*x^2 + b*x + c = 0
          x^2 + B*x + C = 0
    \f]

    The expressions for the roots of \c quadraticEqn are as follows:

    \f[
        x1 = - (b + sign(b) sqrt(b^2 - 4ac)/(2*a))

        x2 = c/(a*x1)
    \f]

    where \c (b^2 - 4ac) is evaluated by fused multiply-adds to avoid
    detrimental cancellation.

    Reference:
    \verbatim
        Cancellation-avoiding quadratic formula (tag:F):
            Ford, W. (2014).
            Numerical linear algebra with applications: Using MATLAB.
            London: Elsevier/Academic Press.
            DOI:10.1016/C2011-0-07533-6

        Kahan's algo. to compute 'b^2-a*c' using fused multiply-adds (tag:JML):
            Jeannerod, C. P., Louvet, N., & Muller, J. M. (2013).
            Further analysis of Kahan's algorithm for the accurate
            computation of 2× 2 determinants.
            Mathematics of Computation, 82(284), 2245-2264.
            DOI:10.1090/S0025-5718-2013-02679-8
    \endverbatim

See also
    Test-quadraticEqn.C

SourceFiles
    quadraticEqnI.H
    quadraticEqn.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_quadraticEqn_H
#define Foam_quadraticEqn_H

#include "Roots.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class quadraticEqn Declaration
\*---------------------------------------------------------------------------*/

class quadraticEqn
:
    public VectorSpace<quadraticEqn, scalar, 3>
{
public:

    //- Component labeling enumeration
    enum components { A, B, C };


    // Constructors

        //- Default construct
        quadraticEqn() = default;

        //- Construct initialized to zero
        inline quadraticEqn(const Foam::zero);

        //- Construct from components
        inline quadraticEqn(const scalar a, const scalar b, const scalar c);


    // Member Functions

    // Access

        scalar a() const noexcept { return this->v_[A]; }
        scalar b() const noexcept { return this->v_[B]; }
        scalar c() const noexcept { return this->v_[C]; }

        scalar& a() noexcept { return this->v_[A]; }
        scalar& b() noexcept { return this->v_[B]; }
        scalar& c() noexcept { return this->v_[C]; }


    // Evaluate

        //- Evaluate the quadratic equation at x
        inline scalar value(const scalar x) const;

        //- Evaluate the derivative of the quadratic equation at x
        inline scalar derivative(const scalar x) const;

        //- Estimate the error of evaluation of the quadratic equation at x
        inline scalar error(const scalar x) const;

        //- Return the roots of the quadratic equation with no particular order
        //  if discriminant > 0: return two distinct real roots
        //  if discriminant < 0: return one of the complex conjugate-pair roots
        //  otherwise          : return two identical real roots
        Roots<2> roots() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "quadraticEqnI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
