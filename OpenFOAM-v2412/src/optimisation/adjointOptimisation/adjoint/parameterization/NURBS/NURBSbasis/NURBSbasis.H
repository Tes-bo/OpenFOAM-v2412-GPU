/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Copyright (C) 2013-2019 FOSS GP
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::NURBSbasis

Description
    NURBSbasis function. Used to construct NURBS curves, surfaces and volumes

SourceFiles
    NURBSbasis.C
    NURBSbasisTemplates.C
    NURBSbasisI.H

\*---------------------------------------------------------------------------*/

#ifndef NURBSbasis_H
#define NURBSbasis_H

#include "dictionary.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class NURBSbasis Declaration
\*---------------------------------------------------------------------------*/

class NURBSbasis
{
    // Private data

        //- Number of control points
        label nCPs_;

        //- Basis order
        label basisDegree_;

        //- Knot values
        scalarField knots_;


    // Private Member Functions

        //- Compute knots
        void computeKnots();


public:

    //- Runtime type information
    TypeName("NURBbasis");

    // Constructors

        //- Construct from number of control points, knot vector and basis order
        NURBSbasis
        (
            const label nCPs,
            const label degree,
            const scalarField& knots
        );

        //- Construct from number of control points and basis order
        //  Assumes uniform knot vector
        NURBSbasis
        (
            const label nCPs,
            const label degree
        );

        //- Construct from dictionary
        //  Assumes uniform knot vector
        NURBSbasis
        (
            const dictionary& dict
        );

        //- Construct as copy
        NURBSbasis(const NURBSbasis&);


    //- Destructor
    virtual ~NURBSbasis() = default;


    // Member Functions

        //- Basis value
        scalar basisValue
        (
            const label iCP,
            const label degree,
            const scalar u
        ) const;

        //- Basis derivative w.r.t u
        scalar basisDerivativeU
        (
            const label iCP,
            const label degree,
            const scalar u
        ) const;

        //- Basis second derivative w.r.t u
        scalar basisDerivativeUU
        (
            const label iCP,
            const label degree,
            const scalar u
        ) const;

        //- Checks to see if given u is affected by given CP.
        bool checkRange
        (
            const scalar u,
            const label CPI,
            const label degree
        ) const;

        //- Adds the new knot's u value, increments the nCPs and returns 'k',
        //- the index of the newly inserted uValue.
        //  uBar should be pre-bounded between 0,1.
        label insertKnot
        (
            const scalar uBar
        );

        inline const scalarField& knots() const;

        inline const label& degree() const;

        inline const label& nCPs() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "NURBSbasisI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
