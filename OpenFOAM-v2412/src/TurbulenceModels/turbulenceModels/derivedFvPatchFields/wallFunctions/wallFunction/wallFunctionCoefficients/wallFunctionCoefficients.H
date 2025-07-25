/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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
    Foam::wallFunctionCoefficients

Description
    Class to host the wall-function coefficients being
    used in the wall function boundary conditions.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Top-level entries
        ...

        // Optional entries
        Cmu             0.09;
        kappa           0.41;
        E               9.8;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                    | Type   | Reqd  | Deflt
      Cmu       | Empirical model coefficient    | scalar | no    | 0.09
      kappa     | von Karman constant            | scalar | no    | 0.41
      E         | Wall roughness parameter       | scalar | no    | 9.8
    \endtable

SourceFiles
    wallFunctionCoefficients.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_wallFunctionCoefficients_H
#define Foam_wallFunctionCoefficients_H

#include "scalarFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class dictionary;
class Ostream;

/*---------------------------------------------------------------------------*\
                  Class wallFunctionCoefficients Declaration
\*---------------------------------------------------------------------------*/

class wallFunctionCoefficients
{
    // Private Data

        //- Empirical model coefficient
        scalar Cmu_;

        //- von Karman constant
        scalar kappa_;

        //- Wall roughness parameter
        scalar E_;

        //- Estimated y+ value at the intersection
        //- of the viscous and inertial sublayers
        scalar yPlusLam_;


public:

    // Constructors

        //- Construct with default coefficients
        wallFunctionCoefficients();

        //- Construct from dictionary
        explicit wallFunctionCoefficients(const dictionary& dict);


    // Member Functions

    // Access

        //- Return the object: Cmu
        scalar Cmu() const noexcept
        {
            return Cmu_;
        }

        //- Return the object: kappa
        scalar kappa() const noexcept
        {
            return kappa_;
        }

        //- Return the object: E
        scalar E() const noexcept
        {
            return E_;
        }

        //- Return the object: yPlusLam
        scalar yPlusLam() const noexcept
        {
            return yPlusLam_;
        }


    // I-O

        //- Write wall-function coefficients as dictionary entries
        void writeEntries(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
