/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::interpolationCell

Description
    Uses the cell value for any location within the cell

\*---------------------------------------------------------------------------*/

#ifndef Foam_interpolationCell_H
#define Foam_interpolationCell_H

#include "interpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                      Class interpolationCell Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class interpolationCell
:
    public interpolation<Type>
{
public:

    //- Runtime type information
    TypeName("cell");


    // Constructors

        //- Construct from components
        explicit interpolationCell
        (
            const GeometricField<Type, fvPatchField, volMesh>& psi
        );


    // Member Functions

        //- Returns the cell value
        inline Type interpolate
        (
            const vector& position,             /* (unused) */
            const label celli,
            const label facei = -1              /* (unused) */
        ) const
        {
            return this->psi_[celli];
        }

        //- Returns the cell value corresponding to the given tetrahedron
        inline Type interpolate
        (
            const barycentric& coordinates,     /* (unused) */
            const tetIndices& tetIs,
            const label facei = -1              /* (unused) */
        ) const
        {
            return this->psi_[tetIs.cell()];
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "interpolationCell.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
