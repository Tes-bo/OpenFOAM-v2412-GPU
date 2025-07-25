/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023-2024 OpenCFD Ltd.
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
    Foam::VF::viewFactorHottel

Description
    Computes view factors according to Hottel's crossed strings method.

    Reference:
    \verbatim
        Hottel, H. C., & Saforim, A. F. (1967).
        Radiative transfer.
        McGraw-Hill Book Company, New York.
    \endverbatim

Usage
    Minimal example in \c <constant>/viewFactorsDict:
    \verbatim
    // Inherited entries
    ...
    \endverbatim

    The inherited entries are elaborated in:
    - \link viewFactorModel.H \endlink

Note
    Only applicable to 2D cases

SourceFiles
    viewFactorModel.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_vf_viewFactorHottel_H
#define Foam_vf_viewFactorHottel_H

#include "viewFactorModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace VF
{

/*---------------------------------------------------------------------------*\
                       Class viewFactorHottel Declaration
\*---------------------------------------------------------------------------*/

class viewFactorHottel
:
    public viewFactorModel
{
    // Private Data

        //- Mesh empty direction
        vector emptyDir_;

        //- Mesh 2D width
        scalar w_;


protected:

    // Protected Member Functions

        //- Calculate view factor using the double-area integral
        static scalar calculateFij
        (
            const point& p0,
            const point& p1,
            const point& p2,
            const point& p3
        );

        //- Calculate
        virtual scalarListList calculate
        (
            const labelListList& visibleFaceFaces,
            const pointField& compactCf,
            const vectorField& compactSf,
            const List<List<vector>>& compactFineSf,
            const List<List<point>>& compactFineCf,
            const DynamicList<List<point>>& compactPoints,
            const DynamicList<label>& compactPatchId
        ) const;


public:

    //- Runtime type information
    TypeName("viewFactorHottel");

    //- Constructor
    viewFactorHottel(const fvMesh& mesh, const dictionary& dict);

    //- Destructor
    virtual ~viewFactorHottel() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace VF
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
