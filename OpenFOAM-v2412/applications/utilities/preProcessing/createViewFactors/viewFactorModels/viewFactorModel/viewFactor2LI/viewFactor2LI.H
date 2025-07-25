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
    Foam::VF::viewFactor2LI

Description
    Computes view factors according to the double line integral (2LI) method.

Usage
    Minimal example in \c <constant>/viewFactorsDict:
    \verbatim
    // Optional entries
    alpha    <scalar>;

    // Inherited entries
    ...
    \endverbatim

    where the entries mean:
    \table
      Property    | Description                       | Type | Reqd | Deflt
      alpha       | Perturbation for common edges     | scalar | no | 0.21
    \endtable

    The inherited entries are elaborated in:
    - \link viewFactorModel.H \endlink

SourceFiles
    viewFactorModel.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_vf_viewFactor2LI_H
#define Foam_vf_viewFactor2LI_H

#include "viewFactorModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace VF
{

/*---------------------------------------------------------------------------*\
                        Class viewFactor2LI Declaration
\*---------------------------------------------------------------------------*/

class viewFactor2LI
:
    public viewFactorModel
{
protected:

    // Protected Data

        //- Perturbation for common edges; default = 0.21
        const scalar alpha_;


    // Protected Member Functions

        //- Calculate view factor using the double-area integral
        static scalar calculateFij
        (
            const List<point>& lPoints,
            const List<point>& rPoints,
            const scalar alpha
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
    TypeName("viewFactor2LI");

    //- Constructor
    viewFactor2LI(const fvMesh& mesh, const dictionary& dict);

    //- Destructor
    virtual ~viewFactor2LI() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace VF
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
