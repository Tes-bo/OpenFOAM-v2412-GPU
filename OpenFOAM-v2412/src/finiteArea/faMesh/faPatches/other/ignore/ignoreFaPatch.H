/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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
    Foam::ignoreFaPatch

Description
    A patch which will not exist in the faMesh. Typical example is a front and
    back plane of a 2-D geometry

SourceFiles
    ignoreFaPatch.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_ignoreFaPatch_H
#define Foam_ignoreFaPatch_H

#include "faPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class ignoreFaPatch Declaration
\*---------------------------------------------------------------------------*/

class ignoreFaPatch
:
    public faPatch
{
public:

    //- Runtime type information
    TypeName("ignore");


    // Constructors

        //- Minimal construct from components
        ignoreFaPatch
        (
            const word& name,
            const label index,
            const faBoundaryMesh& bm,
            const label nbrPolyPatchi = -1,
            const word& patchType = typeName
        );

        //- Construct from components
        ignoreFaPatch
        (
            const word& name,
            const labelUList& edgeLabels,
            const label index,
            const faBoundaryMesh& bm,
            const label nbrPolyPatchi = -1,
            const word& patchType = typeName
        );

        //- Construct from dictionary
        ignoreFaPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const faBoundaryMesh& bm,
            const word& patchType
        );

        //- Copy construct, resetting the boundary mesh
        ignoreFaPatch(const ignoreFaPatch& p, const faBoundaryMesh& bm);

        //- Copy construct, resetting boundary mesh and addressing
        ignoreFaPatch
        (
            const ignoreFaPatch& p,
            const faBoundaryMesh& bm,
            const label index,
            const labelUList& edgeLabels,
            const label nbrPolyPatchi = -1
        );


        //- Construct and return a clone, resetting the boundary mesh
        virtual autoPtr<faPatch> clone(const faBoundaryMesh& bm) const
        {
            return autoPtr<faPatch>(new ignoreFaPatch(*this, bm));
        }

        //- Construct and return a clone, resetting the edge list
        //  and boundary mesh
        virtual autoPtr<faPatch> clone
        (
            const faBoundaryMesh& bm,
            const labelUList& edgeLabels,
            const label index,
            const label nbrPolyPatchi = -1
        ) const
        {
            return autoPtr<faPatch>
            (
                new ignoreFaPatch(*this, bm, index, edgeLabels, nbrPolyPatchi)
            );
        }


    // Member Functions
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
