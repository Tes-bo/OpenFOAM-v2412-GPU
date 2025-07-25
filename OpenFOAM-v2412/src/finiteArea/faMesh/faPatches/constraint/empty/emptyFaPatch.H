/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2021-2022 OpenCFD Ltd.
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
    Foam::emptyFaPatch

Description
    A patch which will not exist in the faMesh. Typical example is a front and
    back plane of a 2-D geometry

Author
    Zeljko Tukovic, FMENA
    Hrvoje Jasak, Wikki Ltd.

SourceFiles
    emptyFaPatch.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_emptyFaPatch_H
#define Foam_emptyFaPatch_H

#include "faPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class emptyFaPatch Declaration
\*---------------------------------------------------------------------------*/

class emptyFaPatch
:
    public faPatch
{
public:

    //- Runtime type information
    TypeName("empty");


    // Constructors

        //- Minimal construct from components
        emptyFaPatch
        (
            const word& name,
            const label index,
            const faBoundaryMesh& bm,
            const label nbrPolyPatchi = -1,
            const word& patchType = typeName
        );

        //- Construct from components
        emptyFaPatch
        (
            const word& name,
            const labelUList& edgeLabels,
            const label index,
            const faBoundaryMesh& bm,
            const label nbrPolyPatchi = -1,
            const word& patchType = typeName
        );

        //- Construct from dictionary
        emptyFaPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const faBoundaryMesh& bm,
            const word& patchType
        );

        //- Copy construct, resetting the boundary mesh
        emptyFaPatch(const emptyFaPatch& p, const faBoundaryMesh& bm);

        //- Copy construct, resetting boundary mesh and addressing
        emptyFaPatch
        (
            const emptyFaPatch& p,
            const faBoundaryMesh& bm,
            const label index,
            const labelUList& edgeLabels,
            const label nbrPolyPatchi = -1
        );


        //- Construct and return a clone, resetting the boundary mesh
        virtual autoPtr<faPatch> clone(const faBoundaryMesh& bm) const
        {
            return autoPtr<faPatch>(new emptyFaPatch(*this, bm));
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
                new emptyFaPatch(*this, bm, index, edgeLabels, nbrPolyPatchi)
            );
        }


    // Member Functions

        virtual label size() const
        {
            return 0;
        }

        //- Return face normals. Over-riding base class return to get zero size
        //
//         virtual const vectorField& edgeNormals() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
