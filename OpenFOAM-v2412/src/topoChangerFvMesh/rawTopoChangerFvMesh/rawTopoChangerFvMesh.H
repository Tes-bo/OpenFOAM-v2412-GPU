/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::rawTopoChangerFvMesh

Description
    topoChangerFvMesh without any added functionality.

    Note: run without FOAM_SETNAN. Temporary has unitialised patch
    data when faces get created out of nothing.

SourceFiles
    rawTopoChangerFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_rawTopoChangerFvMesh_H
#define Foam_rawTopoChangerFvMesh_H

#include "topoChangerFvMesh.H"
#include "bitSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class rawTopoChangerFvMesh Declaration
\*---------------------------------------------------------------------------*/

class rawTopoChangerFvMesh
:
    public topoChangerFvMesh
{
    // Private Member Functions

        //- Set unmapped values
        template<class Type, template<class> class PatchField, class GeoMesh>
        static void setUnmappedValues
        (
            GeometricField<Type, PatchField, GeoMesh>& fld,
            const bitSet& mappedFace,
            const GeometricField<Type, PatchField, GeoMesh>& baseFld
        );

        template<class Type, template<class> class PatchField, class GeoMesh>
        void zeroUnmappedValues(const bitSet& mappedFace) const;

        //- No copy construct
        rawTopoChangerFvMesh(const rawTopoChangerFvMesh&) = delete;

        //- No copy assignment
        void operator=(const rawTopoChangerFvMesh&) = delete;

public:

    //- Runtime type information
    TypeName("rawTopoChangerFvMesh");


    // Constructors

        //- Construct from database
        explicit rawTopoChangerFvMesh
        (
            const IOobject& io,
            const bool doInit=true
        );


    //- Destructor
    virtual ~rawTopoChangerFvMesh() = default;


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "rawTopoChangerFvMeshTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
