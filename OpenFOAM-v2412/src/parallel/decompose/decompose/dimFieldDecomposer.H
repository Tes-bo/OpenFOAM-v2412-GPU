/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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
    Foam::dimFieldDecomposer

Description
    Dimensioned field decomposer.

SourceFiles
    dimFieldDecomposer.C
    dimFieldDecomposerTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_dimFieldDecomposer_H
#define Foam_dimFieldDecomposer_H

#include "fvMesh.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class fvFieldDecomposer Declaration
\*---------------------------------------------------------------------------*/

class dimFieldDecomposer
{
    // Private Data

        //- Reference to processor mesh
        const fvMesh& procMesh_;

        //- Reference to face addressing
        //UNUSED: const labelList& faceAddressing_;

        //- Reference to cell addressing
        const labelList& cellAddressing_;


public:

    //- No copy construct
    dimFieldDecomposer(const dimFieldDecomposer&) = delete;

    //- No copy assignment
    void operator=(const dimFieldDecomposer&) = delete;


    // Constructors

        //- Construct from minimal components
        dimFieldDecomposer
        (
            const fvMesh& procMesh,
            const labelList& cellAddressing
        );

        //- Construct from components with API as per fvFieldDecomposer
        dimFieldDecomposer
        (
            const fvMesh& completeMesh,         //!< unused
            const fvMesh& procMesh,
            const labelList& faceAddressing,    //!< unused
            const labelList& cellAddressing
        );


    //- Destructor
    ~dimFieldDecomposer() = default;


    // Member Functions

        //- Decompose field
        template<class Type>
        tmp<DimensionedField<Type, volMesh>>
        decomposeField
        (
            const DimensionedField<Type, volMesh>& field
        ) const;


        //- Decompose list of fields
        template<class GeoField>
        void decomposeFields(const PtrList<GeoField>& fields) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "dimFieldDecomposerTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
