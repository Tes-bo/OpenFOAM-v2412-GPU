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
    Foam::parFaFieldDistributorCache

Description
    Simple container to manage read/write, redistribute finiteArea fields.

SourceFiles
    parFaFieldDistributorCache.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_parFaFieldDistributorCache_H
#define Foam_parFaFieldDistributorCache_H

#include "faMesh.H"
#include "faMeshDistributor.H"
#include "areaFieldsFwd.H"
#include "edgeFieldsFwd.H"
#include "PtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class parFaFieldDistributorCache Declaration
\*---------------------------------------------------------------------------*/

class parFaFieldDistributorCache
{
    // Private Data

    #undef  declareField
    #define declareField(Type)                                                \
    PtrList<GeometricField<Type, faPatchField, areaMesh>> Type##AreaFields_;  \
    PtrList<GeometricField<Type, faePatchField, edgeMesh>> Type##EdgeFields_;

    declareField(scalar);
    declareField(vector);
    declareField(sphericalTensor);
    declareField(symmTensor);
    declareField(tensor);
    #undef declareField


    // Private Member Functions

        //- Redistribute and write fields of given type (file-scope use)
        template<class GeoField>
        static void redistributeAndWrite
        (
            const faMeshDistributor& distributor,
            PtrList<GeoField>& fields,
            const bool isWriteProc
        );


public:

    // Constructors

        //- No copy construct
        parFaFieldDistributorCache(const parFaFieldDistributorCache&) = delete;

        //- No copy assignment
        void operator=(const parFaFieldDistributorCache&) = delete;

        //- Default construct
        parFaFieldDistributorCache() = default;


    // Member Functions

        //- Read distributed fields
        void read
        (
            const Time& baseRunTime,
            const fileName& proc0CaseName,
            const bool decompose,  // i.e. read from undecomposed case

            const boolList& areaMeshOnProc,
            refPtr<fileOperation>& readHandler,
            const fileName& areaMeshInstance,
            faMesh& mesh
        );

        void redistributeAndWrite
        (
            const faMeshDistributor& distributor,
            const bool isWriteProc
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
