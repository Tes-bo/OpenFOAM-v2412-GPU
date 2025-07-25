/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022-2024 OpenCFD Ltd.
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
    Foam::readFieldsHandler

Description
    A simple field-loader, as per the readFields function object

\*---------------------------------------------------------------------------*/

#ifndef readFieldsHander_H
#define readFieldsHander_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "messageStream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class readFieldsHandler Declaration
\*---------------------------------------------------------------------------*/

class readFieldsHandler
{
    // Private Data

        //- Mesh reference
        fvMesh& mesh_;

        //- Output logging (verbosity)
        bool log;


    // Private Member Functions

        //- Attempt load from io, store on database if successful
        template<class FieldType>
        bool loadAndStore(const IOobject& io)
        {
            if (io.isHeaderClass<FieldType>())
            {
                // Store field on mesh database
                Log << "    Reading " << io.name()
                    << " (" << FieldType::typeName << ')' << endl;

                regIOobject::store(new FieldType(io, mesh_));
                return true;
            }

            return false;
        }

        //- Forward to loadAndStore for supported types
        template<class Type>
        bool loadField(const IOobject& io)
        {
            typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;
            typedef typename VolFieldType::Internal IntVolFieldType;
            typedef GeometricField<Type, fvsPatchField, surfaceMesh>
                SurfaceFieldType;

            return
            (
                loadAndStore<VolFieldType>(io)
             || loadAndStore<IntVolFieldType>(io)
             || loadAndStore<SurfaceFieldType>(io)
            );
        }


        //- Load all fields
        label loadFields(const UList<word>& fieldSet_)
        {
            label nLoaded = 0;

            for (const word& fieldName : fieldSet_)
            {
                // Already loaded?
                const auto* ptr = mesh_.cfindObject<regIOobject>(fieldName);

                if (ptr)
                {
                    ++nLoaded;
                    DebugInfo
                        << "readFields : "
                        << ptr->name() << " (" << ptr->type()
                            << ") already in database" << endl;
                    continue;
                }

                // Load field as necessary
                IOobject io
                (
                    fieldName,
                    mesh_.time().timeName(),
                    mesh_.thisDb(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                );

                const bool ok =
                (
                    io.typeHeaderOk<regIOobject>(false)
                 &&
                    (
                        loadField<scalar>(io)
                     || loadField<vector>(io)
                     || loadField<sphericalTensor>(io)
                     || loadField<symmTensor>(io)
                     || loadField<tensor>(io)
                    )
                );

                if (ok)
                {
                    ++nLoaded;
                }
                else
                {
                    DebugInfo
                        << "readFields : failed to load " << fieldName
                        << endl;
                }
        }

        return nLoaded;
    }


public:

    static const bool debug = false;


    // Constructors

        //- Construct
        explicit readFieldsHandler(fvMesh& mesh, bool verbose=true)
        :
            mesh_(mesh),
            log(verbose)
        {}


    // Member Functions

        bool execute(const UList<word>& fieldNames)
        {
            loadFields(fieldNames);
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
