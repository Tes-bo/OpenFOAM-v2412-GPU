/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2023 OpenCFD Ltd.
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

InNamespace
    Foam

Description
    Gets the indices of (source)particles that have been appended to the
    target cloud and maps the lagrangian fields accordingly.

\*---------------------------------------------------------------------------*/

#ifndef Foam_MapLagrangianFields_H
#define Foam_MapLagrangianFields_H

#include "cloud.H"
#include "GeometricField.H"
#include "meshToMesh0.H"
#include "IOobjectList.H"
#include "CompactIOField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Gets the indices of (source)particles that have been appended to the
//  target cloud and maps the lagrangian fields accordingly.
template<class SourceIOFieldType, class TargetIOFieldType>
void MapLagrangianFields
(
    const string& cloudName,
    const IOobjectList& objects,
    const meshToMesh0& meshToMesh0Interp,
    const labelList& addParticles,
    const char* msg
)
{
    const fvMesh& meshTarget = meshToMesh0Interp.toMesh();

    for (const IOobject& io : objects.csorted<SourceIOFieldType>())
    {
        Info<< "    mapping lagrangian " << msg << ' ' << io.name() << endl;

        // Read field (does not need mesh)
        // Note: some fieldFields are 0 size (e.g. collision records)
        //   if not used, so catch any of those for Field as well

        SourceIOFieldType fieldSource(io);

        // Map
        TargetIOFieldType fieldTarget
        (
            IOobject
            (
                io.name(),
                meshTarget.time().timeName(),
                cloud::prefix/cloudName,
                meshTarget,
                IOobjectOption::NO_READ,
                IOobjectOption::NO_WRITE,
                IOobjectOption::NO_REGISTER
            ),
            min(fieldSource.size(), addParticles.size())  // handle 0 size
        );

        if (!fieldSource.empty())
        {
            forAll(addParticles, i)
            {
                fieldTarget[i] = fieldSource[addParticles[i]];
            }
        }

        fieldTarget.write();
    }
}


//- Gets the indices of (source)particles that have been appended to the
//  target cloud and maps the lagrangian fields accordingly.
template<class Type>
void MapLagrangianFields
(
    const string& cloudName,
    const IOobjectList& objects,
    const meshToMesh0& meshToMesh0Interp,
    const labelList& addParticles
)
{
    MapLagrangianFields
    <
        IOField<Type>,
        IOField<Type>
    >
    (
        cloudName,
        objects,
        meshToMesh0Interp,
        addParticles,
        "Field"
    );

    // Target is CompactIOField to automatically write in
    // compact form for binary format.
    MapLagrangianFields
    <
        IOField<Field<Type>>,
        CompactIOField<Field<Type>, Type>
    >
    (
        cloudName,
        objects,
        meshToMesh0Interp,
        addParticles,
        "FieldField"
    );

    MapLagrangianFields
    <
        CompactIOField<Field<Type>, Type>,
        CompactIOField<Field<Type>, Type>
    >
    (
        cloudName,
        objects,
        meshToMesh0Interp,
        addParticles,
        "FieldField"
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
