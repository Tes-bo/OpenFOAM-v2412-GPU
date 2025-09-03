/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2018-2023 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "parLagrangianDistributor.H"
#include "Time.H"
#include "IOobjectList.H"
#include "mapDistributePolyMesh.H"
#include "cloud.H"
#include "CompactIOField.H"
#include "DynamicList.H"
#include "passivePositionParticleCloud.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Container>
Foam::wordList Foam::parLagrangianDistributor::filterObjects
(
    const IOobjectList& objects,
    const wordRes& selectedFields
)
{
    wordList fieldNames =
    (
        selectedFields.empty()
      ? objects.names<Container>()
      : objects.names<Container>(selectedFields)
    );

    // Parallel synchronise - combine names from all processors
    Pstream::combineReduce(fieldNames, ListOps::uniqueEqOp<word>());
    Foam::sort(fieldNames);  // Consistent order

    return fieldNames;
}


template<class Type>
Foam::label Foam::parLagrangianDistributor::distributeFields
(
    const mapDistributeBase& map,
    const word& cloudName,
    const bool haveCloud,
    const IOobjectList& objects,
    const wordRes& selectedFields
) const
{
    typedef IOField<Type> Container;

    const wordList fieldNames
    (
        filterObjects<IOField<Type>>
        (
            objects,
            selectedFields
        )
    );

    // Read if present
    IOobject srcIOobject
    (
        "none",
        srcMesh_.time().timeName(),
        cloud::prefix/cloudName,
        srcMesh_,
        IOobjectOption::LAZY_READ,
        IOobjectOption::NO_WRITE,
        IOobjectOption::NO_REGISTER
    );

    IOobject tgtIOobject
    (
        "none",
        tgtMesh_.time().timeName(),
        cloud::prefix/cloudName,
        tgtMesh_,
        IOobjectOption::NO_READ,
        IOobjectOption::NO_WRITE,
        IOobjectOption::NO_REGISTER
    );

    //bool reconstruct = false;

    if (verbose_ && fieldNames.size())
    {
        Info<< "    Distributing lagrangian "
            << Container::typeName << "s\n" << nl;
    }

    for (const word& objectName : fieldNames)
    {
        if (verbose_)
        {
            Info<< "        " <<  objectName << nl;
        }

        srcIOobject.resetHeader(objectName);
        tgtIOobject.resetHeader(objectName);

        // Read if present (ie, haveCloud means readOnProc)
        Container field(srcIOobject, haveCloud);

        // Distribute
        map.distribute(field);

        const bool writeOnProc = field.size();

        // Write
        Container(tgtIOobject, std::move(field)).write(writeOnProc);

        //if (!writeOnProc && !reconstruct)
        //{
        //    // When running with -overwrite it should also delete the old
        //    // files. Below works but is not optimal.
        //
        //    Foam::rm(tgtIOobject.objectPath());
        //}
    }

    if (verbose_ && fieldNames.size()) Info<< endl;
    return fieldNames.size();
}


template<class Type>
Foam::label Foam::parLagrangianDistributor::distributeFieldFields
(
    const mapDistributeBase& map,
    const word& cloudName,
    const bool haveCloud,
    const IOobjectList& objects,
    const wordRes& selectedFields
) const
{
    typedef CompactIOField<Field<Type>, Type> Container;

    DynamicList<word> fieldNames;

    // CompactIOField Field names
    fieldNames.push_back
    (
        filterObjects<CompactIOField<Field<Type>, Type>>
        (
            objects,
            selectedFields
        )
    );

    // IOField Field names
    fieldNames.push_back
    (
        filterObjects<IOField<Field<Type>>>
        (
            objects,
            selectedFields
        )
    );

    // Read if present
    IOobject srcIOobject
    (
        "none",
        srcMesh_.time().timeName(),
        cloud::prefix/cloudName,
        srcMesh_,
        IOobjectOption::LAZY_READ,
        IOobjectOption::NO_WRITE,
        IOobjectOption::NO_REGISTER
    );

    IOobject tgtIOobject
    (
        "none",
        tgtMesh_.time().timeName(),
        cloud::prefix/cloudName,
        tgtMesh_,
        IOobjectOption::NO_READ,
        IOobjectOption::NO_WRITE,
        IOobjectOption::NO_REGISTER
    );

    //bool reconstruct = false;

    if (verbose_ && fieldNames.size())
    {
        Info<< "    Distributing lagrangian "
            << Container::typeName << "s\n" << nl;
    }

    for (const word& objectName : fieldNames)
    {
        if (verbose_)
        {
            Info<< "        " <<  objectName << nl;
        }

        srcIOobject.resetHeader(objectName);
        tgtIOobject.resetHeader(objectName);

        // Read if present (ie, haveCloud means readOnProc)
        Container field(srcIOobject, haveCloud);

        // Distribute
        map.distribute(field);

        const bool writeOnProc = field.size();

        // Write
        Container(tgtIOobject, std::move(field)).write(writeOnProc);

        //if (!writeOnProc && !reconstruct)
        //{
        //    // When running with -overwrite it should also delete the old
        //    // files. Below works but is not optimal.
        //
        //    Foam::rm(tgtIOobject.objectPath());
        //}
    }

    if (verbose_ && fieldNames.size()) Info<< endl;
    return fieldNames.size();
}


template<class Container>
Foam::label Foam::parLagrangianDistributor::readFields
(
    const passivePositionParticleCloud& cloud,
    const bool haveCloud,
    const IOobjectList& objects,
    const wordRes& selectedFields
)
{
    const wordList fieldNames
    (
        filterObjects<Container>
        (
            objects,
            selectedFields
        )
    );

    // Read if present
    IOobject readIO
    (
        "none",
        cloud.time().timeName(),
        cloud,
        IOobjectOption::LAZY_READ,
        IOobjectOption::NO_WRITE,
        IOobjectOption::REGISTER
    );

    if (verbose_ && fieldNames.size())
    {
        Info<< "    Reading lagrangian "
            << Container::typeName << "s\n" << nl;
    }

    for (const word& objectName : fieldNames)
    {
        if (verbose_)
        {
            Info<< "        " <<  objectName << nl;
        }

        readIO.resetHeader(objectName);

        // Read if present (ie, haveCloud means readOnProc)
        Container* fieldPtr = new Container(readIO, haveCloud);

        fieldPtr->store();
    }

    if (verbose_ && fieldNames.size()) Info<< endl;
    return fieldNames.size();
}


template<class Container>
Foam::label Foam::parLagrangianDistributor::distributeStoredFields
(
    const mapDistributeBase& map,
    passivePositionParticleCloud& cloud
) const
{
    // Parallel-consistent
    UPtrList<Container> fields(cloud.sorted<Container>());

    IOobject writeIO
    (
        "none",
        tgtMesh_.time().timeName(),
        cloud::prefix/cloud.name(),
        tgtMesh_,
        IOobjectOption::NO_READ,
        IOobjectOption::NO_WRITE,
        IOobjectOption::NO_REGISTER
    );

    //bool reconstruct = false;

    if (verbose_ && fields.size())
    {
        Info<< "    Distributing lagrangian "
            << Container::typeName << "s\n" << nl;
    }

    for (Container& field : fields)
    {
        if (verbose_)
        {
            Info<< "        " <<  field.name() << nl;
        }

        map.distribute(field);

        writeIO.resetHeader(field.name());

        const bool writeOnProc = field.size();

        // Write
        Container(writeIO, std::move(field)).write(writeOnProc);

        //if (!writeOnProc && !reconstruct)
        //{
        //    // When running with -overwrite it should also delete the old
        //    // files. Below works but is not optimal.
        //
        //    Foam::rm(writeIO.objectPath());
        //}
    }

    if (verbose_ && fields.size()) Info<< endl;
    return fields.size();
}


// ************************************************************************* //
