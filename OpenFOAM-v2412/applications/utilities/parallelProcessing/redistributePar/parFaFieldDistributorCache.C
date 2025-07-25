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

\*---------------------------------------------------------------------------*/

#include "parFaFieldDistributorCache.H"

#include "areaFields.H"
#include "edgeFields.H"
#include "fieldsDistributor.H"
#include "faMeshDistributor.H"
#include "faMeshSubset.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class GeoField>
void Foam::parFaFieldDistributorCache::redistributeAndWrite
(
    const faMeshDistributor& distributor,
    PtrList<GeoField>& fields,
    const bool isWriteProc
)
{
    for (GeoField& fld : fields)
    {
        tmp<GeoField> tfld = distributor.distributeField(fld);

        if (isWriteProc)
        {
            tfld().write();
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::parFaFieldDistributorCache::read
(
    const Time& baseRunTime,
    const fileName& proc0CaseName,
    const bool decompose,  // i.e. read from undecomposed case

    const boolList& areaMeshOnProc,
    refPtr<fileOperation>& readHandler,
    const fileName& areaMeshInstance,
    faMesh& mesh
)
{
    Time& runTime = const_cast<Time&>(mesh.time());
    const bool oldProcCase = runTime.processorCase();

    autoPtr<faMeshSubset> subsetterPtr;

    // Missing an area mesh somewhere?
    if (areaMeshOnProc.found(false))
    {
        const bool oldParRun = UPstream::parRun(false);
        const int oldCache = fileOperation::cacheLevel(0);

        // A zero-sized mesh with boundaries.
        // This is used to create zero-sized fields.
        subsetterPtr.reset(new faMeshSubset(mesh, zero{}));

        fileOperation::cacheLevel(oldCache);
        UPstream::parRun(oldParRun);  // Restore parallel state

        // Deregister from polyMesh ...
        auto& obr = const_cast<objectRegistry&>
        (
            subsetterPtr->subMesh().thisDb()
        );

        obr.checkOut(faMesh::typeName);
        obr.checkOut("faBoundaryMesh");
        obr.checkOut("faSchemes");
        obr.checkOut("faSolution");
    }

    // Get original objects (before incrementing time!)
    if (UPstream::master() && decompose)
    {
        runTime.caseName() = baseRunTime.caseName();
        runTime.processorCase(false);
    }

    IOobjectList objects;

    if (readHandler)
    {
        auto oldHandler = fileOperation::fileHandler(readHandler);
        const auto oldComm = UPstream::commWorld(fileHandler().comm());

        objects = IOobjectList(mesh.thisDb(), runTime.timeName());
        readHandler = fileOperation::fileHandler(oldHandler);
        UPstream::commWorld(oldComm);
    }


    if (UPstream::master() && decompose)
    {
        runTime.caseName() = proc0CaseName;
        runTime.processorCase(oldProcCase);
    }

    Info<< "From time " << runTime.timeName()
        << " mesh:" << mesh.thisDb().objectRelPath()
        << " have objects:" << objects.names() << endl;

    if (UPstream::master() && decompose)
    {
        runTime.caseName() = baseRunTime.caseName();
        runTime.processorCase(false);
    }


    #undef  doFieldReading
    #define doFieldReading(Storage)                                   \
    fieldsDistributor::readFields                                     \
    (                                                                 \
        areaMeshOnProc, readHandler, mesh, subsetterPtr, objects,     \
        Storage,                                                      \
        true  /* (deregister field) */                                \
    );

    // areaFields
    doFieldReading(scalarAreaFields_);
    doFieldReading(vectorAreaFields_);
    doFieldReading(sphericalTensorAreaFields_);
    doFieldReading(symmTensorAreaFields_);
    doFieldReading(tensorAreaFields_);

    // edgeFields
    doFieldReading(scalarEdgeFields_);
    doFieldReading(vectorEdgeFields_);
    doFieldReading(tensorEdgeFields_);
    doFieldReading(sphericalTensorEdgeFields_);
    doFieldReading(symmTensorEdgeFields_);
    #undef doFieldReading
}


void Foam::parFaFieldDistributorCache::redistributeAndWrite
(
    const faMeshDistributor& distributor,
    const bool isWriteProc
)
{
    redistributeAndWrite(distributor, scalarAreaFields_, isWriteProc);
    redistributeAndWrite(distributor, vectorAreaFields_, isWriteProc);
    redistributeAndWrite(distributor, sphericalTensorAreaFields_, isWriteProc);
    redistributeAndWrite(distributor, symmTensorAreaFields_, isWriteProc);
    redistributeAndWrite(distributor, tensorAreaFields_, isWriteProc);

    redistributeAndWrite(distributor, scalarEdgeFields_, isWriteProc);
    redistributeAndWrite(distributor, vectorEdgeFields_, isWriteProc);
    redistributeAndWrite(distributor, sphericalTensorEdgeFields_, isWriteProc);
    redistributeAndWrite(distributor, symmTensorEdgeFields_, isWriteProc);
    redistributeAndWrite(distributor, tensorEdgeFields_, isWriteProc);
}


// ************************************************************************* //
