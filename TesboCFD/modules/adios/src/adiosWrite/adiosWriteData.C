/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2023 OpenCFD Ltd.
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

#include "adiosWrite.H"
#include "FlatOutput.H"
#include "cloud.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::functionObjects::adiosWrite::writeData()
{
    // Remove all old ADIOS variables/attributes since new variables may
    // have appeared (eg, via function objects), fields sizes changed ...
    // -> update everything for simplicity/consistency
    adiosCoreWrite::reset();

    // Output file for this time step
    if (!open(time()))
    {
        // Report problem ???
        // WarningInFunction<< Cannot ... file at ...
        return;
    }

    adiosCoreWrite::beginWrite();       // Begin step

    if (UPstream::master())
    {
        adiosCoreWrite::putBaseAttributes();
        adiosCoreWrite::putTimeAttributes(time());
    }

    DynamicList<word>   regNames(regions_.size());
    DynamicList<double> faceTimes(regions_.size());
    DynamicList<double> pointTimes(regions_.size());

    for (const adiosFoam::regionControl& regCtrl : regions_)
    {
        const fvMesh& mesh = time().lookupObject<fvMesh>(regCtrl.name());

        putPatchAttributes(mesh);

        const enum polyMesh::readUpdateState state = regCtrl.updateTimes(mesh);
        if (state != polyMesh::UNCHANGED)
        {
            // Write mesh if dynamic or first time

            Info<< "adiosWrite: write mesh " << mesh.name()
                << " at time " << mesh.time().timeName()
                << (state == polyMesh::POINTS_MOVED ? " (points only)" : "")
                << nl;

            writeMeshPoints(mesh);

            if (state == polyMesh::TOPO_CHANGE)
            {
                writeMeshFaces(mesh);
            }
        }

        writeFields(regCtrl, true);   // verbose
        writeClouds(regCtrl);

        // Information for region attributes
        regNames.append(regCtrl.name());
        faceTimes.append(regCtrl.topoTime());
        pointTimes.append(regCtrl.pointTime());
    }


    if (UPstream::master())
    {
        // Region attributes
        putIntAttribute(adiosFoam::foamAttribute / "nRegions", regions_.size());
        putListAttribute(adiosFoam::foamAttribute / "regions", regNames);
        putListAttribute(adiosFoam::timeAttribute / "faces", faceTimes);
        putListAttribute(adiosFoam::timeAttribute / "points", pointTimes);
    }

    adiosCoreWrite::endWrite();         // End step
    adiosCoreWrite::close();            // close/flush at every timestep
}


Foam::label Foam::functionObjects::adiosWrite::writeFields
(
    const adiosFoam::regionControl& regCtrl,
    bool verbose
)
{
    if (verbose)
    {
        Info<< "  " << regCtrl.info() << endl;
    }

    const fvMesh& mesh = time().lookupObject<fvMesh>(regCtrl.name());

    HashTable<word> wroteFields;
    HashTable<word> badFields;

    const wordList sortedNames = mesh.sortedNames();
    for (const word& objName : sortedNames)
    {
        const regIOobject& obj = *(mesh.find(objName).val());

        const word& name = obj.name();
        const word& type = obj.type();

        if (regCtrl.acceptFieldObject(obj))
        {
            if (writeFieldObject(obj))
            {
                wroteFields.insert(name, type);
            }
            else
            {
                badFields.set(name, type);
            }
        }
    }

    if (verbose)
    {
        for (const word& name : wroteFields.sortedToc())
        {
            Info<< "    " << name << "  (" << wroteFields[name] << ")" << nl;
        }
    }

    if (!badFields.empty())
    {
        const wordList names = badFields.sortedToc();
        wordList types(names.size());

        label nBad = 0;
        for (const word& name : names)
        {
            types[nBad++] = badFields[name];
        }

        WarningInFunction
            << nl
            << badFields.size() << " fields not handled by adiosWrite" << nl
            << "  names: " << flatOutput(names) << nl
            << "  types: " << flatOutput(types) << nl << endl;
    }


    return wroteFields.size();
}


Foam::label Foam::functionObjects::adiosWrite::writeClouds
(
    const adiosFoam::regionControl& regCtrl
)
{
    const fvMesh& mesh = time().lookupObject<fvMesh>(regCtrl.name());
    auto allClouds = mesh.lookupClass<cloud>();

    HashTable<const cloud*> cloudObjects(mesh.lookupClass<cloud>());

    if (!allClouds.empty())
    {
        Info<< type() << ": writeClouds: " << regCtrl.info() << endl;
    }

    wordHashSet pruned;

    forAllConstIters(cloudObjects, iter)
    {
        const cloud* cldPtr = iter.val();
        if (!regCtrl.acceptCloudObject(*cldPtr))
        {
            pruned.insert(cldPtr->name());
        }
    }

    cloudObjects.erase(pruned);

    if (cloudObjects.empty())
    {
        return 0;
    }

    objectRegistry obr
    (
        IOobject
        (
            "adios::writeCloud::",
            mesh.time().constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    const wordList cloudNames(cloudObjects.sortedToc());

    for (const word& cloudName : cloudNames)
    {
        const cloud* cldPtr = *(cloudObjects.cfind(cloudName));

        writeCloudRegistry(cldPtr, obr);

        obr.clear();
    }

    if (UPstream::master())
    {
        const fileName varPath = regCtrl.regionPath();

        // Number of active clouds as region attribute
        putIntAttribute(varPath / "nClouds", cloudNames.size());
        putListAttribute(varPath / "clouds", cloudNames);
    }

    return cloudNames.size();
}


// ************************************************************************* //
