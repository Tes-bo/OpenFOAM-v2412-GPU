/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2020 OpenCFD Ltd.
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
#include "adiosReader.H"
#include "cloud.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::adiosFoam::adiosTime
Foam::functionObjects::adiosWrite::readData(const fileName& bpFile)
{
    Info<< " Read adios step from file: " << bpFile << endl;

    adiosFoam::adiosReader reader(*adiosPtr_, bpFile);

    if (!reader.isGood())
    {
        return adiosFoam::adiosTime();
    }

    if (adiosCore::debug & 2)
    {
        Info<< "Get adiosTime: " << bpFile << endl;
    }
    adiosFoam::adiosTime timeInfo(reader);


    // Direct lookup via time_.lookupClass<fvMesh>() would be nice,
    // but have to trick the compiler not to get the const-version.
    // - instead, get the names and do the lookup ourselves

    const wordList meshNames(time().sortedNames<fvMesh>());

    DynamicList<adiosFoam::fieldInfo> regFields;

    // This is somewhat like demand-driven loading
    for (const word& regName : meshNames)
    {
        fvMesh& mesh = time().lookupObjectRef<fvMesh>(regName);

        if (adiosCore::debug & 2)
        {
            Info<<"region " << regName << " has reader:" << nl;
        }

        // Fields
        const HashSet<adiosFoam::fieldInfo>* fields =
            reader.getFieldInfo(regName);

        if (fields)
        {
            Info<< "  fields: " << flatOutput(fields->sortedToc()) << nl;

            for (const adiosFoam::fieldInfo& fldInfo : *fields)
            {
                auto iter =
                    static_cast<objectRegistry&>(mesh).find(fldInfo.name());

                if (iter.found())
                {
                    reader.readFieldObject
                    (
                        *(iter.val()),
                        fldInfo
                    );
                }
            }
        }


        // Clouds
        const HashSet<adiosFoam::cloudInfo>* clouds =
            reader.getCloudInfo(regName);

        if (clouds)
        {
            Info<< "  clouds: " << flatOutput(clouds->sortedToc()) << nl;

            for (const adiosFoam::cloudInfo& cldInfo : *clouds)
            {
                auto* cloudPtr = mesh.getObjectPtr<cloud>(cldInfo.name());

                if (!cloudPtr)
                {
                    continue;
                }

                if (adiosCore::debug & 2)
                {
                    Info<< "Cloud " << cldInfo
                        << " type: " << cloudPtr->type() << nl;

                    Pout<< "Cloud " << cldInfo
                        << " size: " << cloudPtr->nParcels() << nl;
                }


                auto obrPtr = reader.getCloud(cldInfo);
                if (!obrPtr)
                {
                    continue;
                }
                const objectRegistry& obr = *obrPtr;

                if (adiosCore::debug & 2)
                {
                    Info<<"Cloud registry: " << obr.name()
                        << " has " << obr.classes() << nl;

                    forAllConstIters(obr, iter)
                    {
                        Info<< iter.key() << " " << (*iter)->type() << nl;
                    }
                }

                const label oldSize = cloudPtr->nParcels();

                cloudPtr->readObjects(*obrPtr);

                if (adiosCore::debug & 2)
                {
                    Pout<< "Old size: " << oldSize
                        << " new: " << cloudPtr->nParcels() << nl;
                }
            }
        }
    }


    reader.close();

    return timeInfo;
}


Foam::adiosFoam::adiosTime
Foam::functionObjects::adiosWrite::readData(const instant& when)
{
    return readData(dataDir_/when.name());
}


Foam::adiosFoam::adiosTime
Foam::functionObjects::adiosWrite::readData()
{
    return readData(dataDir_/time().timeName());
}


// ************************************************************************* //
