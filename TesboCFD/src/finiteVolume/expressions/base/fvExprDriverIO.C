/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2010-2018 Bernhard Gschaider
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "fvExprDriver.H"
#include "fvExprDriverWriter.H"
#include "cellSet.H"
#include "faceSet.H"
#include "pointSet.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::refPtr<Foam::labelList>
Foam::expressions::fvExprDriver::getTopoSetLabels
(
    const word& name,
    enum topoSetSource::sourceType setType
) const
{
    refPtr<labelList> selected;

    // Zones first
    // - cheap to handle (no IO) and can simply reference their labels

    switch (setType)
    {
        case topoSetSource::sourceType::CELLZONE_SOURCE:
        {
            const auto& zones = mesh().cellZones();
            const word& zoneTypeName = cellZone::typeName;

            const label zoneID = zones.findZoneID(name);
            if (zoneID < 0)
            {
                FatalErrorInFunction
                    << "No " << zoneTypeName << " named "
                    << name << "found. Has zones: " << zones.names() << endl
                    << exit(FatalError);
            }

            selected.cref(zones[zoneID]);
            break;
        }

        case topoSetSource::sourceType::FACEZONE_SOURCE:
        {
            const auto& zones = mesh().faceZones();
            const word& zoneTypeName = faceZone::typeName;

            const label zoneID = zones.findZoneID(name);
            if (zoneID < 0)
            {
                FatalErrorInFunction
                    << "No " << zoneTypeName << " named "
                    << name << "found. Has zones: " << zones.names() << endl
                    << exit(FatalError);
            }

            selected.cref(zones[zoneID]);
            break;
        }

        case topoSetSource::sourceType::POINTZONE_SOURCE:
        {
            const auto& zones = mesh().pointZones();
            const word& zoneTypeName = pointZone::typeName;

            const label zoneID = zones.findZoneID(name);
            if (zoneID < 0)
            {
                FatalErrorInFunction
                    << "No " << zoneTypeName << " named "
                    << name << "found. Has zones: " << zones.names() << endl
                    << exit(FatalError);
            }

            selected.cref(zones[zoneID]);
            break;
        }

        default:
            break;
    }


    if (selected.good())
    {
        return selected;
    }


    IOobject io(topoSet::findIOobject(mesh(), name));

    switch (setType)
    {
        case topoSetSource::sourceType::CELLSET_SOURCE:
        {
            typedef cellSet classType;

            if (!io.isHeaderClass<classType>())
            {
                FatalErrorInFunction
                    << "Error reading " << classType::typeName
                    << " <" << name << "> : found "
                    << io.headerClassName() << nl
                    << exit(FatalError);
            }

            classType set(io);
            selected.reset(refPtr<labelList>::New(set.sortedToc()));
            break;
        }

        case topoSetSource::sourceType::FACESET_SOURCE:
        {
            typedef faceSet classType;

            if (!io.isHeaderClass<classType>())
            {
                FatalErrorInFunction
                    << "Error reading " << classType::typeName
                    << " <" << name << "> : found "
                    << io.headerClassName() << nl
                    << exit(FatalError);
            }

            classType set(io);
            selected.reset(refPtr<labelList>::New(set.sortedToc()));
            break;
        }

        case topoSetSource::sourceType::POINTSET_SOURCE:
        {
            typedef pointSet classType;

            if (!io.isHeaderClass<classType>())
            {
                FatalErrorInFunction
                    << "Error reading " << classType::typeName
                    << " <" << name << "> : found "
                    << io.headerClassName() << nl
                    << exit(FatalError);
            }

            classType set(io);
            selected.reset(refPtr<labelList>::New(set.sortedToc()));
            break;
        }

        default:
        {
            FatalErrorInFunction
                << "Unexpected sourceType: " << int(setType) << nl
                << " for set <" << name << ">" << nl
                << exit(FatalError);
            break;
        }
    }

    return selected;
}


Foam::word Foam::expressions::fvExprDriver::getHeaderClassName
(
    const polyMesh& mesh,
    const word& name
)
{
    IOobject io
    (
        name,
        mesh.time().timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );
    io.typeHeaderOk<regIOobject>(false);

    DebugInfo
        << "Registry: " << mesh.path()
        << " Name: " << name
        << " Time: " << mesh.time().timeName()
        << " Path: " << io.localFilePath(io.headerClassName())
        << " Class: " << io.headerClassName() << endl;

    return io.headerClassName();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::Ostream& Foam::expressions::fvExprDriver::writeCommon
(
    Ostream& os,
    bool debug
) const
{
    // Write "variables", even if empty
    writeVariableStrings(os, "variables");

    if (debug)
    {
        os.writeEntry("variableValues", variables_);
    }

    if (!storedVariables_.empty() || !delayedVariables_.empty())
    {
        const_cast<fvExprDriver&>
        (
            *this
        ).updateSpecialVariables(true);
    }

    if (!storedVariables_.empty())
    {
        os.writeEntry("storedVariables", storedVariables_);
    }

    if (!delayedVariables_.empty())
    {
        List<exprResultDelayed> list(delayedVariables_.size());

        auto outIter = list.begin();

        forAllConstIters(delayedVariables_, iter)
        {
            *outIter = *iter;
            ++outIter;
        }

        os.writeEntry("delayedVariables", list);
    }

    if (!globalScopes_.empty())
    {
        os.writeEntry("globalScopes", globalScopes_);
    }

    // Write "functions<scalar>" ...
    writeFunctions(os);

    return os;
}


void Foam::expressions::fvExprDriver::createWriterAndRead(const word& name)
{
    if (!writer_ && hasDataToWrite())
    {
        writer_.reset(new fvExprDriverWriter(name + "_" + this->type(), *this));
    }
}


void Foam::expressions::fvExprDriver::tryWrite() const
{
    if (writer_ && mesh().time().writeTime())
    {
        writer_->write();
    }
}


// ************************************************************************* //
