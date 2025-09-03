/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2022 OpenCFD Ltd.
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

#include "meshReader.H"
#include "IOMap.H"
#include "OFstream.H"
#include "Time.H"

// * * * * * * * * * * * * * * * Static Functions  * * * * * * * * * * * * * //

void Foam::meshReader::warnDuplicates
(
    const word& context,
    const wordList& list
)
{
    HashTable<label> hashed(list.size());
    bool duplicates = false;

    for (const word& w : list)
    {
        // Check duplicate name
        auto iter = hashed.find(w);
        if (iter.good())
        {
            ++(*iter);
            duplicates = true;
        }
        else
        {
            hashed.insert(w, 1);
        }
    }

    // Warn about duplicate names
    if (duplicates)
    {
        Info<< nl << "WARNING: " << context << " with identical names:";
        forAllConstIters(hashed, iter)
        {
            if (*iter > 1)
            {
                Info<< "  " << iter.key();
            }
        }
        Info<< nl << endl;
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::meshReader::writeInterfaces(const objectRegistry& registry) const
{
    // write constant/polyMesh/interface
    IOList<labelList> ioObj
    (
        IOobject
        (
            "interfaces",
            registry.time().constant(),
            polyMesh::meshSubDir,
            registry,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );

    ioObj.note() = "as yet unsupported interfaces (baffles)";

    Info<< "Writing " << ioObj.name() << " to " << ioObj.objectPath() << endl;

    OFstream os(ioObj.objectPath());
    ioObj.writeHeader(os);

    os << interfaces_;
    IOobject::writeEndDivider(os);
}


void Foam::meshReader::writeMeshLabelList
(
    const objectRegistry& registry,
    const word& propertyName,
    const labelList& list,
    IOstreamOption streamOpt
) const
{
    // write constant/polyMesh/propertyName
    IOListRef<label> ioObj
    (
        IOobject
        (
            propertyName,
            registry.time().constant(),
            polyMesh::meshSubDir,
            registry,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        list
    );


    ioObj.note() = "persistent data for STARCD <-> OPENFOAM translation";
    Info<< "Writing " << ioObj.name() << " to " << ioObj.objectPath() << endl;

    // NOTE:
    // the cellTableId is an integer and almost always < 1000, thus ASCII
    // will be compacter than binary and makes external scripting easier

    ioObj.writeObject(streamOpt, true);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::meshReader::writeAux(const objectRegistry& registry) const
{
    cellTable_.writeDict(registry);
    writeInterfaces(registry);

    // write origCellId as List<label>
    writeMeshLabelList
    (
        registry,
        "origCellId",
        origCellId_,
        IOstreamOption(IOstreamOption::BINARY)
    );

    // write cellTableId as List<label>
    // this is crucial for later conversion back to ccm/starcd
    writeMeshLabelList
    (
        registry,
        "cellTableId",
        cellTableId_,
        IOstreamOption(IOstreamOption::ASCII)
    );
}


// ************************************************************************* //
