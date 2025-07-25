/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
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

Application
    Test-IOField

Description
    Test the processor-local reading of IOField (used in the lagrangian libs)

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IOField.H"
#include "IOList.H"
#include "primitiveFields.H"
#include "polyMesh.H"
#include "Time.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void doWrite(const IOobject& io, const label sz)
{
    const bool writeOnProc = (sz > 0);

    IOField<Type> fld(io, sz);
    forAll(fld, i)
    {
        fld[i] = i + 1000.25 + (0.25 * i);
    }
    Pout<< "writing:" << fld << endl;
    fld.write(writeOnProc);
}


template<>
void doWrite<bool>(const IOobject& io, const label sz)
{
    const bool writeOnProc = (sz > 0);

    IOField<bool> fld(io, sz);
    forAll(fld, i)
    {
        fld[i] = i % 2;
    }
    Pout<< "writing:" << fld << endl;
    fld.write(writeOnProc);
}


template<class Type>
void doRead(const IOobject& io, const label sz)
{
    const bool readOnProc = (sz > 0);
    Pout<< "    readOnProc:" << readOnProc << endl;
    IOField<Type> fld(io, readOnProc);
    Pout<< "    wanted:" << sz << " actually read:" << fld.size() << endl;

    if (fld.size() != sz)
    {
        FatalErrorInFunction<< "io:" << io.objectPath() << exit(FatalError);
    }
}


template<class Type>
void writeAndRead
(
    const IOobject& io,
    const label sz,
    const word& writeType,
    IOobjectOption::readOption rOpt,
    const word& readType
)
{
    Pout<< "** Writing:" << writeType
        << " Reading:" << readType << endl;

    // The write handler
    fileHandler(fileOperation::New(writeType, true));

    // Delete
    Pout<< "Deleting:" << fileHandler().filePath(io.objectPath()) << endl;
    fileHandler().rm(fileHandler().filePath(io.objectPath()));

    // Write
    Pout<< "Writing:" << fileHandler().objectPath(io, io.name()) << endl;
    doWrite<Type>(io, sz);

    // The read handler
    fileHandler(fileOperation::New(readType, true));

    // Read
    IOobject readIO(io);
    readIO.readOpt(rOpt);
    Pout<< "Reading:"
        << fileHandler().filePath(readIO.objectPath()) << endl;
    doRead<Type>(readIO, sz);

    Pout<< "** Done writing:" << writeType
        << " Reading:" << readType << nl << nl << endl;
}


template<class Type>
void readIfPresent
(
    IOobject& io,
    const label sz,
    const word& readType
)
{
    fileHandler(fileOperation::New(readType, true));

    // Read
    Pout<< "Reading:" << fileHandler().filePath(io.objectPath()) << endl;
    io.readOpt(IOobject::READ_IF_PRESENT);
    doRead<Type>(io, sz);
}



template<class Type>
void doTests(IOobject& io, const label sz)
{
    const wordList handlers
    (
        Foam::fileOperation::wordConstructorTablePtr_->sortedToc()
    );

    Info<< "Found handlers: " << flatOutput(handlers) << nl
        << "Running tests with " << pTraits<Type>::typeName << nl << nl;

    // for (const word& readHandler : handlers)
    // {
    //     readIfPresent<Type>(io, sz, readHandler);
    // }

    for (const word& writeHandler : handlers)
    {
        for (const word& readHandler : handlers)
        {
            writeAndRead<Type>
            (
                io,
                sz,
                writeHandler,
                IOobject::READ_IF_PRESENT,
                readHandler
            );
        }
    }
}


// Main program

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::addBoolOption("bool", "Use bool for tests");
    argList::addBoolOption("scalar", "Use scalar for tests");
    argList::addBoolOption("label", "Use label for tests (default)");
    argList::addBoolOption("ref", "Test writing by ref");

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createPolyMesh.H"

    label sz = 0;
    if (Pstream::myProcNo() % 2)
    {
        sz = 1;
    }

    if (!Pstream::parRun())
    {
        sz = 10;
        Info<< "Serial: using " << sz << nl;
    }

    {
        IOobject io
        (
            "points",
            mesh.time().constant(),
            polyMesh::meshSubDir,
            mesh
        );

        Info<< "points path: " << io.typeFilePath<labelIOList>() << nl;
        Info<< "points path: " << io.typeFilePath<void>() << nl;

        io.resetHeader("bad-points");
        Info<< "bad path: " << io.typeFilePath<void>() << nl;
    }

    IOobject io
    (
        "bla",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER
    );

    {
        dictionary headerDict;
        io.writeHeader(headerDict, "anything", IOstreamOption());
        Info<< "IOobjectHeader" << headerDict << nl;
    }

    label tested = 0;

    if (args.found("bool"))
    {
        doTests<bool>(io, sz);
        ++tested;
    }
    if (args.found("scalar"))
    {
        doTests<scalar>(io, sz);
        ++tested;
    }
    if (!tested || args.found("label"))
    {
        doTests<label>(io, sz);
    }

    if (args.found("ref"))
    {
        Info<< nl << "Testing writing referenced external data" << nl << endl;

        IOobject ioOutput
        (
            args.executable(),
            "constant",
            runTime,
            IOobject::NO_REGISTER  // implicit convert to IOobjectOption
        );

        labelList ints(identity(200));

        ioOutput.rename(args.executable() + "-labels");
        Info<< "write " << ioOutput.objectRelPath() << endl;
        {
            IOListRef<label>(ioOutput, ints).write();
        }

        ioOutput.rename(args.executable() + "-points");
        Info<< "write " << ioOutput.objectRelPath() << endl;
        {
            IOFieldRef<vector>(ioOutput, mesh.points()).write();
        }
    }


    Pout<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
