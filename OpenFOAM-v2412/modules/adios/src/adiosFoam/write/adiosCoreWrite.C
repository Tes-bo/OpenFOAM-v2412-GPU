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

#include "adiosCoreWrite.H"
#include "adiosTime.H"
#include "cloud.H"
#include "foamVersion.H"
#include "endian.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

void Foam::adiosFoam::adiosCoreWrite::read(const dictionary& dict)
{
    // Default values prior to lookup in dictionary
    readMethod_  = "BP";
    writeMethod_ = "MPI";
    writeParams_ = "";

    dict.readIfPresent("debug", adiosCore::debug);

    dict.readIfPresent("readMethod",   readMethod_);
    dict.readIfPresent("writeMethod",  writeMethod_);
    dict.readIfPresent("writeOptions", writeParams_);

    Info<< "  ADIOS writeMethod: " << writeMethod_ << nl
        << "        writeParams: " << writeParams_ << nl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adiosFoam::adiosCoreWrite::adiosCoreWrite(const dictionary& dict)
:
    adiosCore(),
    adiosPtr_(nullptr),
    writeIOPtr_(nullptr),
    writeFilePtr_(nullptr),
    readMethod_(),
    writeMethod_(),
    writeParams_()
{
    Info<< "adiosWrite (" << UPstream::nProcs() << " procs)"  << endl;

    // Read dictionary
    read(dict);

    if (UPstream::parRun())
    {
        adiosPtr_.reset(new adios2::ADIOS(MPI_COMM_WORLD));
    }
    else
    {
        adiosPtr_.reset(new adios2::ADIOS());
    }

    // Use default Engine "BPFile"

    try
    {
        (void) adiosPtr_->DeclareIO("write");
    }
    catch (...)
    {
        // Need catch when in debug mode?
    }

    writeIOPtr_.reset(new adios2::IO(adiosPtr_->AtIO("write")));

    // Adjust parameters
    // Set the actual output method here for the ADIOS group
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::adiosFoam::adiosCoreWrite::~adiosCoreWrite()
{
    // adios_finalize(Pstream::myProcNo());

    // MPI_Barrier(comm_);
    // MPI_Comm_free(&comm_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::adiosFoam::adiosCoreWrite::open(const fileName& dataFile)
{
    DebugInFunction<< "Open write: " << dataFile << endl;

    const fileName path = dataFile.path();
    if (UPstream::master() && !isDir(path))
    {
        Foam::mkDir(path);
    }

    if (writeFilePtr_)
    {
        writeFilePtr_->Close();
        writeFilePtr_ = nullptr;
    }

    writeFilePtr_.reset
    (
        new adios2::Engine(writeIOPtr_->Open(dataFile, adios2::Mode::Write))
    );

    // Any error handling?

    return bool(writeFilePtr_);
}


void Foam::adiosFoam::adiosCoreWrite::close()
{
    // MPI_Barrier(MPI_COMM_WORLD);

    // Close the file
    if (writeFilePtr_)
    {
        DebugInFunction<< "Close write" << endl;

        writeFilePtr_->Close();
        writeFilePtr_ = nullptr;
    }
}


bool Foam::adiosFoam::adiosCoreWrite::beginWrite()
{
    return
    (
        writeFilePtr_
    &&  (adios2::StepStatus::OK == writeFilePtr_->BeginStep())
    );
}


bool Foam::adiosFoam::adiosCoreWrite::endWrite()
{
    if (writeFilePtr_)
    {
        writeFilePtr_->EndStep();
        return true;
    }

    return false;
}


void Foam::adiosFoam::adiosCoreWrite::reset()
{
    // Remove all variable definitions and attributes
    // to enable a new list of definitions in case the mesh changes

    if (writeIOPtr_)
    {
        writeIOPtr_->RemoveAllAttributes();
        writeIOPtr_->RemoveAllVariables();
    }
}


void Foam::adiosFoam::adiosCoreWrite::putBaseAttributes()
{
    putIntAttribute
    (
        adiosFoam::foamAttribute / "api",
        foamVersion::api
    );

    putAttribute
    (
        adiosFoam::foamAttribute / "patch",
        foamVersion::patch
    );

    // Endian type (LSB|MSB)
    putAttribute
    (
        adiosFoam::foamAttribute / "endian",
#ifdef WM_LITTLE_ENDIAN
        "LSB"
#else
        "MSB"
#endif
    );

    // OpenFOAM label size (32|64)
    putIntAttribute
    (
        adiosFoam::foamAttribute / "label",
        int(8 * sizeof(Foam::label))
    );

    // OpenFOAM scalar size (32=single, 64=double)
    putIntAttribute
    (
        adiosFoam::foamAttribute / "scalar",
        int(8 * sizeof(Foam::scalar))
    );

    // Other general information
    putIntAttribute
    (
        adiosFoam::foamAttribute / "nProcs",
        UPstream::nProcs()
    );
}


void Foam::adiosFoam::adiosCoreWrite::putTimeAttributes(const TimeState& t)
{
    putIntAttribute
    (
        adiosFoam::timeAttribute / adiosTime::attrNames[adiosTime::INDEX],
        t.timeIndex()
    );
    putDoubleAttribute
    (
        adiosFoam::timeAttribute / adiosTime::attrNames[adiosTime::VALUE],
        t.timeOutputValue()
    );
    putDoubleAttribute
    (
        adiosFoam::timeAttribute / adiosTime::attrNames[adiosTime::DT],
        t.deltaTValue()
    );
    putDoubleAttribute
    (
        adiosFoam::timeAttribute / adiosTime::attrNames[adiosTime::DT0],
        t.deltaT0Value()
    );
}


// ************************************************************************* //
