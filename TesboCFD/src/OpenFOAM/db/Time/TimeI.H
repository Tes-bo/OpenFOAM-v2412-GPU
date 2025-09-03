/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2023 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::Time::Time
(
    const word& ctrlDictName,
    const argList& args,
    const bool enableFunctionObjects,
    const bool enableLibs,
    IOobjectOption::readOption rOpt  // (default: READ_MODIFIED)
)
:
    Foam::Time
    (
        ctrlDictName,
        args,
        "system",
        "constant",
        enableFunctionObjects,
        enableLibs,
        rOpt
    )
{}


inline Foam::Time::Time
(
    const word& ctrlDictName,
    const fileName& rootPath,
    const fileName& caseName,
    const bool enableFunctionObjects,
    const bool enableLibs,
    IOobjectOption::readOption rOpt  // (default: READ_MODIFIED)
)
:
    Foam::Time
    (
        ctrlDictName,
        rootPath,
        caseName,
        "system",
        "constant",
        enableFunctionObjects,
        enableLibs,
        rOpt
)
{}


inline Foam::Time::Time
(
    const dictionary& dict,
    const fileName& rootPath,
    const fileName& caseName,
    const bool enableFunctionObjects,
    const bool enableLibs,
    IOobjectOption::readOption rOpt  // (default: READ_MODIFIED)
)
:
    Foam::Time
    (
        dict,
        rootPath,
        caseName,
        "system",
        "constant",
        enableFunctionObjects,
        enableLibs,
        rOpt
    )
{}


inline Foam::Time::Time
(
    const fileName& rootPath,
    const fileName& caseName,
    const bool enableFunctionObjects,
    const bool enableLibs
)
:
    Foam::Time
    (
        rootPath,
        caseName,
        "system",
        "constant",
        enableFunctionObjects,
        enableLibs
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::IOstreamOption
Foam::Time::writeStreamOption() const noexcept
{
    return writeStreamOption_;
}


inline Foam::IOstreamOption::streamFormat
Foam::Time::writeFormat() const noexcept
{
    return writeStreamOption_.format();
}


inline Foam::IOstreamOption::streamFormat
Foam::Time::writeFormat(IOstreamOption::streamFormat fmt) noexcept
{
    auto old(writeStreamOption_.format(fmt));
    if (writeStreamOption_.format() != IOstreamOption::ASCII)
    {
        // Disable output compression for non-ASCII
        writeStreamOption_.compression(IOstreamOption::UNCOMPRESSED);
    }
    return old;
}


inline Foam::IOstreamOption::compressionType
Foam::Time::writeCompression() const noexcept
{
    return writeStreamOption_.compression();
}


inline Foam::IOstreamOption::versionNumber
Foam::Time::writeVersion() const noexcept
{
    return writeStreamOption_.version();
}


// ************************************************************************* //
