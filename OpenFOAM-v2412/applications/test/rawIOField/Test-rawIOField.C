/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022-2023 OpenCFD Ltd.
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
    Test-rawIOField

Description
    Reading rawIOField from disk

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "Switch.H"
#include "primitiveFields.H"
#include "pointField.H"
#include "rawIOField.H"
#include "exprTraits.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

#undef  USE_ROOT_CASE
//#define USE_ROOT_CASE


template<class Type>
tmp<Field<Type>> readRawField
(
    const IOobject& io,
    IOobjectOption::readOption withAverage
)
{
    rawIOField<Type> raw(io, withAverage);

    Info<< "File: " << io.objectPath() << nl
        << "Read: " << raw.size()
        << ' ' << pTraits<Type>::typeName << " entries" << nl
        << "Average: " << Switch::name(raw.hasAverage())
        << " = " << raw.average() << endl;

    return tmp<Field<Type>>::New(std::move(static_cast<Field<Type>&>(raw)));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Test behaviour of rawIOField reading (writing?)"
    );

    argList::noCheckProcessorDirectories();

    argList::addBoolOption("scalar", "Read scalar field");
    argList::addBoolOption("vector", "Read vector field");
    argList::addBoolOption("point", "Read point field");
    argList::addBoolOption("average", "Require averaged value entry");
    argList::addBoolOption("try-average", "Optional averaged value entry");

    argList::addArgument("fileName");

    #ifdef USE_ROOT_CASE
        #include "setRootCase.H"
        #include "createTime.H"
    #else
    // Without root case, or time
    argList args(argc, argv);
    #endif

    fileName inputName = args.get<fileName>(1);

    IOobjectOption::readOption withAverage = IOobjectOption::NO_READ;

    if (args.found("average"))
    {
        withAverage = IOobjectOption::MUST_READ;
    }
    else if (args.found("try-average"))
    {
        withAverage = IOobjectOption::READ_IF_PRESENT;
    }

    Info<< "Using case: " << argList::envGlobalPath() << nl
        << "Read file:  " << inputName << nl
        << "with average: " << int(withAverage) << nl
        << endl;


    refPtr<Time> timePtr;

    #ifdef USE_ROOT_CASE
    timePtr.cref(runTime);
    #endif

    // Fallback (eg, no runTime)
    if (!timePtr.good())
    {
        timePtr.reset(Time::NewGlobalTime());
    }

    const auto& tm = *timePtr;

    fileName resolvedName(inputName);
    resolvedName.toAbsolute();

    IOobject io
    (
        resolvedName,   // absolute path
        tm,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER,
        true            // is global object (currently not used)
    );


    if (args.found("scalar"))
    {
        auto tfield = readRawField<scalar>(io, withAverage);
    }
    else if (args.found("point"))
    {
        auto tfield = readRawField<point>(io, withAverage);
    }
    else if (args.found("vector"))
    {
        auto tfield = readRawField<vector>(io, withAverage);
    }
    else
    {
        Info<< "no data type specified!\n";
    }

    Info<< nl << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
