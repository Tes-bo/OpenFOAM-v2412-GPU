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

Namespace
    Foam::adiosFoam

Description
    Core routines and definitions for OpenFOAM ADIOS files.

SourceFiles
    adiosFoam.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_adiosFoam_H
#define Foam_adiosFoam_H

#include "fileName.H"
#include "instantList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace adiosFoam
{

/*---------------------------------------------------------------------------*\
                        Namespace adiosFoam Declaration
\*---------------------------------------------------------------------------*/

//- The default directory name for ADIOS files ("adiosData")
extern const word defaultDirectory;

//- OpenFOAM global attributes within ADIOS file ("/openfoam")
extern const string foamAttribute;

//- OpenFOAM time attributes within ADIOS file ("/time")
extern const string timeAttribute;


//- Search a given directory for ADIOS time directories,
//- but first checks for .bp files (ADIOS-2.4.0 and earlier)
instantList findTimes
(
    const fileName& directory,
    const word& constantName  = "constant"
);


//- Path name for a region
inline const word& regionPath
(
    const word& regionName
)
{
    return regionName;
}


//- Extract region name from a variable name (somewhat fragile)
inline word regionOf
(
    const fileName& varName
)
{
    const auto i = varName.find('/');

    return varName.substr(0, i);
}


//- Path name for fields
inline fileName fieldPath
(
    const word& regionName,
    const fileName& var = fileName::null
)
{
    return regionPath(regionName) / "field" / var;
}


//- Path name for meshes
inline fileName meshPath
(
    const word& regionName,
    const fileName& var = fileName::null
)
{
    return regionPath(regionName) / "polyMesh" / var;
}


//- Path name for clouds
inline fileName cloudPath
(
    const word& regionName,
    const word& cloudName = word::null
)
{
    return regionPath(regionName) / "cloud" / cloudName;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace adiosFoam
} // End namespace Foam

#endif

// ************************************************************************* //
