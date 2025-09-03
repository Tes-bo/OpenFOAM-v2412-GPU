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

Namespace
    Foam::ensightOutput

Description
    A collection of functions for writing ensight file content.

\*---------------------------------------------------------------------------*/

#ifndef Foam_ensightOutputFwd_H
#define Foam_ensightOutputFwd_H

#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class ensightFile;
class ensightGeoFile;


/*---------------------------------------------------------------------------*\
                        Namespace ensightOutput
\*---------------------------------------------------------------------------*/

namespace ensightOutput
{

//- The list type used for component-wise buffering.
//  Always has value_type = \c float (what ensight uses internally)
typedef DynamicList<float> floatBufferType;

//- Static debugging option
extern int debug;

//- Upper limit on number of items for bundled off-processor field transfers.
//- The component-wise transfer uses float (4 bytes).
//  Eg, 5M for 50 ranks of 100k cells each
extern int maxChunk_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace ensightOutput
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
