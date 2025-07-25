/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2024 OpenCFD Ltd.
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

InNamespace
    Foam::ensightOutput

Description
    A collection of functions for writing clouds as ensight file content.

SourceFiles
    ensightOutputCloud.C
    ensightOutputCloudTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_ensightOutputCloud_H
#define Foam_ensightOutputCloud_H

#include "ensightFile.H"
#include "IOField.H"
#include "DynamicList.H"
#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class fvMesh;
class globalIndex;

namespace ensightOutput
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


//- Write cloud positions
bool writeCloudPositions
(
    //! Output file (must be valid on master)
    ensightFile& os,

    //! The positions (measured data) to write.
    //! Also used as intermediate buffer (on master)
    DynamicList<floatVector>& positions,

    //! The global sizes of \p positions (must be valid on master)
    //! and consistent with \p positions dimensions
    const globalIndex& procAddr
);


//- Write cloud positions
bool writeCloudPositions
(
    //! Output file (must be valid on master)
    ensightFile& os,

    //! The positions (measured data) to write.
    //! Also used as intermediate buffer (on master)
    DynamicList<floatVector>& positions
);


//- Write cloud positions
bool writeCloudPositions
(
    //! Output file (must be valid on master)
    ensightFile& os,

    //! The associated volume mesh
    const fvMesh& mesh,

    //! The name of the cloud
    const word& cloudName,

    //! The processor-local existence
    bool exists
);


//- Write cloud field, returning true if the field is non-empty.
template<class Type>
bool writeCloudField
(
    //! Output file (must be valid on master)
    ensightFile& os,

    //! The cloud field
    const UList<Type>& field,

    //! The global sizes of \p field (must be valid on master)
    //! and consistent with \p field dimensions
    const globalIndex& procAddr
);


//- Write cloud field, returning true if the field is non-empty.
template<class Type>
bool writeCloudField
(
    //! Output file (must be valid on master)
    ensightFile& os,

    //! The cloud field
    const UList<Type>& field
);


//- Read cloud field from IOobject (if exists == true) and write,
//- always returning true.
template<class Type>
bool readWriteCloudField
(
    //! Output file (must be valid on master)
    ensightFile& os,

    //! The field object name/location for reading
    const IOobject& fieldObject,

    //! Exists anywhere? May still be missing on a local processor
    bool existsAny
);


/*---------------------------------------------------------------------------*\
                        Namespace ensightOutput::Detail
\*---------------------------------------------------------------------------*/

namespace Detail
{

//- Write cloud field data (serial) with rounding and newlines.
//  \return the current output count
template<class Type>
label writeCloudFieldContent
(
    //! Output file (must be valid on master)
    ensightFile& os,

    const UList<Type>& fld,

    //! The current output count
    label count = 0
);


} // End namespace Detail


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace ensightOutput
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ensightOutputCloudTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
