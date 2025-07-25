/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2022 OpenCFD Ltd.
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
    Foam::vtk

Description
    Namespace for handling VTK output.
    Contains classes and functions for writing VTK file content.

Namespace
    Foam::vtk::legacy

Description
    Namespace for legacy VTK output constants and functions.

SourceFiles
    foamVtkOutput.C
    foamVtkOutputTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_vtk_output_H
#define Foam_vtk_output_H

#include "autoPtr.H"
#include "bitSet.H"
#include "Enum.H"
#include "foamVtkCore.H"
#include "foamVtkFormatter.H"
#include "floatScalar.H"
#include "symmTensor.H"
#include "IOstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class instant;
class globalIndex;

namespace vtk
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// General Functions

    //- Return a default asciiFormatter
    autoPtr<vtk::formatter> newFormatter
    (
        std::ostream& os,
        unsigned prec = IOstream::defaultPrecision()
    );

    //- Return a new formatter based on the specified format type
    autoPtr<vtk::formatter> newFormatter
    (
        std::ostream& os,
        const enum formatType fmtType,
        unsigned prec = IOstream::defaultPrecision()
    );


    //- Write an identity list of labels.
    //  The output does not include the payload size.
    void writeIdentity(vtk::formatter& fmt, const label len, label start=0);

    //- Write a list of uint8_t values.
    //  The output does not include the payload size.
    void writeList(vtk::formatter& fmt, const UList<uint8_t>& values);

    //- Component-wise write of a value (N times)
    template<class Type>
    inline void write(vtk::formatter& fmt, const Type& val, const label n=1);

    //- Component-wise write of a value (N times) in parallel
    //  The value and count may differ on each processor
    template<class Type>
    inline void writeValueParallel
    (
        vtk::formatter& fmt,
        const Type& val,
        const label count=1
    );


    //- Write a list of values.
    //  The output does not include the payload size.
    template<class Type>
    void writeList(vtk::formatter& fmt, const UList<Type>& values);

    //- Write a list of values.
    //  The output does not include the payload size.
    template<class Type, unsigned N>
    void writeList(vtk::formatter& fmt, const FixedList<Type, N>& values);


    //- Write a list of values via indirect addressing.
    //  The output does not include the payload size.
    template<class Type>
    void writeList
    (
        vtk::formatter& fmt,
        const UList<Type>& values,
        const labelUList& addressing
    );

    //- Write a list of values via indirect addressing.
    //  The output does not include the payload size.
    template<class Type>
    void writeList
    (
        vtk::formatter& fmt,
        const UList<Type>& values,
        const bitSet& selected
    );

    //- Write a list of values and a list of values via indirect addressing.
    //  The output does not include the payload size.
    template<class Type>
    void writeLists
    (
        vtk::formatter& fmt,
        const UList<Type>& values1,
        const UList<Type>& values2,
        const labelUList& addressing
    );


    //- Write a list of values.
    //  The output does not include the payload size.
    template<class Type>
    void writeListParallel
    (
        vtk::formatter& fmt,
        const UList<Type>& values
    );

    //- Write a list of values, with constant per-processor offset
    //  The output does not include the payload size.
    void writeListParallel
    (
        vtk::formatter& fmt,
        const UList<label>& values,
        const globalIndex& procOffset
    );

    //- Write a list of values via indirect addressing.
    //  The output does not include the payload size.
    template<class Type>
    void writeListParallel
    (
        vtk::formatter& fmt,
        const UList<Type>& values,
        const labelUList& addressing
    );

    //- Write a list of values via indirect addressing.
    //  The output does not include the payload size.
    template<class Type>
    void writeListParallel
    (
        vtk::formatter& fmt,
        const UList<Type>& values,
        const bitSet& selected
    );

    //- Write a list of values and another list of values.
    //  The output does not include the payload size.
    template<class Type>
    void writeListsParallel
    (
        vtk::formatter& fmt,
        const UList<Type>& values1,
        const UList<Type>& values2
    );

    //- Write a list of values and a list of values via indirect addressing.
    //  The output does not include the payload size.
    template<class Type>
    void writeListsParallel
    (
        vtk::formatter& fmt,
        const UList<Type>& values1,
        const UList<Type>& values2,
        const labelUList& addressing
    );


/*---------------------------------------------------------------------------*\
                               Namespace legacy
\*---------------------------------------------------------------------------*/

namespace legacy
{

// Functions

    //- Emit header for legacy file (vtk DataFile Version 2.0)
    //  Writes "ASCII" or "BINARY" depending on specified type.
    void fileHeader(std::ostream& os, const std::string& title, bool binary);

    //- Emit header for legacy file, with "ASCII" or "BINARY" depending on
    //- the formatter type.
    //  If the contentType is non-empty, it is used for "DATASET" line.
    void fileHeader
    (
        vtk::formatter& fmt,
        const std::string& title,
        const std::string& contentType
    );

    //- Emit header for legacy file, with "ASCII" or "BINARY" depending on
    //- the formatter type.
    //  Includes "DATASET" with the specified dataset type.
    inline void fileHeader
    (
        vtk::formatter& fmt,
        const std::string& title,
        vtk::fileTag contentType
    );

    //- Emit header for legacy file, with "ASCII" or "BINARY" depending on
    //- the formatter type.
    //  Includes "DATASET" of the templated dataset type.
    template<vtk::fileTag ContentType>
    inline void fileHeader(vtk::formatter& fmt, const std::string& title);

    //- Emit header for POINTS (with trailing newline).
    inline void beginPoints(std::ostream& os, label nPoints);

    //- Emit header for VERTICES (with trailing newline).
    //  The nConnectivity is the sum of all connectivity points used,
    //  but \b without additional space for the size prefixes.
    //  The additional prefix sizes are added internally.
    //  \note With nConnectivity == 0, assume one point per element
    inline void beginVerts
    (
        std::ostream& os,
        label nVerts,
        label nConnectivity = 0
    );

    //- Emit header for LINES (with trailing newline).
    //  The nConnectivity is the sum of all connectivity points used,
    //  but \b without additional space for the size prefixes.
    //  The additional prefix sizes are added internally.
    //  \note With nConnectivity == 0, assume two points per element
    inline void beginLines
    (
        std::ostream& os,
        label nLines,
        label nConnectivity = 0
    );

    //- Emit header for POLYGONS (with trailing newline).
    //  The nConnectivity is the sum of all connectivity points used,
    //  but \b without additional space for the size prefixes.
    //  The additional prefix sizes are added internally.
    inline void beginPolys(std::ostream& os, label nPolys, label nConnectivity);


    //- Emit "FIELD FieldData <n>"
    inline void fieldData(vtk::formatter& fmt, label nFields);

    //- Emit legacy FIELD FieldData nFields.
    inline void beginFieldData(vtk::formatter& fmt, label nFields);

    //- Emit legacy CELL_DATA nCells, FIELD FieldData nFields.
    inline void beginCellData
    (
        vtk::formatter& fmt,
        label nCells,
        label nFields
    );

    //- Emit legacy POINT_DATA nPoints, FIELD FieldData nFields.
    inline void beginPointData
    (
        vtk::formatter& fmt,
        label nPoints,
        label nFields
    );


    //- Emit "TimeValue" for a FIELD entry (name as per Catalyst output)
    inline void writeTimeValue(vtk::formatter& fmt, scalar timeValue);

    //- Start output of float field with the specified name.
    template<direction nComp>
    inline void floatField
    (
        vtk::formatter& fmt,
        const word& name,
        const label nEntries
    );

    //- Start output of double field with the specified name.
    template<direction nComp>
    inline void doubleField
    (
        vtk::formatter& fmt,
        const word& name,
        const label nEntries
    );

    //- Start output of int field with the specified name.
    template<direction nComp>
    inline void intField
    (
        vtk::formatter& fmt,
        const word& name,
        const label nEntries
    );

} // End namespace legacy


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace vtk
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "foamVtkOutputI.H"

#ifdef NoRepository
    #include "foamVtkOutputTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
