/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2022 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * Specializations * * * * * * * * * * * * * * * //

namespace Foam
{
namespace vtk
{

//- Template specialization for label
template<>
inline void write<label>(vtk::formatter& fmt, const label& val, const label n)
{
    for (label i=0; i < n; ++i)
    {
        fmt.write(val);
    }
}


//- Template specialization for float
template<>
inline void write<float>(vtk::formatter& fmt, const float& val, const label n)
{
    for (label i=0; i < n; ++i)
    {
        fmt.write(val);
    }
}


//- Template specialization for double
template<>
inline void write<double>(vtk::formatter& fmt, const double& val, const label n)
{
    for (label i=0; i < n; ++i)
    {
        fmt.write(val);
    }
}


//- Template specialization for symmTensor ordering
//  VTK order is (XX, YY, ZZ, XY, YZ, XZ)
template<>
inline void write(vtk::formatter& fmt, const symmTensor& val, const label n)
{
    for (label i=0; i < n; ++i)
    {
        fmt.write(component(val, symmTensor::XX));
        fmt.write(component(val, symmTensor::YY));
        fmt.write(component(val, symmTensor::ZZ));
        fmt.write(component(val, symmTensor::XY));
        fmt.write(component(val, symmTensor::YZ));
        fmt.write(component(val, symmTensor::XZ));
    }
}


} // End namespace vtk
} // End namespace Foam


// * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * * //

inline void Foam::vtk::legacy::fileHeader
(
    vtk::formatter& fmt,
    const std::string& title,
    vtk::fileTag contentType
)
{
    legacy::fileHeader(fmt, title, legacy::contentNames[contentType]);
}


template<Foam::vtk::fileTag ContentType>
inline void Foam::vtk::legacy::fileHeader
(
    vtk::formatter& fmt,
    const std::string& title
)
{
    legacy::fileHeader(fmt, title, legacy::contentNames[ContentType]);
}


inline void Foam::vtk::legacy::beginPoints(std::ostream& os, label nPoints)
{
    os  << nl
        << legacy::fileTagNames[vtk::fileTag::POINTS]
        << ' ' << nPoints
        << " float" << nl;
}


inline void Foam::vtk::legacy::beginVerts
(
    std::ostream& os,
    label nVerts,
    label nConnectivity
)
{
    if (!nConnectivity)
    {
        nConnectivity = nVerts;
    }
    os  << nl
        << legacy::fileTagNames[vtk::fileTag::VERTS]
        << ' ' << nVerts
        << ' ' << (nVerts + nConnectivity) << nl;
}


inline void Foam::vtk::legacy::beginLines
(
    std::ostream& os,
    label nLines,
    label nConnectivity
)
{
    if (!nConnectivity)
    {
        nConnectivity = 2*nLines;
    }
    os  << nl
        << legacy::fileTagNames[vtk::fileTag::LINES]
        << ' ' << nLines
        << ' ' << (nLines + nConnectivity) << nl;
}


inline void Foam::vtk::legacy::beginPolys
(
    std::ostream& os,
    label nPolys,
    label nConnectivity
)
{
    os  << nl
        << legacy::fileTagNames[vtk::fileTag::POLYS]
        << ' ' << nPolys
        << ' ' << (nPolys + nConnectivity) << nl;
}


inline void Foam::vtk::legacy::fieldData
(
    vtk::formatter& fmt,
    label nFields
)
{
    fmt.os()
        << "FIELD FieldData " << nFields << nl;
}


inline void Foam::vtk::legacy::beginFieldData
(
    vtk::formatter& fmt,
    label nFields
)
{
    legacy::fieldData(fmt, nFields);
}


inline void Foam::vtk::legacy::beginCellData
(
    vtk::formatter& fmt,
    label nCells,
    label nFields
)
{
    fmt.os()
        << nl
        << legacy::fileTagNames[vtk::fileTag::CELL_DATA]
        << ' ' << nCells << nl;
    legacy::fieldData(fmt, nFields);
}


inline void Foam::vtk::legacy::beginPointData
(
    vtk::formatter& fmt,
    label nPoints,
    label nFields
)
{
    fmt.os()
        << nl
        << legacy::fileTagNames[vtk::fileTag::POINT_DATA]
        << ' ' << nPoints << nl;
    legacy::fieldData(fmt, nFields);
}


inline void Foam::vtk::legacy::writeTimeValue
(
    vtk::formatter& fmt,
    scalar timeValue
)
{
    legacy::floatField<1>(fmt, "TimeValue", 1);
    fmt.write(timeValue);
    fmt.flush();
}


template<Foam::direction nComp>
inline void Foam::vtk::legacy::doubleField
(
    vtk::formatter& fmt,
    const word& fieldName,
    label nEntries
)
{
    fmt.os()
        << fieldName << ' '
        << int(nComp) << ' ' << nEntries << " double" << nl;
}


template<Foam::direction nComp>
inline void Foam::vtk::legacy::floatField
(
    vtk::formatter& fmt,
    const word& fieldName,
    label nEntries
)
{
    fmt.os()
        << fieldName << ' '
        << int(nComp) << ' ' << nEntries << " float" << nl;
}


template<Foam::direction nComp>
inline void Foam::vtk::legacy::intField
(
    vtk::formatter& fmt,
    const word& fieldName,
    label nEntries
)
{
    fmt.os()
        << fieldName << ' '
        << int(nComp) << ' ' << nEntries << " int" << nl;
}


// ************************************************************************* //
