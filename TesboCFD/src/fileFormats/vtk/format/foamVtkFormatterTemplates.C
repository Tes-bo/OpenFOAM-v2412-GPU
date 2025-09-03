/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2018 OpenCFD Ltd.
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

#include "foamVtkPTraits.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, Foam::direction nComp, int nTuple>
Foam::vtk::formatter& Foam::vtk::formatter::beginDataArray
(
    const word& dataName,
    uint64_t payLoad,
    bool leaveOpen
)
{
    openTag(vtk::fileTag::DATA_ARRAY);
    xmlAttr("type", vtkPTraits<Type>::typeName);
    xmlAttr("Name", dataName);

    if (nComp > 1)
    {
        xmlAttr(fileAttr::NUMBER_OF_COMPONENTS, int(nComp));
    }
    if (nTuple > 0)
    {
        xmlAttr(fileAttr::NUMBER_OF_TUPLES, nTuple);
    }
    xmlAttr("format", name());

    if (formatter::npos != payLoad)
    {
        uint64_t off = offset(payLoad);
        if (formatter::npos != off)
        {
            xmlAttr("offset", off);
        }
    }

    if (!leaveOpen)
    {
        closeTag();
    }

    return *this;
}


template<class Type, Foam::direction nComp, int nTuple>
Foam::vtk::formatter& Foam::vtk::formatter::PDataArray
(
    const word& dataName
)
{
    openTag("PDataArray");
    xmlAttr("type", vtkPTraits<Type>::typeName);

    if (dataName.size())
    {
        xmlAttr("Name", dataName);
    }
    if (nComp > 1)
    {
        xmlAttr(fileAttr::NUMBER_OF_COMPONENTS, int(nComp));
    }
    if (nTuple > 0)
    {
        xmlAttr(fileAttr::NUMBER_OF_TUPLES, nTuple);
    }

    closeTag(true);

    return *this;
}


// ************************************************************************* //
