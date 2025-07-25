/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "ensightOutputSurface.H"
#include "ensightOutput.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::ensightOutputSurface::writeData
(
    ensightFile& os,
    const Field<Type>& fld,
    const bool isPointData
) const
{
    if (isPointData)
    {
        this->writePointData(os, fld);
    }
    else
    {
        this->writeFaceData(os, fld);
    }
}


template<class Type>
void Foam::ensightOutputSurface::writeFaceData
(
    ensightFile& os,
    const Field<Type>& fld
) const
{
    ensightOutput::writeField
    (
        os,
        fld,
        *this,
        false  /* serial only! */
    );
}


template<class Type>
void Foam::ensightOutputSurface::writePointData
(
    ensightFile& os,
    const Field<Type>& fld
) const
{
    const ensightOutputSurface& part = *this;

    // No geometry or field
    if (part.empty() || fld.empty())
    {
        return;
    }


    os.beginPart(part.index());

    ensightOutput::Detail::writeFieldComponents
    (
        os,
        ensightFile::coordinates,
        fld,
        false  /* serial only! */
    );
}


// ************************************************************************* //
