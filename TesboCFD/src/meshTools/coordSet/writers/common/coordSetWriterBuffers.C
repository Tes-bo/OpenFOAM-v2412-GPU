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

\*---------------------------------------------------------------------------*/

#include "coordSetWriter.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::coordSetWriter::writeLine
(
    Ostream& os,
    const UList<word>& values,
    const char* sep
)
{
    if (!values.empty())
    {
        forAll(values, coli)
        {
            if (coli && sep) os << sep;
            os  << values[coli];
        }
        os << nl;
    }
}


void Foam::coordSetWriter::writeLine
(
    Ostream& os,
    const UList<scalar>& values,
    const char* sep
)
{
    if (!values.empty())
    {
        forAll(values, coli)
        {
            if (coli && sep) os << sep;
            os  << values[coli];
        }
        os << nl;
    }
}


void Foam::coordSetWriter::clearBuffers()
{
    #undef doLocalCode
    #define doLocalCode(Type)                                                 \
    {                                                                         \
        Type##Names_.clear();                                                 \
        Type##Fields_.clear();                                                \
    }

    doLocalCode(label);
    doLocalCode(scalar);
    doLocalCode(vector);
    doLocalCode(sphericalTensor);
    doLocalCode(symmTensor);
    doLocalCode(tensor);
    #undef doLocalCode
}


Foam::label Foam::coordSetWriter::nDataColumns() const
{
    label ncol = 0;

    #undef doLocalCode
    #define doLocalCode(Type)                                                 \
        ncol += (Type##Fields_.size() * pTraits<Type>::nComponents);

    doLocalCode(label);
    doLocalCode(scalar);
    doLocalCode(vector);
    doLocalCode(sphericalTensor);
    doLocalCode(symmTensor);
    doLocalCode(tensor);
    #undef doLocalCode

    return ncol;
}


void Foam::coordSetWriter::getBufferLine
(
    DynamicList<scalar>& buf,
    const coordSet& coords,
    const label pointi
) const
{
    buf.clear();

    if (coords.hasVectorAxis())
    {
        const vector& p = coords.vectorCoord(pointi);
        buf.append(p.x());
        buf.append(p.y());
        buf.append(p.z());
    }
    else
    {
        buf.append(coords.scalarCoord(pointi));
    }

    do
    {
        #undef doLocalCode
        #define doLocalCode(Type)                                             \
                                                                              \
            for (const auto& fld : Type##Fields_)                             \
            {                                                                 \
                const auto& val = fld[pointi];                                \
                for (direction d=0; d < pTraits<Type>::nComponents; ++d)      \
                {                                                             \
                    buf.append(component(val, d));                            \
                }                                                             \
            }

        doLocalCode(label);
        doLocalCode(scalar);
        doLocalCode(vector);
        doLocalCode(sphericalTensor);
        doLocalCode(symmTensor);
        doLocalCode(tensor);
        #undef doLocalCode
    }
    while (false);
}


bool Foam::coordSetWriter::writeBuffered()
{
    return false;
}


void Foam::coordSetWriter::writeBufferContents
(
    Ostream& os,
    const coordSet& coords,
    const char* sep
) const
{
    const label npts = coords.size();
    const label ncomp = nDataColumns();

    DynamicList<scalar> compCols(3 + ncomp);

    for (label pointi = 0; pointi < npts; ++pointi)
    {
        getBufferLine(compCols, coords, pointi);
        writeLine(os, compCols, sep);
    }
}


// * * * * * * * * * * * * * * * * * Controls  * * * * * * * * * * * * * * * //

bool Foam::coordSetWriter::buffering() const
{
    return buffering_;
}


// disabled
bool Foam::coordSetWriter::buffering(const bool)
{
    return buffering_;
}


// ************************************************************************* //
