/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "DSMCParcel.H"
#include "IOstreams.H"
#include "IOField.H"
#include "Cloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>
const std::size_t Foam::DSMCParcel<ParcelType>::sizeofFields
(
    sizeof(DSMCParcel<ParcelType>) - sizeof(ParcelType)
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::DSMCParcel<ParcelType>::DSMCParcel
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields,
    bool newFormat
)
:
    ParcelType(mesh, is, readFields, newFormat),
    U_(Zero),
    Ei_(0.0),
    typeId_(-1)
{
    if (readFields)
    {
        if (is.format() == IOstreamOption::ASCII)
        {
            is  >> U_ >> Ei_ >> typeId_;
        }
        else if (!is.checkLabelSize<>() || !is.checkScalarSize<>())
        {
            // Non-native label or scalar size

            is.beginRawRead();

            readRawScalar(is, U_.data(), vector::nComponents);
            readRawScalar(is, &Ei_);
            readRawLabel(is, &typeId_);

            is.endRawRead();
        }
        else
        {
            is.read(reinterpret_cast<char*>(&U_), sizeofFields);
        }
    }

    is.check(FUNCTION_NAME);
}


template<class ParcelType>
void Foam::DSMCParcel<ParcelType>::readFields(Cloud<DSMCParcel<ParcelType>>& c)
{
    const bool readOnProc = c.size();

    ParcelType::readFields(c);

    IOField<vector> U(c.newIOobject("U", IOobject::MUST_READ), readOnProc);
    c.checkFieldIOobject(c, U);

    IOField<scalar> Ei(c.newIOobject("Ei", IOobject::MUST_READ), readOnProc);
    c.checkFieldIOobject(c, Ei);

    IOField<label> typeId
    (
        c.newIOobject("typeId", IOobject::MUST_READ),
        readOnProc
    );
    c.checkFieldIOobject(c, typeId);

    label i = 0;
    for (DSMCParcel<ParcelType>& p : c)
    {
        p.U_ = U[i];
        p.Ei_ = Ei[i];
        p.typeId_ = typeId[i];
        ++i;
    }
}


template<class ParcelType>
void Foam::DSMCParcel<ParcelType>::writeFields
(
    const Cloud<DSMCParcel<ParcelType>>& c
)
{
    ParcelType::writeFields(c);

    const label np = c.size();
    const bool writeOnProc = c.size();

    IOField<vector> U(c.newIOobject("U", IOobject::NO_READ), np);
    IOField<scalar> Ei(c.newIOobject("Ei", IOobject::NO_READ), np);
    IOField<label> typeId(c.newIOobject("typeId", IOobject::NO_READ), np);

    label i = 0;
    for (const DSMCParcel<ParcelType>& p : c)
    {
        U[i] = p.U();
        Ei[i] = p.Ei();
        typeId[i] = p.typeId();
        ++i;
    }

    U.write(writeOnProc);
    Ei.write(writeOnProc);
    typeId.write(writeOnProc);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ParcelType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const DSMCParcel<ParcelType>& p
)
{
    if (os.format() == IOstreamOption::ASCII)
    {
        os  << static_cast<const ParcelType& >(p)
            << token::SPACE << p.U()
            << token::SPACE << p.Ei()
            << token::SPACE << p.typeId();
    }
    else
    {
        os  << static_cast<const ParcelType& >(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.U_),
            DSMCParcel<ParcelType>::sizeofFields
        );
    }

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
