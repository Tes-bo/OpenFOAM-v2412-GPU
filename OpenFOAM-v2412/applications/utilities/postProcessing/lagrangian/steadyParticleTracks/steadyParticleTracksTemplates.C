/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2023 OpenCFD Ltd.
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

#include "steadyParticleTracksTemplates.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::IOField<Type>> Foam::readParticleField
(
    const word& fieldName,
    const IOobjectList& cloudObjects
)
{
    const IOobject* io = cloudObjects.cfindObject<IOField<Type>>(fieldName);
    if (io)
    {
        return tmp<IOField<Type>>::New(*io);
    }

    FatalErrorInFunction
        << "Cloud field name " << fieldName
        << " not found or the incorrect type"
        << abort(FatalError);

    return nullptr;
}


template<class Type>
void Foam::writeVTK(OFstream& os, const Type& value)
{
    os  << component(value, 0);
    for (label d=1; d < pTraits<Type>::nComponents; ++d)
    {
        os  << ' ' << component(value, d);
    }
}


template<class Type>
void Foam::writeVTKField
(
    OFstream& os,
    const IOField<Type>& field,
    const List<labelList>& addr
)
{
    const label step = max(1, floor(8/pTraits<Type>::nComponents));

    Info<< "        writing field " << field.name() << endl;
    os  << nl << field.name() << ' '
        << int(pTraits<Type>::nComponents) << ' '
        << field.size() << " float" << nl;

    ///label offset = 0;
    for (const labelList& ids : addr)
    {
        List<Type> data(UIndirectList<Type>(field, ids));
        label nData = data.size() - 1;
        forAll(data, i)
        {
            writeVTK<Type>(os, data[i]);
            if (((i + 1) % step == 0) || (i == nData))
            {
                os  << nl;
            }
            else
            {
                os  << ' ';
            }
        }
        /// offset += ids.size();
    }
}


template<class Type>
void Foam::processFields
(
    OFstream& os,
    const List<labelList>& addr,
    const IOobjectList& cloudObjects
)
{
    for (const IOobject& io : cloudObjects.csorted<IOField<Type>>())
    {
        Info<< "        reading field " << io.name() << endl;
        IOField<Type> field(io);

        writeVTKField<Type>(os, field, addr);
    }
}


// ************************************************************************* //
