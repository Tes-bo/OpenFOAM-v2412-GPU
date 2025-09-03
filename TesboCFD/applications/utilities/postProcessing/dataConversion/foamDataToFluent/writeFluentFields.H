/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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

InClass
    Foam::writeFluentFields

\*---------------------------------------------------------------------------*/

#ifndef Foam_writeFluentFields_H
#define Foam_writeFluentFields_H

#include "volFields.H"
#include "dictionary.H"
#include "IOobjectList.H"
#include "Ostream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

void writeFluentField
(
    const volScalarField& phi,
    const label fluentFieldIdentifier,
    Ostream& os
);

void writeFluentField
(
    const volVectorField& phi,
    const label fluentFieldIdentifier,
    Ostream& os
);


template<class GeoField>
void readFieldsAndWriteFluent
(
    const dictionary& dict,
    const IOobjectList& objects,
    const fvMesh& mesh,
    Ostream& os
)
{
    for (const IOobject& io : objects.csorted<GeoField>())
    {
        // Lookup field from dictionary and convert field
        const word& fieldName = io.name();
        label unitNumber;
        if
        (
            dict.readIfPresent(fieldName, unitNumber)
         && unitNumber > 0
        )
        {
            // Read field
            GeoField field(io, mesh);

            Info<< "    Converting field " << fieldName << nl;
            writeFluentField(field, unitNumber, os);
        }
    }
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
