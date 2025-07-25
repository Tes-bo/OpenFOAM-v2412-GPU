/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "solidMixtureProperties.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidMixtureProperties::solidMixtureProperties(const dictionary& dict)
:
    components_(),
    properties_()
{
    components_ = dict.toc();
    properties_.setSize(components_.size());

    forAll(components_, i)
    {
        // Handle sub-dictionary or primitive entry

        const dictionary* subDictPtr = dict.findDict(components_[i]);

        if (subDictPtr)
        {
            properties_.set
            (
                i,
                solidProperties::New(*subDictPtr)
            );
        }
        else
        {
            properties_.set
            (
                i,
                solidProperties::New(components_[i])
            );
        }
    }
}


Foam::solidMixtureProperties::solidMixtureProperties
(
    const solidMixtureProperties& s
)
:
    components_(s.components_),
    properties_(s.properties_.clone())
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::solidMixtureProperties> Foam::solidMixtureProperties::New
(
    const dictionary& thermophysicalProperties
)
{
    return autoPtr<solidMixtureProperties>::New(thermophysicalProperties);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::solidMixtureProperties::rho(const scalarField& Y) const
{
    scalar rrho = 0;

    forAll(properties_, i)
    {
        rrho += Y[i]/properties_[i].rho();
    }

    return 1/rrho;
}


Foam::scalar Foam::solidMixtureProperties::Cp(const scalarField& Y) const
{
    scalar Cp = 0;

    forAll(properties_, i)
    {
        Cp += Y[i]*properties_[i].Cp();
    }

    return Cp;
}


// ************************************************************************* //
