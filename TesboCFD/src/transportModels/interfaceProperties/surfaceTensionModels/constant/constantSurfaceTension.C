/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
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

#include "constantSurfaceTension.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceTensionModels
{
    defineTypeNameAndDebug(constant, 0);
    addToRunTimeSelectionTable(surfaceTensionModel, constant, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceTensionModels::constant::constant
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    surfaceTensionModel(mesh),
    sigma_("sigma", dimMass/sqr(dimTime), dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfaceTensionModels::constant::~constant()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::surfaceTensionModels::constant::sigma() const
{
    return volScalarField::New
    (
        "sigma",
        IOobject::NO_REGISTER,
        mesh_,
        sigma_
    );
}


bool Foam::surfaceTensionModels::constant::readDict(const dictionary& dict)
{
    // Handle sub-dictionary format as a special case

    const dictionary* subDictPtr = dict.findDict("sigma");

    if (subDictPtr)
    {
        subDictPtr->readEntry("sigma", sigma_);
    }
    else
    {
        dict.readEntry("sigma", sigma_);
    }

    return true;
}


bool Foam::surfaceTensionModels::constant::writeData(Ostream& os) const
{
    if (surfaceTensionModel::writeData(os))
    {
        os  << sigma_;
        os.endEntry();
        return os.good();
    }

    return false;
}


// ************************************************************************* //
