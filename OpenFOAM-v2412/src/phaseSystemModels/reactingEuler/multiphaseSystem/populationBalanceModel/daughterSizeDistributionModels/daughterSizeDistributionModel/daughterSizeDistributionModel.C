/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2018 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "daughterSizeDistributionModel.H"
#include "breakupModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
    defineTypeNameAndDebug(daughterSizeDistributionModel, 0);
    defineRunTimeSelectionTable(daughterSizeDistributionModel, dictionary);
}
}


// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::diameterModels::daughterSizeDistributionModel>
Foam::diameterModels::daughterSizeDistributionModel::New
(
    const breakupModel& breakup,
    const dictionary& dict
)
{
    const word modelType
    (
        dict.get<word>("daughterSizeDistributionModel")
    );

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "daughterSizeDistributionModel",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return ctorPtr(breakup, dict);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::diameterModels::daughterSizeDistributionModel::
daughterSizeDistributionModel
(
    const breakupModel& breakup,
    const dictionary& dict
)
:
    breakup_(breakup)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diameterModels::daughterSizeDistributionModel::
~daughterSizeDistributionModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::dimensionedScalar&
Foam::diameterModels::daughterSizeDistributionModel::
nik
(
    const label i,
    const label k
) const
{
    return nik_[k][i];
}


void Foam::diameterModels::daughterSizeDistributionModel::correct()
{
    if (nik_.size() == 0)
    {
        forAll(breakup_.popBal().sizeGroups(), k)
        {
            nik_.append(new PtrList<dimensionedScalar>());

            for (label i = 0; i <= k; i++)
            {
                nik_[k].append(new dimensionedScalar (this->calcNik(i, k)));
            }
        }
    }
}


// ************************************************************************* //
