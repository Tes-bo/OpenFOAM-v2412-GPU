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

#include "wallFunctionBlenders.H"
#include "dictionary.H"
#include "Enum.H"
#include "MinMax.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::wallFunctionBlenders::blenderType
>
Foam::wallFunctionBlenders::blenderTypeNames
({
    { blenderType::STEPWISE , "stepwise" },
    { blenderType::MAX , "max" },
    { blenderType::BINOMIAL , "binomial" },
    { blenderType::EXPONENTIAL, "exponential" },
    { blenderType::TANH, "tanh" }
});


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallFunctionBlenders::wallFunctionBlenders()
:
    blender_(blenderType::STEPWISE),
    n_(4.0)
{}


Foam::wallFunctionBlenders::wallFunctionBlenders
(
    const dictionary& dict,
    const blenderType blender,
    const scalar n
)
:
    blender_
    (
        blenderTypeNames.getOrDefault
        (
            "blending",
            dict,
            blender
        )
    ),
    n_
    (
        dict.getCheckOrDefault<scalar>
        (
            "n",
            n,
            scalarMinMax::ge(0)
        )
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::wallFunctionBlenders::writeEntries(Ostream& os) const
{
    os.writeEntry("blending", blenderTypeNames[blender_]);

    if (blender_ == blenderType::BINOMIAL)
    {
        os.writeEntry("n", n_);
    }
}


// ************************************************************************* //
