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

#include "makeCombustionTypes.H"

#include "psiReactionThermo.H"
#include "rhoReactionThermo.H"
#include "EDC.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::combustionModels::EDCversions
>
Foam::combustionModels::EDCversionNames
({
    { EDCversions::v1981, "v1981" },
    { EDCversions::v1996, "v1996" },
    { EDCversions::v2005, "v2005" },
    { EDCversions::v2016, "v2016" },
});

const Foam::combustionModels::EDCversions
Foam::combustionModels::EDCdefaultVersion(EDCversions::v2005);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

makeCombustionTypes(EDC, psiReactionThermo);
makeCombustionTypes(EDC, rhoReactionThermo);

}

// ************************************************************************* //
