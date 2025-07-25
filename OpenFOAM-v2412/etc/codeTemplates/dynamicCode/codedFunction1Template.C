/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
    Copyright (C) YEAR AUTHOR,AFFILIATION
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

#include "codedFunction1Template.H"
#include "addToRunTimeSelectionTable.H"
#include "unitConversion.H"

//{{{ begin codeInclude
${codeInclude}
//}}} end codeInclude

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode
${localCode}
//}}} end localCode

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

// dynamicCode:
// SHA1 = ${SHA1sum}
//
// unique function name that can be checked if the correct library version
// has been loaded
extern "C" void ${typeName}_${SHA1sum}(bool load)
{
    if (load)
    {
        // Code that can be explicitly executed after loading
    }
    else
    {
        // Code that can be explicitly executed before unloading
    }
}


namespace Function1Types
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

//makeFunction1(${typeName}Function1, ${TemplateType});
defineTypeNameAndDebug
(
    ${typeName}Function1_${TemplateType},
    0
);
Function1<${TemplateType}>::addRemovabledictionaryConstructorToTable
    <${typeName}Function1_${TemplateType}>
    addRemovable${typeName}Function1_${TemplateType}ConstructorToTable_;

} // namespace Function1Types
} // namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Function1Types::
${typeName}Function1_${TemplateType}::
${typeName}Function1_${TemplateType}
(
    const word& entryName,
    const dictionary& dict,
    const objectRegistry* obrPtr
)
:
    Function1<${TemplateType}>(entryName, dict, obrPtr)
{
    if (${verbose:-false})
    {
        printMessage("Construct ${typeName} Function1 from dictionary");
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::${TemplateType}
Foam::Function1Types::${typeName}Function1_${TemplateType}::value
(
    const scalar x
) const
{
//{{{ begin code
    ${code}
//}}} end code
}


// ************************************************************************* //
