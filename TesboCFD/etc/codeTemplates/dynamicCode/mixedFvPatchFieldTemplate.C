/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
    Copyright (C) YEAR AUTHOR, AFFILIATION
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

#include "mixedFvPatchFieldTemplate.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "unitConversion.H"
#include "PatchFunction1.H"

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

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

makeRemovablePatchTypeField
(
    fvPatch${FieldType},
    ${typeName}MixedValueFvPatch${FieldType}
);

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::
${typeName}MixedValueFvPatch${FieldType}::
${typeName}MixedValueFvPatch${FieldType}
(
    const fvPatch& p,
    const DimensionedField<${TemplateType}, volMesh>& iF
)
:
    parent_bctype(p, iF)
{
    if (${verbose:-false})
    {
        printMessage("Construct ${typeName} : patch/DimensionedField");
    }
}


Foam::
${typeName}MixedValueFvPatch${FieldType}::
${typeName}MixedValueFvPatch${FieldType}
(
    const ${typeName}MixedValueFvPatch${FieldType}& rhs,
    const fvPatch& p,
    const DimensionedField<${TemplateType}, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    parent_bctype(rhs, p, iF, mapper)
{
    if (${verbose:-false})
    {
        printMessage("Construct ${typeName} : patch/DimensionedField/mapper");
    }
}


Foam::
${typeName}MixedValueFvPatch${FieldType}::
${typeName}MixedValueFvPatch${FieldType}
(
    const fvPatch& p,
    const DimensionedField<${TemplateType}, volMesh>& iF,
    const dictionary& dict
)
:
    parent_bctype(p, iF, dict)
{
    if (${verbose:-false})
    {
        printMessage("Construct ${typeName} : patch/dictionary");
    }
}


Foam::
${typeName}MixedValueFvPatch${FieldType}::
${typeName}MixedValueFvPatch${FieldType}
(
    const ${typeName}MixedValueFvPatch${FieldType}& rhs
)
:
    parent_bctype(rhs),
    dictionaryContent(rhs)
{
    if (${verbose:-false})
    {
        printMessage("Copy construct ${typeName}");
    }
}


Foam::
${typeName}MixedValueFvPatch${FieldType}::
${typeName}MixedValueFvPatch${FieldType}
(
    const ${typeName}MixedValueFvPatch${FieldType}& rhs,
    const DimensionedField<${TemplateType}, volMesh>& iF
)
:
    parent_bctype(rhs, iF)
{
    if (${verbose:-false})
    {
        printMessage("Construct ${typeName} : copy/DimensionedField");
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::
${typeName}MixedValueFvPatch${FieldType}::
~${typeName}MixedValueFvPatch${FieldType}()
{
    if (${verbose:-false})
    {
        printMessage("Destroy ${typeName}");
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void
Foam::
${typeName}MixedValueFvPatch${FieldType}::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (${verbose:-false})
    {
        printMessage("updateCoeffs ${typeName}");
    }

//{{{ begin code
    ${code}
//}}} end code

    this->parent_bctype::updateCoeffs();
}


// ************************************************************************* //
