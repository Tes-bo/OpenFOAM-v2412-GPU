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

#include "codedFvOptionTemplate.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "unitConversion.H"
#include "fvMatrix.H"

//{{{ begin codeInclude
${codeInclude}
//}}} end codeInclude


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
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

defineTypeNameAndDebug(${typeName}FvOption${SourceType}, 0);
addRemovableToRunTimeSelectionTable
(
    option,
    ${typeName}FvOption${SourceType},
    dictionary
);

} // End namespace fv
} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::
${typeName}FvOption${SourceType}::
${typeName}FvOption${SourceType}
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::cellSetOption(name, modelType, dict, mesh)
{
    if (${verbose:-false})
    {
        printMessage("Construct ${typeName} fvOption from dictionary");
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::
${typeName}FvOption${SourceType}::
~${typeName}FvOption${SourceType}()
{
    if (${verbose:-false})
    {
        printMessage("Destroy ${typeName}");
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void
Foam::fv::
${typeName}FvOption${SourceType}::correct
(
    GeometricField<${TemplateType}, fvPatchField, volMesh>& fld
)
{
    if (${verbose:-false})
    {
        Info<< "${typeName}FvOption${SourceType}::correct()\n";
    }

//{{{ begin code
    ${codeCorrect}
//}}} end code
}


void
Foam::fv::
${typeName}FvOption${SourceType}::addSup
(
    fvMatrix<${TemplateType}>& eqn,
    const label fieldi
)
{
    if (${verbose:-false})
    {
        Info<< "${typeName}FvOption${SourceType}::addSup()\n";
    }

//{{{ begin code - warn/fatal if not implemented?
    ${codeAddSup}
//}}} end code
}


void
Foam::fv::
${typeName}FvOption${SourceType}::addSup
(
    const volScalarField& rho,
    fvMatrix<${TemplateType}>& eqn,
    const label fieldi
)
{
    if (${verbose:-false})
    {
        Info<< "${typeName}FvOption${SourceType}::addSup(rho)\n";
    }

//{{{ begin code - warn/fatal if not implemented?
    ${codeAddSupRho}
//}}} end code
}


void
Foam::fv::
${typeName}FvOption${SourceType}::constrain
(
    fvMatrix<${TemplateType}>& eqn,
    const label fieldi
)
{
    if (${verbose:-false})
    {
        Info<< "${typeName}FvOption${SourceType}::constrain()\n";
    }

//{{{ begin code
    ${codeConstrain}
//}}} end code
}


// ************************************************************************* //
