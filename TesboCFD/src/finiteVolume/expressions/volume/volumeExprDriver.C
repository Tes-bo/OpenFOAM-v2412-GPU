/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2022 OpenCFD Ltd.
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

#include "volumeExprDriver.H"
#include "volumeExprScanner.H"
#include "error.H"
#include "fvPatch.H"
#include "fvMesh.H"
#include "className.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace expressions
{
namespace volumeExpr
{

defineTypeNameAndDebug(parseDriver, 0);

addNamedToRunTimeSelectionTable
(
    fvExprDriver,
    parseDriver,
    dictionary,
    volume
);

addNamedToRunTimeSelectionTable
(
    fvExprDriver,
    parseDriver,
    idName,
    volume
);

addNamedToRunTimeSelectionTable
(
    fvExprDriver,
    parseDriver,
    dictionary,
    internalField
);

addNamedToRunTimeSelectionTable
(
    fvExprDriver,
    parseDriver,
    idName,
    internalField
);

} // End namespace volumeExpr
} // End namespace expressions
} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::expressions::volumeExpr::parseDriver::parseDriver
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    parsing::genericRagelLemonDriver(),
    expressions::fvExprDriver(dict),
    mesh_(mesh),
    resultType_(),
    isLogical_(false),
    hasDimensions_(false),
    fieldGeoType_(NO_DATA),
    resultDimensions_()
{
    resetTimeReference(nullptr);
    resetDb(mesh_.thisDb());
}


Foam::expressions::volumeExpr::parseDriver::parseDriver
(
    const fvMesh& mesh,
    const parseDriver& rhs,
    const dictionary& dict
)
:
    parsing::genericRagelLemonDriver(),
    expressions::fvExprDriver(rhs, dict),
    mesh_(mesh),
    resultType_(),
    isLogical_(false),
    hasDimensions_(false),
    fieldGeoType_(NO_DATA),
    resultDimensions_()
{
    resetTimeReference(nullptr);
    resetDb(mesh_.thisDb());
}


Foam::expressions::volumeExpr::parseDriver::parseDriver
(
    const word& meshName,
    const fvMesh& mesh
)
:
    parseDriver(mesh)
{
    //?? Info<< "Warn that meshName is ignored?" << nl;
}


Foam::expressions::volumeExpr::parseDriver::parseDriver
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    parsing::genericRagelLemonDriver(),
    expressions::fvExprDriver(dict),
    mesh_(mesh),
    resultType_(),
    isLogical_(false),
    fieldGeoType_(NO_DATA),
    resultDimensions_()
{
    resetTimeReference(nullptr);
    resetDb(mesh_.thisDb());
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::expressions::volumeExpr::parseDriver::readDict
(
    const dictionary& dict
)
{
    expressions::fvExprDriver::readDict(dict);

    resultDimensions_.clear();  // Avoid stickiness

    hasDimensions_ = resultDimensions_.readIfPresent("dimensions", dict);

    return true;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

unsigned Foam::expressions::volumeExpr::parseDriver::parse
(
    const std::string& expr,
    size_t pos,
    size_t len
)
{
    scanner scan(this->debugScanner());

    scan.process(expr, pos, len, *this);

    return 0;
}


void Foam::expressions::volumeExpr::parseDriver::clearField()
{
    resultField_.reset(nullptr);

    // Characteristics
    resultType_.clear();
    isLogical_ = false;
    fieldGeoType_ = NO_DATA;
}


Foam::autoPtr<Foam::regIOobject>
Foam::expressions::volumeExpr::parseDriver::dupZeroField() const
{
    const auto* regIOobjectPtr = resultField_.get();

    if (!regIOobjectPtr)
    {
        return nullptr;
    }

    autoPtr<regIOobject> zField;

    switch (fieldGeoType_)
    {
        #undef  doLocalCode
        #define doLocalCode(GeoField)                                         \
        {                                                                     \
            const auto* ptr = dynamic_cast<const GeoField*>(regIOobjectPtr);  \
            typedef typename GeoField::value_type Type;                       \
                                                                              \
            if (ptr)                                                          \
            {                                                                 \
                zField.reset                                                  \
                (                                                             \
                    GeoField::New                                             \
                    (                                                         \
                        word(pTraits<Type>::typeName) + word("(zero)"),       \
                        (*ptr).mesh(),                                        \
                        dimensioned<Type>(Zero),                              \
                        /* zeroGradient (volume) or calculated (other) */     \
                        defaultBoundaryType(*ptr)                             \
                    ).ptr()                                                   \
                );                                                            \
                break;                                                        \
            }                                                                 \
        }

        case FieldAssociation::VOLUME_DATA:
        {
            doLocalCode(volScalarField);
            doLocalCode(volVectorField);
            doLocalCode(volTensorField);
            doLocalCode(volSymmTensorField);
            doLocalCode(volSphericalTensorField);
            break;
        }
        case FieldAssociation::FACE_DATA:
        {
            doLocalCode(surfaceScalarField);
            doLocalCode(surfaceVectorField);
            doLocalCode(surfaceTensorField);
            doLocalCode(surfaceSymmTensorField);
            doLocalCode(surfaceSphericalTensorField);
            break;
        }
        case FieldAssociation::POINT_DATA:
        {
            doLocalCode(pointScalarField);
            doLocalCode(pointVectorField);
            doLocalCode(pointTensorField);
            doLocalCode(pointSymmTensorField);
            doLocalCode(pointSphericalTensorField);
            break;
        }
        default: break;
        #undef doLocalCode
    }

    // if (!zField)
    // {
    //     // Report
    // }

    return zField;
}


// ************************************************************************* //
