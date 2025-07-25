/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
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

#include "dsmcFields.H"
#include "volFields.H"
#include "dictionary.H"
#include "dsmcCloud.H"
#include "constants.H"
#include "stringListOps.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(dsmcFields, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        dsmcFields,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

static const word& filteredName
(
    const word& baseName,
    const wordList& names,
    const string& scopePrefix
)
{
    label idx = names.find(baseName);

    if (idx < 0 && !scopePrefix.empty())
    {
        // Take the first matching item
        idx = firstMatchingString(regExp(scopePrefix + baseName), names);
    }

    if (idx < 0)
    {
        return word::null;
    }

    return names[idx];
}

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::dsmcFields::dsmcFields
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::dsmcFields::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    return true;
}


bool Foam::functionObjects::dsmcFields::execute()
{
    return true;
}


bool Foam::functionObjects::dsmcFields::write()
{
    // This is fairly horrible with too many hard-coded names...

    // Pre-filter names to obtain 'Mean' vol fields
    const wordList allMeanNames
    (
        obr_.sortedNames
        (
            regExp("vol.*Field"),  // Any vol field type
            regExp(".+Mean")       // Mean field names
        )
    );

    // The separator is often ':', but could be something else.
    // Replace as first char in [..], so that the regex remains valid,
    // even if the separator happens to be '-'.

    string scopePrefix = ".+[_:]";
    scopePrefix[3] = IOobject::scopeSeparator;


    // Find scoped/unscoped field name and do lookup.
    // Short-circuit with message if not found (name or field)

    // Note: currently just find a match without and with a scoping prefix
    // but could refine to pick the longest name etc, or after finding
    // the first matching field, use the same prefix for all subsequent fields

    #undef  doLocalCode
    #define doLocalCode(Name, FieldType, Member)                               \
                                                                               \
        const FieldType* Member##Ptr = nullptr;                                \
        {                                                                      \
            const word& fldName =                                              \
                filteredName(Name, allMeanNames, scopePrefix);                 \
                                                                               \
            if (!fldName.empty())                                              \
            {                                                                  \
                Member##Ptr = obr_.cfindObject<FieldType>(fldName);            \
            }                                                                  \
                                                                               \
            if (returnReduceOr(!Member##Ptr))                                  \
            {                                                                  \
                Log << type() << ' ' << name() << " : no " << Name             \
                    << " field found - not calculating\n";                     \
                return false;                                                  \
            }                                                                  \
        }                                                                      \
        /* Define the const reference */                                       \
        const FieldType& Member = *Member##Ptr;


    // rhoNMean: always required
    doLocalCode("rhoNMean", volScalarField, rhoNMean);

    // Also check for division by zero
    {
        const scalar minval = min(mag(rhoNMean)).value();

        if (minval <= VSMALL)
        {
            Log << type() << ' ' << name()
                << " : Small value (" << minval << ") in rhoNMean field"
                << " - not calculating to avoid division by zero" << nl;
            return false;
        }
    }


    // The other fields

    doLocalCode("rhoMMean", volScalarField, rhoMMean);
    doLocalCode("momentumMean", volVectorField, momentumMean);
    doLocalCode("linearKEMean", volScalarField, linearKEMean);
    doLocalCode("internalEMean", volScalarField, internalEMean);
    doLocalCode("iDofMean", volScalarField, iDofMean);
    doLocalCode("fDMean", volVectorField, fDMean);
    #undef doLocalCode


    //
    // Everything seem to be okay - can execute
    //
    {
        Log << "Calculating dsmcFields." << endl;

        Log << "    Calculating UMean field." << nl;
        volVectorField UMean
        (
            IOobject
            (
                "UMean",
                obr_.time().timeName(),
                obr_,
                IOobject::NO_READ
            ),
            momentumMean/rhoMMean
        );

        Log << "    Calculating translationalT field." << endl;
        volScalarField translationalT
        (
            IOobject
            (
                "translationalT",
                obr_.time().timeName(),
                obr_,
                IOobject::NO_READ
            ),

            2.0/(3.0*physicoChemical::k.value()*rhoNMean)
           *(linearKEMean - 0.5*rhoMMean*(UMean & UMean))
        );

        Log << "    Calculating internalT field." << endl;
        volScalarField internalT
        (
            IOobject
            (
                "internalT",
                obr_.time().timeName(),
                obr_,
                IOobject::NO_READ
            ),
            (2.0/physicoChemical::k.value())*(internalEMean/iDofMean)
        );

        Log << "    Calculating overallT field." << endl;
        volScalarField overallT
        (
            IOobject
            (
                "overallT",
                obr_.time().timeName(),
                obr_,
                IOobject::NO_READ
            ),
            2.0/(physicoChemical::k.value()*(3.0*rhoNMean + iDofMean))
           *(linearKEMean - 0.5*rhoMMean*(UMean & UMean) + internalEMean)
        );

        Log << "    Calculating pressure field." << endl;
        volScalarField p
        (
            IOobject
            (
                "p",
                obr_.time().timeName(),
                obr_,
                IOobject::NO_READ
            ),
            physicoChemical::k.value()*rhoNMean*translationalT
        );

        volScalarField::Boundary& pBf = p.boundaryFieldRef();

        forAll(mesh_.boundaryMesh(), i)
        {
            const polyPatch& patch = mesh_.boundaryMesh()[i];

            if (isA<wallPolyPatch>(patch))
            {
                pBf[i] =
                    fDMean.boundaryField()[i]
                  & (patch.faceAreas()/mag(patch.faceAreas()));
            }
        }


        // Report

        Log << "    mag(UMean) max/min : "
            << max(mag(UMean)).value() << token::SPACE
            << min(mag(UMean)).value() << nl

            << "    translationalT max/min : "
            << max(translationalT).value() << token::SPACE
            << min(translationalT).value() << nl

            << "    internalT max/min : "
            << max(internalT).value() << token::SPACE
            << min(internalT).value() << nl

            << "    overallT max/min : "
            << max(overallT).value() << token::SPACE
            << min(overallT).value() << nl

            << "    p max/min : "
            << max(p).value() << token::SPACE
            << min(p).value() << endl;


        // Write
        UMean.write();

        translationalT.write();

        internalT.write();

        overallT.write();

        p.write();
    }

    Log << "dsmcFields written." << nl << endl;
    return true;
}


// ************************************************************************* //
