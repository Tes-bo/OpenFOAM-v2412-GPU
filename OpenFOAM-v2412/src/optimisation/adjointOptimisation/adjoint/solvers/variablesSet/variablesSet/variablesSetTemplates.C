/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Copyright (C) 2013-2019 FOSS GP
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "localIOdictionary.H"
#include "FieldField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
GeometricField<Type, PatchField, GeoMesh>* variablesSet::allocateNamedField
(
    const fvMesh& mesh,
    const IOobject& io,
    const word& solverName
)
{
    typedef GeometricField<Type, PatchField, GeoMesh> fieldType;

    // Read-in boundary conditions from given IOobject
    localIOdictionary dict
    (
        IOobject
        (
            io.name(),
            io.instance(),
            io.local(),
            io.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        fieldType::typeName
    );
    dictionary& bField(dict.subDict("boundaryField"));

    // Add solverName to all patch entries.
    // Reduntant if not adjoint fields, but overhead should be small
    for (entry& dEntry : bField)
    {
        if (dEntry.isDict())
        {
            dEntry.dict().add<word>("solverName", solverName, true);
        }
    }
    DebugInfo
        << bField << endl;

    return (new fieldType(io, mesh, dict));
}


template<class Type, template<class> class PatchField, class GeoMesh>
bool variablesSet::readFieldOK
(
    autoPtr<GeometricField<Type, PatchField, GeoMesh>>& fieldPtr,
    const fvMesh& mesh,
    const word& baseName,
    const word& solverName,
    const bool useSolverNameForFields
)
{
    typedef GeometricField<Type, PatchField, GeoMesh> fieldType;

    word customName = baseName + solverName;
    IOobject headerCustomName
    (
        customName,
        mesh.time().timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

    IOobject headerBaseName
    (
        baseName,
        mesh.time().timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

    bool fieldFound(false);

    // Read field with full name (i.e. baseName plus solverName) if present
    if
    (
        headerCustomName.typeHeaderOk<fieldType>(false)
     && useSolverNameForFields
    )
    {
        fieldPtr.reset
        (
            allocateNamedField<Type, PatchField, GeoMesh>
            (
                mesh,
                headerCustomName,
                solverName
            )
        );
        fieldFound = true;
    }
    // else, see whether the base field exists
    else if (headerBaseName.typeHeaderOk<fieldType>(false))
    {
        fieldPtr.reset
        (
            allocateNamedField<Type, PatchField, GeoMesh>
            (
                mesh,
                headerBaseName,
                solverName
            )
        );

        // Rename field if necessary
        if (useSolverNameForFields)
        {
            Info<< "Field " << customName << " not found" << endl;
            Info<< "Reading base field " << baseName << " and renaming ... "
                << endl;
            fieldPtr.ref().rename(customName);
        }
        fieldFound = true;
    }

    return fieldFound;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
autoPtr<GeometricField<Type, PatchField, GeoMesh>>
variablesSet::allocateRenamedField
(
    const autoPtr<GeometricField<Type, PatchField, GeoMesh>>& bf
)
{
    typedef GeometricField<Type, PatchField, GeoMesh> fieldType;
    autoPtr<fieldType> returnField(nullptr);
    if (bf)
    {
        const word timeName = bf().mesh().time().timeName();
        returnField.reset
        (
            new fieldType
            (
                bf().name() + timeName,
                bf()
            )
        );
    }
    return returnField;
}


template<class Type, template<class> class PatchField, class GeoMesh>
void variablesSet::swapAndRename
(
    autoPtr<GeometricField<Type, PatchField, GeoMesh>>& p1,
    autoPtr<GeometricField<Type, PatchField, GeoMesh>>& p2
)
{
    // Swapping pointers is OK for the mean flow fields known by the
    // variablesSet (and, in essence, by the solver).
    // The problem is that turbulence models know references to U and phi
    // which cannot be swapped.
    /*
    const word name1 = p1().name();
    const word name2 = p2().name();
    p1.swap(p2);

    p2().rename("temp");
    p1().rename(name1);
    p2().rename(name2);
    */

    // Copy back-up fields to original instead. Slower but there seems to be
    // no other way
    GeometricField<Type, PatchField, GeoMesh> temp("temp", p1());
    p1() == p2();
    p2() == temp;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void variablesSet::setField
(
    autoPtr<GeometricField<Type, fvPatchField, volMesh>>& fieldPtr,
    const fvMesh& mesh,
    const word& baseName,
    const word& solverName,
    const bool useSolverNameForFields
)
{
    // Try to read in field with custom or base name
    bool fieldFound
    (
        readFieldOK
        (
            fieldPtr,
            mesh,
            baseName,
            solverName,
            useSolverNameForFields
        )
    );

    // No base or custom field found. This is fatal
    if (!fieldFound)
    {
        FatalErrorInFunction
            << "Could not read field with custom ("
            << word(baseName + solverName) << ") "
            << "or base (" << baseName << ") name"
            << exit(FatalError);
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>> variablesSet::allocateField
(
    const fvMesh& mesh,
    const word& baseName,
    const word& solverName,
    const bool useSolverNameForFields
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;

    autoPtr<VolFieldType> fieldPtr(nullptr);
    setField(fieldPtr, mesh, baseName, solverName, useSolverNameForFields);

    return tmp<VolFieldType>(fieldPtr.ptr());
}


template<class Type>
void variablesSet::renameTurbulenceField
(
    GeometricField<Type, fvPatchField, volMesh>& baseField,
    const word& solverName
)
{
    // typedefs
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;
    typedef typename VolFieldType::Boundary Boundary;

    // Name of custom field, to be potentially read in
    const word baseName = baseField.name();
    const word customName = baseName + solverName;
    const fvMesh& mesh = baseField.mesh();

    // Renaming of the base field
    baseField.rename(customName);

    // Create field with baseName and write it, to enable continuation
    // Note: gives problems for multi-point runs since we end up with
    // multiple db entries with the same name (one from here and one from
    // the solver that will construct a turbulenceModel).
    // Handled through solver.write() for now
    /*
    if (!mesh.foundObject<VolFieldType>(baseName))
    {
        autoPtr<VolFieldType> baseCopy(new VolFieldType(baseField));
        baseCopy().IOobject::writeOpt(baseField.writeOpt());
        baseCopy().rename(baseName);
        regIOobject::store(baseCopy);
    }
    */

    // Check whether a field with the custom name exists, read it in and
    // set supplied base field to that
    IOobject headerCustomName
    (
        customName,
        mesh.time().timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE,
        IOobject::NO_REGISTER
    );

    if (headerCustomName.typeHeaderOk<VolFieldType>(true))
    {
        Info<< "Reading custom turbulence field " << customName
            << " and replacing " << baseName << nl << endl;
        VolFieldType customField(headerCustomName, mesh);

        // Copy internalfield
        baseField.primitiveFieldRef() = customField.primitiveField();

        // We might apply different boundary conditions per operating point
        // We need to read them from the custom files and substitute the ones
        // known by the turbulence model field
        Boundary& baseBoundary = baseField.boundaryFieldRef();
        Boundary& customBoundary = customField.boundaryFieldRef();
        forAll(baseBoundary, patchI)
        {
            baseBoundary.set
            (
                patchI,
                customBoundary[patchI].clone(baseField.ref())
            );
        }
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void variablesSet::nullifyField
(
    GeometricField<Type, PatchField, GeoMesh>& field
)
{
    typedef GeometricField<Type, PatchField, GeoMesh> fieldType;
    field == dimensioned<Type>(field.dimensions(), Zero);
    if (field.nOldTimes())
    {
        fieldType& oldTime = field.oldTime();
        variablesSet::nullifyField(oldTime);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
