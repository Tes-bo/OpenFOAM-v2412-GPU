/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2010-2018 Bernhard Gschaider
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

#include "surfaceMesh.H"
#include "fvsPatchField.H"
#include "pointPatchField.H"
#include "typeInfo.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
bool Foam::expressions::fvExprDriver::isGlobalVariable
(
    const word& name,
    const bool wantPointData,
    const label expectedSize
) const
{
    DebugInfo
        << "Looking for global" << (wantPointData ? " point" : "")
        << " field name:" << name;

    const expressions::exprResult& result = lookupGlobal(name);

    DebugInfo
        << " - found (" << result.valueType() << ' '
        << result.isPointData() << ')';


    bool good =
        (result.isType<Type>() && result.isPointData(wantPointData));

    // Do size checking if requested
    if (good && expectedSize >= 0)
    {
        good = returnReduceAnd(result.size() == expectedSize);

        if (debug && !good)
        {
            Info<< " size is";
        }
    }

    DebugInfo << (good ? " good" : " bad") << endl;

    return good;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::expressions::fvExprDriver::getVariable
(
    const word& name,
    const label expectedSize,
    const bool mandatory
) const
{
    refPtr<expressions::exprResult> tvar;

    if (hasVariable(name) && variable(name).isType<Type>())
    {
        tvar.cref(variable(name));
    }
    else if (isGlobalVariable<Type>(name))
    {
        tvar.cref(lookupGlobal(name));
    }


    if (tvar.valid())
    {
        const auto& var = tvar.cref();

        const Field<Type>& vals = var.cref<Type>();

        if
        (
            expectedSize < 0
         || returnReduceAnd(vals.size() == expectedSize)
        )
        {
            // Return a copy of the field
            return tmp<Field<Type>>::New(vals);
        }

        if (!var.isUniform())
        {
            WarningInFunction
                << "Variable " << name
                << " is nonuniform and does not fit the size "
                << expectedSize << ". Using average" << endl;
        }

        return tmp<Field<Type>>::New(expectedSize, gAverage(vals));
    }

    if (mandatory)
    {
        FatalErrorInFunction
            << "Variable (" << name << ") not found." << nl
            << exit(FatalError);
    }

    return nullptr;
}


template<class Type>
bool Foam::expressions::fvExprDriver::foundField
(
    const word& name
) const
{
    if (debug)
    {
        Info<< "fvExprDriver::foundField. Name: " << name
            << " Type: " << Type::typeName
            << " registry:" << searchRegistry()
            << " disk:" << searchFiles() << endl;
    }


    for (int checki = 0; checki < 2; ++checki)
    {
        // Check 0: object context (first)
        // Check 1: regular objectRegistry
        const regIOobject* ioptr = nullptr;

        if (checki == 0)
        {
            ioptr = exprDriver::cfindContextIOobject(name);
        }
        else if (searchRegistry())
        {
            ioptr = this->mesh().cfindIOobject(name);
        }
        if (!ioptr) continue;

        const Type* fldPtr = dynamic_cast<const Type*>(ioptr);

        if (fldPtr)
        {
            if (debug)
            {
                if (checki)
                {
                    Info<< "Found registered:";
                }
                else
                {
                    Info<< "Found context object:";
                }
                Info<< name << endl;
            }
            return true;
        }
        else if (ioptr)
        {
            if (debug)
            {
                if (checki)
                {
                    Info<< "Registered:";
                }
                else
                {
                    Info<< "Context object:";
                }
                Info<< name << " type:"
                    << ioptr->headerClassName() << " != type:"
                    << Type::typeName << nl;
            }
        }
    }


    if (searchFiles() && getTypeOfField(name) == Type::typeName)
    {
        if (debug)
        {
            Info<< "Found file: " << name << nl;
        }
        return true;
    }

    if (debug)
    {
        Info<< name << " not found" << endl;
    }
    return false;
}


template<class Type>
bool Foam::expressions::fvExprDriver::isField
(
    const word& name,
    bool wantPointData,
    label
) const
{
    if (debug)
    {
        Info<< "fvExprDriver::isField <" << name << '>' << endl;
    }

    return
    (
        wantPointData
      ? this->foundField<PointField<Type>>(name)
      :
        (
            this->foundField<VolumeField<Type>>(name)
         || this->foundField<SurfaceField<Type>>(name)
        )
    );
}


template<class GeomField, class Mesh>
Foam::tmp<GeomField> Foam::expressions::fvExprDriver::getOrReadFieldImpl
(
    const word& name,
    const Mesh& meshRef,
    bool mandatory,
    bool getOldTime
)
{
    typedef typename GeomField::value_type Type;

    tmp<GeomField> tfield;

    if (debug)
    {
        Info<< "fvExprDriver::getOrReadField <" << name
            << "> Type: " << GeomField::typeName << endl;
    }


    // Handle variables
    // ~~~~~~~~~~~~~~~~

    refPtr<expressions::exprResult> tvar;

    if (hasVariable(name) && variable(name).isType<Type>())
    {
        tvar.cref(variable(name));
    }
    else if (isGlobalVariable<Type>(name))
    {
        tvar.cref(lookupGlobal(name));
    }

    if (tvar.valid())
    {
        const auto& var = tvar.cref();
        const Type deflt(var.getValue<Type>());

        if (debug)
        {
            Info<< "Getting " << name << " from variables. Default: "
                << deflt << endl;
        }

        if (debug)
        {
            Info<< "Creating field " << name << " of type "
                << GeomField::typeName << nl;
        }

        tfield = GeomField::New
        (
            name,
            meshRef,
            dimensioned<Type>(deflt),
            // Patch is zeroGradient (volFields) or calculated (other)
            defaultBoundaryType(GeomField::null())
        );
        auto& fld = tfield.ref();

        if (debug)
        {
            Info<< "New field: " << name << " ownedByRegistry"
                << fld.ownedByRegistry() << endl;
        }

        const Field<Type>& vals = var.cref<Type>();

        if (debug)
        {
            Pout<< "sizes: " << vals.size() << ' ' << fld.size() << endl;
        }

        if (returnReduceAnd(vals.size() == fld.size()))
        {
            fld.primitiveFieldRef() = vals;
        }
        else
        {
            const Type avg = gAverage(vals);

            bool noWarn = false;

            if (!noWarn)
            {
                MinMax<Type> range = gMinMax(vals);

                if (range.mag() > SMALL)
                {
                    WarningInFunction
                        << "The min/max ranges differ " << range
                        << " - using average " << avg << nl;
                }
            }

            fld.primitiveFieldRef() = avg;
        }

        correctField(fld);

        return tfield;
    }


    // Find context or registered field
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const GeomField* origFldPtr = nullptr;

    for (int checki = 0; !origFldPtr && checki < 2; ++checki)
    {
        // Check 0: object context (first)
        // Check 1: regular objectRegistry

        if (checki == 0)
        {
            origFldPtr = exprDriver::cfindContextObject<GeomField>(name);
        }
        else if (searchRegistry())
        {
            origFldPtr =
                meshRef.thisDb().template cfindObject<GeomField>(name);
        }
    }

    if (origFldPtr)
    {
        // Found from context or registry

        if (debug)
        {
            Info<< "Retrieve context/registered:" << name << nl;
        }

        const GeomField& origFld = *origFldPtr;

        // Make a deep copy of the data to return. Avoids shadowing
        // the original object, but most importantly the backend
        // parser (eg, lemon) will be working with combining via plain
        // pointers. We thus lose any of the tmp<> shallow copy semantics
        // anyhow. Additionally, need to disable dimension checking here or
        // elsewhere too.

        tfield = GeomField::New(name + "_exprDriverCopy", origFld);

        if (getOldTime)
        {
            if (debug)
            {
                Info<< "Getting oldTime of " << name << " has "
                    << origFld.nOldTimes() << endl;
            }

            if (!origFld.nOldTimes() && this->prevIterIsOldTime())
            {
                if (debug)
                {
                    Info<< "No oldTime, using previous iteration" << endl;
                }
                tfield.ref().oldTime() = origFld.prevIter();
            }
        }
    }
    else if (searchFiles() && getTypeOfField(name) == GeomField::typeName)
    {
        if (debug)
        {
            Info<< "Reading " << name << " from disc" << endl;
        }

        tfield.reset
        (
            this->readAndRegister<GeomField>(name, meshRef)
        );
        // oldTime automatically read
    }


    if (debug)
    {
        Info<< "field: valid()=" << Switch::name(tfield.valid()) << endl;
    }

    if (tfield.valid())
    {
        GeomField& fld = tfield.ref();

        if (debug)
        {
            Info<< "Valid " << name << " found. Removing dimensions" << nl;
        }

        fld.dimensions().clear();

        if (fld.nOldTimes())
        {
            if (debug)
            {
                Info<< "Removing dimensions of oldTime of " << name
                    << " has " << fld.nOldTimes() << nl;
            }

            // Switch dimension checking off
            const bool oldDimChecking = dimensionSet::checking(false);

            // go through ALL old times
            GeomField* fp = &(fld);

            while (fp->nOldTimes())
            {
                fp = &(fp->oldTime());
                fp->dimensions().clear();
            }

            // Restore old value of dimension checking
            dimensionSet::checking(oldDimChecking);
        }
    }
    else if (mandatory)
    {
        FatalErrorInFunction
            << "Could not find field " << name
            << " in registry or on file-system" << nl
            << exit(FatalError);
    }

    return tfield;
}


template<class T>
Foam::autoPtr<T> Foam::expressions::fvExprDriver::getTopoSet
(
    const fvMesh& mesh,
    const word& name,
    SetOrigin& origin
) const
{
    // Avoid possible name clashes
    const word regName = name + "RegisteredNameFor" + T::typeName;

    if (debug)
    {
        Info<< "Looking for " << T::typeName << " named " << name
            << " or registered as " << regName << " with mesh "
            << "Caching:" << cacheSets()
            << " Found:" << (mesh.foundObject<T>(name))
            << " Found registered:" << mesh.foundObject<T>(regName)
            << endl;
    }


    origin = SetOrigin::INVALID;
    autoPtr<T> setPtr;

    if
    (
        !cacheSets()
     ||
        (
            !mesh.thisDb().foundObject<T>(regName)
         && !mesh.thisDb().foundObject<T>(name)
        )
    )
    {
        if (debug)
        {
            Info<< "Constructing new " << T::typeName << ' ' << name << nl;

            if (debug > 1)
            {
                Pout<< mesh.thisDb().names();
            }
        }

        origin = SetOrigin::FILE;
        setPtr.reset(new T(mesh, name, IOobject::MUST_READ));

        if (cacheSets())
        {
            if (debug)
            {
                Info<< "Registering a copy of " << name << " with mesh" << nl;
            }

            autoPtr<T> toCache(new T(mesh, regName, *setPtr));
            toCache->store(toCache);
        }
    }
    else
    {
        const T* ptr = mesh.thisDb().template cfindObject<T>(name);

        if (ptr)
        {
            if (debug)
            {
                Info<< "Getting existing " << name << endl;
            }

            origin = SetOrigin::MEMORY;
            setPtr.reset(new T(mesh, name, *ptr));
        }
        else
        {
            if (debug)
            {
                Info<< "Getting existing " << regName << endl;
            }

            origin = SetOrigin::CACHE;
            setPtr.reset(new T(mesh, name, mesh.lookupObject<T>(regName)));
        }
    }


    return setPtr;
}


template<class T>
bool Foam::expressions::fvExprDriver::updateSet
(
    autoPtr<T>& setPtr,
    const word& name,
    SetOrigin origin
) const
{
    const label oldSize = setPtr->size();

    bool updated = false;
    const auto& mesh = refCast<const polyMesh>(setPtr->db());

    if (debug)
    {
        Info<< "UpdateSet: " << setPtr->name() << " Id: " << name
            << " Origin: " << int(origin) << endl;
    }

    switch (origin)
    {
        case SetOrigin::FILE:
        {
            IOobject header
            (
                name,
                mesh.time().timeName(),
                polyMesh::meshSubDir/"sets",
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            );

            if (header.typeHeaderOk<T>())
            {
                if (debug)
                {
                    Pout<< "Rereading from "
                        << header.localFilePath(T::typeName) << endl;
                }
                setPtr.reset(new T(header));
                updated = true;
            }
            break;
        }

        case SetOrigin::NEW:
        case SetOrigin::MEMORY:
        case SetOrigin::CACHE:
        {
            if (origin == SetOrigin::NEW)
            {
                WarningInFunction
                    << "State NEW shouldn't exist"
                    << endl;
            }

            word sName = name;

            const T* ptr = mesh.thisDb().template cfindObject<T>(name);

            if (ptr)
            {
                if (debug)
                {
                    Info<< "Found " << name
                        << " and rereading it" << endl;
                }

                setPtr.reset(new T(mesh, name, *ptr));
            }
            else
            {
                FatalErrorInFunction
                    << name << " Not found" << endl
                    << "In registry: " << mesh.thisDb().names() << endl
                    << exit(FatalError);
            }
            updated = true;
            break;
        }

        case INVALID:
        {
            FatalErrorInFunction
                << T::typeName << ' ' << name << " is invalid" << endl
                << exit(FatalError);
            break;
        }

        default:
        {
            if (debug)
            {
                Info<< "Origin " << int(origin) << " not implemented" << endl;
            }
            break;
        }
    }

    if (debug)
    {
        Pout<< name << " old size " << oldSize << " new: "
            << setPtr->size() << endl;
    }

    return updated;
}


// ************************************************************************* //
