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

#include "primitivePatchInterpolation.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::expressions::patchExpr::parseDriver::getVariableIfAvailable
(
    const word& name
) const
{
    bool hasPointData = false;

    refPtr<expressions::exprResult> tvar;

    if (hasVariable(name) && variable(name).isType<Type>())
    {
        tvar.cref(variable(name));
        hasPointData = tvar().isPointData();
    }
    else if (isGlobalVariable<Type>(name))
    {
        tvar.cref(lookupGlobal(name));
    }


    if (tvar.valid())
    {
        const auto& var = tvar.cref();
        const Field<Type>& vals = var.cref<Type>();

        const label len = (hasPointData ? this->pointSize() : this->size());

        if (returnReduceAnd(vals.size() == len))
        {
            // Return a copy of the field
            return tmp<Field<Type>>::New(vals);
        }

        if (!var.isUniform())
        {
            WarningInFunction
                << "Variable " << name
                << " is nonuniform and does not fit the size"
                << ". Using average" << endl;
        }

        return tmp<Field<Type>>::New(this->size(), gAverage(vals));
    }

    return nullptr;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::expressions::patchExpr::parseDriver::getVolField(const word& name)
{
    return getField<Type>(name);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::expressions::patchExpr::parseDriver::getSurfaceField(const word& name)
{
    return getField<Type>(name);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::expressions::patchExpr::parseDriver::getPointField(const word& name)
{
    return getField<Type>(name);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::expressions::patchExpr::parseDriver::getField(const word& name)
{
    tmp<Field<Type>> tfield = getVariableIfAvailable<Type>(name);

    if (tfield.valid())
    {
        return tfield;
    }

    const objectRegistry& obr = this->mesh().thisDb();
    const label patchIndex = patch_.index();


    // Local, temporary storage and/or lookup values
    bool found = false;
    tmp<VolumeField<Type>> vfield;
    tmp<SurfaceField<Type>> sfield;
    tmp<PointField<Type>> pfield;

    for (int checki = 0; !found && checki < 2; ++checki)
    {
        // Check 0: object context (first)
        // Check 1: regular objectRegistry
        const regIOobject* ioptr =
        (
            (checki == 0)
          ? exprDriver::cfindContextIOobject(name)
          : obr.cfindIOobject(name)
        );
        if (!ioptr) continue;

        if (!found)
        {
            vfield.cref(dynamic_cast<const VolumeField<Type>*>(ioptr));
            found = vfield.valid();
        }
        if (!found)
        {
            sfield.cref(dynamic_cast<const SurfaceField<Type>*>(ioptr));
            found = sfield.valid();
        }
        if (!found)
        {
            pfield.cref(dynamic_cast<const PointField<Type>*>(ioptr));
            found = pfield.valid();
        }
    }


    // Finally, search files if necessary (and permitted)
    if (!found && searchFiles())
    {
        const word fldType = this->getTypeOfField(name);

        if (fldType == VolumeField<Type>::typeName)
        {
            vfield = this->readAndRegister<VolumeField<Type>>(name, mesh());
        }
        else if (fldType == SurfaceField<Type>::typeName)
        {
            sfield = this->readAndRegister<SurfaceField<Type>>(name, mesh());
        }
        else if (fldType == PointField<Type>::typeName)
        {
            pfield = this->readAndRegister<PointField<Type>>
            (
                name,
                pointMesh::New(mesh())
            );
        }
    }


    if (vfield.valid())
    {
        return tmp<Field<Type>>::New
        (
            vfield().boundaryField()[patchIndex]
        );
    }
    if (sfield.valid())
    {
        return tmp<Field<Type>>::New
        (
            sfield().boundaryField()[patchIndex]
        );
    }
    if (pfield.valid())
    {
        return pfield().boundaryField()[patchIndex].patchInternalField();
    }


    FatalErrorInFunction
        << "No field '" << name << "' of type "
        << pTraits<Type>::typeName << nl << nl;

    FatalError
        << VolumeField<Type>::typeName << " Fields: "
        << flatOutput(obr.sortedNames<VolumeField<Type>>()) << nl;

    FatalError
        << SurfaceField<Type>::typeName << " Fields: "
        << flatOutput(obr.sortedNames<SurfaceField<Type>>()) << nl;

    FatalError
        << PointField<Type>::typeName << " Fields: "
        << flatOutput(obr.sortedNames<PointField<Type>>()) << nl;

    FatalError
        << exit(FatalError);

    return tmp<Field<Type>>::New();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::expressions::patchExpr::parseDriver::patchInternalField
(
    const word& name
)
{
    tmp<Field<Type>> tfield = getVariableIfAvailable<Type>(name);

    if (tfield.valid())
    {
        return tfield;
    }

    const objectRegistry& obr = this->mesh().thisDb();
    const label patchIndex = patch_.index();


    // Local, temporary storage and/or lookup values
    bool found = false;
    tmp<VolumeField<Type>> vfield;
    tmp<PointField<Type>> pfield;

    for (int checki = 0; !found && checki < 2; ++checki)
    {
        // Check 0: object context (first)
        // Check 1: regular objectRegistry
        const regIOobject* ioptr =
        (
            (checki == 0)
          ? exprDriver::cfindContextIOobject(name)
          : obr.cfindIOobject(name)
        );
        if (!ioptr) continue;

        if (!found)
        {
            vfield.cref(dynamic_cast<const VolumeField<Type>*>(ioptr));
            found = vfield.valid();
        }
        if (!found)
        {
            pfield.cref(dynamic_cast<const PointField<Type>*>(ioptr));
            found = pfield.valid();
        }
    }


    // Finally, search files if necessary (and permitted)
    if (!found && searchFiles())
    {
        const word fldType = this->getTypeOfField(name);

        if (fldType == VolumeField<Type>::typeName)
        {
            vfield = this->readAndRegister<VolumeField<Type>>(name, mesh());
        }
        else if (fldType == PointField<Type>::typeName)
        {
            pfield = this->readAndRegister<PointField<Type>>
            (
                name,
                pointMesh::New(mesh())
            );
        }
    }


    if (vfield.valid())
    {
        return vfield().boundaryField()[patchIndex].patchInternalField();
    }
    if (pfield.valid())
    {
        return pfield().boundaryField()[patchIndex].patchInternalField();
    }


    FatalErrorInFunction
        << "No field '" << name << "' of type "
        << pTraits<Type>::typeName << nl << nl;

    FatalError
        << VolumeField<Type>::typeName << " Fields: "
        << flatOutput(obr.sortedNames<VolumeField<Type>>()) << nl;

    FatalError
        << PointField<Type>::typeName << " Fields: "
        << flatOutput(obr.sortedNames<PointField<Type>>()) << nl;

    FatalError
        << exit(FatalError);

    return tmp<Field<Type>>::New();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::expressions::patchExpr::parseDriver::patchNeighbourField
(
    const word& name
)
{
    tmp<Field<Type>> tfield = getVariableIfAvailable<Type>(name);

    if (tfield.valid())
    {
        return tfield;
    }

    const objectRegistry& obr = this->mesh().thisDb();
    const label patchIndex = patch_.index();


    // Local, temporary storage and/or lookup values
    bool found = false;
    tmp<VolumeField<Type>> vfield;

    for (int checki = 0; !found && checki < 2; ++checki)
    {
        // Check 0: object context (first)
        // Check 1: regular objectRegistry
        const regIOobject* ioptr =
        (
            (checki == 0)
          ? exprDriver::cfindContextIOobject(name)
          : obr.cfindIOobject(name)
        );
        if (!ioptr) continue;

        if (!found)
        {
            vfield.cref(dynamic_cast<const VolumeField<Type>*>(ioptr));
            found = vfield.valid();
        }
    }


    // Finally, search files if necessary (and permitted)
    if (!found && searchFiles())
    {
        const word fldType = this->getTypeOfField(name);

        if (fldType == VolumeField<Type>::typeName)
        {
            vfield = this->readAndRegister<VolumeField<Type>>(name, mesh());
        }
    }


    if (vfield.valid())
    {
        return vfield().boundaryField()[patchIndex].patchNeighbourField();
    }


    FatalErrorInFunction
        << "No field '" << name << "' of type "
        << pTraits<Type>::typeName << nl << nl;

    FatalError
        << VolumeField<Type>::typeName << " Fields: "
        << flatOutput(obr.sortedNames<VolumeField<Type>>()) << nl;

    FatalError
        << exit(FatalError);

    return tmp<Field<Type>>::New();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::expressions::patchExpr::parseDriver::patchNormalField
(
    const word& name
)
{
    tmp<Field<Type>> tfield = getVariableIfAvailable<Type>(name);

    if (tfield.valid())
    {
        return tfield;
    }

    const objectRegistry& obr = this->mesh().thisDb();
    const label patchIndex = patch_.index();


    // Local, temporary storage and/or lookup values
    bool found = false;
    tmp<VolumeField<Type>> vfield;

    for (int checki = 0; !found && checki < 2; ++checki)
    {
        // Check 0: object context (first)
        // Check 1: regular objectRegistry
        const regIOobject* ioptr =
        (
            (checki == 0)
          ? exprDriver::cfindContextIOobject(name)
          : obr.cfindIOobject(name)
        );
        if (!ioptr) continue;

        if (!found)
        {
            vfield.cref(dynamic_cast<const VolumeField<Type>*>(ioptr));
            found = vfield.valid();
        }
    }


    // Finally, search files if necessary (and permitted)
    if (!found && searchFiles())
    {
        const word fldType = this->getTypeOfField(name);

        if (fldType == VolumeField<Type>::typeName)
        {
            vfield = this->readAndRegister<VolumeField<Type>>(name, mesh());
        }
    }


    if (vfield.valid())
    {
        return vfield().boundaryField()[patchIndex].snGrad();
    }


    FatalErrorInFunction
        << "No field '" << name << "' of type "
        << pTraits<Type>::typeName << nl << nl;

    FatalError
        << VolumeField<Type>::typeName << " Fields: "
        << flatOutput(obr.sortedNames<VolumeField<Type>>()) << nl;

    FatalError
        << exit(FatalError);

    return tmp<Field<Type>>::New();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::expressions::patchExpr::parseDriver::faceToPoint
(
    const Field<Type>& field
) const
{
    primitivePatchInterpolation interp(patch_.patch());

    return interp.faceToPointInterpolate(field);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::expressions::patchExpr::parseDriver::pointToFace
(
    const Field<Type>& field
) const
{
    primitivePatchInterpolation interp(patch_.patch());

    return interp.pointToFaceInterpolate(field);
}


// ************************************************************************* //
