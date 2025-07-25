/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2023 PCOpt/NTUA
    Copyright (C) 2013-2023 FOSS GP
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const word& RASModelVariables::TMVar1BaseName() const
{
    return TMVar1BaseName_;
}


inline const word& RASModelVariables::TMVar2BaseName() const
{
    return TMVar2BaseName_;
}


inline const word& RASModelVariables::nutBaseName() const
{
    return nutBaseName_;
}


inline bool RASModelVariables::hasTMVar1() const
{
    return bool(TMVar1Ptr_);
}


inline bool RASModelVariables::hasTMVar2() const
{
    return bool(TMVar2Ptr_);
}


inline bool RASModelVariables::hasNut() const
{
    return bool(nutPtr_);
}


inline bool RASModelVariables::hasDist() const
{
    return bool(distPtr_);
}


inline const volScalarField& RASModelVariables::TMVar1() const
{
    if (solverControl_.useAveragedFields())
    {
        return TMVar1MeanPtr_.cref();
    }

    return TMVar1Ptr_.cref();
}


inline volScalarField& RASModelVariables::TMVar1()
{
    if (solverControl_.useAveragedFields())
    {
        return TMVar1MeanPtr_.ref();
    }

    return TMVar1Ptr_.ref();
}


inline const volScalarField& RASModelVariables::TMVar2() const
{
    if (solverControl_.useAveragedFields())
    {
        return TMVar2MeanPtr_.cref();
    }

    return TMVar2Ptr_.cref();
}

inline volScalarField& RASModelVariables::TMVar2()
{
    if (solverControl_.useAveragedFields())
    {
        return TMVar2MeanPtr_.ref();
    }

    return TMVar2Ptr_.ref();
}

inline const volScalarField& RASModelVariables::nutRef() const
{
    if (solverControl_.useAveragedFields() && nutMeanPtr_)
    {
        return nutMeanPtr_.cref();
    }

    return nutPtr_.cref();
}


inline volScalarField& RASModelVariables::nutRef()
{
    if (solverControl_.useAveragedFields() && nutMeanPtr_)
    {
        return nutMeanPtr_.ref();
    }

    return nutPtr_.ref();
}


inline tmp<volScalarField> RASModelVariables::nut() const
{
    if (hasNut())
    {
        return tmp<volScalarField>(nutRef());
    }

    return volScalarField::New
    (
        "dummylaminarNut",
        IOobject::NO_REGISTER,
        mesh_,
        dimensionedScalar(dimViscosity, Zero)
    );
}


inline const volScalarField& RASModelVariables::d() const
{
    return distPtr_.cref();
}


inline volScalarField& RASModelVariables::d()
{
    return distPtr_.ref();
}


inline tmp<scalarField> RASModelVariables::TMVar1(const label patchi) const
{
    if (hasTMVar1())
    {
        return TMVar1().boundaryField()[patchi];
    }

    return tmp<scalarField>::New(mesh_.boundary()[patchi].size(), Zero);
}


inline tmp<scalarField> RASModelVariables::TMVar2(const label patchi) const
{
    if (hasTMVar2())
    {
        return TMVar2().boundaryField()[patchi];
    }

    return tmp<scalarField>::New(mesh_.boundary()[patchi].size(), Zero);
}


inline tmp<scalarField> RASModelVariables::nut(const label patchi) const
{
    if (hasNut())
    {
        return nutRef().boundaryField()[patchi];
    }

    return tmp<scalarField>::New(mesh_.boundary()[patchi].size(), Zero);
}


inline tmp<fvPatchScalarField>
RASModelVariables::nutPatchField(const label patchi) const
{
    if (hasNut())
    {
        return nutRef().boundaryField()[patchi];
    }

    // Using dummy internalField
    return
        tmp<fvPatchScalarField>::New(mesh_.boundary()[patchi], mesh_.V(), Zero);
}


inline const volScalarField& RASModelVariables::TMVar1Inst() const
{
    return TMVar1Ptr_.cref();
}


inline volScalarField& RASModelVariables::TMVar1Inst()
{
    return TMVar1Ptr_.ref();
}


inline const volScalarField& RASModelVariables::TMVar2Inst() const
{
    return TMVar2Ptr_.cref();
}


inline volScalarField& RASModelVariables::TMVar2Inst()
{
    return TMVar2Ptr_.ref();
}


inline const volScalarField& RASModelVariables::nutRefInst() const
{
    return nutPtr_.cref();
}


inline volScalarField& RASModelVariables::nutRefInst()
{
    return nutPtr_.ref();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
