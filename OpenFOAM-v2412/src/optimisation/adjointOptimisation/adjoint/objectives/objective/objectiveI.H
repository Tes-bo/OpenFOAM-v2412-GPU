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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::word& Foam::objective::objectiveName() const
{
    return objectiveName_;
}


inline bool Foam::objective::shouldWrite() const
{
    return shouldWrite_;
}


inline void Foam::objective::setWrite(const bool shouldWrite)
{
    shouldWrite_ = shouldWrite;
}


inline const Foam::volScalarField& Foam::objective::dJdb() const
{
    return *dJdbPtr_;
}


inline const Foam::scalarField& Foam::objective::dJdbField() const
{
    return *dJdbFieldPtr_;
}


inline const Foam::fvPatchVectorField&
Foam::objective::boundarydJdb(const label patchI) const
{
    return bdJdbPtr_()[patchI];
}


inline const Foam::fvPatchVectorField&
Foam::objective::dSdbMultiplier(const label patchI) const
{
    return bdSdbMultPtr_()[patchI];
}


inline const Foam::fvPatchVectorField&
Foam::objective::dndbMultiplier(const label patchI) const
{
    return bdndbMultPtr_()[patchI];
}


inline const Foam::fvPatchVectorField&
Foam::objective::dxdbMultiplier(const label patchI) const
{
    return bdxdbMultPtr_()[patchI];
}


inline const Foam::fvPatchVectorField&
Foam::objective::dxdbDirectMultiplier(const label patchI) const
{
    return bdxdbDirectMultPtr_()[patchI];
}


inline const Foam::vectorField& Foam::objective::boundaryEdgeMultiplier
(
    const label patchI,
    const label edgeI
) const
{
    if (!bdxdbDirectMultPtr_)
    {
        FatalErrorInFunction
            << "Unallocated boundaryEdgeMultiplier field"
            << exit(FatalError);
    }
    return bEdgeContribution_()[patchI][edgeI];
}


inline const Foam::boundaryVectorField& Foam::objective::boundarydJdb() const
{
    return *bdJdbPtr_;
}


inline const Foam::boundaryVectorField& Foam::objective::dSdbMultiplier() const
{
    return *bdSdbMultPtr_;
}


inline const Foam::boundaryVectorField& Foam::objective::dndbMultiplier() const
{
    return *bdndbMultPtr_;
}


inline const Foam::boundaryVectorField& Foam::objective::dxdbMultiplier() const
{
    return *bdxdbMultPtr_;
}


inline const Foam::boundaryVectorField&
Foam::objective::dxdbDirectMultiplier() const
{
    return *bdxdbDirectMultPtr_;
}


inline const Foam::vectorField3& Foam::objective::boundaryEdgeMultiplier() const
{
    if (!bdxdbDirectMultPtr_)
    {
        FatalErrorInFunction
            << "Unallocated boundaryEdgeMultiplier field"
            << endl << endl
            << exit(FatalError);
    }
    return *bEdgeContribution_;
}


inline const Foam::volScalarField& Foam::objective::divDxDbMultiplier() const
{
    return *divDxDbMultPtr_;
}


inline const Foam::volTensorField& Foam::objective::gradDxDbMultiplier() const
{
    return *gradDxDbMultPtr_;
}


inline bool Foam::objective::hasdJdb() const noexcept
{
    return bool(dJdbPtr_);
}


inline bool Foam::objective::hasdJdbField() const noexcept
{
    return bool(dJdbFieldPtr_);
}


inline bool Foam::objective::hasBoundarydJdb() const noexcept
{
    return bool(bdJdbPtr_);
}


inline bool Foam::objective::hasdSdbMult() const noexcept
{
    return bool(bdSdbMultPtr_);
}


inline bool Foam::objective::hasdndbMult() const noexcept
{
    return bool(bdndbMultPtr_);
}


inline bool Foam::objective::hasdxdbMult() const noexcept
{
    return bool(bdxdbMultPtr_);
}


inline bool Foam::objective::hasdxdbDirectMult() const noexcept
{
    return bool(bdxdbDirectMultPtr_);
}


inline bool Foam::objective::hasBoundaryEdgeContribution() const noexcept
{
    return bool(bEdgeContribution_);
}


inline bool Foam::objective::hasDivDxDbMult() const noexcept
{
    return bool(divDxDbMultPtr_);
}


inline bool Foam::objective::hasGradDxDbMult() const noexcept
{
    return bool(gradDxDbMultPtr_);
}


inline bool Foam::objective::hasIntegrationStartTime() const noexcept
{
    return bool(integrationStartTimePtr_);
}


inline bool Foam::objective::hasIntegrationEndTime() const noexcept
{
    return bool(integrationEndTimePtr_);
}


inline void Foam::objective::setComputed(const bool isComputed) noexcept
{
    computed_ = isComputed;
}


// ************************************************************************* //
