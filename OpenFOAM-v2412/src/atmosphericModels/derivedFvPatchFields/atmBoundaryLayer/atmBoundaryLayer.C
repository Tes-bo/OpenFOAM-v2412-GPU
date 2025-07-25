/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2016 OpenFOAM Foundation
    Copyright (C) 2018-2022 OpenCFD Ltd.
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

#include "atmBoundaryLayer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

atmBoundaryLayer::atmBoundaryLayer(const Time& time, const polyPatch& pp)
:
    initABL_(false),
    kappa_(0.41),
    Cmu_(0.09),
    C1_(0.0),
    C2_(1.0),
    ppMin_((boundBox(pp.localPoints())).min()),
    time_(time),
    patch_(pp),
    flowDir_(nullptr),
    zDir_(nullptr),
    Uref_(nullptr),
    Zref_(nullptr),
    z0_(nullptr),
    d_(nullptr)
{}


atmBoundaryLayer::atmBoundaryLayer
(
    const Time& time,
    const polyPatch& pp,
    const dictionary& dict
)
:
    initABL_(dict.getOrDefault<bool>("initABL", true)),
    kappa_
    (
        dict.getCheckOrDefault<scalar>("kappa", 0.41, scalarMinMax::ge(SMALL))
    ),
    Cmu_(dict.getCheckOrDefault<scalar>("Cmu", 0.09, scalarMinMax::ge(SMALL))),
    C1_(dict.getOrDefault("C1", 0.0)),
    C2_(dict.getOrDefault("C2", 1.0)),
    ppMin_((boundBox(pp.localPoints())).min()),
    time_(time),
    patch_(pp),
    flowDir_(Function1<vector>::New("flowDir", dict, &time)),
    zDir_(Function1<vector>::New("zDir", dict, &time)),
    Uref_(Function1<scalar>::New("Uref", dict, &time)),
    Zref_(Function1<scalar>::New("Zref", dict, &time)),
    z0_(PatchFunction1<scalar>::New(pp, "z0", dict)),
    d_(PatchFunction1<scalar>::New(pp, "d", dict))
{}


atmBoundaryLayer::atmBoundaryLayer
(
    const atmBoundaryLayer& abl,
    const fvPatch& patch,
    const fvPatchFieldMapper& mapper
)
:
    initABL_(abl.initABL_),
    kappa_(abl.kappa_),
    Cmu_(abl.Cmu_),
    C1_(abl.C1_),
    C2_(abl.C2_),
    ppMin_(abl.ppMin_),
    time_(abl.time_),
    patch_(patch.patch()),
    flowDir_(abl.flowDir_.clone()),
    zDir_(abl.zDir_.clone()),
    Uref_(abl.Uref_.clone()),
    Zref_(abl.Zref_.clone()),
    z0_(abl.z0_.clone(patch_)),
    d_(abl.d_.clone(patch_))
{}


atmBoundaryLayer::atmBoundaryLayer(const atmBoundaryLayer& abl)
:
    initABL_(abl.initABL_),
    kappa_(abl.kappa_),
    Cmu_(abl.Cmu_),
    C1_(abl.C1_),
    C2_(abl.C2_),
    ppMin_(abl.ppMin_),
    time_(abl.time_),
    patch_(abl.patch_),
    flowDir_(abl.flowDir_.clone()),
    zDir_(abl.zDir_.clone()),
    Uref_(abl.Uref_.clone()),
    Zref_(abl.Zref_.clone()),
    z0_(abl.z0_.clone(patch_)),
    d_(abl.d_.clone(patch_))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

vector atmBoundaryLayer::flowDir() const
{
    const scalar t = time_.timeOutputValue();
    const vector dir(flowDir_->value(t));
    const scalar magDir = mag(dir);

    if (magDir < SMALL)
    {
        FatalErrorInFunction
            << "magnitude of " << flowDir_->name() << " = " << magDir
            << " vector must be greater than zero"
            << abort(FatalError);
    }

    return dir/magDir;
}


vector atmBoundaryLayer::zDir() const
{
    const scalar t = time_.timeOutputValue();
    const vector dir(zDir_->value(t));
    const scalar magDir = mag(dir);

    if (magDir < SMALL)
    {
        FatalErrorInFunction
            << "magnitude of " << zDir_->name() << " = " << magDir
            << " vector must be greater than zero"
            << abort(FatalError);
    }

    return dir/magDir;
}


tmp<scalarField> atmBoundaryLayer::Ustar(const scalarField& z0) const
{
    const scalar t = time_.timeOutputValue();
    const scalar Uref = Uref_->value(t);
    const scalar Zref = Zref_->value(t);

    if (Zref < 0)
    {
        FatalErrorInFunction
            << "Negative entry in " << Zref_->name() << " = " << Zref
            << abort(FatalError);
    }

    // (derived from RH:Eq. 6, HW:Eq. 5)
    return kappa_*Uref/log((Zref + z0)/z0);
}


void atmBoundaryLayer::autoMap(const fvPatchFieldMapper& mapper)
{
    if (z0_)
    {
        z0_->autoMap(mapper);
    }
    if (d_)
    {
        d_->autoMap(mapper);
    }
}


void atmBoundaryLayer::rmap
(
    const atmBoundaryLayer& abl,
    const labelList& addr
)
{
    if (z0_)
    {
        z0_->rmap(abl.z0_(), addr);
    }
    if (d_)
    {
        d_->rmap(abl.d_(), addr);
    }
}


tmp<vectorField> atmBoundaryLayer::U(const vectorField& pCf) const
{
    const scalar t = time_.timeOutputValue();
    const scalarField d(d_->value(t));
    const scalarField z0(max(z0_->value(t), ROOTVSMALL));
    const scalar groundMin = zDir() & ppMin_;

    // (YGCJ:Table 1, RH:Eq. 6, HW:Eq. 5)
    scalarField Un
    (
        (Ustar(z0)/kappa_)*log(((zDir() & pCf) - groundMin - d + z0)/z0)
    );

    return flowDir()*Un;
}


tmp<scalarField> atmBoundaryLayer::k(const vectorField& pCf) const
{
    const scalar t = time_.timeOutputValue();
    const scalarField d(d_->value(t));
    const scalarField z0(max(z0_->value(t), ROOTVSMALL));
    const scalar groundMin = zDir() & ppMin_;

    // (YGCJ:Eq. 21; RH:Eq. 7, HW:Eq. 6 when C1=0 and C2=1)
    return
        sqr(Ustar(z0))/sqrt(Cmu_)
       *sqrt(C1_*log(((zDir() & pCf) - groundMin - d + z0)/z0) + C2_);
}


tmp<scalarField> atmBoundaryLayer::epsilon(const vectorField& pCf) const
{
    const scalar t = time_.timeOutputValue();
    const scalarField d(d_->value(t));
    const scalarField z0(max(z0_->value(t), ROOTVSMALL));
    const scalar groundMin = zDir() & ppMin_;

    // (YGCJ:Eq. 22; RH:Eq. 8, HW:Eq. 7 when C1=0 and C2=1)
    return
        pow3(Ustar(z0))/(kappa_*((zDir() & pCf) - groundMin - d + z0))
       *sqrt(C1_*log(((zDir() & pCf) - groundMin - d + z0)/z0) + C2_);
}


tmp<scalarField> atmBoundaryLayer::omega(const vectorField& pCf) const
{
    const scalar t = time_.timeOutputValue();
    const scalarField d(d_->value(t));
    const scalarField z0(max(z0_->value(t), ROOTVSMALL));
    const scalar groundMin = zDir() & ppMin_;

    // (YGJ:Eq. 13)
    return Ustar(z0)/(kappa_*sqrt(Cmu_)*((zDir() & pCf) - groundMin - d + z0));
}


void atmBoundaryLayer::write(Ostream& os) const
{
    os.writeEntry("initABL", initABL_);
    os.writeEntry("kappa", kappa_);
    os.writeEntry("Cmu", Cmu_);
    os.writeEntry("C1", C1_);
    os.writeEntry("C2", C2_);
    if (flowDir_)
    {
        flowDir_->writeData(os);
    }
    if (zDir_)
    {
        zDir_->writeData(os);
    }
    if (Uref_)
    {
        Uref_->writeData(os);
    }
    if (Zref_)
    {
        Zref_->writeData(os);
    }
    if (z0_)
    {
        z0_->writeData(os) ;
    }
    if (d_)
    {
        d_->writeData(os);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
