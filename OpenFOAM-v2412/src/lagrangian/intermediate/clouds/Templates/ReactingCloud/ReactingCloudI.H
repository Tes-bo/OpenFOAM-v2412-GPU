/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2023 OpenCFD Ltd.
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

template<class CloudType>
inline const Foam::ReactingCloud<CloudType>&
Foam::ReactingCloud<CloudType>::cloudCopy() const
{
    return *cloudCopyPtr_;
}


template<class CloudType>
inline const typename CloudType::particleType::constantProperties&
Foam::ReactingCloud<CloudType>::constProps() const
{
    return constProps_;
}


template<class CloudType>
inline typename CloudType::particleType::constantProperties&
Foam::ReactingCloud<CloudType>::constProps()
{
    return constProps_;
}


template<class CloudType>
inline const Foam::CompositionModel<Foam::ReactingCloud<CloudType>>&
Foam::ReactingCloud<CloudType>::composition() const
{
    return *compositionModel_;
}


template<class CloudType>
inline const Foam::PhaseChangeModel<Foam::ReactingCloud<CloudType>>&
Foam::ReactingCloud<CloudType>::phaseChange() const
{
    return *phaseChangeModel_;
}


template<class CloudType>
inline Foam::PhaseChangeModel<Foam::ReactingCloud<CloudType>>&
Foam::ReactingCloud<CloudType>::phaseChange()
{
    return *phaseChangeModel_;
}


template<class CloudType>
inline void Foam::ReactingCloud<CloudType>::transferToCarrier
(
    const parcelType& p,
    const typename parcelType::trackingData& td
)
{
    const auto& comp = this->composition();

    const label celli = p.cell();

    const scalar m = p.nParticle()*p.mass();

    this->rhokTrans()[celli] += m;

    this->UTrans()[celli] += m*p.U();

    const scalar pc = td.pc();
    const scalar T = p.T();
    const auto& Y = p.Y();

    forAll(Y, i)
    {
        const scalar dm = m*p.Y[i];
        const label gid = comp.localToCarrierId(0, i);
        const scalar hs = comp.carrier().Hs(gid, pc, T);

        this->rhoTrans(gid)[celli] += dm;
        this->hsTrans()[celli] += dm*hs;
    }
}


template<class CloudType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ReactingCloud<CloudType>::rhoTrans(const label i)
{
    return rhoTrans_[i];
}


template<class CloudType>
inline
const Foam::PtrList<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>&
Foam::ReactingCloud<CloudType>::rhoTrans() const
{
    return rhoTrans_;
}


template<class CloudType>
inline Foam::PtrList<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>&
Foam::ReactingCloud<CloudType>::rhoTrans()
{
    return rhoTrans_;
}


template<class CloudType>
inline Foam::tmp<Foam::fvScalarMatrix> Foam::ReactingCloud<CloudType>::SYi
(
    const label i,
    volScalarField& Yi
) const
{
    if (this->solution().coupled())
    {
        if (this->solution().semiImplicit("Yi"))
        {
            auto trhoTrans = volScalarField::New
            (
                IOobject::scopedName(this->name(), "rhoTrans"),
                IOobject::NO_REGISTER,
                this->mesh(),
                dimensionedScalar(dimMass/dimTime/dimVolume, Zero)
            );
            auto& sourceField = trhoTrans.ref();

            sourceField.primitiveFieldRef() =
                rhoTrans_[i]/(this->db().time().deltaTValue()*this->mesh().V());

            const dimensionedScalar YiSMALL("YiSMALL", dimless, SMALL);

            return
                fvm::Sp(neg(sourceField)*sourceField/(Yi + YiSMALL), Yi)
              + pos0(sourceField)*sourceField;
        }
        else
        {
            auto tfvm = tmp<fvScalarMatrix>::New(Yi, dimMass/dimTime);
            auto& fvm = tfvm.ref();

            fvm.source() = -rhoTrans_[i]/this->db().time().deltaTValue();

            return tfvm;
        }
    }

    return tmp<fvScalarMatrix>::New(Yi, dimMass/dimTime);
}


template<class CloudType>
inline Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::ReactingCloud<CloudType>::Srho(const label i) const
{
    auto tRhoi = volScalarField::Internal::New
    (
        IOobject::scopedName(this->name(), "rhoTrans"),
        IOobject::NO_REGISTER,
        this->mesh(),
        dimensionedScalar
        (
            rhoTrans_[0].dimensions()/dimTime/dimVolume, Zero
        )
    );
    scalarField& rhoi = tRhoi.ref();

    if (this->solution().coupled())
    {
        rhoi = rhoTrans_[i]/(this->db().time().deltaTValue()*this->mesh().V());
    }

    return tRhoi;
}


template<class CloudType>
inline Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::ReactingCloud<CloudType>::Srho() const
{
    auto trhoTrans = volScalarField::Internal::New
    (
        IOobject::scopedName(this->name(), "rhoTrans"),
        IOobject::NO_REGISTER,
        this->mesh(),
        dimensionedScalar
        (
            rhoTrans_[0].dimensions()/dimTime/dimVolume, Zero
        )
    );
    scalarField& sourceField = trhoTrans.ref();

    if (this->solution().coupled())
    {
        forAll(rhoTrans_, i)
        {
            sourceField += rhoTrans_[i];
        }

        sourceField /= this->db().time().deltaTValue()*this->mesh().V();
    }

    return trhoTrans;
}


template<class CloudType>
inline Foam::tmp<Foam::fvScalarMatrix>
Foam::ReactingCloud<CloudType>::Srho(volScalarField& rho) const
{
    if (this->solution().coupled())
    {
        auto trhoTrans = volScalarField::New
        (
            IOobject::scopedName(this->name(), "rhoTrans"),
            IOobject::NO_REGISTER,
            this->mesh(),
            dimensionedScalar(dimMass/dimTime/dimVolume, Zero)
        );
        scalarField& sourceField = trhoTrans.ref();

        if (this->solution().semiImplicit("rho"))
        {
            forAll(rhoTrans_, i)
            {
                sourceField += rhoTrans_[i];
            }
            sourceField /= this->db().time().deltaTValue()*this->mesh().V();

            return fvm::SuSp(trhoTrans()/rho, rho);
        }
        else
        {
            auto tfvm = tmp<fvScalarMatrix>::New(rho, dimMass/dimTime);
            auto& fvm = tfvm.ref();

            forAll(rhoTrans_, i)
            {
                sourceField += rhoTrans_[i];
            }

            fvm.source() = -trhoTrans()/this->db().time().deltaT();

            return tfvm;
        }
    }

    return tmp<fvScalarMatrix>::New(rho, dimMass/dimTime);
}


// ************************************************************************* //
