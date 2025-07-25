/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenCFD Ltd.
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

#include "InterfaceForce.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::InterfaceForce<CloudType>::InterfaceForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    ParticleForce<CloudType>(owner, mesh, dict, typeName, true),
    alphaName_(this->coeffs().lookup("alpha")),
    C_(this->coeffs().getScalar("C")),
    gradInterForceInterpPtr_(nullptr)
{}


template<class CloudType>
Foam::InterfaceForce<CloudType>::InterfaceForce(const InterfaceForce& pf)
:
    ParticleForce<CloudType>(pf),
    alphaName_(pf.alphaName_),
    C_(pf.C_),
    gradInterForceInterpPtr_(pf.gradInterForceInterpPtr_)
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::InterfaceForce<CloudType>::~InterfaceForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::InterfaceForce<CloudType>::cacheFields(const bool store)
{
    static word resultName("gradAlpha");

    volVectorField* resultPtr =
        this->mesh().template getObjectPtr<volVectorField>(resultName);

    if (store)
    {
        if (!resultPtr)
        {
            const volScalarField& alpha = this->mesh().template
                lookupObject<volScalarField>(alphaName_);

            resultPtr = new volVectorField
            (
                resultName,
                fvc::grad(alpha*(1-alpha))
            );

            resultPtr->store();
        }

        const volVectorField& gradInterForce = *resultPtr;

        gradInterForceInterpPtr_.reset
        (
            interpolation<vector>::New
            (
                this->owner().solution().interpolationSchemes(),
                gradInterForce
            ).ptr()
        );
    }
    else
    {
        gradInterForceInterpPtr_.clear();

        if (resultPtr)
        {
            resultPtr->checkOut();
        }
    }
}


template<class CloudType>
Foam::forceSuSp Foam::InterfaceForce<CloudType>::calcNonCoupled
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    forceSuSp value(Zero);

    const interpolation<vector>& interp = gradInterForceInterpPtr_();

    value.Su() =
       C_
      *mass
      *interp.interpolate(p.coordinates(), p.currentTetIndices());

    return value;
}


// ************************************************************************* //
