/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "volFields.H"
#include "surfaceFields.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class GeoType>
bool Foam::ROMmodels::DMD::createAndWriteImpl() const
{
    typedef typename GeoType::value_type Type;

    const wordList modeReNames(modeNames(word("Re")));
    const wordList modeImNames(modeNames(word("Im")));

    for (const label i : modes_)
    {
        const auto* modeRePtr = mesh_.cfindObject<GeoType>(modeReNames[i]);

        if (!modeRePtr) return false;

        const auto* modeImPtr = mesh_.cfindObject<GeoType>(modeImNames[i]);

        if (!modeImPtr) return false;
    }


    forAll(times_, timei)
    {
        runTime_.setTime(times_[timei], timei);

        Info<< "\nTime = " << runTime_.timeName() << endl;

        // Calculate the eigenvalue exponent corresponding to specified time
        const scalar k = (times_[timei].value() - startTime_)/deltaT_;

        GeoType reconstructedFld
        (
            IOobject
            (
                IOobject::scopedName(fieldName_, "reconstructed"),
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh_,
            dimensioned<Type>(dimless, Zero),
            fvPatchFieldBase::zeroGradientType()
        );

        forAll(modes_, i)
        {
            const label j = modes_[i];
            const auto& modeRe = mesh_.lookupObject<GeoType>(modeReNames[j]);
            const auto& modeIm = mesh_.lookupObject<GeoType>(modeImNames[j]);

            const complex evalk(pow(evals_[i], k));

            // (K:Eq. 84)
            reconstructedFld +=
            (
                (modeRe*amps_[i].Re() - modeIm*amps_[i].Im())*evalk.Re()
              - (modeRe*amps_[i].Im() + modeIm*amps_[i].Re())*evalk.Im()
            );
        }

        reconstructedFld.correctBoundaryConditions();

        reconstructedFld.dimensions().reset(dims_);

        reconstructedFld.write();
    }

    return true;
}


// ************************************************************************* //
