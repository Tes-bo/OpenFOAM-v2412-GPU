/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "kinematicThinFilm.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(kinematicThinFilm, 0);
addToRunTimeSelectionTable(liquidFilmBase, kinematicThinFilm, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

kinematicThinFilm::kinematicThinFilm
(
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    liquidFilmModel(modelType, mesh, dict)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void kinematicThinFilm::preEvolveRegion()
{
    rhoSp_.storePrevIter();
    USp_.storePrevIter();
    pnSp_.storePrevIter();

    // Update mass exchange sources
    liquidFilmModel::preEvolveRegion();

    // gas pressure map from primary region
    ppf_ = pg();
}


void kinematicThinFilm::evolveRegion()
{
    DebugInFunction << endl;

    const areaVectorField& ns = regionMesh().faceAreaNormals();

    const areaVectorField gs(g_ - ns*(ns & g_));

    phi2s_ = fac::interpolate(h_)*phif_;

    for (int oCorr=1; oCorr<=nOuterCorr_; ++oCorr)
    {
        pf_.storePrevIter();

        faVectorMatrix UsEqn
        (
            fam::ddt(h_, Uf_)
          + fam::div(phi2s_, Uf_)
          ==
            gs*h_
          + turbulence_->Su(Uf_)
          + faOptions()(h_, Uf_, sqr(dimVelocity))
          + forces_.correct(Uf_)
          + USp_
        );

        UsEqn.relax();

        faOptions().constrain(UsEqn);

        if (momentumPredictor_)
        {
            solve(UsEqn == -fac::grad(pf_*h_)/rho_ + pf_*fac::grad(h_)/rho_);
        }

        for (int corr=1; corr<=nCorr_; ++corr)
        {
            areaScalarField UsA(UsEqn.A());

            Uf_ = UsEqn.H()/UsA;
            Uf_.correctBoundaryConditions();
            faOptions().correct(Uf_);

            phif_ =
                (fac::interpolate(Uf_) & regionMesh().Le())
                - fac::interpolate(1.0/(rho_*UsA))
                * fac::lnGrad(pf_*h_)*regionMesh().magLe()
                + fac::interpolate(pf_/(rho_*UsA))
                * fac::lnGrad(h_)*regionMesh().magLe();

            for (int nFilm=1; nFilm<=nFilmCorr_; ++nFilm)
            {
                faScalarMatrix hEqn
                (
                    fam::ddt(h_)
                  + fam::div(phif_, h_)
                 ==
                    faOptions()(rho_, h_, dimVelocity)
                  + rhoSp_
                );

                hEqn.relax();
                faOptions().constrain(hEqn);
                hEqn.solve();
                faOptions().correct(h_);

                if (nFilm == nFilmCorr_)
                {
                    phi2s_ = hEqn.flux();
                }
            }

            // Bound h_
            h_ = max(h_, h0_);

            pf_ = rho_*gn_*h_ - sigma_*fac::laplacian(h_) + pnSp_ + ppf_;
            pf_.correctBoundaryConditions();
            pf_.relax();

            Uf_ -= (1.0/(rho_*UsA))*fac::grad(pf_*h_)
                - (pf_/(rho_*UsA))*fac::grad(h_);
            Uf_.correctBoundaryConditions();
            faOptions().correct(Uf_);
        }
    }

    Info<< "Film h min/max   = " << gMinMax(h_) << nl
        << "Film mag(U) min/max   = " << gMinMaxMag(Uf_) << endl;
}


void kinematicThinFilm::postEvolveRegion()
{
    // Reset sources
    liquidFilmModel::postEvolveRegion();

    // Correct thermo
    correctThermoFields();

    // Correct turbulence
    turbulence_->correct();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
