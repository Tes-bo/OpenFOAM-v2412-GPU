/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2023 OpenCFD Ltd.
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

#include "contactAngleForce.H"
#include "addToRunTimeSelectionTable.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(contactAngleForce, 0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

contactAngleForce::contactAngleForce
(
    const word& typeName,
    liquidFilmBase& film,
    const dictionary& dict
)
:
    force(typeName, film, dict),
    Ccf_(coeffDict_.get<scalar>("Ccf")),
    hCrit_(coeffDict_.getOrDefault<scalar>("hCrit", GREAT))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<faVectorMatrix> contactAngleForce::correct(areaVectorField& U)
{
    auto tforce = tmp<areaVectorField>::New
    (
        IOobject
        (
            IOobject::scopedName(typeName, "contactForce"),
            film().regionMesh().time().timeName(),
            film().regionMesh().thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        film().regionMesh(),
        dimensionedVector(dimForce/dimDensity/dimArea, Zero)
    );
    vectorField& force = tforce.ref().primitiveFieldRef();

    const labelUList& own = film().regionMesh().owner();
    const labelUList& nei = film().regionMesh().neighbour();

    const DimensionedField<scalar, areaMesh>& magSf = film().regionMesh().S();
    const scalarField& magSff = magSf.field();

    const edgeScalarField& invDx = film().regionMesh().deltaCoeffs();

    const areaScalarField& sigma = film().sigma();
    const areaScalarField& mu = film().mu();
    const areaScalarField& rhof = film().rho();
    const areaVectorField& Uf = film().Uf();
    const areaScalarField& hf = film().h();

    tmp<areaScalarField> talpha = film().alpha();
    const areaScalarField& alpha = talpha();

    tmp<areaScalarField> ttheta = theta();
    const areaScalarField& theta = ttheta();

    tmp<areaVectorField> tgradAlpha = fac::grad(alpha);
    const areaVectorField& gradAlpha = tgradAlpha();

    forAll(nei, edgei)
    {
        const label faceO = own[edgei];
        const label faceN = nei[edgei];

        label facei = -1;
        if ((alpha[faceO] > 0.5) && (alpha[faceN] < 0.5))
        {
            facei = faceO;
        }
        else if ((alpha[faceO] < 0.5) && (alpha[faceN] > 0.5))
        {
            facei = faceN;
        }

        if (facei != -1)
        {
            const vector n(normalised(gradAlpha[facei]));
            const scalar cosTheta = cos(degToRad(theta[facei]));

            // (MHDX:Eq. 13)
            force[facei] +=
                Ccf_*n*sigma[facei]*(1 - cosTheta)
               /invDx[edgei]/rhof[facei]/magSff[facei];

            // (NDPC:Eq. 11)
            if (hf[facei] > hCrit_)
            {
                force[facei] -= mu[facei]*Uf[facei]/hCrit_;
            }
        }
    }

    for (const faPatchScalarField& sigmaBf : sigma.boundaryField())
    {
        const faPatch& p = sigmaBf.patch();

        if (!p.coupled())
        {
            const labelUList& faces = p.edgeFaces();

            forAll(sigmaBf, edgei)
            {
                const label face0 = faces[edgei];
                force[face0] = Zero;
            }
        }
    }

    if (film().regionMesh().time().writeTime())
    {
        tforce().write();
        gradAlpha.write();
    }

    auto tfvm = tmp<faVectorMatrix>::New(U, dimForce/dimDensity);

    tfvm.ref() += tforce;

    return tfvm;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
