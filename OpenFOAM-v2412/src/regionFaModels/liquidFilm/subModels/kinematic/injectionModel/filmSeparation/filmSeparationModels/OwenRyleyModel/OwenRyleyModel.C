/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2021-2024 OpenCFD Ltd.
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

#include "OwenRyleyModel.H"
#include "processorFaPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace filmSeparationModels
{
    defineTypeNameAndDebug(OwenRyleyModel, 0);
    addToRunTimeSelectionTable(filmSeparationModel, OwenRyleyModel, dictionary);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

tmp<areaScalarField> OwenRyleyModel::calcInvR1
(
    const areaVectorField& U
) const
{
    const dimensionedScalar smallU(dimVelocity, ROOTVSMALL);
    const areaVectorField UHat(U/(mag(U) + smallU));

    auto tinvR1 = areaScalarField::New
    (
        "invR1",
        IOobjectOption::NO_REGISTER,
        (UHat & (UHat & -gradNHat_))
    );
    scalarField& invR1 = tinvR1.ref().primitiveFieldRef();

    // Apply defined patch radii
    if (definedPatchRadii_ > 0)
    {
        invR1 = 1.0/max(1e-6, definedPatchRadii_);
    }

    // Filter out large radii
    for (auto& x : invR1)
    {
        if (mag(x) < 1e-6)
        {
            x = -1;
        }
    }

    return tinvR1;
}


tmp<scalarField> OwenRyleyModel::calcCosAngle
(
    const edgeScalarField& phi,
    const scalarField& invR1
) const
{
    const areaVectorField& U = film().Uf();
    const dimensionedScalar smallU(dimVelocity, ROOTVSMALL);
    const areaVectorField UHat(U/(mag(U) + smallU));

    const vector gHat(normalised(film().g().value()));

    const faMesh& mesh = film().regionMesh();
    const labelUList& own = mesh.edgeOwner();
    const labelUList& nbr = mesh.edgeNeighbour();

    scalarField phiMax(mesh.nFaces(), -GREAT);
    scalarField cosAngle(UHat.size(), Zero);

    // Internal edges
    forAll(nbr, edgei)
    {
        const label cellO = own[edgei];
        const label cellN = nbr[edgei];
        const scalar phiEdge = phi[edgei];

        if (phiEdge > phiMax[cellO])
        {
            phiMax[cellO] = phiEdge;
            cosAngle[cellO] = -gHat & UHat[cellN];
        }
        if (-phiEdge > phiMax[cellN])
        {
            phiMax[cellN] = -phiEdge;
            cosAngle[cellN] = -gHat & UHat[cellO];
        }
    }

    // Processor edges
    for (const auto& phip : phi.boundaryField())
    {
        if (isA<processorFaPatch>(phip.patch()))
        {
            const auto& edgeFaces = phip.patch().edgeFaces();
            const auto& UHatp = UHat.boundaryField()[phip.patch().index()];
            forAll(phip, edgei)
            {
                const label cellO = edgeFaces[edgei];
                if (phip[edgei] > phiMax[cellO])
                {
                    phiMax[cellO] = phip[edgei];
                    cosAngle[cellO] = -gHat & UHatp[edgei];
                }
            }
        }
    }

    cosAngle *= pos(invR1);


    if (debug && mesh.time().writeTime())
    {
        {
            areaScalarField areaCosAngle
            (
                mesh.newIOobject("cosAngle"),
                mesh,
                dimensionedScalar(dimless, Zero)
            );
            areaCosAngle.primitiveFieldRef() = cosAngle;
            areaCosAngle.correctBoundaryConditions();
            areaCosAngle.write();
        }

        {
            areaScalarField areaInvR1
            (
                mesh.newIOobject("InvR1"),
                mesh,
                dimensionedScalar(inv(dimLength), Zero)
            );
            areaInvR1.primitiveFieldRef() = invR1;
            areaInvR1.write();
        }
    }

    return clamp(cosAngle, scalarMinMax(-1, 1));
}


tmp<scalarField> OwenRyleyModel::netForce() const
{
    const faMesh& mesh = film().regionMesh();

    const areaScalarField& h = film().h();
    const areaVectorField& U = film().Uf();
    const edgeScalarField& phi = film().phi2s();
    const areaScalarField& rho = film().rho();
    const areaScalarField& sigma = film().sigma();

    const scalarField magSqrU(magSqr(U));
    const scalarField invR1(calcInvR1(U));
    const scalarField cosAngle(calcCosAngle(phi, invR1));

    const scalar magG = mag(film().g().value());

    // Initialize the net-force magnitude
    auto tFnet = tmp<scalarField>::New(mesh.nFaces(), Zero);
    auto& Fnet = tFnet.ref();

    forAll(Fnet, i)
    {
        if ((invR1[i] > minInvR1()) && (h[i]*invR1[i] > minHbyR1()))
        {
            const scalar R1 = 1.0/(invR1[i] + ROOTVSMALL);
            const scalar R2 = R1 + h[i];

            // Inertial force (OR:Eq. 11; '2' is an exponent in the Eq.)
            const scalar Fi = -4.0/3.0*h[i]*rho[i]*magSqrU[i]*invR1[i];

            // Body force (OR:Eq. 11)
            const scalar Fb =
                -0.5*rho[i]*magG*invR1[i]*(sqr(R1) - sqr(R2))*cosAngle[i];

            // Surface force (OR:Eq. 11)
            const scalar Fs = sigma[i]/R2;

            Fnet[i] = Fi + Fb + Fs;
        }
    }

    if (debug && mesh.time().writeTime())
    {
        {
            areaScalarField areaFnet
            (
                mesh.newIOobject("Fnet"),
                mesh,
                dimensionedScalar(dimForce, Zero)
            );
            areaFnet.primitiveFieldRef() = Fnet;
            areaFnet.write();
        }
    }

    return tFnet;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

OwenRyleyModel::OwenRyleyModel
(
    const regionModels::areaSurfaceFilmModels::liquidFilmBase& film,
    const dictionary& dict
)
:
    filmSeparationModel(film, dict),
    fThreshold_(dict.getOrDefault<scalar>("fThreshold", 1e-8)),
    definedPatchRadii_(dict.getOrDefault<scalar>("definedPatchRadii", 0)),
    minHbyR1_(dict.getOrDefault<scalar>("deltaByR1Min", 0)),
    minInvR1_(dict.getOrDefault<scalar>("minInvR1", 5)),
    gradNHat_(fac::grad(film.regionMesh().faceAreaNormals()))
{
    if (mag(film.g().value()) < ROOTVSMALL)
    {
        FatalErrorInFunction
            << "Gravitational acceleration must be non-zero."
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField> OwenRyleyModel::separatedMassRatio() const
{
    const faMesh& mesh = film().regionMesh();

    tmp<scalarField> tFnet = netForce();
    const auto& Fnet = tFnet.cref();

    // Initialize the mass ratio of film separation
    auto tseparated = tmp<scalarField>::New(mesh.nFaces(), Zero);
    auto& separated = tseparated.ref();

    forAll(Fnet, i)
    {
        if ((Fnet[i] + fThreshold_) < 0)
        {
            separated[i] = 1;
        }
    }

    return tseparated;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace filmSeparationModels
} // End namespace Foam


// ************************************************************************* //

