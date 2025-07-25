/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 OpenCFD Ltd.
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

#include "FriedrichModel.H"
#include "processorFaPatch.H"
#include "unitConversion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace filmSeparationModels
{
    defineTypeNameAndDebug(FriedrichModel, 0);
    addToRunTimeSelectionTable(filmSeparationModel, FriedrichModel, dictionary);


const Foam::Enum
<
    FriedrichModel::separationType
>
FriedrichModel::separationTypeNames
({
    { separationType::FULL, "full" },
    { separationType::PARTIAL , "partial" },
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bitSet FriedrichModel::calcCornerEdges() const
{
    bitSet cornerEdges(mesh().nEdges(), false);

    const areaVectorField& faceCentres = mesh().areaCentres();
    const areaVectorField& faceNormals = mesh().faceAreaNormals();

    const labelUList& own = mesh().edgeOwner();
    const labelUList& nbr = mesh().edgeNeighbour();

    // Check if internal face-normal vectors diverge (no separation)
    // or converge (separation may occur)
    forAll(nbr, edgei)
    {
        const label faceO = own[edgei];
        const label faceN = nbr[edgei];

        cornerEdges[edgei] = isCornerEdgeSharp
        (
            faceCentres[faceO],
            faceCentres[faceN],
            faceNormals[faceO],
            faceNormals[faceN]
        );
    }


    // Skip the rest of the routine if the simulation is a serial run
    if (!Pstream::parRun()) return cornerEdges;

    // Check if processor face-normal vectors diverge (no separation)
    // or converge (separation may occur)
    const faBoundaryMesh& patches = mesh().boundary();

    for (const faPatch& fap : patches)
    {
        if (isA<processorFaPatch>(fap))
        {
            const label patchi = fap.index();
            const auto& edgeFaces = fap.edgeFaces();
            const label internalEdgei = fap.start();

            const auto& faceCentresp = faceCentres.boundaryField()[patchi];
            const auto& faceNormalsp = faceNormals.boundaryField()[patchi];

            forAll(faceNormalsp, bndEdgei)
            {
                const label faceO = edgeFaces[bndEdgei];
                const label meshEdgei = internalEdgei + bndEdgei;

                cornerEdges[meshEdgei] = isCornerEdgeSharp
                (
                    faceCentres[faceO],
                    faceCentresp[bndEdgei],
                    faceNormals[faceO],
                    faceNormalsp[bndEdgei]
                );
            }
        }
    }

    return cornerEdges;
}


bool FriedrichModel::isCornerEdgeSharp
(
    const vector& faceCentreO,
    const vector& faceCentreN,
    const vector& faceNormalO,
    const vector& faceNormalN
) const
{
    // Calculate the relative position of centres of faces sharing an edge
    const vector relativePosition(faceCentreN - faceCentreO);

    // Calculate the relative normal of faces sharing an edge
    const vector relativeNormal(faceNormalN - faceNormalO);

    // Return true if the face normals converge, meaning that the edge is sharp
    return ((relativeNormal & relativePosition) < -1e-8);
}


scalarList FriedrichModel::calcCornerAngles() const
{
    scalarList cornerAngles(mesh().nEdges(), Zero);

    const areaVectorField& faceNormals = mesh().faceAreaNormals();

    const labelUList& own = mesh().edgeOwner();
    const labelUList& nbr = mesh().edgeNeighbour();

    // Process internal edges
    forAll(nbr, edgei)
    {
        if (!cornerEdges_[edgei]) continue;

        const label faceO = own[edgei];
        const label faceN = nbr[edgei];

        cornerAngles[edgei] = calcCornerAngle
        (
            faceNormals[faceO],
            faceNormals[faceN]
        );
    }


    // Skip the rest of the routine if the simulation is a serial run
    if (!Pstream::parRun()) return cornerAngles;

    // Process processor edges
    const faBoundaryMesh& patches = mesh().boundary();

    for (const faPatch& fap : patches)
    {
        if (isA<processorFaPatch>(fap))
        {
            const label patchi = fap.index();
            const auto& edgeFaces = fap.edgeFaces();
            const label internalEdgei = fap.start();

            const auto& faceNormalsp = faceNormals.boundaryField()[patchi];

            forAll(faceNormalsp, bndEdgei)
            {
                const label faceO = edgeFaces[bndEdgei];
                const label meshEdgei = internalEdgei + bndEdgei;

                if (!cornerEdges_[meshEdgei]) continue;

                cornerAngles[meshEdgei] = calcCornerAngle
                (
                    faceNormals[faceO],
                    faceNormalsp[bndEdgei]
                );
            }
        }
    }

    return cornerAngles;
}


scalar FriedrichModel::calcCornerAngle
(
    const vector& faceNormalO,
    const vector& faceNormalN
) const
{
    const scalar magFaceNormal = mag(faceNormalO)*mag(faceNormalN);

    // Avoid any potential exceptions during the cosine calculations
    if (magFaceNormal < SMALL) return 0;

    scalar cosAngle = (faceNormalO & faceNormalN)/magFaceNormal;
    cosAngle = clamp(cosAngle, -1, 1);

    return std::acos(cosAngle);
}


bitSet FriedrichModel::calcSeparationFaces() const
{
    bitSet separationFaces(mesh().faces().size(), false);

    const edgeScalarField& phis = film().phi2s();

    const labelUList& own = mesh().edgeOwner();
    const labelUList& nbr = mesh().edgeNeighbour();

    // Process internal faces
    forAll(nbr, edgei)
    {
        if (!cornerEdges_[edgei]) continue;

        const label faceO = own[edgei];
        const label faceN = nbr[edgei];

        isSeparationFace
        (
            separationFaces,
            phis[edgei],
            faceO,
            faceN
        );
    }


    // Skip the rest of the routine if the simulation is a serial run
    if (!Pstream::parRun()) return separationFaces;

    // Process processor faces
    const faBoundaryMesh& patches = mesh().boundary();

    for (const faPatch& fap : patches)
    {
        if (isA<processorFaPatch>(fap))
        {
            const label patchi = fap.index();
            const auto& edgeFaces = fap.edgeFaces();
            const label internalEdgei = fap.start();

            const auto& phisp = phis.boundaryField()[patchi];

            forAll(phisp, bndEdgei)
            {
                const label faceO = edgeFaces[bndEdgei];
                const label meshEdgei(internalEdgei + bndEdgei);

                if (!cornerEdges_[meshEdgei]) continue;

                isSeparationFace
                (
                    separationFaces,
                    phisp[bndEdgei],
                    faceO
                );
            }
        }
    }

    return separationFaces;
}


void FriedrichModel::isSeparationFace
(
    bitSet& separationFaces,
    const scalar phiEdge,
    const label faceO,
    const label faceN
) const
{
    const scalar tol = 1e-8;

    // Assuming there are no sources/sinks at the edge
    if (phiEdge > tol)  // From owner to neighbour
    {
        separationFaces[faceO] = true;
    }
    else if ((phiEdge < -tol) && (faceN != -1))  // From neighbour to owner
    {
        separationFaces[faceN] = true;
    }
}


scalarList FriedrichModel::calcSeparationAngles
(
    const bitSet& separationFaces
) const
{
    scalarList separationAngles(mesh().faces().size(), Zero);

    const labelUList& own = mesh().edgeOwner();
    const labelUList& nbr = mesh().edgeNeighbour();

    // Process internal faces
    forAll(nbr, edgei)
    {
        if (!cornerEdges_[edgei]) continue;

        const label faceO = own[edgei];
        const label faceN = nbr[edgei];

        if (separationFaces[faceO])
        {
            separationAngles[faceO] = cornerAngles_[edgei];
        }

        if (separationFaces[faceN])
        {
            separationAngles[faceN] = cornerAngles_[edgei];
        }
    }


    // Skip the rest of the routine if the simulation is a serial run
    if (!Pstream::parRun()) return separationAngles;

    // Process processor faces
    const edgeScalarField& phis = film().phi2s();
    const faBoundaryMesh& patches = mesh().boundary();

    for (const faPatch& fap : patches)
    {
        if (isA<processorFaPatch>(fap))
        {
            const label patchi = fap.index();
            const auto& edgeFaces = fap.edgeFaces();
            const label internalEdgei = fap.start();

            const auto& phisp = phis.boundaryField()[patchi];

            forAll(phisp, bndEdgei)
            {
                const label faceO = edgeFaces[bndEdgei];
                const label meshEdgei(internalEdgei + bndEdgei);

                if (!cornerEdges_[meshEdgei]) continue;

                if (separationFaces[faceO])
                {
                    separationAngles[faceO] = cornerAngles_[meshEdgei];
                }
            }
        }
    }

    return separationAngles;
}


tmp<scalarField> FriedrichModel::Fratio() const
{
    const areaVectorField Up(film().Up());
    const areaVectorField& Uf = film().Uf();
    const areaScalarField& h = film().h();
    const areaScalarField& rho = film().rho();
    const areaScalarField& mu = film().mu();
    const areaScalarField& sigma = film().sigma();

    // Identify the faces where separation may occur
    const bitSet separationFaces(calcSeparationFaces());

    // Calculate the corner angles corresponding to the separation faces
    const scalarList separationAngles(calcSeparationAngles(separationFaces));

    // Initialize the force ratio
    auto tFratio = tmp<scalarField>::New(mesh().faces().size(), Zero);
    auto& Fratio = tFratio.ref();

    // Process internal faces
    forAll(separationFaces, i)
    {
        // Skip the routine if the face is not a candidate for separation
        if (!separationFaces[i]) continue;

        // Calculate the corner-angle trigonometric values
        const scalar sinAngle = std::sin(separationAngles[i]);
        const scalar cosAngle = std::cos(separationAngles[i]);

        // Reynolds number (FLW:Eq. 16)
        const scalar Re = h[i]*mag(Uf[i])*rho[i]/mu[i];

        // Weber number (FLW:Eq. 17)
        const scalar We =
            h[i]*rhop_*sqr(mag(Up[i]) - mag(Uf[i]))/(2.0*sigma[i]);

        // Characteristic breakup length (FLW:Eq. 15)
        const scalar Lb =
            0.0388*Foam::sqrt(h[i])*Foam::pow(Re, 0.6)*Foam::pow(We, -0.5);

        // Force ratio - denominator (FLW:Eq. 20)
        const scalar den =
            sigma[i]*(sinAngle + 1.0) + rho[i]*magG_*h[i]*Lb*cosAngle;

        if (mag(den) > 0)
        {
            // Force ratio (FLW:Eq. 20)
            Fratio[i] = rho[i]*sqr(mag(Uf[i]))*h[i]*sinAngle/den;
        }
    }


    // Skip the rest of the routine if the simulation is a serial run
    if (!Pstream::parRun()) return tFratio;

    // Process processor faces
    const faBoundaryMesh& patches = mesh().boundary();

    for (const faPatch& fap : patches)
    {
        if (isA<processorFaPatch>(fap))
        {
            const label patchi = fap.index();
            const label internalEdgei = fap.start();

            const auto& hp = h.boundaryField()[patchi];
            const auto& Ufp = Uf.boundaryField()[patchi];
            const auto& Upp = Up.boundaryField()[patchi];
            const auto& rhop = rho.boundaryField()[patchi];
            const auto& sigmap = sigma.boundaryField()[patchi];
            const auto& mup = mu.boundaryField()[patchi];

            forAll(hp, i)
            {
                // Skip the routine if the face is not a candidate for separation
                if (!separationFaces[i]) continue;

                const label meshEdgei = internalEdgei + i;

                // Calculate the corner-angle trigonometric values
                const scalar sinAngle = std::sin(cornerAngles_[meshEdgei]);
                const scalar cosAngle = std::cos(cornerAngles_[meshEdgei]);

                // Reynolds number (FLW:Eq. 16)
                const scalar Re = hp[i]*mag(Ufp[i])*rhop[i]/mup[i];

                // Weber number (FLW:Eq. 17)
                const scalar We =
                    hp[i]*rhop_*sqr(mag(Upp[i]) - mag(Ufp[i]))/(2.0*sigmap[i]);

                // Characteristic breakup length (FLW:Eq. 15)
                const scalar Lb =
                    0.0388*Foam::sqrt(hp[i])
                   *Foam::pow(Re, 0.6)*Foam::pow(We, -0.5);

                // Force ratio - denominator (FLW:Eq. 20)
                const scalar den =
                    sigmap[i]*(sinAngle + 1.0)
                  + rhop[i]*magG_*hp[i]*Lb*cosAngle;

                if (mag(den) > 0)
                {
                    // Force ratio (FLW:Eq. 20)
                    Fratio[i] = rhop[i]*sqr(mag(Ufp[i]))*hp[i]*sinAngle/den;
                }
            }
        }
    }

    return tFratio;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

FriedrichModel::FriedrichModel
(
    const regionModels::areaSurfaceFilmModels::liquidFilmBase& film,
    const dictionary& dict
)
:
    filmSeparationModel(film, dict),
    separation_
    (
        separationTypeNames.getOrDefault
        (
            "separationType",
            dict,
            separationType::FULL
        )
    ),
    rhop_(dict.getScalar("rhop")),
    magG_(mag(film.g().value())),
    C0_(dict.getOrDefault<scalar>("C0", 0.882)),
    C1_(dict.getOrDefault<scalar>("C1", -1.908)),
    C2_(dict.getOrDefault<scalar>("C2", 1.264)),
    cornerEdges_(calcCornerEdges()),
    cornerAngles_(calcCornerAngles())
{
    if (rhop_ < VSMALL)
    {
        FatalIOErrorInFunction(dict)
            << "Primary-phase density, rhop: " << rhop_ << " must be non-zero."
            << abort(FatalIOError);
    }

    if (mag(C2_) < VSMALL)
    {
        FatalIOErrorInFunction(dict)
            << "Empirical constant, C2 = " << C2_ << "cannot be zero."
            << abort(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField> FriedrichModel::separatedMassRatio() const
{
    tmp<scalarField> tFratio = Fratio();
    const auto& Fratio = tFratio.cref();

    // Initialize the mass ratio of film separation
    auto tseparated = tmp<scalarField>::New(mesh().faces().size(), Zero);
    auto& separated = tseparated.ref();


    switch (separation_)
    {
        case separationType::FULL:
        {
            forAll(Fratio, i)
            {
                if (Fratio[i] > 1)
                {
                    separated[i] = 1;
                }
            }
            break;
        }
        case separationType::PARTIAL:
        {
            forAll(Fratio, i)
            {
                if (Fratio[i] > 1)
                {
                    // (ZJD:Eq. 16)
                    separated[i] = C0_ + C1_*Foam::exp(-Fratio[i]/C2_);
                }
            }
            break;
        }
        default:
            break;  // This should not happen.
    }

    if (debug && mesh().time().writeTime())
    {
        {
            areaScalarField areaFratio
            (
                mesh().newIOobject("Fratio"),
                mesh(),
                dimensionedScalar(dimForce, Zero)
            );
            areaFratio.primitiveFieldRef() = Fratio;
            areaFratio.write();
        }
    }


    return tseparated;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace filmSeparationModels
} // End namespace Foam


// ************************************************************************* //

