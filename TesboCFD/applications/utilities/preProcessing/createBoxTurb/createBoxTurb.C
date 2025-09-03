/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenCFD Ltd.
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

Application
    createBoxTurb

Description
    Create a box of isotropic turbulence based on a user-specified
    energy spectrum.

    Based on the reference
    \verbatim
    Saad, T., Cline, D., Stoll, R., Sutherland, J.C.
    "Scalable Tools for Generating Synthetic Isotropic Turbulence with
    Arbitrary Spectra"
    AIAA Journal, Vol. 55, No. 1 (2017), pp. 327-331.
    \endverbatim

    The \c -createBlockMesh option creates a block mesh and exits, which
    can then be decomposed and the utility run in parallel.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "block.H"
#include "mathematicalConstants.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::vector randomUnitVector(Random& rndGen)
{
    // Sample point on a sphere
    scalar t = rndGen.globalPosition<scalar>(-1, 1);
    scalar phim = rndGen.globalSample01<scalar>()*mathematical::twoPi;
    scalar thetam = Foam::acos(t);

    return vector
    (
        Foam::sin(thetam*Foam::cos(phim)),
        Foam::sin(thetam*Foam::sin(phim)),
        Foam::cos(thetam)
    );
}


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Create a box of isotropic turbulence based on a user-specified"
        " energy spectrum."
    );

    argList::addBoolOption
    (
        "createBlockMesh",
        "create the block mesh and exit"
    );

    #include "setRootCase.H"

    #include "createTime.H"
    #include "createFields.H"

    if (args.found("createBlockMesh"))
    {
        // Create a box block mesh with cyclic patches
        #include "createBlockMesh.H"
        Info<< "\nEnd\n" << endl;
        return 0;
    }

    #include "createMesh.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Minimum wave number
    scalar kappa0 = mathematical::twoPi/cmptMin(L);

    // Maximum wave number
    scalar kappaMax = mathematical::pi/cmptMin(delta);

    Info<< "Wave number min/max = " << kappa0 << ", " << kappaMax << endl;

    Info<< "Generating velocity field" << endl;

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector(dimVelocity, Zero)
    );

    vectorField& Uc = U.primitiveFieldRef();
    const scalar deltaKappa = (kappaMax - kappa0)/scalar(nModes - 1);
    const vectorField& C(mesh.C());
    for (label modei = 1; modei <= nModes; ++modei)
    {
        // Equidistant wave mode
        scalar kappaM = kappa0 + deltaKappa*(modei-1);

        Info<< "Processing mode:" << modei << " kappaM:" << kappaM << endl;

        // Energy
        scalar E = Ek->value(kappaM);

        // Wave amplitude
        scalar qm = Foam::sqrt(E*deltaKappa);

        // Wave number unit vector
        const vector kappaHatm(randomUnitVector(rndGen));

        vector kappaTildem(0.5*kappaM*cmptMultiply(kappaHatm, delta));
        for (direction i = 0; i < 3; ++i)
        {
            kappaTildem[i] = 2/delta[i]*Foam::sin(kappaTildem[i]);
        }

        // Intermediate unit vector zeta
        const vector zetaHatm(randomUnitVector(rndGen));

        // Unit vector sigma
        vector sigmaHatm(zetaHatm^kappaTildem);
        sigmaHatm /= mag(kappaTildem);

        // Phase angle
        scalar psim = 0.5*rndGen.position(-mathematical::pi, mathematical::pi);

        // Add the velocity contribution per mode
        Uc += 2*qm*cos(kappaM*(kappaHatm & C) + psim)*sigmaHatm;
    }

    U.write();

    {
        Info<< "Generating kinetic energy field" << endl;
        volScalarField k("k", 0.5*magSqr(U));
        k.write();
        Info<< "min/max/average k = "
            << gMin(k) << ", " << gMax(k) << ", " << gAverage(k)
            << endl;
    }

    {
        Info<< "Generating div(U) field" << endl;
        volScalarField divU(fvc::div(U));
        divU.write();
        Info<< "min/max/average div(U) = "
            << gMin(divU) << ", " << gMax(divU) << ", " << gAverage(divU)
            << endl;
    }

    Info<< nl;
    runTime.printExecutionTime(Info);

    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //
