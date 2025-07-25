/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2023 PCOpt/NTUA
    Copyright (C) 2013-2023 FOSS GP
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

#include "BezierDesignVariables.H"
#include "IOmanip.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(BezierDesignVariables, 0);
    addToRunTimeSelectionTable
    (
        shapeDesignVariables,
        BezierDesignVariables,
        dictionary
    );
}


// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

void Foam::BezierDesignVariables::readBounds
(
    autoPtr<scalar> lowerBoundPtr,
    autoPtr<scalar> upperBoundPtr
)
{
    designVariables::readBounds(lowerBoundPtr, upperBoundPtr);

    if (dict_.found("lowerCPBounds"))
    {
        vector lowerCPBounds(dict_.get<vector>("lowerCPBounds"));
        lowerBounds_.reset(new scalarField(getVars().size(), Zero));
        setBounds(lowerBounds_, lowerCPBounds);
    }

    if (dict_.found("upperCPBounds"))
    {
        vector upperCPBounds(dict_.get<vector>("upperCPBounds"));
        upperBounds_.reset(new scalarField(getVars().size(), Zero));
        setBounds(upperBounds_, upperCPBounds);
    }
}


void Foam::BezierDesignVariables::setBounds
(
    autoPtr<scalarField>& bounds,
    const vector& cpBounds
)
{
    bounds.reset(new scalarField(getVars().size(), Zero));
    const label nCPs(bezier_.nBezier());
    for (label iCP = 0; iCP < nCPs; ++iCP)
    {
        bounds()[iCP] = cpBounds.x();
        bounds()[nCPs + iCP] = cpBounds.y();
        bounds()[2*nCPs + iCP] = cpBounds.z();
    }
}


Foam::tmp<Foam::vectorField>
Foam::BezierDesignVariables::computeBoundaryDisplacement
(
    const scalarField& correction
)
{
    // Reset boundary movement field
    dx_.primitiveFieldRef() = Zero;

    // Compute boundary movement using the derivatives of grid nodes
    // wrt to the Bezier control points and the correction
    const label nCPs(bezier_.nBezier());
    auto tcpMovement(tmp<vectorField>::New(nCPs, Zero));
    vectorField& cpMovement = tcpMovement.ref();
    const boolListList& confineMovement = bezier_.confineMovement();

    forAll(cpMovement, cpI)
    {
        if (!confineMovement[0][cpI])
        {
            cpMovement[cpI].x() = correction[cpI];
        }
        if (!confineMovement[1][cpI])
        {
            cpMovement[cpI].y() = correction[nCPs + cpI];
        }
        if (!confineMovement[2][cpI])
        {
            cpMovement[cpI].z() = correction[2*nCPs + cpI];
        }

        dx_ += (bezier_.dxidXj()[cpI] & cpMovement[cpI]);
    }

    return tcpMovement;
}


void Foam::BezierDesignVariables::decomposeVarID
(
    label& cpI,
    label& dir,
    const label varID
) const
{
    const label nBezier = bezier_.nBezier();
    cpI = varID%nBezier;
    dir = varID/nBezier;
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::BezierDesignVariables::BezierDesignVariables
(
    fvMesh& mesh,
    const dictionary& dict
)
:
    shapeDesignVariables(mesh, dict),
    bezier_
    (
        mesh,
        IOdictionary
        (
            IOobject
            (
                "optimisationDict",
                mesh_.time().globalPath()/"system",
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            )
        )
    ),
    dx_
    (
        IOobject
        (
            "dx",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pointMesh::New(mesh_),
        dimensionedVector(dimless, Zero)
    )
{
    // Set the size of the design variables field
    scalarField::setSize(3*bezier_.nBezier(), Zero);

    // Set the active design variables
    activeDesignVariables_ = bezier_.getActiveDesignVariables();

    // Read bounds
    readBounds();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::BezierDesignVariables::update(scalarField& correction)
{
    // Translate the correction field to control point movements
    computeBoundaryDisplacement(correction);

    // Transfer movement to the displacementMethod
    displMethodPtr_->setMotionField(dx_);

    // Update the design variables
    scalarField::operator+=(correction);

    // Do the actual mesh movement
    moveMesh();
}


Foam::scalar Foam::BezierDesignVariables::computeEta(scalarField& correction)
{
    // Transfer the correction field to control point movement
    computeBoundaryDisplacement(correction);

    const scalar maxDisplacement(max(mag(dx_)).value());

    Info<< "maxAllowedDisplacement/maxDisplacement at the boundary\t"
        << maxInitChange_() << "/" << maxDisplacement << endl;

    const scalar eta = maxInitChange_()/maxDisplacement;
    Info<< "Setting eta value to " << eta << endl;
    correction *= eta;

    return eta;
}


bool Foam::BezierDesignVariables::globalSum() const
{
    return false;
}


Foam::tmp<Foam::vectorField> Foam::BezierDesignVariables::dxdbFace
(
    const label patchI,
    const label varID
) const
{
    label cpI(-1), dir(-1);
    decomposeVarID(cpI, dir, varID);
    return bezier_.dxdbFace(patchI, cpI, dir);
}


Foam::tmp<Foam::vectorField> Foam::BezierDesignVariables::dndb
(
    const label patchI,
    const label varID
) const
{
    label cpI(-1), dir(-1);
    decomposeVarID(cpI, dir, varID);
    return bezier_.dndbBasedSensitivities(patchI, cpI, dir, false);
}


Foam::tmp<Foam::vectorField> Foam::BezierDesignVariables::dSdb
(
    const label patchI,
    const label varID
) const
{
    label cpI(-1), dir(-1);
    decomposeVarID(cpI, dir, varID);
    return bezier_.dndbBasedSensitivities(patchI, cpI, dir, true);
}


Foam::tmp<Foam::volVectorField>
Foam::BezierDesignVariables::dCdb(const label varID) const
{
    label cpI(-1), dir(-1);
    decomposeVarID(cpI, dir, varID);
    label patchI(-1);
    // There is no mechanism in place to identify the parametertised patch.
    // Look over all patches and grab one with a non-zero dxdb
    for (const label pI : parametertisedPatches_)
    {
        tmp<vectorField> dxdbFace = bezier_.dxdbFace(pI, cpI, dir);
        if (gSum(mag(dxdbFace)) > SMALL)
        {
            patchI = pI;
        }
    }
    return solveMeshMovementEqn(patchI, varID);
}


// ************************************************************************* //
