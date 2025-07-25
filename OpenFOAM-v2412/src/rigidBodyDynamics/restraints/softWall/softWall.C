/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "softWall.H"
#include "rigidBodyModel.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace RBD
{
namespace restraints
{
    defineTypeNameAndDebug(softWall, 0);

    addToRunTimeSelectionTable
    (
        restraint,
        softWall,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RBD::restraints::softWall::softWall
(
    const word& name,
    const dictionary& dict,
    const rigidBodyModel& model
)
:
    restraint(name, dict, model)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::RBD::restraints::softWall::~softWall()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::RBD::restraints::softWall::restrain
(
    scalarField& tau,
    Field<spatialVector>& fx,
    const rigidBodyModelState& state
) const
{

    point p = bodyPoint(refAttachmentPt_);

    // Current axis of the spring
    vector r = p - anchor_;

    vector force = vector::zero;
    vector moment = vector::zero;

    vector v = bodyPointVelocity(refAttachmentPt_).l();

    scalar m = model_.bodies()[bodyID_].m();

    scalar d = (wallNormal_/mag(wallNormal_)) & r;

    vector rDir = r/(mag(r) + VSMALL);

    scalar wn =  3.14/C_;
    scalar damping = psi_*2*m*wn;
    scalar stiffness = sqr(wn)*m;

    if (d < 0)
    {
        force =
            (-damping*(rDir & v) + stiffness*d)*rDir;

        moment = (p ^ force);
    }

    if (model_.debug)
    {
        Info<< " stiffness :" << stiffness*d << nl
            << " damping :" << -damping*mag(rDir & v) << nl
            << " force : " << force << nl
            << " d : " << d << nl
            << " r : " << r << nl
            << " p : " << p << nl
            << " velocity : " << v
            << endl;
    }

    // Accumulate the force for the restrained body
    fx[bodyIndex_] += spatialVector(moment, force);
}


bool Foam::RBD::restraints::softWall::read
(
    const dictionary& dict
)
{
    restraint::read(dict);

    coeffs_.readEntry("anchor", anchor_);
    coeffs_.readEntry("refAttachmentPt", refAttachmentPt_);
    coeffs_.readEntry("psi", psi_);
    coeffs_.readEntry("C", C_);
    coeffs_.readEntry("wallNormal", wallNormal_);

    return true;
}


void Foam::RBD::restraints::softWall::write
(
    Ostream& os
) const
{
    restraint::write(os);

    os.writeEntry("anchor", anchor_);
    os.writeEntry("refAttachmentPt", refAttachmentPt_);
    os.writeEntry("psi", psi_);
    os.writeEntry("C", C_);
    os.writeEntry("wallNormal", wallNormal_);

}


// ************************************************************************* //
