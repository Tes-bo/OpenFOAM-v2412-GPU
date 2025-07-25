/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
    Copyright (C) 2020-2021 OpenCFD Ltd.
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

#include "rigidBodyMotion.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::RBD::rigidBodyMotion::read(const dictionary& dict)
{
    if (rigidBodyModel::read(dict))
    {
        aRelax_ =
            Function1<scalar>::NewIfPresent
            (
                "accelerationRelaxation",
                dict,
                word::null,
                &time()
            );
        aDamp_ = dict.getOrDefault<scalar>("accelerationDamping", 1);
        report_ = dict.getOrDefault<Switch>("report", false);

        return true;
    }

    return false;
}


void Foam::RBD::rigidBodyMotion::write(Ostream& os) const
{
    rigidBodyModel::write(os);

    if (aRelax_)
    {
        aRelax_->writeData(os);
    }
    os.writeEntry("accelerationDamping", aDamp_);
    os.writeEntry("report", report_);
}


// ************************************************************************* //
