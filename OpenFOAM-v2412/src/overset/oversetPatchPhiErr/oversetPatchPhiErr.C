/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "oversetPatchPhiErr.H"
#include "oversetFvPatchField.H"
#include "volFields.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void Foam::oversetPatchPhiErr
(
    const fvScalarMatrix& m,
    const surfaceScalarField& phi
)
{
    const volScalarField::Boundary& bm = m.psi().boundaryField();

    forAll(bm, patchi)
    {
        const auto& fvp = bm[patchi];

        if (isA<oversetFvPatchField<scalar>>(fvp))
        {
            const oversetFvPatchField<scalar>& op =
                refCast<const oversetFvPatchField<scalar>>(fvp);

            op.fringeFlux(m, phi);
        }
    }
}


// ************************************************************************* //
