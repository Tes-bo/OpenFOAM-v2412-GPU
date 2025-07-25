/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
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

#include "fvMotionSolver.H"
#include "fixedValuePointPatchFields.H"
#include "cellMotionFvPatchFields.H"
#include "facePointPatch.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::wordList Foam::fvMotionSolver::cellMotionBoundaryTypes
(
    const typename GeometricField<Type, pointPatchField, pointMesh>::
    Boundary& pmUbf
) const
{
    wordList cmUbf(fvMesh_.boundary().size());

    forAll(pmUbf, patchi)
    {
        const auto& pfld = pmUbf[patchi];
        const auto* fppPtr = isA<facePointPatch>(pfld.patch());
        if (fppPtr)
        {
            const auto& fpp = *fppPtr;
            const label polyPatchi = fpp.patch().index();

            if (isA<fixedValuePointPatchField<Type>>(pfld))
            {
                cmUbf[polyPatchi] = cellMotionFvPatchField<Type>::typeName;
            }
            else
            {
                // Take over pointPatch type
                cmUbf[polyPatchi] = pfld.type();
            }

            if (debug)
            {
                Pout<< "Patch:" << fvMesh_.boundary()[patchi].patch().name()
                    << " pointType:" << pfld.type()
                    << " cellType:" << cmUbf[patchi] << endl;
            }
        }
    }

    return cmUbf;
}


// ************************************************************************* //
