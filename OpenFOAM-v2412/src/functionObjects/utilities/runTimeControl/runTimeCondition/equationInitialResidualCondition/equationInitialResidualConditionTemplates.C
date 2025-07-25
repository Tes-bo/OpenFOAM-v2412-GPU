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

#include "fvMesh.H"

template<class Type>
void Foam::functionObjects::runTimeControls::
equationInitialResidualCondition::setResidual
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& fieldName,
    const label componenti,
    bool& canSet,
    scalar& residual
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> volFieldType;

    if (canSet && mesh.foundObject<volFieldType>(fieldName))
    {
        const List<SolverPerformance<Type>> sp(dict.lookup(fieldName));
        const Type& allComponents = sp.first().initialResidual();

        if (componenti != -1)
        {
            if (componenti > pTraits<Type>::nComponents - 1)
            {
                FatalErrorInFunction
                    << "Requested component [" << componenti
                    <<  "] for field " << fieldName
                    << " is out of range 0.."
                    << pTraits<Type>::nComponents - 1
                    << exit(FatalError);
            }

            residual = component(allComponents, componenti);
        }
        else
        {
            residual = cmptMax(allComponents);
        }

        canSet = false;
    }
}
