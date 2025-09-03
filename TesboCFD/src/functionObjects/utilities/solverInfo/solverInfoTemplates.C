/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2016 OpenFOAM Foundation
    Copyright (C) 2015-2021 OpenCFD Ltd.
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

#include "solverInfo.H"
#include "volFields.H"
#include "ListOps.H"
#include "zeroGradientFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::functionObjects::solverInfo::writeFileHeader
(
    Ostream& os,
    const word& fieldName
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    if (foundObject<fieldType>(fieldName))
    {
        writeTabbed(os, fieldName + "_solver");

        typename pTraits<Type>::labelType validComponents
        (
            mesh_.validComponents<Type>()
        );

        for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; ++cmpt)
        {
            if (component(validComponents, cmpt) != -1)
            {
                const word cmptName(pTraits<Type>::componentNames[cmpt]);
                const word fieldBase(fieldName + cmptName);

                writeTabbed(os, fieldBase + "_initial");
                writeTabbed(os, fieldBase + "_final");
                writeTabbed(os, fieldBase + "_iters");
            }
        }

        writeTabbed(os, fieldName + "_converged");
    }
}


template<class Type>
void Foam::functionObjects::solverInfo::initialiseResidualField
(
    const word& fieldName
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> volFieldType;

    if (foundObject<volFieldType>(fieldName))
    {
        const Foam::dictionary& solverDict =
            mesh_.data().solverPerformanceDict();

        if (solverDict.found(fieldName))
        {
            typename pTraits<Type>::labelType validComponents
            (
                mesh_.validComponents<Type>()
            );

            for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; ++cmpt)
            {
                if (component(validComponents, cmpt) != -1)
                {
                    const word resultName
                    (
                        fieldName + word(pTraits<Type>::componentNames[cmpt])
                    );

                    createResidualField(resultName);
                }
            }
        }
    }
}


template<class Type>
void Foam::functionObjects::solverInfo::updateSolverInfo(const word& fieldName)
{
    typedef GeometricField<Type, fvPatchField, volMesh> volFieldType;
    typedef typename pTraits<Type>::labelType labelType;

    if (foundObject<volFieldType>(fieldName))
    {
        const Foam::dictionary& solverDict =
            mesh_.data().solverPerformanceDict();

        if (solverDict.found(fieldName))
        {
            const List<SolverPerformance<Type>> sp
            (
                solverDict.lookup(fieldName)
            );

            const SolverPerformance<Type>& sp0 = sp.first();
            const word& solverName = sp0.solverName();
            const Type& initialResidual = sp0.initialResidual();
            const Type& finalResidual = sp0.finalResidual();
            const labelType nIterations = sp0.nIterations();
            const Switch converged(sp0.converged());

            const labelType validComponents(mesh_.validComponents<Type>());

            file() << token::TAB << solverName;

            for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; ++cmpt)
            {
                if (component(validComponents, cmpt) != -1)
                {
                    const scalar ri = component(initialResidual, cmpt);
                    const scalar rf = component(finalResidual, cmpt);
                    const label n = component(nIterations, cmpt);

                    file()
                        << token::TAB << ri
                        << token::TAB << rf
                        << token::TAB << n;

                    const word cmptName(pTraits<Type>::componentNames[cmpt]);
                    const word resultName(fieldName + cmptName);
                    setResult(resultName + "_initial", ri);
                    setResult(resultName + "_final", rf);
                    setResult(resultName + "_iters", n);
                }
            }

            file() << token::TAB << converged;
        }
    }
}


// ************************************************************************* //
