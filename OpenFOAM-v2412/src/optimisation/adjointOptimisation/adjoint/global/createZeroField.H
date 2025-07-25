/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Copyright (C) 2013-2019 FOSS GP
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifndef createZeroField_H
#define createZeroField_H

namespace Foam
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
autoPtr<GeometricField<Type, fvPatchField, volMesh>>
createZeroFieldPtr
(
    const fvMesh& mesh,
    const word& name,
    const dimensionSet dims,
    bool printAllocation = false
)
{
    if (printAllocation)
    {
        Info<< "Allocating new volField " << name << nl << endl;
    }

    return
    (
        autoPtr<GeometricField<Type, fvPatchField, volMesh>>::New
        (
            IOobject
            (
               name,
               mesh.time().timeName(),
               mesh,
               IOobject::NO_READ,
               IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<Type>(dims, Zero)
        )
    );
}


template<class Type>
autoPtr<typename GeometricField<Type, fvPatchField, volMesh>::Boundary>
createZeroBoundaryPtr
(
    const fvMesh& mesh,
    bool printAllocation = false
)
{
    if (printAllocation)
    {
        Info<< "Allocating new boundaryField " << nl << endl;
    }

    typedef typename GeometricField<Type, fvPatchField, volMesh>::Boundary
        Boundary;

    // Make sure that the patchFields to be generated will be of type
    // calculated, even if they are of constraint type.
    // Necessary to avoid unexpected behaviour when computing sensitivities
    // on symmetry patches (not a good practice either way)
    const fvBoundaryMesh& bm = mesh.boundary();
    wordList actualPatchTypes(bm.size(), word::null);
    forAll(actualPatchTypes, patchi)
    {
        if
        (
            fvPatchField<Type>::patchConstructorTablePtr_->contains
            (
                bm[patchi].type()
            )
        )
        {
            actualPatchTypes[patchi] = bm[patchi].type();
        }
    }

    autoPtr<Boundary> bPtr
    (
        new Boundary
        (
            mesh.boundary(),
            DimensionedField<Type, volMesh>::null(),
            wordList(bm.size(), fvPatchFieldBase::calculatedType()),
            actualPatchTypes
        )
    );

    // Values are not assigned! Assign manually
    Boundary& bRef = bPtr();
    forAll(bRef, pI)
    {
        bRef[pI] = pTraits<Type>::zero;
    }

    return (bPtr);
}


template<class Type>
autoPtr<List<Field<Type>>>
createZeroBoundaryPointFieldPtr
(
    const fvMesh& mesh,
    bool printAllocation = false
)
{
    if (printAllocation)
    {
        Info<< "Allocating new point boundaryField " << nl << endl;
    }

    autoPtr<List<Field<Type>>> bPtr
    (
        new List<Field<Type>>
        (
            mesh.boundary().size()
        )
    );

    List<Field<Type>>& bRef = bPtr();
    forAll(bRef, pI)
    {
        bRef[pI] =
            Field<Type>
            (
                mesh.boundaryMesh()[pI].nPoints(),
                pTraits<Type>::zero
            );
    }

    return (bPtr);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
