/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "singleCellFvMesh.H"
#include "calculatedFvPatchFields.H"
#include "directFvPatchFieldMapper.H"
#include "Time.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::singleCellFvMesh::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    // 1. Create the complete field with dummy patch fields

    tmp<GeometricField<Type, fvPatchField, volMesh>> tresult
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                vf.name(),
                time().timeName(),
                *this,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            *this,
            Type(gAverage(vf.primitiveField())),
            vf.dimensions(),
            fvPatchFieldBase::calculatedType()
        )
    );
    auto& result = tresult.ref();


    // 2. Change the fvPatchFields to the correct type using a mapper
    //  constructor (with reference to the now correct internal field)

    auto& bf = result.boundaryFieldRef();

    if (agglomerate())
    {
        forAll(vf.boundaryField(), patchi)
        {
            const labelList& agglom = patchFaceAgglomeration_[patchi];
            label nAgglom = max(agglom)+1;

            // Use inverse of agglomeration. This is from agglomeration to
            // original (fine) mesh patch face.
            labelListList coarseToFine(invertOneToMany(nAgglom, agglom));
            inplaceReorder(patchFaceMap_[patchi], coarseToFine);
            scalarListList coarseWeights(nAgglom);
            forAll(coarseToFine, coarseI)
            {
                const labelList& fineFaces = coarseToFine[coarseI];
                coarseWeights[coarseI] = scalarList
                (
                    fineFaces.size(),
                    1.0/fineFaces.size()
                );
            }

            bf.set
            (
                patchi,
                fvPatchField<Type>::New
                (
                    vf.boundaryField()[patchi],
                    boundary()[patchi],
                    result.internalField(),
                    agglomPatchFieldMapper(coarseToFine, coarseWeights)
                )
            );
        }
    }
    else
    {
        forAll(vf.boundaryField(), patchi)
        {
            labelList map(identity(vf.boundaryField()[patchi].size()));

            bf.set
            (
                patchi,
                fvPatchField<Type>::New
                (
                    vf.boundaryField()[patchi],
                    boundary()[patchi],
                    result.internalField(),
                    directFvPatchFieldMapper(map)
                )
            );
        }
    }

    return tresult;
}


// ************************************************************************* //
