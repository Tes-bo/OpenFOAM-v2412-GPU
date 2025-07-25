/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Container, class Type>
Foam::tmp<Foam::Field<Type>>
Foam::wallDistAddressing::collectPatchFields
(
    const Container& bfld
) const
{
    label n = 0;
    for (const auto& patchi : patchIDs_)
    {
        n += bfld[patchi].size();
    }

    auto tresult(tmp<Field<Type>>::New(n));
    auto& result = tresult.ref();

    n = 0;
    for (const auto& patchi : patchIDs_)
    {
        const auto& pfld = bfld[patchi];

        SubList<Type>(result, pfld.size(), n) = pfld;
        n += pfld.size();
    }
    return tresult;
}


template<class VolField>
void Foam::wallDistAddressing::extract
(
    const Field<typename VolField::value_type>& wallFld,
    VolField& fld
) const
{
    {
        const label nUntrafoCells =
        (
            untransformedPatchStarts_.size()
          ? untransformedPatchStarts_[0]
          : untransformedItems_.size()
        );
        for (label i = 0; i < nUntrafoCells; i++)
        {
            const label celli = untransformedItems_[i];
            const label sloti = untransformedSlots_[i];
            fld[celli] = wallFld[sloti];
        }
    }
    {
        const label nTrafoCells =
        (
            transformedPatchStarts_.size()
          ? transformedPatchStarts_[0]
          : transformedItems_.size()
        );
        for (label i = 0; i < nTrafoCells; i++)
        {
            const label celli = transformedItems_[i];
            const label sloti = transformedSlots_[i];
            fld[celli] = wallFld[sloti];
        }
    }

    forAll(fld.boundaryField(), patchi)
    {
        const auto& pfld = fld.boundaryField()[patchi];
        Field<typename VolField::value_type> patchField(pfld.size());

        {
            const label start = untransformedPatchStarts_[patchi];
            const label end = untransformedPatchStarts_[patchi+1];

            for (label i = start; i < end; i++)
            {
                const label facei = untransformedItems_[i];
                const label sloti = untransformedSlots_[i];
                patchField[facei-pfld.patch().start()] = wallFld[sloti];
            }
        }

        {
            const label start = transformedPatchStarts_[patchi];
            const label end = transformedPatchStarts_[patchi+1];

            for (label i = start; i < end; i++)
            {
                const label facei = transformedItems_[i];
                const label sloti = transformedSlots_[i];
                patchField[facei-pfld.patch().start()] = wallFld[sloti];
            }
        }

        fld.boundaryFieldRef()[patchi] = patchField;
    }
    fld.correctBoundaryConditions();
}


template<class VolField, class TransformOp>
const VolField& Foam::wallDistAddressing::map
(
    VolField& fld,
    const TransformOp& top
) const
{
    const auto patchFld
    (
        collectPatchFields
        <
            typename VolField::Boundary,
            typename VolField::value_type
        >
        (
            fld.boundaryField()
        )
    );

    // Distribute and transform
    mapPtr_().distribute
    (
        mesh_.globalData().globalTransforms(),
        patchFld.ref(),
        top
    );

    // Extract into volField
    extract(patchFld(), fld);

    return fld;
}


// ************************************************************************* //
