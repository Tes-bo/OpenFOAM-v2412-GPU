/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023-2024 OpenCFD Ltd.
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

template<class Type>
void Foam::VF::raySearchEngine::interpolate
(
    GeometricField<Type, fvPatchField, volMesh>& fld,
    const List<List<Type>>& values
) const
{
    label compacti = 0;

    auto& vfbf = fld.boundaryFieldRef();

    if (agglomMeshPtr_)
    {
        const auto& coarseMesh = agglomMeshPtr_();
        const auto& finalAgglom = coarseMesh.patchFaceAgglomeration();

        for (const label patchi : patchIDs_)
        {
            const labelList& agglom = finalAgglom[patchi];

            if (agglom.empty()) continue;

            label nAgglom = max(agglom) + 1;
            const labelListList coarseToFine(invertOneToMany(nAgglom, agglom));
            const labelList& coarsePatchFace =
                coarseMesh.patchFaceMap()[patchi];

            forAll(coarseToFine, i)
            {
                const label coarseFacei = coarsePatchFace[i];
                const labelList& fineFaces = coarseToFine[coarseFacei];
                const Type sumValues = sum(values[compacti]);

                for (const label fineFacei : fineFaces)
                {
                    vfbf[patchi][fineFacei] = sumValues;
                }

                ++compacti;
            }
        }
    }
    else
    {
        label compacti = 0;
        for (const label patchi : patchIDs_)
        {
            auto& vfp = vfbf[patchi];

            for (Type& vfi : vfp)
            {
                vfi = sum(values[compacti++]);
            }
        }
    }
}


// ************************************************************************* //
