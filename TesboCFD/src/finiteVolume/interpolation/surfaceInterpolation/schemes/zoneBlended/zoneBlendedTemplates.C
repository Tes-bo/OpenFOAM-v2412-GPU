/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::zoneBlended<Type>::setSchemes
(
    const fvMesh& mesh,
    const dictionary& dict
)
{
    zoneNames_.resize(dict.size());
    schemePtrs_.resize(dict.size());

    schemePtrs_[0] =
        surfaceInterpolationScheme<Type>::New(mesh, dict.lookup("default"));

    zoneNames_[0] = "default";

    corrected_ = scheme(0).corrected();

    label schemei = 1;

    for (const auto& e : dict)
    {
        if (e.isDict())
        {
            FatalIOErrorInFunction(dict)
                << "Entries must be given in zoneName-scheme pairs"
                << abort(FatalIOError);
        }

        const word& key = e.keyword();

        if (key == "default")
        {
            // Background scheme - already handled at index 0
        }
        else
        {
            zoneNames_[schemei] = key;

            schemePtrs_[schemei] =
                surfaceInterpolationScheme<Type>::New(mesh, e.stream());

            corrected_ = corrected_ || scheme(schemei).corrected();

            ++schemei;
        }
    }
}


template<class Type>
void Foam::zoneBlended<Type>::setSchemes
(
    const fvMesh& mesh,
    const surfaceScalarField& faceFlux,
    const dictionary& dict
)
{
    zoneNames_.resize(dict.size());
    schemePtrs_.resize(dict.size());

    schemePtrs_[0] =
        surfaceInterpolationScheme<Type>::New
        (
            mesh,
            faceFlux,
            dict.lookup("default")
        );

    zoneNames_[0] = "default";

    corrected_ = scheme(0).corrected();

    label schemei = 1;

    for (const auto& e : dict)
    {
        if (e.isDict())
        {
            FatalIOErrorInFunction(dict)
                << "Entries must be given in faceZoneName-scheme pairs"
                << abort(FatalIOError);
        }

        const word& key = e.keyword();

        if (key == "default")
        {
            // Background scheme - already handled at index 0
        }
        else
        {
            zoneNames_[schemei] = key;

            schemePtrs_[schemei] =
                surfaceInterpolationScheme<Type>::New
                (
                    mesh,
                    faceFlux,
                    e.stream()
                );

            corrected_ = corrected_ || scheme(schemei).corrected();

            ++schemei;
        }
    }
}


template<class Type>
template<class FieldType>
void Foam::zoneBlended<Type>::setFaceZoneValues
(
    FieldType& dest,
    const FieldType& src,
    const faceZone& fz
) const
{
    const auto& mesh = dest.mesh();
    const auto& pbm = mesh.boundaryMesh();
    const auto& srcBf = src.boundaryField();
    auto& destBf = dest.boundaryFieldRef();

    for (const label facei : fz)
    {
        if (mesh.isInternalFace(facei))
        {
            dest[facei] = src[facei];
        }
        else
        {
            const labelPair pf = pbm.whichPatchFace(facei);
            auto& pdest = destBf[pf.first()];
            if (pdest.size())
            {
                pdest[pf.second()] = srcBf[pf.first()][pf.second()];
            }
        }
    }
}


template<class Type>
template<class FieldType>
void Foam::zoneBlended<Type>::zeroFaceZoneValues
(
    FieldType& dest,
    const faceZone& fz
) const
{
    const auto& mesh = dest.mesh();
    const auto& pbm = mesh.boundaryMesh();
    auto& destBf = dest.boundaryFieldRef();

    for (const label facei : fz)
    {
        if (mesh.isInternalFace(facei))
        {
            dest[facei] = pTraits<Type>::zero;
        }
        else
        {
            const labelPair pf = pbm.whichPatchFace(facei);
            auto& pdest = destBf[pf.first()];
            if (pdest.size())
            {
                pdest[pf.second()] = pTraits<Type>::zero;
            }
        }
    }
}


// ************************************************************************* //
