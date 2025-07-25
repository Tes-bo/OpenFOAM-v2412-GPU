/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2023 OpenCFD Ltd.
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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class GeoField>
Foam::label Foam::functionObjects::vtkWrite::writeVolFieldsImpl
(
    autoPtr<vtk::internalWriter>& internalWriter,
    UPtrList<vtk::patchWriter>& patchWriters,
    const fvMeshSubset& proxy,
    const wordHashSet& candidateNames
) const
{
    const fvMesh& baseMesh = proxy.baseMesh();

    label count = 0;

    for
    (
        const GeoField& origField
      : baseMesh.csorted<GeoField>(candidateNames)
    )
    {
        bool ok = false;
        auto tfield = fvMeshSubsetProxy::interpolate(proxy, origField);
        const auto& field = tfield();

        // Internal
        if (internalWriter)
        {
            ok = true;
            internalWriter->write(field);
        }

        // Boundary
        for (vtk::patchWriter& writer : patchWriters)
        {
            ok = true;
            writer.write(field);
        }

        if (ok)
        {
            ++count;

            if (verbose_)
            {
                if (count == 1)
                {
                    Log << "    " << GeoField::typeName << '(';
                }
                else
                {
                    Log << ' ';
                }
                Log << origField.name();
            }
        }
    }

    if (verbose_ && count)
    {
        Log << ')' << endl;
    }

    return count;
}


template<class GeoField>
Foam::label Foam::functionObjects::vtkWrite::writeVolFieldsImpl
(
    autoPtr<vtk::internalWriter>& internalWriter,
    const autoPtr<volPointInterpolation>& pInterp,
    UPtrList<vtk::patchWriter>& patchWriters,
    const UPtrList<PrimitivePatchInterpolation<primitivePatch>>& patchInterps,
    const fvMeshSubset& proxy,
    const wordHashSet& candidateNames
) const
{
    const fvMesh& baseMesh = proxy.baseMesh();

    label count = 0;

    for
    (
        const GeoField& origField
      : baseMesh.csorted<GeoField>(candidateNames)
    )
    {
        bool ok = false;
        auto tfield = fvMeshSubsetProxy::interpolate(proxy, origField);
        const auto& field = tfield();

        // Internal
        if (internalWriter && pInterp)
        {
            ok = true;
            internalWriter->write(field, *pInterp);
        }

        // Boundary
        label writeri = 0;
        for (vtk::patchWriter& writer : patchWriters)
        {
            if (patchInterps.test(writeri))
            {
                ok = true;
                writer.write(field, patchInterps[writeri]);
            }
            ++writeri;
        }

        if (ok)
        {
            ++count;

            if (verbose_)
            {
                if (count == 1)
                {
                    Log << "    " << GeoField::typeName << "->point(";
                }
                else
                {
                    Log << ' ';
                }
                Log << origField.name();
            }
        }
    }

    if (verbose_ && count)
    {
        Log << ')' << endl;
    }

    return count;
}


// ************************************************************************* //
