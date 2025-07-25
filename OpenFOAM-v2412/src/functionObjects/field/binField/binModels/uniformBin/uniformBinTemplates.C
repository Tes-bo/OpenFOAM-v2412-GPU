/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2022 OpenCFD Ltd.
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

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::binModels::uniformBin::writeFileHeader
(
    OFstream& os
) const
{
    writeHeader(os, "bins");

    const tensor& R = coordSysPtr_->R();
    for (direction i = 0; i < vector::nComponents; ++i)
    {
        writeHeaderValue(os, "e" + Foam::name(i) + " bins", nBins_[i]);
        writeHeaderValue(os, "    start", binLimits_.min()[i]);
        writeHeaderValue(os, "    end", binLimits_.max()[i]);
        writeHeaderValue(os, "    delta", binWidth_[i]);
        writeHeaderValue(os, "    direction", R.col(i));
    }
    writeCommented(os, "bin end co-ordinates:");
    os  << nl;

    // Compute and print bin end points in binning directions
    for (direction i = 0; i < vector::nComponents; ++i)
    {
        scalar binEnd = binLimits_.min()[i];

        writeCommented(os, "e"+Foam::name(i)+" co-ords   :");
        for (label j = 0; j < nBins_[i]; ++j)
        {
            binEnd += binWidth_[i];
            os  << tab << binEnd;
        }
        os  << nl;
    }

    writeHeader(os, "");
    writeCommented(os, "Time");

    for (label i = 0; i < nBin_; ++i)
    {
        const word ibin(Foam::name(i) + ':');
        writeTabbed(os, writeComponents<Type>("total" + ibin));
        writeTabbed(os, writeComponents<Type>("internal" + ibin));

        if (decomposePatchValues_)
        {
            writeTabbed(os, writeComponents<Type>("normal" + ibin));
            writeTabbed(os, writeComponents<Type>("tangential" + ibin));
        }
        else
        {
            writeTabbed(os, writeComponents<Type>("patch" + ibin));
        }
    }

    os  << endl;
}

template<class Type>
bool Foam::binModels::uniformBin::processField(const label fieldi)
{
    const word& fieldName = fieldNames_[fieldi];

    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;

    const VolFieldType* fieldPtr = mesh_.findObject<VolFieldType>(fieldName);

    if (!fieldPtr)
    {
        return false;
    }

    if (writeToFile() && !writtenHeader_)
    {
        writeFileHeader<Type>(filePtrs_[fieldi]);
    }

    const VolFieldType& fld = *fieldPtr;

    // Total number of fields
    //
    // 0: internal
    // 1: patch total
    //
    // OR
    //
    // 0: internal
    // 1: patch normal
    // 2: patch tangential
    label nField = 2;
    if (decomposePatchValues_)
    {
        nField += 1;
    }

    List<List<Type>> data(nField);
    for (auto& binList : data)
    {
        binList.resize(nBin_, Zero);
    }

    for (const label zonei : cellZoneIDs_)
    {
        const cellZone& cZone = mesh_.cellZones()[zonei];

        for (const label celli : cZone)
        {
            const label bini = cellToBin_[celli];

            if (bini >= 0)
            {
                data[0][bini] += fld[celli];
            }
        }
    }

    for (const label patchi : patchIDs_)
    {
        const polyPatch& pp = mesh_.boundaryMesh()[patchi];
        const vectorField np(mesh_.boundary()[patchi].nf());

        forAll(pp, facei)
        {
            const label localFacei =
                pp.start() - mesh_.nInternalFaces() + facei;
            const label bini = faceToBin_[localFacei];

            if (bini >= 0)
            {
                const Type& v = fld.boundaryField()[patchi][facei];

                if (!decomposePatchValues(data, bini, v, np[facei]))
                {
                    data[1][bini] += v;
                }
            }
        }
    }

    for (auto& binList : data)
    {
        reduce(binList, sumOp<List<Type>>());
    }

    if (writeToFile())
    {
        writeBinnedData(data, filePtrs_[fieldi]);
    }

    return true;
}


// ************************************************************************* //
