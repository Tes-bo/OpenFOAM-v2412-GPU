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

#include "fieldExtents.H"
#include "processorPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(fieldExtents, 0);
    addToRunTimeSelectionTable(functionObject, fieldExtents, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::fieldExtents::writeFileHeader(Ostream& os)
{
    if (!fieldSet_.updateSelection())
    {
        return;
    }

    if (writtenHeader_)
    {
        writeBreak(os);
    }
    else
    {
        writeHeader(os, "Field extents");
        writeHeaderValue(os, "Reference position", C0_);
        writeHeaderValue(os, "Threshold", threshold_);
    }

    writeCommented(os, "Time");

    for (const word& fieldName : fieldSet_.selectionNames())
    {
        if (internalField_)
        {
            writeTabbed(os, fieldName + "_internal");
        }
        for (const label patchi : patchIDs_)
        {
            const word& patchName = mesh_.boundaryMesh()[patchi].name();
            writeTabbed(os, fieldName + "_" + patchName);
        }
    }

    os  << endl;

    writtenHeader_ = true;
}


template<>
Foam::tmp<Foam::volScalarField> Foam::functionObjects::fieldExtents::calcMask
(
    const GeometricField<scalar, fvPatchField, volMesh>& field
) const
{
    return
        pos
        (
            field
          - dimensionedScalar("t", field.dimensions(), threshold_)
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldExtents::fieldExtents
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name, typeName, dict),
    internalField_(true),
    threshold_(0),
    C0_(Zero),
    fieldSet_(mesh_)
{
    read(dict);

    // Note: delay creating the output file header to handle field names
    // specified using regular expressions
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fieldExtents::read(const dictionary& dict)
{
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    if (fvMeshFunctionObject::read(dict) && writeFile::read(dict))
    {
        dict.readIfPresent<bool>("internalField", internalField_);

        threshold_ = dict.get<scalar>("threshold");

        dict.readIfPresent<vector>("referencePosition", C0_);

        wordRes patchNames;
        if (dict.readIfPresent("patches", patchNames))
        {
            patchIDs_ = pbm.indices(patchNames);
        }
        else
        {
            labelHashSet patchSet(2*pbm.size());

            // All non-processor and non-empty patches
            forAll(pbm, patchi)
            {
                const polyPatch& pp = pbm[patchi];

                if
                (
                    !isA<processorPolyPatch>(pp)
                 && !isA<emptyPolyPatch>(pp)
                )
                {
                    patchSet.insert(patchi);
                }
            }

            patchIDs_ = patchSet.sortedToc();
        }

        if (!internalField_ && patchIDs_.empty())
        {
            IOWarningInFunction(dict)
                << "No internal field or patches selected - no field extent "
                << "information will be generated" << endl;
        }

        fieldSet_.read(dict);

        return true;
    }

    return false;
}


bool Foam::functionObjects::fieldExtents::execute()
{
    return true;
}


bool Foam::functionObjects::fieldExtents::write()
{
    writeFileHeader(file());

    Log << type() << " " << name() <<  " write:" << nl;

    for (const word& fieldName : fieldSet_.selectionNames())
    {
        calcFieldExtents<scalar>(fieldName, true);
        calcFieldExtents<vector>(fieldName);
        calcFieldExtents<sphericalTensor>(fieldName);
        calcFieldExtents<symmTensor>(fieldName);
        calcFieldExtents<tensor>(fieldName);
    }

    Log << endl;

    return true;
}


// ************************************************************************* //
