/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2024 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify i
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

#include "externalCoupled.H"
#include "OSspecific.H"
#include "Fstream.H"
#include "volFields.H"
#include "externalCoupledMixedFvPatchFields.H"
#include "mixedFvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "SpanStream.H"
#include "StringStream.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
bool Foam::functionObjects::externalCoupled::readData
(
    const UPtrList<const fvMesh>& meshes,
    const wordRe& groupName,
    const word& fieldName
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> volFieldType;
    typedef externalCoupledMixedFvPatchField<Type> patchFieldType;

    wordList regionNames(meshes.size());
    forAll(meshes, i)
    {
        regionNames[i] = meshes[i].dbDir();
    }

    // File only opened on master; contains data for all processors, for all
    // patchIDs.
    autoPtr<IFstream> masterFilePtr;
    if (UPstream::master())
    {
        const fileName transferFile
        (
            groupDir(commDirectory(), compositeName(regionNames), groupName)
          / fieldName + ".in"
        );

        Log << type() << ": reading data from " << transferFile << endl;

        masterFilePtr.reset(new IFstream(transferFile));

        if (!masterFilePtr().good())
        {
            FatalIOErrorInFunction(masterFilePtr())
                << "Cannot open file for region " << compositeName(regionNames)
                << ", field " << fieldName
                << exit(FatalIOError);
        }
    }


    const wordRes patchSelection(Foam::one{}, groupName);

    label nFound = 0;

    for (const fvMesh& mesh : meshes)
    {
        auto* vfptr = mesh.getObjectPtr<volFieldType>(fieldName);

        if (!vfptr)
        {
            continue;
        }
        ++nFound;

        auto& bf = vfptr->boundaryFieldRef();

        // Get the patches
        const labelList patchIDs
        (
            mesh.boundaryMesh().patchSet(patchSelection).sortedToc()
        );

        // Handle column-wise reading of patch data. Supports most easy types
        for (const label patchi : patchIDs)
        {
            if (isA<patchFieldType>(bf[patchi]))
            {
                // Explicit handling of externalCoupledMixed bcs - they
                // have specialised reading routines.

                patchFieldType& pf = refCast<patchFieldType>
                (
                    bf[patchi]
                );

                // Read from master into local stream
                std::string lines;
                readLines
                (
                    bf[patchi].size(),  // number of lines to read
                    masterFilePtr,
                    lines
                );

                ISpanStream isstr(lines);

                // Pass responsibility for all reading over to bc
                pf.readData(isstr);

                // Update the value from the read coefficient. Bypass any
                // additional processing by derived type.
                pf.patchFieldType::evaluate();
            }
            else if (isA<mixedFvPatchField<Type>>(bf[patchi]))
            {
                mixedFvPatchField<Type>& pf = refCast<mixedFvPatchField<Type>>
                (
                    bf[patchi]
                );

                // Read columns from file for
                // value, snGrad, refValue, refGrad, valueFraction
                List<scalarField> data;
                readColumns
                (
                    bf[patchi].size(),              // number of lines to read
                    4*pTraits<Type>::nComponents+1, // nColumns: 4*Type+1*scalar
                    masterFilePtr,
                    data
                );

                // Transfer read data to bc.
                // Skip value, snGrad
                direction columni = 2*pTraits<Type>::nComponents;

                Field<Type>& refValue = pf.refValue();
                for
                (
                    direction cmpt = 0;
                    cmpt < pTraits<Type>::nComponents;
                    cmpt++
                )
                {
                    refValue.replace(cmpt, data[columni++]);
                }
                Field<Type>& refGrad = pf.refGrad();
                for
                (
                    direction cmpt = 0;
                    cmpt < pTraits<Type>::nComponents;
                    cmpt++
                )
                {
                    refGrad.replace(cmpt, data[columni++]);
                }
                pf.valueFraction() = data[columni];

                // Update the value from the read coefficient. Bypass any
                // additional processing by derived type.
                pf.mixedFvPatchField<Type>::evaluate();
            }
            else if (isA<fixedGradientFvPatchField<Type>>(bf[patchi]))
            {
                fixedGradientFvPatchField<Type>& pf =
                    refCast<fixedGradientFvPatchField<Type>>(bf[patchi]);

                // Read columns for value and gradient
                List<scalarField> data;
                readColumns
                (
                    bf[patchi].size(),              // number of lines to read
                    2*pTraits<Type>::nComponents,   // nColumns: Type
                    masterFilePtr,
                    data
                );

                // Transfer gradient to bc
                Field<Type>& gradient = pf.gradient();
                for
                (
                    direction cmpt = 0;
                    cmpt < pTraits<Type>::nComponents;
                    cmpt++
                )
                {
                    gradient.replace
                    (
                        cmpt,
                        data[pTraits<Type>::nComponents+cmpt]
                    );
                }

                // Update the value from the read coefficient. Bypass any
                // additional processing by derived type.
                pf.fixedGradientFvPatchField<Type>::evaluate();
            }
            else if (isA<fixedValueFvPatchField<Type>>(bf[patchi]))
            {
                fixedValueFvPatchField<Type>& pf =
                    refCast<fixedValueFvPatchField<Type>>(bf[patchi]);

                // Read columns for value only
                List<scalarField> data;
                readColumns
                (
                    bf[patchi].size(),              // number of lines to read
                    pTraits<Type>::nComponents,     // number of columns to read
                    masterFilePtr,
                    data
                );

                // Transfer read value to bc
                Field<Type> value(bf[patchi].size());
                for
                (
                    direction cmpt = 0;
                    cmpt < pTraits<Type>::nComponents;
                    cmpt++
                )
                {
                    value.replace(cmpt, data[cmpt]);
                }

                pf == value;

                // Update the value from the read coefficient. Bypass any
                // additional processing by derived type.
                pf.fixedValueFvPatchField<Type>::evaluate();
            }
            else
            {
                FatalErrorInFunction
                    << "Unsupported boundary condition " << bf[patchi].type()
                    << " for patch " << bf[patchi].patch().name()
                    << " in region " << mesh.name()
                    << exit(FatalError);
            }

            initialisedCoupling_ = true;
        }
    }

    return nFound;
}


template<class Type>
bool Foam::functionObjects::externalCoupled::writeData
(
    const UPtrList<const fvMesh>& meshes,
    const wordRe& groupName,
    const word& fieldName
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> volFieldType;
    typedef externalCoupledMixedFvPatchField<Type> patchFieldType;

    wordList regionNames(meshes.size());
    forAll(meshes, i)
    {
        regionNames[i] = meshes[i].dbDir();
    }

    // File only opened on master; contains data for all processors, for all
    // patchIDs
    autoPtr<OFstream> masterFilePtr;
    if (Pstream::master())
    {
        const fileName transferFile
        (
            groupDir(commDirectory(), compositeName(regionNames), groupName)
          / fieldName + ".out"
        );

        Log << type() << ": writing data to " << transferFile << endl;

        masterFilePtr.reset(new OFstream(transferFile));

        if (!masterFilePtr().good())
        {
            FatalIOErrorInFunction(masterFilePtr())
                << "Cannot open file for region " << compositeName(regionNames)
                << ", field " << fieldName
                << exit(FatalIOError);
        }
    }


    const wordRes patchSelection(Foam::one{}, groupName);

    bool headerDone = false;
    label nFound = 0;

    for (const fvMesh& mesh : meshes)
    {
        const auto* vfptr = mesh.getObjectPtr<volFieldType>(fieldName);

        if (!vfptr)
        {
            continue;
        }
        ++nFound;

        const auto& bf = vfptr->boundaryField();

        // Get the patches
        const labelList patchIDs
        (
            mesh.boundaryMesh().patchSet(patchSelection).sortedToc()
        );

        // Handle column-wise writing of patch data. Supports most easy types
        for (const label patchi : patchIDs)
        {
            // const globalIndex globalFaces(bf[patchi].size());

            if (isA<patchFieldType>(bf[patchi]))
            {
                // Explicit handling of externalCoupledMixed bcs - they
                // have specialised writing routines

                const patchFieldType& pf = refCast<const patchFieldType>
                (
                    bf[patchi]
                );
                OStringStream os;

                // Pass responsibility for all writing over to bc
                pf.writeData(os);

                // Collect contributions from all processors and output them on
                // master
                if (UPstream::master())
                {
                    // Output master data first
                    if (!headerDone)
                    {
                        pf.writeHeader(masterFilePtr());
                        headerDone = true;
                    }
                    masterFilePtr() << os.str().c_str();

                    for (const int proci : UPstream::subProcs())
                    {
                        string str;
                        IPstream::recv(str, proci);

                        masterFilePtr() << str.c_str();
                    }
                }
                else
                {
                    OPstream::send(os.str(), UPstream::masterNo());
                }
            }
            else if (isA<mixedFvPatchField<Type>>(bf[patchi]))
            {
                const mixedFvPatchField<Type>& pf =
                    refCast<const mixedFvPatchField<Type>>(bf[patchi]);

                const globalIndex glob
                (
                    globalIndex::gatherOnly{},
                    pf.size()
                );

                Field<Type> value(glob.gather(pf));
                Field<Type> snGrad(glob.gather(pf.snGrad()()));
                Field<Type> refValue(glob.gather(pf.refValue()));
                Field<Type> refGrad(glob.gather(pf.refGrad()));
                scalarField valueFraction(glob.gather(pf.valueFraction()));

                if (UPstream::master())
                {
                    forAll(value, facei)
                    {
                        masterFilePtr()
                            << value[facei] << token::SPACE
                            << snGrad[facei] << token::SPACE
                            << refValue[facei] << token::SPACE
                            << refGrad[facei] << token::SPACE
                            << valueFraction[facei] << nl;
                    }
                }
            }
            else
            {
                // Output the value and snGrad

                const globalIndex glob
                (
                    globalIndex::gatherOnly{},
                    bf[patchi].size()
                );

                Field<Type> value(glob.gather(bf[patchi]));
                Field<Type> snGrad(glob.gather(bf[patchi].snGrad()()));

                if (UPstream::master())
                {
                    forAll(value, facei)
                    {
                        masterFilePtr()
                            << value[facei] << token::SPACE
                            << snGrad[facei] << nl;
                    }
                }
            }
        }
    }

    return nFound;
}


// ************************************************************************* //
