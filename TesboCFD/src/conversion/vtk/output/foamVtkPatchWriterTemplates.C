/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2021 OpenCFD Ltd.
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

#include "foamVtkPatchWriter.H"
#include "foamVtkOutput.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField>
void Foam::vtk::patchWriter::write
(
    const GeometricField<Type, PatchField, pointMesh>& field
)
{
    if (isState(outputState::POINT_DATA))
    {
        ++nPointData_;
    }
    else
    {
        reportBadState(FatalErrorInFunction, outputState::POINT_DATA)
            << " for field " << field.name() << nl << endl
            << exit(FatalError);
    }


    label nPoints = nLocalPoints_;

    if (parallel_)
    {
        reduce(nPoints, sumOp<label>());
    }


    this->beginDataArray<Type>(field.name(), nPoints);

    if (parallel_ ? Pstream::master() : true)
    {
        for (const label patchId : patchIDs_)
        {
            const auto& pfld = field.boundaryField()[patchId];

            // Only valuePointPatchField is actually derived from Field
            const auto* vpp = isA<Field<Type>>(pfld);
            if (vpp)
            {
                vtk::writeList(format(), *vpp);
            }
            else
            {
                vtk::writeList(format(), pfld.patchInternalField()());
            }
        }
    }


    if (parallel_)
    {
        // Patch Ids are identical across all processes
        const label nPatches = patchIDs_.size();

        if (Pstream::master())
        {
            Field<Type> recv;

            // Receive each patch field and write
            for (const int subproci : Pstream::subProcs())
            {
                IPstream fromProc(UPstream::commsTypes::scheduled, subproci);

                for (label i=0; i < nPatches; ++i)
                {
                    fromProc >> recv;

                    vtk::writeList(format(), recv);
                }
            }
        }
        else
        {
            // Send each patch field
            OPstream toProc
            (
                UPstream::commsTypes::scheduled,
                UPstream::masterNo()
            );

            for (const label patchId : patchIDs_)
            {
                const auto& pfld = field.boundaryField()[patchId];

                // Only valuePointPatchField is actually derived from Field
                const auto* vpp = isA<Field<Type>>(pfld);
                if (vpp)
                {
                    toProc << *vpp;
                }
                else
                {
                    toProc << pfld.patchInternalField();
                }
            }
        }
    }


    this->endDataArray();
}


template<class Type, template<class> class PatchField>
void Foam::vtk::patchWriter::write
(
    const GeometricField<Type, PatchField, volMesh>& field
)
{
    if (isState(outputState::CELL_DATA))
    {
        ++nCellData_;
    }
    else
    {
        reportBadState(FatalErrorInFunction, outputState::CELL_DATA)
            << " for field " << field.name() << nl << endl
            << exit(FatalError);
    }

    label nFaces = nLocalPolys_;

    if (parallel_)
    {
        reduce(nFaces, sumOp<label>());
    }


    this->beginDataArray<Type>(field.name(), nFaces);

    if (parallel_ ? Pstream::master() : true)
    {
        for (const label patchId : patchIDs_)
        {
            const auto& pfld = field.boundaryField()[patchId];

            if (useNearCellValue_)
            {
                vtk::writeList(format(), pfld.patchInternalField()());
            }
            else
            {
                vtk::writeList(format(), pfld);
            }
        }
    }

    if (parallel_)
    {
        // Patch Ids are identical across all processes
        const label nPatches = patchIDs_.size();

        if (Pstream::master())
        {
            Field<Type> recv;

            // Receive each patch field and write
            for (const int subproci : Pstream::subProcs())
            {
                IPstream fromProc(UPstream::commsTypes::scheduled, subproci);

                for (label i=0; i < nPatches; ++i)
                {
                    fromProc >> recv;

                    vtk::writeList(format(), recv);
                }
            }
        }
        else
        {
            // Send each patch field
            OPstream toProc
            (
                UPstream::commsTypes::scheduled,
                UPstream::masterNo()
            );

            for (const label patchId : patchIDs_)
            {
                const auto& pfld = field.boundaryField()[patchId];

                if (useNearCellValue_)
                {
                    toProc << pfld.patchInternalField()();
                }
                else
                {
                    toProc << static_cast<const Field<Type>&>(pfld);
                }
            }
        }
    }


    this->endDataArray();
}


template<class Type>
void Foam::vtk::patchWriter::write
(
    const GeometricField<Type, fvPatchField, volMesh>& field,
    const PrimitivePatchInterpolation<primitivePatch>& pInter
)
{
    if (isState(outputState::POINT_DATA))
    {
        ++nPointData_;
    }
    else
    {
        reportBadState(FatalErrorInFunction, outputState::POINT_DATA)
            << " for field " << field.name() << nl << endl
            << exit(FatalError);
    }

    label nPoints = nLocalPoints_;

    if (parallel_)
    {
        reduce(nPoints, sumOp<label>());
    }


    this->beginDataArray<Type>(field.name(), nPoints);

    if (parallel_ ? Pstream::master() : true)
    {
        for (const label patchId : patchIDs_)
        {
            const auto& pfld = field.boundaryField()[patchId];

            if (useNearCellValue_)
            {
                auto tfield =
                    pInter.faceToPointInterpolate
                    (
                        pfld.patchInternalField()()
                    );

                vtk::writeList(format(), tfield());
            }
            else
            {
                auto tfield = pInter.faceToPointInterpolate(pfld);

                vtk::writeList(format(), tfield());
            }
        }
    }


    if (parallel_)
    {
        // Patch Ids are identical across all processes
        const label nPatches = patchIDs_.size();

        if (UPstream::master())
        {
            Field<Type> recv;

            // Receive each patch field and write
            for (const int subproci : UPstream::subProcs())
            {
                IPstream fromProc(UPstream::commsTypes::scheduled, subproci);

                for (label i=0; i < nPatches; ++i)
                {
                    fromProc >> recv;

                    vtk::writeList(format(), recv);
                }
            }
        }
        else
        {
            // Send each patch field
            OPstream toProc
            (
                UPstream::commsTypes::scheduled,
                UPstream::masterNo()
            );

            for (const label patchId : patchIDs_)
            {
                const auto& pfld = field.boundaryField()[patchId];

                if (useNearCellValue_)
                {
                    auto tfield =
                        pInter.faceToPointInterpolate
                        (
                            pfld.patchInternalField()()
                        );

                    toProc << tfield();
                }
                else
                {
                    auto tfield = pInter.faceToPointInterpolate(pfld);

                    toProc << tfield();
                }
            }
        }
    }


    this->endDataArray();
}


// ************************************************************************* //
