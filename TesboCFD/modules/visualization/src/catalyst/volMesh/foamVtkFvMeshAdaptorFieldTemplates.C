/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2023 OpenCFD Ltd.
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

#ifndef foamVtkFvMeshAdaptorFieldTemplates_C
#define foamVtkFvMeshAdaptorFieldTemplates_C

// OpenFOAM includes
#include "error.H"
#include "emptyFvPatchField.H"
#include "wallPolyPatch.H"
#include "volPointInterpolation.H"
#include "zeroGradientFvPatchField.H"
#include "interpolatePointToCell.H"

// vtk includes
#include "vtkFloatArray.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
#include "vtkSmartPointer.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//
// volume-fields
//

template<class Type>
void Foam::vtk::fvMeshAdaptor::convertVolField
(
    const PtrList<patchInterpolator>& patchInterpList,
    const GeometricField<Type, fvPatchField, volMesh>& fld
)
{
    const fvMesh& mesh = fld.mesh();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // Interpolated field (demand driven)
    autoPtr<GeometricField<Type, pointPatchField, pointMesh>> ptfPtr;
    if (interpFields_)
    {
        ptfPtr.reset
        (
            volPointInterpolation::New(mesh).interpolate(fld).ptr()
        );
    }

    // INTERNAL
    convertVolFieldInternal(fld, ptfPtr);

    // BOUNDARY
    for (const label patchId : patchIds_)
    {
        const polyPatch& pp = patches[patchId];
        const word& longName = pp.name();

        auto iter = cachedVtp_.find(longName);
        if (!iter.good() || !iter.val().dataset)
        {
            // Should not happen, but for safety require a vtk geometry
            Pout<<"Cache miss for VTP patch " << longName << endl;
            continue;
        }

        foamVtpData& vtpData = iter.val();
        auto dataset = vtpData.dataset;

        // This is slightly roundabout, but we deal with groups and with
        // single patches.
        // For groups (spanning several patches) it is fairly messy to
        // get interpolated point fields. We would need to create a indirect
        // patch each time to obtain the mesh points. We thus skip that.
        //
        // To improve code reuse, we allocate the CellData as a zeroed-field
        // ahead of time.

        vtkSmartPointer<vtkFloatArray> cdata = vtk::Tools::zeroField<Type>
        (
            fld.name(),
            dataset->GetNumberOfPolys()
        );

        vtkSmartPointer<vtkFloatArray> pdata;

        const fvPatchField<Type>& ptf = fld.boundaryField()[patchId];

        if
        (
            isType<emptyFvPatchField<Type>>(ptf)
         ||
            (
                extrapPatches_
             && !polyPatch::constraintType(patches[patchId].type())
            )
        )
        {
            fvPatch p(ptf.patch().patch(), mesh.boundary());

            tmp<Field<Type>> tpptf
            (
                fvPatchField<Type>(p, fld).patchInternalField()
            );

            vtk::Tools::transcribeFloatData(cdata, tpptf());

            if
            (
                patchId < patchInterpList.size()
             && patchInterpList.set(patchId)
            )
            {
                pdata = vtk::Tools::convertFieldToVTK
                (
                    fld.name(),
                    patchInterpList[patchId].faceToPointInterpolate(tpptf)()
                );
            }
        }
        else
        {
            vtk::Tools::transcribeFloatData(cdata, ptf);

            if
            (
                patchId < patchInterpList.size()
             && patchInterpList.set(patchId)
            )
            {
                pdata = vtk::Tools::convertFieldToVTK
                (
                    fld.name(),
                    patchInterpList[patchId].faceToPointInterpolate(ptf)()
                );
            }
        }

        if (cdata)
        {
            dataset->GetCellData()->AddArray(cdata);
        }
        if (pdata)
        {
            dataset->GetPointData()->AddArray(pdata);
        }
    }
}


template<class Type>
void Foam::vtk::fvMeshAdaptor::convertVolFields
(
    const PtrList<patchInterpolator>& patchInterpList,
    const wordRes& selectFields
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    // Restrict to GeometricField<Type, ...>

    for (const auto& fld : mesh_.thisDb().csorted<fieldType>(selectFields))
    {
        convertVolField(patchInterpList, fld);
    }
}


template<class Type>
void Foam::vtk::fvMeshAdaptor::convertVolFieldInternal
(
    const GeometricField<Type, fvPatchField, volMesh>& fld,
    autoPtr<GeometricField<Type, pointPatchField, pointMesh>>& ptfPtr
)
{
    if (!usingInternal())
    {
        return;
    }

    const auto& longName = internalName();

    auto iter = cachedVtu_.find(longName);
    if (!iter.good() || !iter.val().dataset)
    {
        // Should not happen, but for safety require a vtk geometry
        Pout<<"Cache miss for VTU " << longName << endl;
        return;
    }
    foamVtuData& vtuData = iter.val();
    auto dataset = vtuData.dataset;

    dataset->GetCellData()->AddArray
    (
        vtuData.convertField(fld)
    );

    if (ptfPtr)
    {
        dataset->GetPointData()->AddArray
        (
            convertPointField(ptfPtr(), fld, vtuData)
        );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
vtkSmartPointer<vtkFloatArray> Foam::vtk::fvMeshAdaptor::convertPointField
(
    const GeometricField<Type, pointPatchField, pointMesh>& pfld,
    const GeometricField<Type, fvPatchField, volMesh>& vfld,
    const foamVtuData& vtuData
)
{
    const int nComp(pTraits<Type>::nComponents);
    const labelUList& addPointCellLabels = vtuData.additionalIds();
    const labelUList& pointMap = vtuData.pointMap();

    // Use a pointMap or address directly into mesh
    const label nPoints = (pointMap.size() ? pointMap.size() : pfld.size());

    auto data = vtkSmartPointer<vtkFloatArray>::New();
    data->SetNumberOfComponents(nComp);
    data->SetNumberOfTuples(nPoints + addPointCellLabels.size());

    // Note: using the name of the original volField
    // not the name generated by the interpolation "volPointInterpolate(<name>)"

    if (notNull(vfld))
    {
        data->SetName(vfld.name().c_str());
    }
    else
    {
        data->SetName(pfld.name().c_str());
    }

    DebugInfo
        << "convert Point field: " << pfld.name()
        << " size="  << (nPoints + addPointCellLabels.size())
        << " (" << nPoints << " + " << addPointCellLabels.size()
        << ") nComp=" << nComp << endl;


    float scratch[pTraits<Type>::nComponents];

    vtkIdType pointi = 0;
    if (pointMap.size())
    {
        for (const label meshPointi : pointMap)
        {
            vtk::Tools::foamToVtkTuple(scratch, pfld[meshPointi]);
            data->SetTuple(pointi++, scratch);
        }
    }
    else
    {
        for (const Type& val : pfld)
        {
            vtk::Tools::foamToVtkTuple(scratch, val);
            data->SetTuple(pointi++, scratch);
        }
    }

    // Continue with additional points
    // - correspond to cell centres

    if (notNull(vfld))
    {
        for (const label meshCelli : addPointCellLabels)
        {
            vtk::Tools::foamToVtkTuple(scratch, vfld[meshCelli]);
            data->SetTuple(pointi++, scratch);
        }
    }
    else
    {
        for (const label meshCelli : addPointCellLabels)
        {
            vtk::Tools::foamToVtkTuple
            (
                scratch, interpolatePointToCell(pfld, meshCelli)
            );
            data->SetTuple(pointi++, scratch);
        }
    }

    return data;
}


#endif

// ************************************************************************* //
