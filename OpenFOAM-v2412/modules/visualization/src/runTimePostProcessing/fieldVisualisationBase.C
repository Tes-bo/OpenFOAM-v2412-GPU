/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2021 OpenCFD Ltd.
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

// OpenFOAM includes
#include "fieldVisualisationBase.H"
#include "runTimePostProcessing.H"

#include "doubleVector.H"
#include "foamVtkTools.H"

// VTK includes
#include "vtkArrowSource.h"
#include "vtkCellDataToPointData.h"
#include "vtkCellData.h"
#include "vtkColorTransferFunction.h"
#include "vtkCompositeDataSet.h"
#include "vtkDataObjectTreeIterator.h"
#include "vtkFieldData.h"
#include "vtkGlyph3D.h"
#include "vtkLookupTable.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkRenderer.h"
#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::functionObjects::runTimePostPro::fieldVisualisationBase::
    colourByType
>
Foam::functionObjects::runTimePostPro::fieldVisualisationBase::
colourByTypeNames
({
    { colourByType::cbColour, "colour" },
    { colourByType::cbField, "field" },
});

const Foam::Enum
<
    Foam::functionObjects::runTimePostPro::fieldVisualisationBase::
    colourMapType
>
Foam::functionObjects::runTimePostPro::fieldVisualisationBase::
colourMapTypeNames
({
    { colourMapType::cmCoolToWarm, "coolToWarm" },
    { colourMapType::cmCoolToWarm, "blueWhiteRed" },
    { colourMapType::cmColdAndHot, "coldAndHot" },
    { colourMapType::cmFire, "fire" },
    { colourMapType::cmRainbow, "rainbow" },
    { colourMapType::cmGreyscale, "greyscale" },
    { colourMapType::cmGreyscale, "grayscale" },
    { colourMapType::cmXray, "xray" },
});


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::functionObjects::runTimePostPro::fieldVisualisationBase::fieldSummary
Foam::functionObjects::runTimePostPro::fieldVisualisationBase::
queryFieldSummary
(
    const word& fieldName,
    vtkDataSet* dataset
)
{
    fieldSummary queried;

    if (dataset)
    {
        auto* array = vtkDataArray::SafeDownCast
        (
            dataset->GetCellData()->GetAbstractArray(fieldName.c_str())
        );

        if (array)
        {
            queried.nComponents_ = array->GetNumberOfComponents();
            queried.association_ |= FieldAssociation::CELL_DATA;
            queried.range_ += vtk::Tools::rangeOf(array);
        }

        array = vtkDataArray::SafeDownCast
        (
            dataset->GetPointData()->GetAbstractArray(fieldName.c_str())
        );

        if (array)
        {
            queried.nComponents_ = array->GetNumberOfComponents();
            queried.association_ |= FieldAssociation::POINT_DATA;
            queried.range_ += vtk::Tools::rangeOf(array);
        }
    }

    return queried;
}


Foam::functionObjects::runTimePostPro::fieldVisualisationBase::fieldSummary
Foam::functionObjects::runTimePostPro::fieldVisualisationBase::
queryFieldSummary
(
    const word& fieldName,
    vtkCompositeDataSet* data
)
{
    fieldSummary queried;

    auto iter = vtkSmartPointer<vtkDataObjectTreeIterator>::New();

    iter->SetDataSet(data);
    iter->VisitOnlyLeavesOn();
    iter->SkipEmptyNodesOn();

    for
    (
        iter->InitTraversal();
        !iter->IsDoneWithTraversal();
        iter->GoToNextItem()
    )
    {
        auto* dataset = vtkDataSet::SafeDownCast
        (
            iter->GetCurrentDataObject()
        );

        if (dataset)
        {
            fieldSummary local(queryFieldSummary(fieldName, dataset));

            if (!queried.nComponents_)
            {
                queried.nComponents_ = local.nComponents_;
            }

            queried.association_ |= local.association_;
            queried.range_ += local.range_;
        }
    }

    return queried;
}


Foam::functionObjects::runTimePostPro::fieldVisualisationBase::FieldAssociation
Foam::functionObjects::runTimePostPro::fieldVisualisationBase::
queryFieldAssociation
(
    const word& fieldName,
    vtkDataSet* dataset
)
{
    unsigned where(FieldAssociation::NO_DATA);

    if (dataset)
    {
        if (dataset->GetCellData()->HasArray(fieldName.c_str()))
        {
            where |= FieldAssociation::CELL_DATA;
        }
        if (dataset->GetPointData()->HasArray(fieldName.c_str()))
        {
            where |= FieldAssociation::POINT_DATA;
        }
    }

    return FieldAssociation(where);
}


Foam::functionObjects::runTimePostPro::fieldVisualisationBase::FieldAssociation
Foam::functionObjects::runTimePostPro::fieldVisualisationBase::
queryFieldAssociation
(
    const word& fieldName,
    vtkCompositeDataSet* data
)
{
    unsigned where(FieldAssociation::NO_DATA);

    auto iter = vtkSmartPointer<vtkDataObjectTreeIterator>::New();

    iter->SetDataSet(data);
    iter->VisitOnlyLeavesOn();
    iter->SkipEmptyNodesOn();

    for
    (
        iter->InitTraversal();
        !iter->IsDoneWithTraversal();
        iter->GoToNextItem()
    )
    {
        auto* dataset = vtkDataSet::SafeDownCast
        (
            iter->GetCurrentDataObject()
        );

        where |= queryFieldAssociation(fieldName, dataset);
    }

    return FieldAssociation(where);
}


void Foam::functionObjects::runTimePostPro::fieldVisualisationBase::addMagField
(
    const word& fieldName,
    vtkFieldData* fieldData
)
{
    if (!fieldData)
    {
        return;
    }

    auto* input = vtkDataArray::SafeDownCast
    (
        fieldData->GetAbstractArray(fieldName.c_str())
    );

    if (!input)
    {
        return;
    }

    const word magFieldName = "mag(" + fieldName + ")";

    auto* output = vtkDataArray::SafeDownCast
    (
        fieldData->GetAbstractArray(magFieldName.c_str())
    );

    if (output)
    {
        return;
    }


    // Simplfy and only handle scalar/vector input

    const int nCmpt = input->GetNumberOfComponents();
    const vtkIdType len = input->GetNumberOfTuples();

    if (nCmpt == 1)
    {
        auto data = vtkSmartPointer<vtkFloatArray>::New();

        data->SetName(magFieldName.c_str());
        data->SetNumberOfComponents(1);
        data->SetNumberOfTuples(len);

        double scratch;
        for (vtkIdType i=0; i < len; ++i)
        {
            input->GetTuple(i, &scratch);

            scratch = Foam::mag(scratch);
            data->SetTuple(i, &scratch);
        }

        fieldData->AddArray(data);
    }
    else if (nCmpt == 3)
    {
        auto data = vtkSmartPointer<vtkFloatArray>::New();

        data->SetName(magFieldName.c_str());
        data->SetNumberOfComponents(1);
        data->SetNumberOfTuples(len);

        doubleVector scratch;
        for (vtkIdType i=0; i < len; ++i)
        {
            input->GetTuple(i, scratch.v_);

            scratch.x() = Foam::mag(scratch);

            data->SetTuple(i, scratch.v_);
        }

        fieldData->AddArray(data);
    }
}


void Foam::functionObjects::runTimePostPro::fieldVisualisationBase::addMagField
(
    const word& fieldName,
    vtkDataSet* dataset
)
{
    if (dataset)
    {
        addMagField(fieldName, dataset->GetCellData());
        addMagField(fieldName, dataset->GetPointData());
    }
}


void Foam::functionObjects::runTimePostPro::fieldVisualisationBase::addMagField
(
    const word& fieldName,
    vtkCompositeDataSet* data
)
{
    auto iter = vtkSmartPointer<vtkDataObjectTreeIterator>::New();

    iter->SetDataSet(data);
    iter->VisitOnlyLeavesOn();
    iter->SkipEmptyNodesOn();

    for
    (
        iter->InitTraversal();
        !iter->IsDoneWithTraversal();
        iter->GoToNextItem()
    )
    {
        auto* dataset = vtkDataSet::SafeDownCast
        (
            iter->GetCurrentDataObject()
        );
        addMagField(fieldName, dataset);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::runTimePostPro::fieldVisualisationBase::
fieldSummary::reduce()
{
    if (Pstream::parRun())
    {
        Foam::reduce(nComponents_, maxOp<int>());
        Foam::reduce(association_, bitOrOp<unsigned>());
        Foam::reduce(range_, minMaxOp<scalar>());
    }
}


Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const InfoProxy
    <
        functionObjects::runTimePostPro::fieldVisualisationBase::fieldSummary
    >& proxy
)
{
    os  << "nComponents:" << proxy.t_.nComponents_
        << " association:" << label(proxy.t_.association_)
        << " min/max:" << proxy.t_.range_;

    return os;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::runTimePostPro::fieldVisualisationBase::
setColourMap
(
    vtkLookupTable* lut
) const
{
    lut->SetNumberOfColors(nColours_);

    auto ctf = vtkSmartPointer<vtkColorTransferFunction>::New();

    switch (colourMap_)
    {
        case cmCoolToWarm:  // ParaView: "Cool To Warm"
        {
            ctf->SetColorSpaceToDiverging();
            ctf->AddRGBPoint(0.0, 0.231372, 0.298039, 0.752941);
            ctf->AddRGBPoint(0.5, 0.865003, 0.865003, 0.865003);
            ctf->AddRGBPoint(1.0, 0.705882, 0.0156863, 0.14902);
            // ctf->SetNanColor(1, 1, 0);
            break;
        }

        case cmColdAndHot:  // ParaView : "Cold and Hot"
        {
            ctf->SetColorSpaceToRGB();
            ctf->AddRGBPoint(0, 0, 1, 1);
            ctf->AddRGBPoint(0.45, 0, 0, 1);
            ctf->AddRGBPoint(0.5, 0, 0, 0.5019608);
            ctf->AddRGBPoint(0.55, 1, 0, 0);
            ctf->AddRGBPoint(1, 1, 1, 0);
            break;
        }

        case cmFire:  // ParaView: Black-Body Radiation
        {
            ctf->SetColorSpaceToRGB();
            ctf->AddRGBPoint(0, 0, 0, 0);
            ctf->AddRGBPoint(0.4, 0.901961, 0, 0);
            ctf->AddRGBPoint(0.8, 0.901961, 0.901961, 0);
            ctf->AddRGBPoint(1, 1, 1, 1);
            // ctf->SetNanColor(0, 0.49804, 1);
            break;
        }

        case cmRainbow:
        {
            ctf->SetColorSpaceToHSV();
            ctf->AddRGBPoint(0, 0, 0, 1);
            ctf->AddRGBPoint(0.5, 0, 1, 0);
            ctf->AddRGBPoint(1, 1, 0, 0);
            // ctf->SetNanColor(0.498039, 0.498039, 0.498039);
            break;
        }

        case cmGreyscale: // ParaView: grayscale
        {
            ctf->SetColorSpaceToRGB();
            ctf->AddRGBPoint(0, 0, 0, 0);
            ctf->AddRGBPoint(1, 1, 1, 1);
            // ctf->SetNanColor(1, 0, 0);
            break;
        }

        case cmXray: // ParaView: "X ray"
        {
            ctf->SetColorSpaceToRGB();
            ctf->AddRGBPoint(0, 1, 1, 1);
            ctf->AddRGBPoint(1, 0, 0, 0);
            // ctf->SetNanColor(1, 0, 0);
            break;
        }
    }


    double rgba[4] = { 0, 0, 0, 1 };
    for (int i = 0; i < nColours_; ++i)
    {
        ctf->GetColor(double(i)/double(nColours_), rgba);
        lut->SetTableValue(i, rgba);
    }
}


void Foam::functionObjects::runTimePostPro::fieldVisualisationBase::
addScalarBar
(
    const scalar position,
    vtkRenderer* renderer,
    vtkLookupTable* lut
) const
{
    // Add the scalar bar - only once!
    if (renderer && Pstream::master())
    {
        scalarBar_.add(colours_["text"]->value(position), renderer, lut);
    }
}


void Foam::functionObjects::runTimePostPro::fieldVisualisationBase::
setField
(
    const scalar position,
    const word& colourFieldName,
    const FieldAssociation fieldAssociation,
    vtkMapper* mapper,
    vtkRenderer* renderer
) const
{
    mapper->InterpolateScalarsBeforeMappingOn();

    switch (colourBy_)
    {
        case cbColour:
        {
            mapper->ScalarVisibilityOff();
            break;
        }

        case cbField:
        {
            // Create look-up table for colours
            auto lut = vtkSmartPointer<vtkLookupTable>::New();
            setColourMap(lut);
            lut->SetVectorMode(vtkScalarsToColors::MAGNITUDE);
            lut->SetTableRange(range_.first(), range_.second());

            // Configure the mapper
            const char* fieldName = colourFieldName.c_str();
            mapper->SelectColorArray(fieldName);

            // Use either point or cell data
            // - if both point and cell data exists, preferentially choose
            //   point data.  This is often the case when using glyphs.

            if (fieldAssociation & FieldAssociation::POINT_DATA)
            {
                mapper->SetScalarModeToUsePointFieldData();
            }
            else if (fieldAssociation & FieldAssociation::CELL_DATA)
            {
                mapper->SetScalarModeToUseCellFieldData();
            }
            else
            {
                WarningInFunction
                    << "Unable to determine cell or point data type "
                    << "- assuming point data";
                mapper->SetScalarModeToUsePointFieldData();
            }
            mapper->SetScalarRange(range_.first(), range_.second());
            mapper->SetColorModeToMapScalars();
            mapper->SetLookupTable(lut);
            mapper->ScalarVisibilityOn();

            // Add the scalar bar
            addScalarBar(position, renderer, lut);
            break;
        }
    }

    mapper->Modified();
}


void Foam::functionObjects::runTimePostPro::fieldVisualisationBase::
addGlyphs
(
    const scalar position,
    const word& scaleFieldName,
    const fieldSummary& scaleFieldInfo,
    const word& colourFieldName,
    const fieldSummary& colourFieldInfo,
    const scalar maxGlyphLength,

    vtkPolyData* data,
    vtkActor* actor,
    vtkRenderer* renderer
) const
{
    // Determine whether we have CellData/PointData and (scalar/vector)
    // or if we need to a cell->point data filter.

    if (!scaleFieldInfo.exists())
    {
        WarningInFunction
            << "Cannot add glyphs. No such cell or point field: "
            << scaleFieldName << endl;
        return;
    }

    if (!scaleFieldInfo.isScalar() && !scaleFieldInfo.isVector())
    {
        WarningInFunction
            << "Glyphs can only be added to scalar or vector data. "
            << "Unable to process field " << scaleFieldName << endl;
        return;
    }


    // Setup glyphs

    // The min/max data range for the input data (cell or point),
    // which will be slightly less after using a cell->point filter
    // (since it averages), but is still essentially OK.


    auto glyph = vtkSmartPointer<vtkGlyph3D>::New();
    glyph->ScalingOn();

    auto glyphMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    glyphMapper->SetInputConnection(glyph->GetOutputPort());

    vtkSmartPointer<vtkCellDataToPointData> cellToPoint;

    // The data source is filtered or original (PointData)
    if (!scaleFieldInfo.hasPointData() || !colourFieldInfo.hasPointData())
    {
        // CellData - Need a cell->point filter
        cellToPoint = vtkSmartPointer<vtkCellDataToPointData>::New();
        cellToPoint->SetInputData(data);

        glyph->SetInputConnection(cellToPoint->GetOutputPort());
    }
    else
    {
        glyph->SetInputData(data);
    }


    if (scaleFieldInfo.nComponents_ == 1)
    {
        auto sphere = vtkSmartPointer<vtkSphereSource>::New();
        sphere->SetCenter(0, 0, 0);
        sphere->SetRadius(0.5);

        // Setting higher resolution slows the rendering significantly
        // sphere->SetPhiResolution(20);
        // sphere->SetThetaResolution(20);

        glyph->SetSourceConnection(sphere->GetOutputPort());

        if (maxGlyphLength > 0)
        {
            // Using range from the data:
            // glyph->SetRange
            // (
            //     scaleFieldInfo.range_.first(),
            //     scaleFieldInfo.range_.second()
            // );

            // Set range according to user-supplied limits
            glyph->ClampingOn();
            glyph->SetRange(range_.first(), range_.second());

            // If range[0] != min(value), maxGlyphLength behaviour will not
            // be correct...
            glyph->SetScaleFactor(maxGlyphLength);
        }
        else
        {
            glyph->SetScaleFactor(1);
        }
        glyph->SetScaleModeToScaleByScalar();
        glyph->OrientOff();
        glyph->SetColorModeToColorByScalar();
        glyph->SetInputArrayToProcess
        (
            0, // index (0) = scalars
            0, // port
            0, // connection
            vtkDataObject::FIELD_ASSOCIATION_POINTS,
            scaleFieldName.c_str()
        );
    }
    else if (scaleFieldInfo.nComponents_ == 3)
    {
        auto arrow = vtkSmartPointer<vtkArrowSource>::New();
        arrow->SetTipResolution(10);
        arrow->SetTipRadius(0.1);
        arrow->SetTipLength(0.35);
        arrow->SetShaftResolution(10);
        arrow->SetShaftRadius(0.03);

        glyph->SetSourceConnection(arrow->GetOutputPort());

        if (maxGlyphLength > 0)
        {
            // Set range according data limits
            glyph->ClampingOn();
            glyph->SetRange
            (
                scaleFieldInfo.range_.first(),
                scaleFieldInfo.range_.second()
            );
            glyph->SetScaleFactor(maxGlyphLength);
        }
        else
        {
            glyph->SetScaleFactor(1);
        }
        glyph->SetScaleModeToScaleByVector();
        glyph->OrientOn();
        glyph->SetVectorModeToUseVector();
        glyph->SetColorModeToColorByVector();
        glyph->SetInputArrayToProcess
        (
            1, // index (1) = vectors
            0, // port
            0, // connection
            vtkDataObject::FIELD_ASSOCIATION_POINTS,
            scaleFieldName.c_str()
        );
    }


    // Apply colouring etc.
    // We already established PointData, which as either in the original,
    // or generated with vtkCellDataToPointData filter.

    {
        glyph->Update();

        setField
        (
            position,
            colourFieldName,
            FieldAssociation::POINT_DATA,  // Original or after filter
            glyphMapper,
            renderer
        );

        glyphMapper->Update();

        actor->SetMapper(glyphMapper);

        renderer->AddActor(actor);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::runTimePostPro::fieldVisualisationBase::
fieldVisualisationBase
(
    const dictionary& dict,
    const HashPtrTable<Function1<vector>>& colours
)
:
    colours_(colours),
    fieldName_(dict.get<word>("field")),
    smooth_(dict.getOrDefault("smooth", false)),
    colourBy_(cbColour),
    colourMap_(cmRainbow),
    nColours_(256),
    range_(),
    scalarBar_()
{
    colourByTypeNames.readEntry("colourBy", dict, colourBy_);

    switch (colourBy_)
    {
        case cbColour:
        {
            scalarBar_.hide();
            break;
        }

        case cbField:
        {
            dict.readEntry("range", range_);

            colourMapTypeNames.readIfPresent("colourMap", dict, colourMap_);

            dict.readIfPresent("nColours", nColours_);
            if (nColours_ < 2)
            {
                nColours_ = 2;
            }

            const dictionary* sbar = dict.findDict("scalarBar");
            if (sbar)
            {
                scalarBar_.read(*sbar);
            }
            else
            {
                scalarBar_.hide();
            }
            break;
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::runTimePostPro::fieldVisualisationBase::
~fieldVisualisationBase()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::HashPtrTable<Foam::Function1<Foam::vector>, Foam::word>&
Foam::functionObjects::runTimePostPro::fieldVisualisationBase::colours() const
{
    return colours_;
}


const Foam::word&
Foam::functionObjects::runTimePostPro::fieldVisualisationBase::fieldName()
const
{
    return fieldName_;
}


// ************************************************************************* //
