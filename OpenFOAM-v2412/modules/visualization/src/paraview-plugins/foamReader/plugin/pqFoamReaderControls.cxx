/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenCFD Ltd.
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

#include "pqFoamReaderControls.h"

#include <QCheckBox>
#include <QComboBox>
#include <QFrame>
#include <QGridLayout>
#include <QPushButton>

#include "pqPVApplicationCore.h"
#include "pqPipelineRepresentation.h"
#include "pqView.h"
#include "vtkSMDocumentation.h"
#include "vtkSMIntVectorProperty.h"
#include "vtkSMPropertyGroup.h"
#include "vtkSMSourceProxy.h"
#include "vtkSMEnumerationDomain.h"


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

// file-scope
// Add horizontal divider to layout
static void addHline(QGridLayout* layout, int row, int nCols)
{
    QFrame* hline = new QFrame(layout->parentWidget());
    hline->setFrameStyle(QFrame::HLine | QFrame::Sunken);

    layout->addWidget(hline, row, 0, 1, nCols);
}


// file-scope
// Widget properties
static QWidget* setWidgetProperties
(
    QWidget* widget,
    vtkSMProperty* prop
)
{
    widget->setFocusPolicy(Qt::NoFocus); // avoid dotted border

    vtkSMDocumentation* doc = prop->GetDocumentation();
    if (doc)
    {
        const char* txt = doc->GetDescription();
        if (txt)
        {
            QString tip = QString(txt).simplified();
            if (tip.size())
            {
               widget->setToolTip(tip);
            }
        }
    }

    return widget;
}


// file-scope
// Button properties
static QAbstractButton* setButtonProperties
(
    QAbstractButton* b,
    vtkSMProperty* prop
)
{
    setWidgetProperties(b, prop);
    b->setText(prop->GetXMLLabel());

    vtkSMIntVectorProperty* intProp =
       vtkSMIntVectorProperty::SafeDownCast(prop);

    // Initial checked state for integer (bool) properties
    if (intProp)
    {
        b->setChecked(intProp->GetElement(0));
    }

    return b;
}


// file-scope
// Fill combo-box from XML enumeration
static QComboBox* setComboBoxContent
(
    QComboBox* b,
    vtkSMIntVectorProperty* prop
)
{
    vtkSMEnumerationDomain* propEnum =
        vtkSMEnumerationDomain::SafeDownCast
        (
            prop->FindDomain("vtkSMEnumerationDomain")
        );

    if (propEnum)
    {
        unsigned int n = propEnum->GetNumberOfEntries();
        for (unsigned int idx=0; idx < n; ++idx)
        {
            const int val   = propEnum->GetEntryValue(idx);
            const char* txt = propEnum->GetEntryText(idx);

            b->insertItem(val, txt);
        }

        // Set default
        const int val = prop->GetElement(0);
        unsigned int idx = 0;
        if (!propEnum->IsInDomain(val, idx))
        {
            idx = 0;
        }
        b->setCurrentIndex(idx);
    }

    return b;
}


// file-scope
// Translate a combo-box index to a lookup value
static int comboBoxValue(vtkSMIntVectorProperty* prop, int idx)
{
    vtkSMEnumerationDomain* propEnum =
        vtkSMEnumerationDomain::SafeDownCast
        (
            prop->FindDomain("vtkSMEnumerationDomain")
        );

    if (propEnum)
    {
        return propEnum->GetEntryValue(idx);
    }
    else
    {
        return idx;
    }
}


static vtkSMIntVectorProperty* lookupIntProp
(
    vtkSMPropertyGroup* group,
    const char* name
)
{
    vtkSMProperty* prop = group->GetProperty(name);

    if (prop)
    {
        return vtkSMIntVectorProperty::SafeDownCast(prop);
    }

    return nullptr;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void pqFoamReaderControls::fireCommand
(
    vtkSMIntVectorProperty* prop,
    int val
)
{
    vtkSMProxy* pxy = this->proxy();

    prop->SetElement(0, val); // Set int value, toggle bool, etc

    // Fire off command
    prop->Modified();
    pxy->UpdateProperty(pxy->GetPropertyName(prop));
}


void pqFoamReaderControls::updateParts()
{
    vtkSMProxy* pxy = this->proxy();

    pxy->UpdatePropertyInformation(pxy->GetProperty("PartArrayStatus"));
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void pqFoamReaderControls::refreshPressed()
{
    vtkSMProxy* pxy = this->proxy();

    // Fire off command
    refresh_->Modified();
    pxy->UpdateProperty(pxy->GetPropertyName(refresh_));

    vtkSMSourceProxy::SafeDownCast(pxy)->UpdatePipeline();

    // Trigger a rendering (all views)
    pqPVApplicationCore::instance()->render();
}


void pqFoamReaderControls::cacheMesh(int idx)
{
    fireCommand(meshCaching_, comboBoxValue(meshCaching_, idx));
}


void pqFoamReaderControls::showPatchNames(bool checked)
{
    fireCommand(showPatchNames_, checked);

    // update the active view
    if (this->view())
    {
        this->view()->render();
    }
    // OR: update all views
    // pqPVApplicationCore::instance()->render();
}


void pqFoamReaderControls::showGroupsOnly(bool checked)
{
    fireCommand(showGroupsOnly_, checked);
    updateParts();
}


void pqFoamReaderControls::includeSets(bool checked)
{
    fireCommand(includeSets_, checked);
    updateParts();
}


void pqFoamReaderControls::includeZones(bool checked)
{
    fireCommand(includeZones_, checked);
    updateParts();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pqFoamReaderControls::pqFoamReaderControls
(
    vtkSMProxy* proxy,
    vtkSMPropertyGroup* group,
    QWidget* parent
)
:
    Superclass(proxy, parent),
    refresh_(group->GetProperty("Refresh")),
    showPatchNames_(lookupIntProp(group, "ShowPatchNames")),
    showGroupsOnly_(lookupIntProp(group, "ShowGroupsOnly")),
    includeSets_(lookupIntProp(group, "IncludeSets")),
    includeZones_(lookupIntProp(group, "IncludeZones")),
    meshCaching_(lookupIntProp(group, "MeshCaching"))
{
    typedef vtkSMIntVectorProperty intProp;

    QGridLayout* form = new QGridLayout(this);

    const int nCols = 3;

    // ROW
    // ~~~
    int row = 0;

    if (refresh_)
    {
        QPushButton* b = new QPushButton(this);
        setButtonProperties(b, refresh_);
        form->addWidget(b, row, 0, Qt::AlignLeft);

        connect(b, SIGNAL(clicked()), this, SLOT(refreshPressed()));
    }

    intProp* zeroTime = lookupIntProp(group, "ZeroTime");
    if (zeroTime)
    {
        QCheckBox* b = new QCheckBox(this);
        setButtonProperties(b, zeroTime);
        form->addWidget(b, row, 1, Qt::AlignLeft);

        addPropertyLink
        (
            b, "checked", SIGNAL(toggled(bool)), zeroTime
        );
    }

    // LINE
    addHline(form, ++row, nCols);

    // ROW
    // ~~~
    ++row;

    if (includeSets_)
    {
        QCheckBox* b = new QCheckBox(this);
        setButtonProperties(b, includeSets_);
        form->addWidget(b, row, 0, Qt::AlignLeft);

        addPropertyLink
        (
            b, "checked", SIGNAL(toggled(bool)), includeSets_
        );
        connect
        (
            b, SIGNAL(toggled(bool)), this, SLOT(includeSets(bool))
        );
    }

    if (showGroupsOnly_)
    {
        QCheckBox* b = new QCheckBox(this);
        setButtonProperties(b, showGroupsOnly_);
        form->addWidget(b, row, 1, Qt::AlignLeft);

        addPropertyLink
        (
            b, "checked", SIGNAL(toggled(bool)), showGroupsOnly_
        );
        connect
        (
            b, SIGNAL(toggled(bool)), this, SLOT(showGroupsOnly(bool))
        );
    }


    // ROW
    // ~~~
    ++row;

    if (includeZones_)
    {
        QCheckBox* b = new QCheckBox(this);
        setButtonProperties(b, includeZones_);
        form->addWidget(b, row, 0, Qt::AlignLeft);

        addPropertyLink
        (
            b, "checked", SIGNAL(toggled(bool)), includeZones_
        );
        connect
        (
            b, SIGNAL(toggled(bool)), this, SLOT(includeZones(bool))
        );
    }

    if (showPatchNames_)
    {
        QCheckBox* b = new QCheckBox(this);
        setButtonProperties(b, showPatchNames_);
        form->addWidget(b, row, 1, Qt::AlignLeft);

        connect
        (
            b, SIGNAL(toggled(bool)), this, SLOT(showPatchNames(bool))
        );
    }

    // LINE
    addHline(form, ++row, nCols);

    // ROW
    // ~~~
    ++row;

    intProp* interpolate = lookupIntProp(group, "InterpolateFields");
    if (interpolate)
    {
        QCheckBox* b = new QCheckBox(this);
        setButtonProperties(b, interpolate);
        form->addWidget(b, row, 0, Qt::AlignLeft);

        addPropertyLink
        (
            b, "checked", SIGNAL(toggled(bool)), interpolate
        );
    }

    intProp* extrapolate = lookupIntProp(group, "ExtrapolatePatches");
    if (extrapolate)
    {
        QCheckBox* b = new QCheckBox(this);
        setButtonProperties(b, extrapolate);
        form->addWidget(b, row, 1, Qt::AlignLeft);

        addPropertyLink
        (
            b, "checked", SIGNAL(toggled(bool)), extrapolate
        );
    }

    // LINE
    addHline(form, ++row, nCols);

    // ROW
    // ~~~
    ++row;

    intProp* updateGui = lookupIntProp(group, "UpdateGUI");
    if (updateGui)
    {
        QPushButton* b = new QPushButton(this);
        setButtonProperties(b, updateGui);
        form->addWidget(b, row, 0, Qt::AlignLeft);

        addPropertyLink
        (
            b, "checked", SIGNAL(clicked()), updateGui
        );
    }

    intProp* usePolyhedron = lookupIntProp(group, "UseVTKPolyhedron");
    if (usePolyhedron)
    {
        QCheckBox* b = new QCheckBox(this);
        setButtonProperties(b, usePolyhedron);
        form->addWidget(b, row, 1, Qt::AlignLeft);

        addPropertyLink
        (
            b, "checked", SIGNAL(toggled(bool)), usePolyhedron
        );
    }

    if (meshCaching_)
    {
        QComboBox* b = new QComboBox(this);
        form->addWidget(b, row, 2, Qt::AlignLeft);

        setWidgetProperties(b, meshCaching_);
        setComboBoxContent(b, meshCaching_);

        addPropertyLink
        (
            b, "indexChanged", SIGNAL(currentIndexChanged(int)), meshCaching_
        );

        connect
        (
            b, SIGNAL(currentIndexChanged(int)), this, SLOT(cacheMesh(int))
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

pqFoamReaderControls::~pqFoamReaderControls()
{}


// ************************************************************************* //
