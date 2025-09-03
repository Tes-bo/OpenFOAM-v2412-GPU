/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2021 OpenCFD Ltd.
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

#include "pqFoamBlockMeshControls.h"

#include <QCheckBox>
#include <QGridLayout>
#include <QPushButton>

#include "pqPVApplicationCore.h"
#include "pqView.h"
#include "vtkSMDocumentation.h"
#include "vtkSMIntVectorProperty.h"
#include "vtkSMPropertyGroup.h"
#include "vtkSMSourceProxy.h"


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

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

    auto* intProp = vtkSMIntVectorProperty::SafeDownCast(prop);

    // Initial checked state for integer (bool) properties
    if (intProp)
    {
        b->setChecked(intProp->GetElement(0));
    }

    return b;
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

void pqFoamBlockMeshControls::fireCommand(vtkSMProperty* prop)
{
    vtkSMProxy* pxy = this->proxy();

    // Fire off command
    prop->Modified();
    pxy->UpdateProperty(pxy->GetPropertyName(prop));
}


void pqFoamBlockMeshControls::fireCommand
(
    vtkSMIntVectorProperty* prop,
    int val
)
{
    vtkSMProxy* pxy = this->proxy();

    prop->SetElement(0, val); // Set int value, toogle bool, etc

    // Fire off command
    prop->Modified();
    pxy->UpdateProperty(pxy->GetPropertyName(prop));
}


void pqFoamBlockMeshControls::updateParts()
{
    vtkSMProxy* pxy = this->proxy();

    pxy->UpdatePropertyInformation(pxy->GetProperty("BlockArrayStatus"));
    pxy->UpdatePropertyInformation(pxy->GetProperty("CurvedEdgesArrayStatus"));
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void pqFoamBlockMeshControls::refreshPressed()
{
    fireCommand(refresh_);

    vtkSMSourceProxy::SafeDownCast(this->proxy())->UpdatePipeline();

    // Trigger a rendering (all views)
    pqPVApplicationCore::instance()->render();

    updateParts();
}


void pqFoamBlockMeshControls::showPatchNames(bool checked)
{
    fireCommand(showPatchNames_, checked);

    // Update the active view
    if (this->view())
    {
        this->view()->render();
    }

    // OR: update all views
    // pqPVApplicationCore::instance()->render();
}


void pqFoamBlockMeshControls::showPointNumbers(bool checked)
{
    fireCommand(showPointNumbers_, checked);

    // Update the active view
    if (this->view())
    {
        this->view()->render();
    }

    // OR: update all views
    // pqPVApplicationCore::instance()->render();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pqFoamBlockMeshControls::pqFoamBlockMeshControls
(
    vtkSMProxy* proxy,
    vtkSMPropertyGroup* group,
    QWidget* parent
)
:
    Superclass(proxy, parent),
    refresh_(group->GetProperty("Refresh")),
    showPatchNames_(lookupIntProp(group, "ShowPatchNames")),
    showPointNumbers_(lookupIntProp(group, "ShowPointNumbers"))
{
    QGridLayout* form = new QGridLayout(this);

    if (refresh_)
    {
        QPushButton* b = new QPushButton(this);
        setButtonProperties(b, refresh_);
        form->addWidget(b, 0, 0, Qt::AlignLeft);

        connect
        (
            b, SIGNAL(clicked()), this, SLOT(refreshPressed())
        );
    }

    if (showPatchNames_)
    {
        QCheckBox* b = new QCheckBox(this);
        setButtonProperties(b, showPatchNames_);
        form->addWidget(b, 0, 1, Qt::AlignLeft);

        connect
        (
            b, SIGNAL(toggled(bool)), this, SLOT(showPatchNames(bool))
        );
    }

    if (showPointNumbers_)
    {
        QCheckBox* b = new QCheckBox(this);
        setButtonProperties(b, showPointNumbers_);
        form->addWidget(b, 0, 2, Qt::AlignLeft);

        connect
        (
            b, SIGNAL(toggled(bool)), this, SLOT(showPointNumbers(bool))
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

pqFoamBlockMeshControls::~pqFoamBlockMeshControls()
{}


// ************************************************************************* //
