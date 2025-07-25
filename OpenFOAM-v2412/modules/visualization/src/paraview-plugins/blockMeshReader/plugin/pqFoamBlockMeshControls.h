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

Class
    pqFoamBlockMeshControls

Description
    Customized property controls for the ParaView blockMesh reader.

    Refresh, ShowPatchNames, ShowPointNumbers.

SourceFiles
    pqFoamBlockMeshControls.cxx

\*---------------------------------------------------------------------------*/
#ifndef pqFoamBlockMeshControls_h
#define pqFoamBlockMeshControls_h

#include "pqPropertyWidget.h"

// Forward Declarations (ParaView)
class vtkSMProperty;
class vtkSMIntVectorProperty;


/*---------------------------------------------------------------------------*\
                      Class pqFoamBlockMeshControls Declaration
\*---------------------------------------------------------------------------*/

class pqFoamBlockMeshControls
:
    public pqPropertyWidget
{
    Q_OBJECT;
    typedef pqPropertyWidget Superclass;

    // Private Data

        //- Refresh (push button)
        vtkSMProperty* refresh_;

        //- Show Patch Names (bool property)
        vtkSMIntVectorProperty* showPatchNames_;

        //- Show Point Numbers (bool property)
        vtkSMIntVectorProperty* showPointNumbers_;


    // Private Member Functions

    //- Update property
    void fireCommand(vtkSMProperty* prop);

    //- Update int property or toggle bool property
    void fireCommand(vtkSMIntVectorProperty* prop, int val);

    //- Update "BlockArrayStatus", "CurvedEdgesArrayStatus" information
    void updateParts();


    //- No copy construct
    pqFoamBlockMeshControls(const pqFoamBlockMeshControls&) = delete;

    //- No copy assignment
    void operator=(const pqFoamBlockMeshControls&) = delete;


protected slots:

    // Protected Member Functions

    //- Trigger refresh
    void refreshPressed();

    //- Sync property with changed checkbox state, update rendered view(s)
    void showPatchNames(bool checked);

    //- Sync property with changed checkbox state, update rendered view(s)
    void showPointNumbers(bool checked);


public:

    //- Construct from components
    pqFoamBlockMeshControls
    (
        vtkSMProxy* proxy,
        vtkSMPropertyGroup* group,
        QWidget* parent = nullptr
    );


    //- Destructor
    virtual ~pqFoamBlockMeshControls();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
