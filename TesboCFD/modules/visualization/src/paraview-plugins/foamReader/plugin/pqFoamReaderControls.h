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
    pqFoamReaderControls

Description
    A custom property group widget for the PVFoamReader.

SourceFiles
    pqFoamReaderControls.cxx

\*---------------------------------------------------------------------------*/
#ifndef pqFoamReaderControls_h
#define pqFoamReaderControls_h

#include "pqPropertyWidget.h"

// Forward Declarations
class vtkSMProperty;
class vtkSMIntVectorProperty;


/*---------------------------------------------------------------------------*\
                    Class pqFoamReaderControls Declaration
\*---------------------------------------------------------------------------*/

class pqFoamReaderControls
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

        //- Show Groups Only (bool property)
        vtkSMIntVectorProperty* showGroupsOnly_;

        //- IncludeSets (bool property)
        vtkSMIntVectorProperty* includeSets_;

        //- IncludeZones (bool property)
        vtkSMIntVectorProperty* includeZones_;

        //- MeshCaching (enum property)
        vtkSMIntVectorProperty* meshCaching_;


    // Private Member Functions

    //- Update property
    void fireCommand(vtkSMProperty* prop);

    //- Update int property or toggle bool property
    void fireCommand(vtkSMIntVectorProperty* prop, int val);


    //- No copy construct
    pqFoamReaderControls(const pqFoamReaderControls&) = delete;

    //- No copy assignment
    void operator=(const pqFoamReaderControls&) = delete;


private slots:

    // Private Member Functions

    //- Update "PartArrayStatus" property information
    void updateParts();


protected slots:

    // Protected Member Functions

    void refreshPressed();
    void cacheMesh(int idx);
    void showPatchNames(bool checked);
    void showGroupsOnly(bool checked);
    void includeSets(bool checked);
    void includeZones(bool checked);


public:

    //- Construct from components
    pqFoamReaderControls
    (
        vtkSMProxy* proxy,
        vtkSMPropertyGroup* group,
        QWidget* parent = nullptr
    );

    //- Destructor
    virtual ~pqFoamReaderControls();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
