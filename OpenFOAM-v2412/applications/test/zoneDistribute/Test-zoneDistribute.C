/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
    Copyright (C) 2019-2020 DLR
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

Application
    Test-zoneDistribute

Description
    Test of zoneDistribute validated with mapDistribute

    Original code supplied by Henning Scheufler, DLR (2019)

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "centredCPCCellToCellStencilObject.H"
#include "zoneDistribute.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    const extendedCentredCellToCellStencil& stencil =
        centredCPCCellToCellStencilObject::New(mesh);

    const List<List<label>>& stencilAddr = stencil.stencil();

    List<List<vector>> stencilCentre(mesh.nCells());

    stencil.collectData
    (
       mesh.C(),
       stencilCentre
    );

    zoneDistribute exchangeFields(mesh);


    boolList interfaceCell(mesh.nCells(),true);
    exchangeFields.setUpCommforZone(interfaceCell);


    const labelListList& stencil_zoneDist = exchangeFields.getStencil();
    const globalIndex& gblNumbering = exchangeFields.globalNumbering();

    Map<vectorField> mapCC(exchangeFields.getFields(interfaceCell,mesh.C()));

    // compare stencils
    Pout<< "size of the stencil match "
        << (stencilAddr.size() == stencil_zoneDist.size()) << endl;

    label stencilMatch = 0;

    forAll(stencilAddr,celli)
    {
        const vectorField& neiCC = mapCC[celli];
        if (neiCC.size() != stencilCentre[celli].size())
        {
            continue;
        }

        bool foundAllLabel = true;
        for (const vector& cc : neiCC)
        {
            if (!stencilCentre[celli].found(cc))
            {
                foundAllLabel = false;
                break;
            }
        }

        if (foundAllLabel)
        {
            ++stencilMatch;
        }

    }

    if (stencilMatch == mesh.nCells())
    {
        Pout << "all Values are identical "  << endl;
    }
    else
    {
        Pout << "values did not match in : " << stencilMatch << " of "
             << mesh.nCells() << " cases" << endl;
    }

    return 0;
}


// ************************************************************************* //
