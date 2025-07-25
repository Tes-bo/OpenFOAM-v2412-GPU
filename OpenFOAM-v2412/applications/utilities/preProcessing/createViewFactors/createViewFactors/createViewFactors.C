/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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
    createViewFactors

Group
    grpPreProcessingUtilities

Description
    Creates view factors to be used in the view-factor radiation model.

    Operands:
    \table
      Operand | Type                | Location
      input   | dictionary          | \<constant\>/viewFactorsDict
      input   | dictionary          | \<constant\>/finalAgglom
      output  | scalarListList      | \<constant\>/F
      output  | mapDistribute       | \<constant\>/mapDist
      output  | labelListList       | \<constant\>/globalFaceFaces
      output  | volScalarField      | \<time\>/viewVectorField
      output  | OBJ                 | allVisibleFaces.obj
    \endtable

    where the dictionaries mean:
    \table
      Dictionary      | Description
      viewFactorsDict | Main-control dictionary
      finalAgglom     | (Optional) Agglomeration addressing (from faceAgglomerate)
      F               | View factors (matrix)
      mapDist         | Map used for parallel running
      globalFaceFaces | Face addressing
      viewVectorField | View factors as a volume field
      allVisibleFaces.obj | The visualisation of the rays
    \endtable

Usage
    Minimal example in \c <constant>/viewFactorsDict:
    \verbatim
    // Inherited entries
    raySearchEngine     <word>;
    agglomerate         <bool>;
    nRayPerFace         <label>;
    writeViewFactors    <bool>;
    writeRays           <bool>;
    ...
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Reqd | Deflt
      raySearchEngine  | Ray search engine type         | word | yes  | -
      agglomerate  | Flag to agglomeration              | bool | yes  | -
      nRayPerFace  | Number of rays issued per face     | label | yes | -
      writeViewFactors | Flag to write the view factor field | bool | yes |-
      writeRays    | Flag to write the ray geometry     | bool | no | false
    \endtable

    Options for the \c raySearchEngine entry:
    \verbatim
      voxel    | Ray search engine discretising space into uniform voxels
    \endverbatim

    The inherited entries are elaborated in:
    - \link viewFactorModel.H \endlink
    - \link raySearchEngine.H \endlink

Note

  - Participating patches must be in the \c vewFactorWall group, i.e. using the
    \c inGroups entry of the "\<case\>/polyMesh/boundary" file.

    \verbatim
    myPatch
    {
        type            wall;
        inGroups        2(wall viewFactorWall);
        ...
    }
    \endverbatim

    Reads:

    - <constant>/viewFactorsDict : main controls
    - <constant>/finalAgglom : agglomeration addressing (from faceAgglomerate)


    Generates:

    - <constant>/F : view factors (matrix)
    - <constant>/mapDist : map used for parallel running
    - <constant>/globalFaceFaces : face addressing

SeeAlso
- Foam::VF::raySearchEngine
- Foam::VE::viewFactorModel

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "viewFactorModel.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    IOdictionary dict
    (
        IOobject
        (
            "viewFactorsDict",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ
        )
    );


    // Calculate the view factors
    auto modelPtr = VF::viewFactorModel::New(mesh, dict);

    modelPtr->calculate();

    Info<< nl;

    runTime.printExecutionTime(Info);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //