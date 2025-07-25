/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenCFD Ltd.
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
    test

Description
    Finite volume method test code for 2-D space.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "vector2D.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef GeometricField<vector2D, fvPatchField, volMesh> volVector2DField;
    typedef fvPatchField<vector2D> fvPatchVector2DField;

    defineTemplateTypeNameAndDebug(volVector2DField::Internal, 0);
    defineTemplateTypeNameAndDebug(volVector2DField, 0);

    defineTemplateRunTimeSelectionTable(fvPatchVector2DField, patch);
    defineTemplateRunTimeSelectionTable(fvPatchVector2DField, patchMapper);
    defineTemplateRunTimeSelectionTable(fvPatchVector2DField, dictionary);
}


int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"

    GeometricField<vector2D, fvPatchField, volMesh> fld
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
