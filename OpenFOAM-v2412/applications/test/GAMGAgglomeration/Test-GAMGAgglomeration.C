/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Test-GAMGAgglomeration

Description
    Test application for GAMG agglomeration. Hardcoded to expect GAMG on p.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "GAMGAgglomeration.H"
#include "OFstream.H"
#include "meshTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::addBoolOption
    (
        "writeObj",
        "write obj files of agglomeration"
    );
    argList::addBoolOption
    (
        "normalise",
        "normalise agglomeration (0..1)"
    );

    #include "setRootCase.H"
    #include "createTime.H"

    bool writeObj = args.found("writeObj");
    bool normalise = args.found("normalise");

    #include "createMesh.H"

    const fvSolution& sol = static_cast<const fvSolution&>(mesh);
    const dictionary& pDict = sol.subDict("solvers").subDict("p");

    const GAMGAgglomeration& agglom = GAMGAgglomeration::New
    (
        mesh,
        pDict
    );

    labelList cellToCoarse(identity(mesh.nCells()));
    CompactListList<label> coarseToCell
    (
        invertOneToManyCompact(mesh.nCells(), cellToCoarse)
    );

    ++runTime;

    // Write initial agglomeration
    {
        volScalarField scalarAgglomeration
        (
            mesh.thisDb().newIOobject("agglomeration"),
            mesh,
            Foam::zero{},
            dimless,
            fvPatchFieldBase::zeroGradientType()
        );
        scalarField& fld = scalarAgglomeration.primitiveFieldRef();
        forAll(fld, celli)
        {
            fld[celli] = cellToCoarse[celli];
        }
        if (normalise)
        {
            fld /= max(fld);
        }
        scalarAgglomeration.correctBoundaryConditions();
        scalarAgglomeration.write();

        Info<< "Writing initial cell distribution to "
            << runTime.timeName() << endl;
    }


    for (label level = 0; level < agglom.size(); level++)
    {
        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        const labelList& addr = agglom.restrictAddressing(level);
        label coarseSize = max(addr)+1;

        Info<< "Level : " << level << endl
            << "    current size      : "
            << returnReduce(addr.size(), sumOp<label>()) << endl
            << "    agglomerated size : "
            << returnReduce(coarseSize, sumOp<label>()) << endl;

        labelList newAddr;
        label newCoarseSize = 0;
        bool ok = GAMGAgglomeration::checkRestriction
        (
            newAddr,
            newCoarseSize,

            agglom.meshLevel(level).lduAddr(),
            addr,
            coarseSize
        );
        if (!ok)
        {
            WarningInFunction
                << "At level " << level
                << " there are " << coarseSize
                << " agglomerated cells but " << newCoarseSize
                << " disconnected regions" << endl
                << "    This means that some agglomerations (coarse cells)"
                << "    consist of multiple disconnected regions."
                << endl;
        }


        forAll(addr, finei)
        {
            labelUIndList(cellToCoarse, coarseToCell[finei]) = addr[finei];
        }
        coarseToCell = invertOneToManyCompact(coarseSize, cellToCoarse);

        // Write agglomeration
        {
            volScalarField scalarAgglomeration
            (
                mesh.thisDb().newIOobject("agglomeration"),
                mesh,
                Foam::zero{},
                dimless,
                fvPatchFieldBase::zeroGradientType()
            );

            scalarField& fld = scalarAgglomeration.primitiveFieldRef();
            forAll(fld, celli)
            {
                fld[celli] = cellToCoarse[celli];
            }
            if (normalise)
            {
                fld /= max(fld);
            }
            scalarAgglomeration.correctBoundaryConditions();
            scalarAgglomeration.write();
        }

        if (writeObj)
        {
            OFstream str(runTime.path()/runTime.timeName()/"aggomeration.obj");
            label vertI = 0;

            // Write all mesh cc
            forAll(mesh.cellCentres(), celli)
            {
                meshTools::writeOBJ(str, mesh.cellCentres()[celli]);
                vertI++;
            }

            // Determine coarse cc
            forAll(coarseToCell, coarsei)
            {
                const auto& cellLabels = coarseToCell[coarsei];

                point coarseCc = average
                (
                    pointField(mesh.cellCentres(), cellLabels)
                );
                meshTools::writeOBJ(str, coarseCc);
                vertI++;

                for (label celli : cellLabels)
                {
                    str << "l " << celli+1 << ' ' << vertI << nl;
                }
            }
        }

        Info<< endl;
    }


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
