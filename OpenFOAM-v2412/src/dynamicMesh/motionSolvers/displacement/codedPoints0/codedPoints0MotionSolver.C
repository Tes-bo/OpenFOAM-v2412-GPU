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

#include "codedPoints0MotionSolver.H"
#include "dictionary.H"
#include "Time.H"
#include "dynamicCode.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(codedPoints0MotionSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        codedPoints0MotionSolver,
        dictionary
    );
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::dlLibraryTable& Foam::codedPoints0MotionSolver::libs() const
{
    return mesh().time().libs();
}


Foam::string Foam::codedPoints0MotionSolver::description() const
{
    return "points0MotionSolver " + name();
}


void Foam::codedPoints0MotionSolver::clearRedirect() const
{
    redirectMotionSolverPtr_.reset(nullptr);
}


const Foam::dictionary&
Foam::codedPoints0MotionSolver::codeDict() const
{
    return motionSolver::coeffDict();
}


void Foam::codedPoints0MotionSolver::prepare
(
    dynamicCode& dynCode,
    const dynamicCodeContext& context
) const
{
    // Set additional rewrite rules
    dynCode.setFilterVariable("typeName", name_);

    // Compile filtered C template
    dynCode.addCompileFile(codeTemplateC);

    // Copy filtered H template
    dynCode.addCopyFile(codeTemplateH);

    #ifdef FULLDEBUG
    dynCode.setFilterVariable("verbose", "true");
    DetailInfo
        <<"compile " << name_ << " sha1: " << context.sha1() << endl;
    #endif

    // Define Make/options
    dynCode.setMakeOptions
    (
        "EXE_INC = -g \\\n"
        "-I$(LIB_SRC)/finiteVolume/lnInclude \\\n"
        "-I$(LIB_SRC)/meshTools/lnInclude \\\n"
        "-I$(LIB_SRC)/dynamicMesh/lnInclude \\\n"
        "-I$(LIB_SRC)/fvMotionSolvers/lnInclude \\\n"
      + context.options()
      + "\n\nLIB_LIBS = \\\n"
        "    -lfiniteVolume \\\n"
        "    -lmeshTools \\\n"
        "    -ldynamicMesh \\\n"
        "    -lfvMotionSolvers \\\n"
      + context.libs()
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::codedPoints0MotionSolver::codedPoints0MotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict
)
:
    motionSolver(mesh, dict, typeName),
    codedBase(),
    name_(dict.getCompat<word>("name", {{"redirectType", 1706}})),
    redirectMotionSolverPtr_(nullptr)
{
    updateLibrary(name_);
    redirectMotionSolver();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::motionSolver&
Foam::codedPoints0MotionSolver::redirectMotionSolver() const
{
    if (!redirectMotionSolverPtr_)
    {
        // Get the dictionary for the solver and override the
        // solver name (in case it is not a subdictionary and contains
        // the 'coded' as the motionSolver)
        dictionary constructDict(motionSolver::coeffDict());
        constructDict.set("solver", name_);
        constructDict.set("motionSolver", name_);

        IOobject io(*this);
        io.readOpt(IOobject::NO_READ);

        redirectMotionSolverPtr_ = motionSolver::New
        (
            mesh(),
            IOdictionary(io, constructDict)
        );
    }

    return *redirectMotionSolverPtr_;
}


Foam::tmp<Foam::pointField> Foam::codedPoints0MotionSolver::curPoints() const
{
    updateLibrary(name_);
    return redirectMotionSolver().curPoints();
}


void Foam::codedPoints0MotionSolver::solve()
{
    updateLibrary(name_);
    redirectMotionSolver().solve();
}


void Foam::codedPoints0MotionSolver::movePoints(const pointField& fld)
{
    updateLibrary(name_);
    return redirectMotionSolver().movePoints(fld);
}


void Foam::codedPoints0MotionSolver::updateMesh(const mapPolyMesh& mpm)
{
    updateLibrary(name_);
    return redirectMotionSolver().updateMesh(mpm);
}


// ************************************************************************* //
