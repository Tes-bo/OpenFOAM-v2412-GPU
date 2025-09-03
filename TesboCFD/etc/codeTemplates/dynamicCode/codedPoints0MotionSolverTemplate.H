/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2021 OpenCFD Ltd.
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

Description
    Template for use with dynamic code generation of a
    points0MotionSolver.

SourceFiles
    codedPoints0MotionSolverTemplate.C

\*---------------------------------------------------------------------------*/

#ifndef codedPoints0MotionSolverTemplate_H
#define codedPoints0MotionSolverTemplate_H

#include "points0MotionSolver.H"
#include "dictionaryContent.H"

//{{{ begin codeInclude
${codeInclude}
//}}} end codeInclude

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     A templated codedPoints0MotionSolver
\*---------------------------------------------------------------------------*/

class ${typeName}Points0MotionSolver
:
    public points0MotionSolver,
    public dictionaryContent
{
    // Private Member Functions

        //- Report a message with the SHA1sum
        inline static void printMessage(const char* message)
        {
            Info<< message << " sha1: " << SHA1sum << '\n';
        }

        //- No copy construct
        ${typeName}Points0MotionSolver
        (
            const ${typeName}Points0MotionSolver&
        ) = delete;

        //- No copy assignment
        void operator=
        (
            const ${typeName}Points0MotionSolver&
        ) = delete;


public:

    //- SHA1 representation of the code content
    static constexpr const char* const SHA1sum = "${SHA1sum}";

    //- Runtime type information
    TypeName("${typeName}");


    // Constructors

        //- Construct from mesh and dictionary
        ${typeName}Points0MotionSolver
        (
            const polyMesh& mesh,
            const IOdictionary&
        );


    //- Destructor
    virtual ~${typeName}Points0MotionSolver();


    // Member Functions

        //- Code context as a dictionary
        const dictionary& codeContext() const noexcept
        {
            return dictionaryContent::dict();
        }

        //- Provide current points for motion. Uses current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        virtual void solve()
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
