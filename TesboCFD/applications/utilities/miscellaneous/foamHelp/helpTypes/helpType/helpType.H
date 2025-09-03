/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2014 OpenFOAM Foundation
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
    Foam::helpType

Description
    Base class for foam help classes

SourceFiles
    helpType.C
    helpTypeNew.C

\*---------------------------------------------------------------------------*/

#ifndef helpType_H
#define helpType_H

#include "runTimeSelectionTables.H"
#include "autoPtr.H"
#include "argList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class fvMesh;

/*---------------------------------------------------------------------------*\
                          Class helpType Declaration
\*---------------------------------------------------------------------------*/

class helpType
{
protected:

    //- Return file path to the Doxygen sources (if available)
    fileName doxygenPath() const;

    //- Display the list of documentation options
    void displayDocOptions
    (
        const string& searchStr,
        const bool exactMatch,
        const word& ext
    ) const;

    //- Display the help documentation in a browser
    void displayDoc
    (
        const word& className,
        const string& searchStr,
        const bool exactMatch,
        const word& ext
    ) const;


public:

    //- Runtime type information
    TypeName("helpType");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        helpType,
        dictionary,
        (),
        ()
    );


    //- Default construct
    helpType() = default;

    //- Selector
    static autoPtr<helpType> New(const word& helpTypeName);

    //- Destructor
    virtual ~helpType() = default;


    // Member Functions

        //- Initialise - typically setting static variables,
        //  e.g. command line arguments
        virtual void init();

        //- Execute the help
        virtual void execute(const argList& args, const fvMesh& mesh) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
