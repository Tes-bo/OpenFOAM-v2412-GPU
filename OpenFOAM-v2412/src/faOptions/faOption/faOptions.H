/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2022 OpenCFD Ltd.
------------------------------------------------------------------------------
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
    Foam::fa::options

Description
    Finite-area options

SourceFiles
    faOptions.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_fa_options_H
#define Foam_fa_options_H

#include "faOptionList.H"
#include "IOdictionary.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fa
{

/*---------------------------------------------------------------------------*\
                        Class options Declaration
\*---------------------------------------------------------------------------*/

class options
:
    public IOdictionary,
    public optionList
{
    // Private Member Functions

        //- Create IO object if dictionary is present
        IOobject createIOobject(const fvMesh& mesh) const;

        //- No copy construct
        options(const options&) = delete;

        //- No copy assignment
        void operator=(const options&) = delete;


public:

    // Declare name of the class and its debug switch
    ClassName("faOptions");


    // Constructors

        //- Construct from components with list of field names
        explicit options(const fvMesh& mesh);

        //- Construct faOptions and register to database if not present
        //- otherwise lookup and return
        static options& New(const fvMesh& mesh);


    //- Destructor
    virtual ~options() = default;


    // Member Functions

        //- Inherit read from optionList
        using optionList::read;

        //- Read dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
