/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 Bernhard Gschaider
    Copyright (C) 2019-2021 OpenCFD Ltd.
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
    Foam::expressions::fvExprDriverWriter

Description
    Registered input/output for an expressions::fvExprDriver

SourceFiles
    fvExprDriverWriter.C

\*---------------------------------------------------------------------------*/

#ifndef expressions_fvExprDriverWriter_H
#define expressions_fvExprDriverWriter_H

#include "fvExprDriver.H"
#include "regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace expressions
{

/*---------------------------------------------------------------------------*\
                     Class fvExprDriverWriter Declaration
\*---------------------------------------------------------------------------*/

class fvExprDriverWriter
:
    public regIOobject
{
    // Private Data

        //- The driver to read/write
        fvExprDriver& driver_;


    // Private Member Functions

        //- No default construct
        fvExprDriverWriter() = delete;

        //- No copy construct
        fvExprDriverWriter(const fvExprDriverWriter&) = delete;

        //- No copy assignment
        void operator=(const fvExprDriverWriter&) = delete;


public:

    //- Runtime type information
    TypeNameNoDebug("fvExprDriverWriter");


    // Constructors

        //- Construct for named driver
        fvExprDriverWriter(const word& name, fvExprDriver& driver);


    //- Destructor
    virtual ~fvExprDriverWriter() = default;


    // Member Functions

        virtual bool readData(Istream& is);
        virtual bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace expressions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
