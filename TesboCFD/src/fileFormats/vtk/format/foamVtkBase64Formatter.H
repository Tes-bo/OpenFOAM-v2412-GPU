/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2022 OpenCFD Ltd.
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
    Foam::vtk::base64Formatter

Description
    Inline base-64 encoded binary output.
    Uses an output filter layer to write base-64 encoded content.

\*---------------------------------------------------------------------------*/

#ifndef Foam_vtk_base64Formatter_H
#define Foam_vtk_base64Formatter_H

#include "foamVtkFormatter.H"
#include "foamVtkBase64Layer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace vtk
{

/*---------------------------------------------------------------------------*\
                    Class vtk::base64Formatter Declaration
\*---------------------------------------------------------------------------*/

class base64Formatter
:
    public vtk::foamVtkBase64Layer
{
    // Private Data Members

        static const char* name_;
        static const outputOptions opts_;


    // Private Member Functions

        //- No copy construct
        base64Formatter(const base64Formatter&) = delete;

        //- No copy assignment
        void operator=(const base64Formatter&) = delete;


public:

    // Constructors

        //- Construct and attach to an output stream
        explicit base64Formatter(std::ostream& os);


    //- Destructor. Closes/flushes the underlying layer.
    virtual ~base64Formatter();


    // Member Functions

        //- The output is INLINE_BASE64.
        virtual const vtk::outputOptions& opts() const;

        //- Name for the XML output type ("binary")
        virtual const char* name() const;


        //- End the encoding sequence (padding the final characters with '=')
        //  and write a newline to the output if any encoding was done.
        virtual void flush();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace vtk
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
