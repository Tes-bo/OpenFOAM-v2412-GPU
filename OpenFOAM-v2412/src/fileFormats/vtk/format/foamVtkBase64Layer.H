/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenCFD Ltd.
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
    Foam::vtk::foamVtkBase64Layer

Description
    Base-64 encoded output layer - normally only used indirectly by formatters.

\*---------------------------------------------------------------------------*/

#ifndef Foam_vtk_base64Layer_H
#define Foam_vtk_base64Layer_H

#include "foamVtkFormatter.H"
#include "base64Layer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace vtk
{

/*---------------------------------------------------------------------------*\
                   Class vtk::foamVtkBase64Layer Declaration
\*---------------------------------------------------------------------------*/

class foamVtkBase64Layer
:
    public vtk::formatter,
    protected Foam::base64Layer
{
    // Private Data Members

        static const char* encoding_;


protected:

    // Protected Member Functions

        //- Write
        void write(const char* s, std::streamsize n);

        //- No copy construct
        foamVtkBase64Layer(const foamVtkBase64Layer&) = delete;

        //- No copy assignment
        void operator=(const foamVtkBase64Layer&) = delete;


    // Constructors

        //- Construct and attach to an output stream
        explicit foamVtkBase64Layer(std::ostream& os);

public:

    //- Destructor. Closes/flushes the underlying layer
    virtual ~foamVtkBase64Layer();


    // Member Functions

        //- Name for the XML append encoding ("base64").
        virtual const char* encoding() const;

        //- Write leading size for binary output
        //  \return True - format uses this information
        virtual bool writeSize(const uint64_t numbytes);

        virtual void write(const uint8_t val);
        virtual void write(const label val);
        virtual void write(const float val);
        virtual void write(const double val);

        //- End the encoding sequence (padding the final characters with '=')
        virtual void flush();

        //- The encoded length for base64 encoded output.
        virtual std::size_t encodedLength(std::size_t n) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace vtk
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
