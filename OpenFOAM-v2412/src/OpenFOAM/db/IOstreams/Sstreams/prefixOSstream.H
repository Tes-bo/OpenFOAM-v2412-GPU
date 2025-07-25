/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2014 OpenFOAM Foundation
    Copyright (C) 2020-2023 OpenCFD Ltd.
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
    Foam::prefixOSstream

Description
    Version of OSstream that prints a prefix on each line.

    This is useful for running in parallel as it allows the processor number
    to be automatically prepended to each message line.

SourceFiles
    prefixOSstream.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_prefixOSstream_H
#define Foam_prefixOSstream_H

#include "OSstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class prefixOSstream Declaration
\*---------------------------------------------------------------------------*/

class prefixOSstream
:
    public OSstream
{
    // Private Data

        bool printPrefix_;

        string prefix_;


    // Private Member Functions

        //- Emit pending prefix before any output
        inline void checkWritePrefix();


public:

    // Constructors

        //- Construct and set stream status
        prefixOSstream
        (
            std::ostream& os,
            const string& streamName,
            IOstreamOption streamOpt = IOstreamOption()
        );


    // Member Functions

    // Enquiry

        //- Return the stream prefix
        const string& prefix() const noexcept { return prefix_; }

        //- Return non-const access to the stream prefix
        string& prefix() noexcept { return prefix_; }


    // Write Functions

        //- Inherit write methods from OSstream
        using OSstream::writeQuoted;

        //- Write token to stream or otherwise handle it.
        //  \return false if the token type was not handled by this method
        virtual bool write(const token& tok) override;

        //- Write character
        virtual Ostream& write(const char c) override;

        //- Write character/string content, with/without surrounding quotes
        virtual Ostream& writeQuoted
        (
            const char* str,
            std::streamsize len,
            const bool quoted=true
        ) override;

        //- Write character string
        virtual Ostream& write(const char* str) override;

        //- Write word
        virtual Ostream& write(const word& val) override;

        //- Write string
        virtual Ostream& write(const std::string& val) override;

        //- Write int32_t
        virtual Ostream& write(const int32_t val) override;

        //- Write int64_t
        virtual Ostream& write(const int64_t val) override;

        //- Write float
        virtual Ostream& write(const float val) override;

        //- Write double
        virtual Ostream& write(const double val) override;

        //- Write binary block
        virtual Ostream& write(const char* buf, std::streamsize count) override;

        //- Add indentation characters
        virtual void indent() override;


    // Print

        //- Print stream description to Ostream
        virtual void print(Ostream& os) const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
