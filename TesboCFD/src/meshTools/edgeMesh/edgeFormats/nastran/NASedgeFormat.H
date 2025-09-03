/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017-2024 OpenCFD Ltd.
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
    Foam::fileFormats::NASedgeFormat

Description
    Nastran edge reader.

    - Interprets "CBEAM", "CROD" and "PLOTEL" entries as edges.
    - Handles Nastran short, long formats and comma-separated free format.
    - Properly handles the Nastran compact floating point notation: \n
    \verbatim
        GRID          28        10.20269-.030265-2.358-8
    \endverbatim

SourceFiles
    NASedgeFormat.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_NASedgeFormat_H
#define Foam_NASedgeFormat_H

#include "edgeMesh.H"
#include "NASCore.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                      Class NASedgeFormat Declaration
\*---------------------------------------------------------------------------*/

class NASedgeFormat
:
    public edgeMesh,
    public NASCore
{
public:

    // Constructors

        //- Construct from file name
        explicit NASedgeFormat(const fileName& filename);


    // Selectors

        //- Read file and return edge mesh
        static autoPtr<edgeMesh> New(const fileName& name)
        {
            return autoPtr<edgeMesh>::NewFrom<NASedgeFormat>(name);
        }


    //- Destructor
    virtual ~NASedgeFormat() = default;


    // Member Functions

        //- Read from a file
        virtual bool read(const fileName& filename) override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
