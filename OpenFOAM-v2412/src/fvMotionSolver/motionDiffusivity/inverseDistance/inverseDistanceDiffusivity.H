/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2012 OpenFOAM Foundation
    Copyright (C) 2017-2018 OpenCFD Ltd.
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
    Foam::inverseDistanceDiffusivity

Description
    Inverse distance to the given patches motion diffusivity.

SourceFiles
    inverseDistanceDiffusivity.C

\*---------------------------------------------------------------------------*/

#ifndef inverseDistanceDiffusivity_H
#define inverseDistanceDiffusivity_H

#include "uniformDiffusivity.H"
#include "wordRes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class inverseDistanceDiffusivity Declaration
\*---------------------------------------------------------------------------*/

class inverseDistanceDiffusivity
:
    public uniformDiffusivity
{
    // Private data

        //- Patches selected to base the distance on
        //  These can contain patch names or regular expressions to search for.
        wordRes patchNames_;


    // Private Member Functions


        //- No copy construct
        inverseDistanceDiffusivity(const inverseDistanceDiffusivity&) = delete;

        //- No copy assignment
        void operator=(const inverseDistanceDiffusivity&) = delete;


public:

    //- Runtime type information
    TypeName("inverseDistance");


    // Constructors

        //- Construct for the given fvMesh and data Istream
        inverseDistanceDiffusivity(const fvMesh& mesh, Istream& mdData);


    //- Destructor
    virtual ~inverseDistanceDiffusivity() = default;


    // Member Functions

        //- Correct the motion diffusivity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
