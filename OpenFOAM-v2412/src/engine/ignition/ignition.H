/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::ignition

Description
    Foam::ignition

SourceFiles
    ignition.C
    ignitionIO.C

\*---------------------------------------------------------------------------*/

#ifndef ignition_H
#define ignition_H

#include "ignitionSite.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                           Class ignition Declaration
\*---------------------------------------------------------------------------*/

class ignition
{
    // Private data

        const fvMesh& mesh_;

        bool ignite_;

        PtrList<ignitionSite> ignSites_;


    // Private Member Functions

        //- No copy construct
        ignition(const ignition&) = delete;

        //- No copy assignment
        void operator=(const ignition&) = delete;


public:

    // Constructors

        //- Construct from Istream and database
        ignition(const dictionary&, const Time&, const fvMesh&);

        //- Construct from Istream and engineTime
        ignition(const dictionary&, const engineTime&, const fvMesh&);


    // Member Functions

        //- Return the set of ignition sites
        const PtrList<ignitionSite>& sites()
        {
            return ignSites_;
        }


        // Check

            bool ignite() const
            {
                return ignite_;
            }

            //- Are any of the ignition site currently igniting
            bool igniting() const;

            //- Has the mixture been ignited?
            bool ignited() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
