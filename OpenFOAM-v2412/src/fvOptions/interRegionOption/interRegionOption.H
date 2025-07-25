/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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
    Foam::fv::interRegionOption

Description
    Intermediate class for handling inter-region exchanges.

Usage
    Minimal example by using \c constant/fvOptions:
    \verbatim
    <userDefinedName1>
    {
        // Mandatory/Optional (inherited) entries
        ...

        // Mandatory entries (unmodifiable)
        nbrRegion         <nbrRegionName>;

        // Optional entries (unmodifiable)
        master            true;

        // Mandatory/Optional (derived) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                         | Type  | Reqd | Dflt
      nbrRegion  | Name of the neighbour region to map | word  | yes  | -
      master     | Flag to determine whether master or <!--
                 --> secondary region                  | bool  | no   | true
    \endtable

    The inherited entries are elaborated in:
     - \link fvOption.H \endlink
     - \link meshToMesh.H \endlink

SourceFiles
    interRegionOption.C

\*---------------------------------------------------------------------------*/

#ifndef interRegionOption_H
#define interRegionOption_H

#include "fvOption.H"
#include "volFields.H"
#include "autoPtr.H"
#include "meshToMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                Class interRegionOption Declaration
\*---------------------------------------------------------------------------*/

class interRegionOption
:
    public fv::option
{
protected:

    // Protected Data

        //- Master or slave region
        bool master_;

        //- Name of the neighbour region to map
        word nbrRegionName_;

        //- Mesh to mesh interpolation object
        autoPtr<meshToMesh> meshInterpPtr_;


    // Protected Member Functions

        //- Set the mesh to mesh interpolation object
        void setMapper();


public:

    //- Runtime type information
    TypeName("interRegionOption");


    // Constructors

        //- Construct from dictionary
        interRegionOption
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        interRegionOption(const interRegionOption&) = delete;

        //- No copy assignment
        void operator=(const interRegionOption&) = delete;


    //- Destructor
    virtual ~interRegionOption() = default;


    // Member Functions

        // Access

            //- Return const access to the neighbour region name
            inline const word& nbrRegionName() const;

            //- Return const access to the mapToMap pointer
            inline const meshToMesh& meshInterp() const;


        // IO

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "interRegionOptionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
