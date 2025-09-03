/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2024 OpenCFD Ltd.
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
    Foam::boundaryRegion

Description
    The boundaryRegion persistent data saved as a Map<dictionary>.

    The meshReader supports boundaryRegion information.

    The <tt>constant/boundaryRegion</tt> file is an \c IOMap<dictionary>
    that is used to save the information persistently.
    It contains the boundaryRegion information of the following form:

    \verbatim
        (
            INT
            {
                BoundaryType    WORD;
                Label           WORD;
            }
            ...
        )
    \endverbatim

SourceFiles
    boundaryRegion.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_boundaryRegion_H
#define Foam_boundaryRegion_H

#include "Map.H"
#include "dictionary.H"
#include "labelList.H"
#include "wordRes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class objectRegistry;

/*---------------------------------------------------------------------------*\
                      Class boundaryRegion Declaration
\*---------------------------------------------------------------------------*/

class boundaryRegion
:
    public Map<dictionary>
{
public:

    //- No copy construct
    boundaryRegion(const boundaryRegion&) = delete;


    // Constructors

        //- Default construct
        boundaryRegion() noexcept = default;

        //- Read construct from registry, name, instance
        explicit boundaryRegion
        (
            const objectRegistry& obr,
            const word& name = "boundaryRegion",
            const fileName& instance = "constant"
        );


    //- Destructor
    ~boundaryRegion() = default;


    // Member Functions

        //- Add to the end, return index
        label push_back(const dictionary& dict);

        //- The max table index, -1 if empty
        label maxIndex() const;

        //- The index corresponding to entry with 'Label' of given name,
        //- or -1 if not found
        label findIndex(const word& name) const;

        //- The 'Label' name corresponding to id,
        //- or boundaryRegion_ID if not otherwise defined
        word name(const label id) const;

        //- Return the extracted Map of (id => name)
        Map<word> names() const;

        //- Return the extracted Map of (id => names) selected by patterns
        Map<word> names(const wordRes& patterns) const;

        //- Return the extracted Map of (id => type)
        Map<word> boundaryTypes() const;

        //- Return BoundaryType corresponding to patch 'name',
        //- "patch" if not found
        word boundaryType(const word& name) const;


        //- Read constant/boundaryRegion
        void readDict
        (
            const objectRegistry& obr,
            const word& name = "boundaryRegion",
            const fileName& instance = "constant"
        );

        //- Write constant/boundaryRegion for later reuse
        void writeDict
        (
            const objectRegistry& obr,
            const word& name = "boundaryRegion",
            const fileName& instance = "constant"
        ) const;


    // Member Operators

        //- Copy assignment
        void operator=(const boundaryRegion&);

        //- Assign from Map<dictionary>
        void operator=(const Map<dictionary>&);


    // Friend Functions

        //- Rename regions
        //  each dictionary entry is a single word:
        //  \verbatim
        //      newPatchName    originalName;
        //  \endverbatim
        void rename(const dictionary&);


    // Housekeeping

        //- Add to the end, return index
        label append(const dictionary& dict) { return push_back(dict); }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
