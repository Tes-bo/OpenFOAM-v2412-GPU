/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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
    Foam::cellTable

Description
    The cellTable persistent data saved as a Map<dictionary>.

    The meshReader supports cellTable information.

    The <tt>constant/cellTable</tt> file is an \c IOMap<dictionary> that is
    used to save the information persistently. It contains the cellTable
    information of the following form:

    \verbatim
        (
            ID
            {
                Label           WORD;
                MaterialType    WORD;
                MaterialId      INT;
                PorosityId      INT;
                ColorIdx        INT;
                ...
            }
        ...
        )
    \endverbatim

    If the \a Label is missing, a value <tt>cellTable_{ID}</tt> will be
    inferred. If the \a MaterialType is missing, the value \a fluid will
    be inferred.

SourceFiles
    cellTable.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_cellTable_H
#define Foam_cellTable_H

#include "Map.H"
#include "dictionary.H"
#include "wordRes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class polyMesh;
class objectRegistry;

/*---------------------------------------------------------------------------*\
                         Class cellTable Declaration
\*---------------------------------------------------------------------------*/

class cellTable
:
    public Map<dictionary>
{
    // Private Member Functions

        //- Add required entries - MaterialType
        void addDefaults();

        void setEntry(const label id, const word& key, const word& value);


public:

    //- No copy construct
    cellTable(const cellTable&) = delete;


    // Constructors

        //- Default construct
        cellTable() noexcept = default;

        //- Read construct from registry, name, instance
        explicit cellTable
        (
            const objectRegistry& obr,
            const word& name = "cellTable",
            const fileName& instance = "constant"
        );


    //- Destructor
    ~cellTable() = default;


    // Member Functions

        //- Add to the end, return index
        label push_back(const dictionary& dict);

        //- The max table index, -1 if empty
        label maxIndex() const;

        //- The index corresponding to entry with 'Label' of given name,
        //- or -1 if not found
        label findIndex(const word& name) const;

        //- The 'Label' name corresponding to id,
        //- or cellTable_ID if not otherwise defined
        word name(const label id) const;

        //- Return the extracted Map of (id => name)
        Map<word> names() const;

        //- Return the extracted Map of (id => names) selected by patterns
        Map<word> names(const wordRes& patterns) const;

        //- Return the extracted Map of (id => name) for materialType
        //- (fluid | solid | shell)
        Map<word> selectType(const word& materialType) const;

        //- Return a Map of (id => name) for fluids
        Map<word> fluids() const;

        //- Return a Map of (id => name) for shells
        Map<word> shells() const;

        //- Return a Map of (id => name) for solids
        Map<word> solids() const;

        //- Return a Map of (id => fluid|solid|shell)
        Map<word> materialTypes() const;

        //- Assign material Type
        void setMaterial(const label, const word&);

        //- Assign name
        void setName(const label, const word&);

        //- Assign default name if not already set
        void setName(const label);

        //- Read constant/cellTable
        void readDict
        (
            const objectRegistry&,
            const word& name = "cellTable",
            const fileName& instance = "constant"
        );

        //- Write constant/cellTable for later reuse
        void writeDict
        (
            const objectRegistry&,
            const word& name = "cellTable",
            const fileName& instance = "constant"
        ) const;


    // Member Operators

        //- Assignment
        void operator=(const cellTable&);

        //- Assign from Map<dictionary>
        void operator=(const Map<dictionary>&);

        //- Assign from cellZones
        void operator=(const polyMesh&);


    // Friend Functions

        //- Classify tableIds into cellZones according to the cellTable
        void addCellZones(polyMesh&, const labelList& tableIds) const;

        //- Combine tableIds together
        //  each dictionary entry is a wordList
        void combine(const dictionary& mapDict, labelList& tableIds);


    // Housekeeping

        //- Add to the end, return index
        label append(const dictionary& dict) { return push_back(dict); }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
