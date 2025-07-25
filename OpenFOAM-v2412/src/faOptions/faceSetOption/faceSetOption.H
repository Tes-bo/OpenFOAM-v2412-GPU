/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2022 OpenCFD Ltd.
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
    Foam::fa::faceSetOption

Description
    Intermediate abstract class for handling
    face-set options for the derived faOptions.

Usage
    Minimal example by using \c constant/faOptions:
    \verbatim
    <userDefinedName1>
    {
        // Mandatory/Optional (inherited) entries
        ...

        // Mandatory entries (unmodifiable)
        selectionMode   all;

        // Optional entries (runtime modifiable)
        timeStart       1.0;

        // Conditional mandatory entries (runtime modifiable)

        // when timeStart entry is present
        duration        1.4;

        // when selectionMode=faceZone
        faceZones       (<name> ...);
        //or: faceZone  <name>;

        // when selectionMode=patch
        patches         (<name> ...)
        //or: patch     <name>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property      | Description                        | Type  | Reqd | Dflt
      selectionMode | Mode of face selection - see below | word  | yes  | -
      timeStart     | Start time of faOption             | scalar | no  | -1
      duration      | Duration of faOption execution <!--
                          --> starting from timeStart    | scalar | cndtnl | 0
      faceSet       | Name of operand faceSet            | word   | cndtnl | -
      faceZone      | Name of operand faceZone(s)        | wordRe | cndtnl | -
      faceZones     | Names of operand faceZones         | wordRes| cndtnl | -
      patch         | Name of operand poly patch(s)      | wordRe | cndtnl | -
      patches       | Names of operand poly patches      | wordRes| cndtnl | -
    \endtable

    Options for the \c selectionMode entry:
    \verbatim
      all       | Use all faces in the computational domain
      faceSet   | Use subset corresponding to specified (volume) faceSet
      faceZones | Use subset corresponding to specified (volume) faceZones
      faceZone  | Use subset corresponding to specified (volume) faceZones
      patch     | Use subset corresponding of specified volume patches
      patches   | Use subset corresponding of specified volume patches
    \endverbatim

    The inherited entries are elaborated in:
     - \link faOption.H \endlink

Note
  - Source/sink options are to be added to the right-hand side of equations.

SourceFiles
    faceSetOption.C
    faceSetOptionI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_fa_faceSetOption_H
#define Foam_fa_faceSetOption_H

#include "faOption.H"
#include "faceSet.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fa
{

/*---------------------------------------------------------------------------*\
                         Class faceSetOption Declaration
\*---------------------------------------------------------------------------*/

class faceSetOption
:
    public fa::option
{
public:

    // Public Enumeration

        //- Enumeration for selection mode types
        enum selectionModeType
        {
            smAll,          //!< "all" finite-area faces
            smFaceSet,      //!< "faceSet" : subset with (volume) face set
            smFaceZone,     //!< "faceZone" : subset with (volume) zone faces
            smPatch,        //!< "patch" : subset with (volume) patches
        };

        //- List of selection mode type names
        static const Enum<selectionModeType> selectionModeTypeNames_;


protected:

    // Protected Data

        //- Time start
        scalar timeStart_;

        //- Duration
        scalar duration_;

        //- Face selection mode
        selectionModeType selectionMode_;

        //- Face selection names (for set, zone or patch selections)
        wordRes selectionNames_;

        //- Set of faces to apply source to
        labelList faces_;

        //- Sum of face area
        scalar A_;


    // Protected Functions

        //- Set face selection name from dictionary input
        void setSelection(const dictionary& dict);

        //- Set face selection based on user input selection mode
        void setFaceSelection();

        //- Recalculate the area
        void setArea();

        //- Zero all non-selected locations within field
        template<class Type>
        inline void subsetFilter(List<Type>& field) const;


public:

    //- Runtime type information
    TypeName("faceSetOption");


    // Constructors

        //- Construct from components
        faceSetOption
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~faceSetOption() = default;


    // Member Functions

        // Access

            //- Return const access to the time start
            inline scalar timeStart() const noexcept;

            //- Return const access to the duration
            inline scalar duration() const noexcept;

            //- Return true if within time limits
            inline bool inTimeLimits(const scalar timeValue) const;

            //- True if sub-selection should be used
            inline bool useSubMesh() const noexcept;

            //- Return the face selection mode
            inline selectionModeType selectionMode() const noexcept;

            //- Return const access to the selection names
            //- (set, zone or patch selection)
            inline const wordRes& selectionNames() const noexcept;

            //- Return const access to the first set/zone/patch name
            inline const wordRe& zoneName() const;

            //- Return const access to the total face area
            inline scalar A() const noexcept;

            //- Return const access to the local finite-area face selection
            inline const labelList& faces() const noexcept;


        // Edit

            //- Adjust the time start, return the old value
            inline scalar timeStart(scalar val) noexcept;

            //- Adjust the duration, return the old value
            inline scalar duration(scalar val) noexcept;


        // Checks

            //- Is the source active?
            virtual bool isActive();


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "faceSetOptionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
