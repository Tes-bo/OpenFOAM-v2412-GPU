/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017-2023 OpenCFD Ltd.
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
    Foam::fv::cellSetOption

Description
    Intermediate abstract class for handling
    cell-set options for the derived fvOptions.

Usage
    Minimal example by using \c constant/fvOptions:
    \verbatim
    fvOption1
    {
        // Mandatory entries
        selectionMode     <word>;

        // Optional entries
        timeStart         <scalar>;
        updateSelection   <bool>;

        // Conditional entries

            // when timeStart entry is present
            duration          <scalar>;

            // when selectionMode=cellSet
            cellSet           <word>;

            // when selectionMode=cellZone
            cellZone          <word>;

            //OR: cellZones   (<word> ...);

            // when selectionMode=points
            points            (<point1> <point2> ... <pointN>);

            // when selectionMode=movingPoints
            movingPoints
            (
                <word>        <Function1<vector>>;
                // e.g.
                point1        <Function1<vector>>;
                pointN        <Function1<vector>>;
                ...
            );

            // when selectionMode=geometric
            selection
            {
                topoSet1      <dictionary>;

                box1
                {
                    action  use;
                    source  box;
                    min     (-0.1 -0.01 -0.1);
                    max     (0.1 0.30 0.1);
                }
                ball
                {
                    action  add;
                    source  sphere;
                    origin  (-0.1 -0.01 -0.1);
                    radius  0.25;
                }
                ...
            }

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property      | Description                        | Type  | Reqd | Deflt
      selectionMode | Mode of cell selection - see below | word  | yes  | -
      timeStart     | Start time of fvOption            | scalar | no   | -1
      updateSelection | Flag to enable selection updates | bool  | no   | false
      duration      | Duration of fvOption execution <!--
                    --> starting from timeStart       | scalar | choice | 0
      cellSet       | Name of operand cellSet         | word   | choice | -
      cellZone      | Name of operand cellZone        | wordRe | choice | -
      cellZones     | Name of operand cellZones      | wordRes | choice | -
      points        | Set of points in global coordinate <!--
                    --> system                    | vectorList | choice | -
      movingPoints  | Set of moving points in global coordinate system <!--
                    -->                  | Function1\<vector\> | choice | -
      selection     | Dictionary of geometric selections | dict | choice | -
    \endtable

    Options for the \c selectionMode entry:
    \verbatim
      all       | Use all cells in the computational domain
      cellZone  | Use specified cellZone
      cellSet   | Use specified cellSet
      points    | Use cells containing a given set of points
      movingPoints | Use cells containing a given set of moving points
      geometric | Select cells based on topoSetCellSource actions
      cellType  | For overset only : use cells of type 'porous
    \endverbatim

    The inherited entries are elaborated in:
     - \link fvOption.H \endlink
     - \link Function1.H \endlink

Note
  - Source/sink options are to be added to the right-hand side of equations.
  - The geometric selection uses \c topoSetCellSource to select cells.
    Any \c searchableSurface selections must describe a closed volume.
    Ie, its \c hasVolumeType() method must be \c true.

See also
    Foam::cellBitSet::select

SourceFiles
    cellSetOption.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_fv_cellSetOption_H
#define Foam_fv_cellSetOption_H

#include "fvOption.H"
#include "fvMesh.H"
#include "dictionary.H"
#include "Time.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                         Class cellSetOption Declaration
\*---------------------------------------------------------------------------*/

class cellSetOption
:
    public fv::option
{
public:

    // Public Data

        //- Enumeration for selection mode types
        enum selectionModeType : char
        {
            smAll,          //!< "all" cells
            smCellSet,      //!< "cellSet"
            smCellZone,     //!< "cellZone"
            smPoints,       //!< "points"
            smMovingPoints, //!< "movingPoints"
            smGeometric,    //!< "geometric"
            smCellType      //!< "overset type cells"
        };

        //- List of selection mode type names
        static const Enum<selectionModeType> selectionModeTypeNames_;


protected:

    // Protected Data

        //- Cell selection mode
        selectionModeType selectionMode_;

        //- Flag to enable dictionary-based updates of selections
        bool updateSelection_;

        //- Start time of fvOption
        scalar timeStart_;

        //- Duration of fvOption execution starting from timeStart
        scalar duration_;

        //- Face selection names (for set or zone selections)
        wordRes selectionNames_;

        //- List of points for "points" selectionMode
        List<point> points_;

        //- List of points for "movingPoints" selectionMode
        PtrList<Function1<point>> movingPoints_;

        //- Dictionary entries for "geometric" (topoSetCellSource) selection
        dictionary geometricSelection_;

        //- Set of cells to apply source to
        labelList cells_;

        //- Sum of cell volumes
        scalar V_;


    // Protected Functions

        //- Set cell selection name or points selection from dictionary input
        void setSelection(const dictionary& dict);

        //- Set the cell selection based on user input selection mode
        void setCellSelection();

        //- Recalculate the volume
        void setVol();


public:

    //- Runtime type information
    TypeName("cellSetOption");


    // Constructors

        //- Construct from components
        cellSetOption
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~cellSetOption() = default;


    // Member Functions

        // Access

            //- Return const access to the time start
            inline scalar timeStart() const noexcept;

            //- Return const access to the duration
            inline scalar duration() const noexcept;

            //- True if within time limits
            inline bool inTimeLimits(const scalar timeValue) const;

            //- True if sub-selection should be used
            inline bool useSubMesh() const noexcept;

            //- Return the cell selection mode
            inline selectionModeType selectionMode() const noexcept;

            //- Return const access to the selection names
            //- (set or zone selection)
            inline const wordRes& selectionNames() const noexcept;

            //- Return const access to the first set/zone name
            inline const wordRe& zoneName() const;

            //- Return const access to the total cell volume
            inline scalar V() const noexcept;

            //- Return const access to the cell selection
            inline const labelList& cells() const noexcept;

            //- Return flag for selection updates
            bool isSelectionUpdated() const noexcept
            {
                return updateSelection_;
            }


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


    // Housekeeping

        //- The name of the cell set/zone [as a word]
        //- for "cellSet" / "cellZone" selection modes)
        const word& cellSetName() const { return zoneName(); }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cellSetOptionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
