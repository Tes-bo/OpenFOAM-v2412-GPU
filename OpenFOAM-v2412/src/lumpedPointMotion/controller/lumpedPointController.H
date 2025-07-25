/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 OpenCFD Ltd.
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
    Foam::lumpedPointController

Description
    Simple connectivity of point labels to specify a controller for lumped
    point movement.

    \heading Dictionary parameters
    \table
        Property    | Description                           | Required | Default
        pointLabels | List of point labels                  | yes |
    \endtable

Note
    If the calling program itself specified a point-label mapping
    (eg, original ids from FEA), these can be used initially and remapped.

SourceFiles
    lumpedPointController.C

\*---------------------------------------------------------------------------*/

#ifndef lumpedPointController_H
#define lumpedPointController_H

#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class dictionary;
template<class T> class Map;

/*---------------------------------------------------------------------------*\
                    Class lumpedPointController Declaration
\*---------------------------------------------------------------------------*/

class lumpedPointController
{
    // Private Data

        //- The lumped points associated with the controller
        labelList pointLabels_;


public:

    // Constructors

        //- Default construct
        lumpedPointController() noexcept;

        //- Copy construct from point ids
        explicit lumpedPointController(const labelUList& pointLabels);

        //- Move construct from point ids
        explicit lumpedPointController(labelList&& pointLabels);

        //- Construct from dictionary
        explicit lumpedPointController(const dictionary& dict);


        //- Factory method
        static autoPtr<lumpedPointController> New(const dictionary& dict)
        {
            return autoPtr<lumpedPointController>::New(dict);
        }


    //- Destructor
    ~lumpedPointController() = default;


    // Member Functions

        //- No controller points specified?
        bool empty() const
        {
            return pointLabels_.empty();
        }

        //- Number of controller points specified
        label size() const
        {
            return pointLabels_.size();
        }

        //- The controller points
        const labelList& pointLabels() const
        {
            return pointLabels_;
        }

        //- Check point mapping (the count) or remap the point labels
        void remapPointLabels
        (
            const label nPoints,
            const Map<label>& originalIds
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
