/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2022 OpenCFD Ltd.
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
    Foam::decompositionConstraints::geometric

Description
    Keep faces together based on geometric considerations from a
    list of topoSetFaceSource.
    The faces selected (inside) of each source are to be kept together
    during the decomposition.

    \verbatim
    constraints
    {
        geom1
        {
            type    geometric;

            grow    false;

            selection
            {
                box1
                {
                    source  box;
                    min     (-0.1 -0.01 -0.1);
                    max     (0.1 0.30 0.1);
                }
                ball
                {
                    source  sphere;
                    origin  (-0.1 -0.01 -0.1);
                    radius  0.25;
                }
                blob
                {
                    source      surface;
                    surfaceType triSurfaceMesh;
                    surfaceName blob.obj;
                }
            }
        }
    }
    \endverbatim

    \heading Dictionary parameters
    \table
        Property    | Description                           | Required | Default
        type        | geometric                             | yes |
        grow        | Grow cells with partial connectivity  | no  | false
        selection   | Dictionary of geometric selections    | yes |
    \endtable

Note
    Any searchableSurface selections must describe a closed volume.
    Ie, its hasVolumeType() method must be true.

    The combined face selections are treated as an OR operation.

    The "grow" mode includes an additional check to include cell faces
    for any cell that already has two or more of its faces "unblocked".
    This could indicate a connection over a corner, but does not distinguish
    between connectivity introduced by the constraint and the connectivity
    defined by other constraints.

SourceFiles
    geometricConstraint.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_geometricConstraint_H
#define Foam_geometricConstraint_H

#include "decompositionConstraint.H"
#include "PtrList.H"
#include "topoSetFaceSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace decompositionConstraints
{

/*---------------------------------------------------------------------------*\
                          Class geometric Declaration
\*---------------------------------------------------------------------------*/

class geometric
:
    public decompositionConstraint
{
    // Private Data

        //- Face selections as topoSetFaceSource
        PtrList<topoSetFaceSource> sources_;

        //- Dictionary of face selections for topoSetFaceSource
        dictionary selection_;

        //- Additional check of cell connection (via corners?)
        bool grow_;


public:

    //- Runtime type information
    TypeName("geometric");


    // Constructors

        //- Construct with constraint dictionary
        explicit geometric(const dictionary& dict);

        //- Move construct from components.
        //  (topoSetFaceSource cannot be cloned)
        explicit geometric(PtrList<topoSetFaceSource>&& selections);


    //- Destructor
    virtual ~geometric() = default;


    // Member Functions

        //- Add this constraint to list of constraints
        virtual void add
        (
            const polyMesh& mesh,
            boolList& blockedFace,
            PtrList<labelList>& specifiedProcessorFaces,
            labelList& specifiedProcessor,
            List<labelPair>& explicitConnections
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace decompositionConstraints
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
