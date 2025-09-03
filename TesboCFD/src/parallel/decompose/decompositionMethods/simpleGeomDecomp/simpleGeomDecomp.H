/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2018-2023 OpenCFD Ltd.
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
    Foam::simpleGeomDecomp

Description
    Simple geometric decomposition, selectable as \c simple

    Method coefficients:
    \table
        Property  | Description                             | Required | Default
        n         | (nx ny nz)                              | yes |
        order     | order of operation (unused)             | no  | xyz
        delta     | delta (jitter) for rotation matrix      | no  | 0.001
        transform | cartesian coordinate transformation     | no  |
    \endtable

SourceFiles
    simpleGeomDecomp.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_simpleGeomDecomp_H
#define Foam_simpleGeomDecomp_H

#include "geomDecomp.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class simpleGeomDecomp Declaration
\*---------------------------------------------------------------------------*/

class simpleGeomDecomp
:
    public geomDecomp
{
    // Private Member Functions

        labelList decomposeOneProc(const pointField& points) const;

        labelList decomposeOneProc
        (
            const pointField& points,
            const scalarField& weights
        ) const;


public:

    // Generated Methods

        //- No copy construct
        simpleGeomDecomp(const simpleGeomDecomp&) = delete;

        //- No copy assignment
        void operator=(const simpleGeomDecomp&) = delete;


    //- Runtime type information
    TypeName("simple");


    // Constructors

        //- Construct with number of x/y/z division
        //- (no coefficients or constraints)
        explicit simpleGeomDecomp(const Vector<label>& divisions);

        //- Construct given decomposition dictionary and optional region name
        explicit simpleGeomDecomp
        (
            const dictionary& decompDict,
            const word& regionName = ""
        );


    //- Destructor
    virtual ~simpleGeomDecomp() = default;


    // Member Functions

        //- Sends all points to the master for decomposition.
        virtual bool parallelAware() const
        {
            return true;
        }

        //- Return for every coordinate the wanted processor number.
        //- using uniform or specified point weights.
        virtual labelList decompose
        (
            const pointField& points,
            const scalarField& weights = scalarField::null()
        ) const;

        //- Decompose with uniform or specified point weights.
        virtual labelList decompose
        (
            const polyMesh& mesh,  //!< To check mesh dimensions
            const pointField& points,
            const scalarField& pointWeights = scalarField::null()
        ) const
        {
            checkDecompositionDirections(mesh.geometricD());
            return decompose(points, pointWeights);
        }

        //- Explicitly provided connectivity
        virtual labelList decompose
        (
            const CompactListList<label>& globalCellCells,  //!< unused
            const pointField& cc,
            const scalarField& cWeights = scalarField::null()
        ) const
        {
            return decompose(cc, cWeights);
        }

        //- Explicitly provided connectivity
        virtual labelList decompose
        (
            const labelListList& globalCellCells,  //!< unused
            const pointField& cc,
            const scalarField& cWeights = scalarField::null()
        ) const
        {
            return decompose(cc, cWeights);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
