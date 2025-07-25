/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
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
    Foam::leastSquaresFaVectors

Description
    Least-squares gradient scheme vectors for the Finite Area method

Author
    Hrvoje Jasak, Wikki Ltd.

SourceFiles
    leastSquaresFaVectors.C

\*---------------------------------------------------------------------------*/

#ifndef leastSquaresFaVectors_H
#define leastSquaresFaVectors_H

#include "MeshObject.H"
#include "faMesh.H"
#include "edgeFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class leastSquaresFaVectors Declaration
\*---------------------------------------------------------------------------*/

class leastSquaresFaVectors
:
    public MeshObject<faMesh, MoveableMeshObject, leastSquaresFaVectors>
{
    // Private Typedefs

        typedef MeshObject
        <
            faMesh,
            MoveableMeshObject,
            leastSquaresFaVectors
        > MeshObject_type;


    // Private Data

        //- Least-squares gradient vectors
        mutable std::unique_ptr<edgeVectorField> pVectorsPtr_;
        mutable std::unique_ptr<edgeVectorField> nVectorsPtr_;


    // Private Member Functions

        //- Construct Least-squares gradient vectors
        void makeLeastSquaresVectors() const;


public:

    // Declare name of the class and its debug switch
    TypeName("leastSquaresFaVectors");


    // Constructors

        //- Construct given an faMesh
        explicit leastSquaresFaVectors(const faMesh&);


    //- Destructor
    virtual ~leastSquaresFaVectors() = default;


    // Member functions

        //- Return reference to owner least square vectors
        const edgeVectorField& pVectors() const;

        //- Return reference to neighbour least square vectors
        const edgeVectorField& nVectors() const;


        //- Delete the least square vectors when the mesh moves
        virtual bool movePoints();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
