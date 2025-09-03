/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 M. Janssens
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
    Foam::processorColour

Description
    Colouring processors such that no neighbours have the same colour

SourceFiles
    processorColour.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_processorColour_H
#define Foam_processorColour_H

#include "MeshObject.H"
#include "lduMesh.H"
#include "lduInterfacePtrsList.H"
#include "primitiveFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class lduMatrix;

/*---------------------------------------------------------------------------*\
                      Class processorColour Declaration
\*---------------------------------------------------------------------------*/

class processorColour
:
    public MeshObject<lduMesh, MoveableMeshObject, processorColour>,
    public labelList
{
    // Private Typedefs

        typedef MeshObject
        <
            lduMesh,
            MoveableMeshObject,
            processorColour
        > MeshObject_type;


    // Private Member Functions

        //- No copy construct
        processorColour(const processorColour&) = delete;

        //- No copy assignment
        void operator=(const processorColour&) = delete;


protected:

    // Protected data

        //- Max number of colours
        label nColours_;


    // Protected Member Fucntions

        static void walkFront
        (
            const lduMesh& mesh,
            DynamicList<label>& front,
            labelList& cellColour
        );

public:

    //- Runtime type information
    TypeName("processorColour");


    // Constructors

        //- Construct given mesh
        explicit processorColour(const lduMesh& mesh);


    // Selectors

        //- Should use the MeshObject provided one but that needs a
        //  mesh.name() ...
        static const processorColour& New(const lduMesh& mesh);


    //- Destructor
    virtual ~processorColour() = default;


    // Member Functions

        // Access

            label nColours() const
            {
                return nColours_;
            }

            label myColour() const
            {
                return operator[](Pstream::myProcNo());
            }


        //- Calculate processor colouring from processor connectivity. Sets
        //- colour per processor such that no neighbouring processors have the
        //- same colour. Returns number of colours used.
        static label colour(const lduMesh& mesh, labelList& procColour);

        //- Calculate (locally) per cell colour according to walking from
        //- global patches. Returns number of colours used.
        //- Note: asssumes all non-processor interfaces are global.
        static label cellColour
        (
            const lduMesh& mesh,
            labelList& cellColour,
            labelList& patchToColour
        );

        virtual bool movePoints()
        {
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
