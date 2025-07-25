/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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
    Foam::fv::MapFieldConstraint

Group
    grpFvOptionsConstraints

Description
    The \c MapFieldConstraint constrains values of given fields of \c Type
    with a source field from an external mesh, where
    \c \<Type\>={scalar,vector,sphericalTensor,symmTensor,tensor}.
    Optionally, the source field can be translated and/or rotated as a function
    of time.

Usage
    Minimal example by using \c constant/fvOptions:
    \verbatim
    \<Type\>MapFieldConstraint1
    {
        // Mandatory entries
        type                \<Type\>MapFieldConstraint;
        field               <word>;
        srcMesh             <fileName>;
        mapMethod           <word>;

        // Optional entries
        consistent          <bool>;
        patchMapMethod      <word>;
        transform
        {
            // Optional entries
            position    <Function1<vector>>;
            origin      <vector>;

            direction   <Function1<vector>>;
            normal      <vector>;
        }

        // Conditional entries

            // when consistent=false
            patchMap        <HashTable<word>>;  // (<patchSrc> <patchTgt>);
            cuttingPatches  <wordList>;  // (<patchTgt1> ... <patchTgtN>);

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                           | Type | Reqd | Deflt
      type      | Type name: \<Type\>MapFieldConstraint | word | yes  | -
      field     | Name of operand field                 | word | yes  | -
      srcMesh  | Directory path to mesh to map from  | fileName | yes | -
      mapMethod | Mapping method                        | word | yes  | -
      consistent | Flag to determine if meshes have consistent boundaries <!--
                -->                                     | bool | no   | false
      patchMapMethod | Name of patch-map method         | word | no   | -
      patchMap     | Coincident source/target patches in two cases   <!--
                   -->                           | wordHashTable | no | -
      cuttingPatches | Target patches cutting the source domain      <!--
                   -->                                | wordList | no | -
      transform    | Transform settings for source mesh points       <!--
                   -->                                  | dict | no   | -
      position     | Position of source mesh as a function of time   <!--
                   -->                     | Function1\<vector\> | no | -
      direction    | Direction of source mesh as a function of time  <!--
                   -->                     | Function1\<vector\> | no | -
      origin       | Origin of source mesh              | vector | no | -
      normal       | Normal of reference plane representing source mesh <!--
                   -->                                  | vector | no | -
    \endtable

    The inherited entries are elaborated in:
      - \link fvOption.H \endlink
      - \link Function1.H \endlink

SourceFiles
    MapFieldConstraint.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_fv_MapFieldConstraint_H
#define Foam_fv_MapFieldConstraint_H

#include "fvOption.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class meshToMesh;
template<class Type> class Function1;

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class MapFieldConstraint Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class MapFieldConstraint
:
    public fv::option
{
    // Private Classes

        class transform
        {
            // Private Data

                //- Position of source mesh as a function of time
                autoPtr<Function1<point>> positionPtr_;

                //- Direction of source mesh as a function of time
                autoPtr<Function1<point>> directionPtr_;

                //- Cached points of source mesh
                pointField points_;

                //- Origin of source mesh
                point origin_;

                //- Normal of reference plane representing source mesh
                vector normal_;

                //- Flag to deduce if transformation is active
                bool active_;


        public:

            // Constructors

                //- Default construct
                transform();

                //- No copy construct
                transform(const transform&) = delete;

                //- No copy assignment
                void operator=(const transform&) = delete;


            // Member Functions

            // Access

                //- Return flag to deduce if transformation is active
                bool isActive() const noexcept { return active_; }


            // Evaluation

                //- Translate source mesh as a function of time
                void translate(refPtr<fvMesh>& srcMeshPtr, const scalar time);

                //- Rotate source mesh as a function of time
                void rotate(refPtr<fvMesh>& srcMeshPtr, const scalar time);


            // I-O

                //- Initialize the class members
                bool initialize(const fvMesh& srcMesh, const dictionary& dict);
        };


    // Private Data

        //- Transformation settings for source mesh
        transform transform_;

        //- Time database for source mesh to map from
        autoPtr<Time> srcTimePtr_;

        //- Source mesh to map from
        refPtr<fvMesh> srcMeshPtr_;

        //- Mesh-to-mesh interpolation from source mesh to target mesh
        autoPtr<meshToMesh> interpPtr_;

        //- List of coincident source/target patches in two cases
        HashTable<word> patchMap_;

        //- Set of cells to apply source to
        labelList cells_;

        //- List of names of target patches cutting the source domain
        wordList cuttingPatches_;

        //- Name of map method
        word mapMethodName_;

        //- Name of patch-map method
        word patchMapMethodName_;

        //- Flag to determine if meshes have consistent boundaries
        bool consistent_;


    // Private Member Functions

        //- Helper function to set source mesh
        //  Fetch fvMesh from a given Time database
        //  Otherwise, load it from disk and cache it to the database
        void setSourceMesh
        (
            refPtr<fvMesh>& meshRef,
            const autoPtr<Time>& runTimePtr
        );

        //- Helper function to create the mesh-to-mesh interpolation
        void createInterpolation
        (
            const fvMesh& srcMesh,
            const fvMesh& tgtMesh
        );

        //- Return requested field from object registry
        //- otherwise read it from disk and register it to the object registry
        template<class VolFieldType>
        VolFieldType& getOrReadField
        (
            const fvMesh& thisMesh,
            const word& fieldName
        ) const;

        //- Return the local cell indices of the target mesh
        labelList tgtCellIDs() const;


public:

    //- Runtime type information
    TypeName("MapFieldConstraint");


    // Constructors

        //- Construct from components
        MapFieldConstraint
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        MapFieldConstraint(const MapFieldConstraint&) = delete;

        //- No copy assignment
        void operator=(const MapFieldConstraint&) = delete;


    //- Destructor
    virtual ~MapFieldConstraint() = default;


    // Member Functions

        //- Read source dictionary
        virtual bool read(const dictionary& dict);

        //- Set value on field
        virtual void constrain(fvMatrix<Type>& eqn, const label);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "MapFieldConstraint.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
