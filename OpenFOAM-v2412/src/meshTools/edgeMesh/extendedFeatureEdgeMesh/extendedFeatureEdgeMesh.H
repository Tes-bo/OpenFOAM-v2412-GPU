/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2014 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::extendedFeatureEdgeMesh

Description

    extendedEdgeMesh + IO.

SourceFiles
    extendedFeatureEdgeMesh.C

\*---------------------------------------------------------------------------*/

#ifndef extendedFeatureEdgeMesh_H
#define extendedFeatureEdgeMesh_H

#include "extendedEdgeMesh.H"
#include "regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class objectRegistry;

/*---------------------------------------------------------------------------*\
                       Class extendedFeatureEdgeMesh Declaration
\*---------------------------------------------------------------------------*/

class extendedFeatureEdgeMesh
:
    public regIOobject,
    public extendedEdgeMesh
{

public:

    //- Runtime type information
    TypeName("extendedFeatureEdgeMesh");


    // Constructors

        //- Construct (read) given an IOobject
        explicit extendedFeatureEdgeMesh(const IOobject& io);

        //- Copy construct with IOobject
        extendedFeatureEdgeMesh
        (
            const IOobject& io,
            const extendedEdgeMesh&
        );

        //- Construct given a surface with selected edges,point
        //  (surfaceFeatures), an objectRegistry and a
        //  fileName to write to.
        //  Extracts, classifies and reorders the data from surfaceFeatures.
        extendedFeatureEdgeMesh
        (
            const surfaceFeatures& sFeat,
            const objectRegistry& obr,
            const fileName& sFeatFileName,
            const boolList& surfBaffleRegions
        );

        //- Construct from PrimitivePatch
        extendedFeatureEdgeMesh
        (
            const IOobject& io,
            const PrimitivePatch<faceList, pointField>& surf,
            const labelUList& featureEdges,
            const labelUList& regionFeatureEdges,
            const labelUList& featurePoints
        );

        //- Construct from all components
        extendedFeatureEdgeMesh
        (
            const IOobject& io,
            const pointField& pts,
            const edgeList& eds,
            label concaveStart,
            label mixedStart,
            label nonFeatureStart,
            label internalStart,
            label flatStart,
            label openStart,
            label multipleStart,
            const vectorField& normals,
            const List<sideVolumeType>& normalVolumeTypes,
            const vectorField& edgeDirections,
            const labelListList& normalDirections,
            const labelListList& edgeNormals,
            const labelListList& featurePointNormals,
            const labelListList& featurePointEdges,
            const labelList& regionEdges
        );


    //- Destructor
    virtual ~extendedFeatureEdgeMesh() = default;


    // IO

        //- Give precedence to the regIOobject write
        using regIOobject::write;

        //- ReadData function required for regIOobject read operation
        virtual bool readData(Istream&);

        //- WriteData function required for regIOobject write operation
        virtual bool writeData(Ostream&) const;

        //- Is object global
        virtual bool global() const
        {
            return true;
        }

        //- Return complete path + object name if the file exists
        //- either in the case/processor or case otherwise null
        virtual fileName filePath() const
        {
            return globalFilePath(type());
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Global file type for extendedFeatureEdgeMesh
template<>
struct is_globalIOobject<extendedFeatureEdgeMesh> : std::true_type {};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
