/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2017-2020 OpenCFD Ltd.
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
    Foam::cellModel

Description
    Maps a geometry to a set of cell primitives.

    This enables geometric cell data to be calculated without access
    to the primitive geometric level. This means mapping a 3D
    geometry to a set of pyramids which are each described by a cell
    face and the cell centre point.

    Also includes a static collection of cell models (normally loaded from
    etc/cellModels), and a means of looking them up.

SourceFiles
    cellModelI.H
    cellModel.C
    cellModels.C
    cellModelIO.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_cellModel_H
#define Foam_cellModel_H

#include "pointField.H"
#include "edgeList.H"
#include "faceList.H"
#include "InfoProxy.H"
#include "autoPtr.H"
#include "PtrList.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class cellModel;
Ostream& operator<<(Ostream& os, const cellModel& cm);

template<>
Ostream& operator<<(Ostream&, const InfoProxy<cellModel>&);


/*---------------------------------------------------------------------------*\
                          Class cellModel Declaration
\*---------------------------------------------------------------------------*/

class cellModel
{
public:

        //- Enumeration of commonly used cellModel types.
        //  The indices must match those in "etc/cellModels"
        enum modelType
        {
            UNKNOWN = 0,     //!< unknown
            HEX = 3,         //!< hex
            WEDGE = 4,       //!< wedge
            PRISM = 5,       //!< prism
            PYR = 6,         //!< pyr
            TET = 7,         //!< tet
            SPLITHEX = 8,    //!< splitHex
            TETWEDGE = 9,    //!< tetWedge
        };

        //- Names of commonly used cellModels corresponding to modelType.
        //  The names must match those in "etc/cellModels"
        static const Enum<modelType> modelNames;


    // Lookup Methods

        //- Look up pointer to cellModel by enumeration, or nullptr on failure.
        static const cellModel* ptr(const modelType model);

        //- Look up pointer to cellModel by name, or nullptr on failure.
        static const cellModel* ptr(const word& modelName);

        //- Look up pointer to cellModel by index, or nullptr on failure
        static const cellModel* ptr(const label modelIndex);


        //- Look up reference to cellModel by enumeration. Fatal on failure
        static const cellModel& ref(const modelType model);

        //- Look up reference to cellModel by name. Fatal on failure
        static const cellModel& ref(const word& modelName);

        //- Look up reference to cellModel by index. Fatal on failure
        static const cellModel& ref(const label modelIndex);


private:

    // Private Static Data

        //- PtrList of predefined models
        static PtrList<cellModel> models_;

        //- Lookup of model pointers (in models_) by index
        static List<const cellModel*> modelPtrs_;


    // Private Data

        //- (Unique) model name
        word name_;

        //- Index in the model list
        label index_;

        //- Number of points in the model which determines the geometry
        label nPoints_;

        //- Faces of the model
        faceList faces_;

        //- Edges of the model
        edgeList edges_;


    // Private Member Functions

        //- Construct from central "etc/cellModels" file.
        static void constructModels();


public:

    // Constructors

        //- Construct from Istream
        explicit cellModel(Istream& is);

        //- Return a new cellModel created from Istream
        static autoPtr<cellModel> New(Istream& is)
        {
            return autoPtr<cellModel>::New(is);
        }

        //- Return clone
        autoPtr<cellModel> clone() const
        {
            return autoPtr<cellModel>::New(*this);
        }


    // Member Functions

        //- Return model name
        inline const word& name() const noexcept;

        //- Return index of model in the model list
        inline label index() const noexcept;

        //- Return number of points
        inline label nPoints() const noexcept;

        //- Return number of edges
        inline label nEdges() const noexcept;

        //- Return number of faces
        inline label nFaces() const noexcept;

        //- Return a raw list of model edges
        inline const edgeList& modelEdges() const noexcept;

        //- Return a raw list of model faces
        inline const faceList& modelFaces() const noexcept;

        //- Return list of cell edges
        inline edgeList edges(const labelUList& pointLabels) const;

        //- Return list of cell faces
        inline faceList faces(const labelUList& pointLabels) const;

        //- Return the cell edge for specified model edge
        inline Foam::edge edge
        (
            const label modelEdgei,
            const labelUList& pointLabels
        ) const;

        //- Return the cell face for specified model face
        inline Foam::face face
        (
            const label modelFacei,
            const labelUList& pointLabels
        ) const;


        //- Centroid of the cell
        vector centre
        (
            const labelList& pointLabels,
            const UList<point>& points
        ) const;

        //- Cell volume
        scalar mag
        (
            const labelList& pointLabels,
            const UList<point>& points
        ) const;


        //- Return info proxy,
        //- used to print information to a stream
        InfoProxy<cellModel> info() const noexcept { return *this; }

        //- The writeData member function required by regIOobject
        bool writeData(Ostream& os) const
        {
            os << *this;
            return os.good();
        }


    // Ostream Operator

        friend Ostream& operator<<(Ostream& os, const cellModel& cm);
};


// Global Operators

//- Equality: true when model pointers are identical
inline bool operator==(const cellModel& lhs, const cellModel& rhs) noexcept;

//- Inequality: true when model pointers are not identical
inline bool operator!=(const cellModel& lhs, const cellModel& rhs) noexcept;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cellModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
