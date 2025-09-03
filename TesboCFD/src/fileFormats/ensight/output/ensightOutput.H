/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2022 OpenCFD Ltd.
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

Namespace
    Foam::ensightOutput

Description
    A collection of functions for writing ensight file content.

SourceFiles
    ensightOutput.C
    ensightOutputTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_ensightOutput_H
#define Foam_ensightOutput_H

#include "ensightOutputFwd.H"
#include "ensightFile.H"
#include "ensightGeoFile.H"
#include "ensightCells.H"
#include "ensightFaces.H"
#include "ensightPTraits.H"

#include "faceListFwd.H"
#include "cellListFwd.H"

#include "ListOps.H"
#include "ListListOps.H"
#include "IndirectList.H"
#include "CompactListList.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Local definitions, could be relocated to ListListOps directly

namespace Foam
{
namespace ListListOps
{

//- Return the sizes of the sub-lists
template<class T, class Addr, class AccessOp>
labelList subSizes
(
    const IndirectListBase<T, Addr>& lists,
    AccessOp aop
)
{
    labelList output(lists.size());
    auto out = output.begin();

    for (const T& sub : lists)
    {
        *out = aop(sub).size();
        ++out;
    }

    return output;
}


//- Inplace renumber the values (not the indices) of a list of lists.
//  Negative IntListType elements are left untouched.
template<class IntListType>
void inplaceRenumber
(
    const labelUList& oldToNew,
    IntListType& lists
)
{
    for (auto& sub : lists)
    {
        for (auto& item : sub)
        {
            if (item >= 0)
            {
                item = oldToNew[item];
            }
        }
    }
}

} // End namespace ListListOps
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class cellShape;
class polyMesh;

namespace ensightOutput
{

/*---------------------------------------------------------------------------*\
                                Geometry Output
\*---------------------------------------------------------------------------*/

//- Write list of faces
void writeFaceList
(
    ensightGeoFile& os,
    const UList<face>& faces,
    const label pointOffset = 0  //!< Additional point offset for each face
);

//- Write list of faces (indirect addressing)
void writeFaceList
(
    ensightGeoFile& os,
    const UIndirectList<face>& faces,
    const label pointOffset = 0  //!< Additional point offset for each face
);

//- Write list of faces (stored in compact form)
void writeFaceList
(
    ensightGeoFile& os,
    const CompactListList<label>& faces,
    const label pointOffset = 0  //!< Additional point offset for each face
);

//- Write cell connectivity via cell shapes
void writeCellShapes
(
    ensightGeoFile& os,
    const UList<cellShape>& shapes,
    const label pointOffset = 0  //!< Additional point offset
);


//- Write the point ids per poly element.
//  Points have been already renumbered
void writePolysPoints
(
    ensightGeoFile& os,
    const cellUList& meshCells,
    const labelUList& addr,     //!< Cell ids to write
    const faceUList& meshFaces,
    const labelUList& faceOwner
);

//- Write the point ids per poly element, with point renumbering
void writePolysPoints
(
    ensightGeoFile& os,
    const polyMesh& mesh,
    const labelUList& addr,     //!< Cell ids to write
    const labelList& pointMap   //!< Point map to use
);


//- Write the regular face connectivity for specified type and
//- and specified faces
void writeFaceConnectivity
(
    ensightGeoFile& os,
    const ensightFaces::elemType etype,
    const label nTotal,
    const UIndirectList<face>& faces,
    //! Prefer collective write?
    bool parallel
);


//- Write the regular face connectivity for specified type
void writeFaceConnectivity
(
    ensightGeoFile& os,
    const ensightFaces::elemType etype,
    const label nTotal,
    const faceUList& faces,
    //! Prefer collective write?
    bool parallel
);


//- Write the face connectivity for the part
void writeFaceConnectivity
(
    ensightGeoFile& os,
    const ensightFaces& part,
    const faceUList& faces,
    //! Prefer collective write?
    bool parallel
);


//- Write the \b presorted face connectivity for the part
//  This is a special case when the list of faces is already in
//  ensight sorted order
void writeFaceConnectivityPresorted
(
    ensightGeoFile& os,
    const ensightFaces& part,
    const faceUList& faces,
    //! Prefer collective write?
    bool parallel
);


/*---------------------------------------------------------------------------*\
                                Field Output
\*---------------------------------------------------------------------------*/

//- Write a field of cell values as an indirect list,
//- using the cell ids from ensightCells
template<class Type>
bool writeField
(
    //! Component scratch buffer
    ensightOutput::floatBufferType& scratch,

    //! Output file (must be valid on master)
    ensightFile& os,

    //! The field content to be written
    const Field<Type>& fld,

    //! The addressing (element-wise) into the field
    const ensightCells& part,

    //! Prefer collective write?
    bool parallel
);

//- Write a field of faces values as an indirect list,
//- using the face ids from ensightFaces
template<class Type>
bool writeField
(
    //! Component scratch buffer
    ensightOutput::floatBufferType& scratch,

    //! Output file (must be valid on master)
    ensightFile& os,

    //! The field content to be written
    const Field<Type>& fld,

    //! The addressing (element-wise) into the field
    const ensightFaces& part,

    //! Prefer collective write?
    bool parallel
);


//- Write a field of cell values as an indirect list,
//- using the cell ids from ensightCells
template<class Type>
bool writeField
(
    //! Output file (must be valid on master)
    ensightFile& os,

    //! The field content to be written
    const Field<Type>& fld,

    //! The addressing (element-wise) into the field
    const ensightCells& part,

    //! Prefer collective write?
    bool parallel
)
{
    ensightOutput::floatBufferType scratch;
    return ensightOutput::writeField(scratch, os, fld, part, parallel);
}


//- Write a field of faces values as an indirect list,
//- using the face ids from ensightFaces
template<class Type>
bool writeField
(
    //! Output file (must be valid on master)
    ensightFile& os,

    //! The field content to be written
    const Field<Type>& fld,

    //! The addressing (element-wise) into the field
    const ensightFaces& part,

    //! Prefer collective write?
    bool parallel
)
{
    ensightOutput::floatBufferType scratch;
    return ensightOutput::writeField(scratch, os, fld, part, parallel);
}


/*---------------------------------------------------------------------------*\
                        Namespace ensightOutput::Detail
\*---------------------------------------------------------------------------*/

//- \brief Implementation details and output backends that would not normally
//- be called directly by a user.

namespace Detail
{

//- Return sizes of faces in the list
labelList getFaceSizes(const UList<face>& faces);

//- Return sizes of faces in the list
labelList getFaceSizes(const UIndirectList<face>& faces);

//- The number of faces per poly element
labelList getPolysNFaces(const polyMesh& mesh, const labelUList& addr);

//- The number of points for each face of the poly elements
labelList getPolysNPointsPerFace(const polyMesh& mesh, const labelUList& addr);

//- Generate 0-based point ids for each poly element face
//  The returned CompactListList is divided per output face
CompactListList<label> getPolysFacePoints
(
    const polyMesh& mesh,
    const labelUList& addr,     //!< Cell ids to write
    const labelList& pointMap   //!< Point map to use
);


//- Write CompactListList<label> by components
void writeLabelListList
(
    ensightGeoFile& os,
    const labelUList& offsets,
    const labelUList& values,
    const label pointOffset
);


//- Write a list of faces or cell shapes with one-entity per line
template<class LabelListListType>
void writeLabelListList
(
    ensightGeoFile& os,
    const LabelListListType& listOfLists,
    const label pointOffset
);


//- Copy specified field component into a scalar buffer.
//- Works for various lists types. Must be adequately sized before calling
template<template<typename> class FieldContainer, class Type>
void copyComponent
(
    //! Input field data
    const FieldContainer<Type>& input,

    //! Component to be extracted
    const direction cmpt,

    //! [out] Component scratch buffer
    UList<float>& cmptBuffer
);


//- Write coordinates (component-wise) for the given part.
//
//  Has internal check for (nPoints != 0)
template<template<typename> class FieldContainer>
bool writeCoordinates
(
    //! Output file (must be valid on master)
    ensightGeoFile& os,

    //! The ensight part id
    const label partId,

    //! The ensight part description
    const word& partName,

    //! The total number of points
    const label nPoints,

    //! The point field to be written
    const FieldContainer<Foam::point>& fld,

    //! Prefer collective write?
    bool parallel
);


//- Write field content (component-wise) for the given ensight element type.
template<template<typename> class FieldContainer, class Type>
void writeFieldComponents
(
    //! Component scratch buffer
    ensightOutput::floatBufferType& scratch,

    //! Output file (must be valid on master)
    ensightFile& os,

    //! The ensight element type (ignored if nullptr)
    const char* key,

    //! The field content to be written for this element type
    const FieldContainer<Type>& fld,

    //! Prefer collective write?
    bool parallel
);

//- Write field content (component-wise) for the given ensight element type.
template<template<typename> class FieldContainer, class Type>
void writeFieldComponents
(
    //! Output file (must be valid on master)
    ensightFile& os,

    //! The ensight element type (can be nullptr)
    const char* key,

    //! The field content to be written for this element type
    const FieldContainer<Type>& fld,

    //!< Prefer collective write?
    bool parallel
)
{
    ensightOutput::floatBufferType scratch;
    Detail::writeFieldComponents(scratch, os, key, fld, parallel);
}


//- Write a sub-field of faces values as an indirect list,
//- using the sub-list sizing information from ensightFaces
template<class Type>
bool writeFaceSubField
(
    //! Component scratch buffer
    ensightOutput::floatBufferType& scratch,

    //! Output file (must be valid on master)
    ensightFile& os,

    //! The field content to be written
    const Field<Type>& fld,

    //! The addressing (element-wise) into the field
    const ensightFaces& part,

    //! Prefer collective write?
    bool parallel
);


//- Write a sub-field of faces values as an indirect list,
//- using the sub-list sizing information from ensightFaces
template<class Type>
bool writeFaceSubField
(
    //! Output file (must be valid on master)
    ensightFile& os,

    //! The field content to be written
    const Field<Type>& fld,

    //! The addressing (element-wise) into the field
    const ensightFaces& part,

    //! Prefer collective write?
    bool parallel
)
{
    ensightOutput::floatBufferType scratch;
    return Detail::writeFaceSubField(scratch, os, fld, part, parallel);
}


//- Write a field of faces values as an indirect list,
//- using the face order from ensightFaces
template<class Type>
bool writeFaceLocalField
(
    //! Component scratch buffer
    ensightOutput::floatBufferType& scratch,

    //! Output file (must be valid on master)
    ensightFile& os,

    //! The field content to be written
    const Field<Type>& fld,

    //! The addressing (element-wise) into the field
    const ensightFaces& part,

    //! Prefer collective write?
    bool parallel
);

//- Write a field of faces values as an indirect list,
//- using the face order from ensightFaces
template<class Type>
bool writeFaceLocalField
(
    //! Output file (must be valid on master)
    ensightFile& os,

    //! The field content to be written
    const Field<Type>& fld,

    //! The addressing (element-wise) into the field
    const ensightFaces& part,

    //! Prefer collective write?
    bool parallel
)
{
    ensightOutput::floatBufferType scratch;
    return Detail::writeFaceLocalField(scratch, os, fld, part, parallel);
}


} // End namespace Detail


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace ensightOutput
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ensightOutputTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
