/*---------------------------------------------------------------------------* \
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2023 OpenCFD Ltd.
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
    Foam::PDRblock

Description
    A single block x-y-z rectilinear mesh addressable as i,j,k with
    simplified creation. Some of the input is similar to blockMeshDict,
    but since this specialization is for a single-block that is aligned
    with the x-y-z directions, it provides a different means of specifying
    the mesh.

    Dictionary controls
    \table
        Property    | Description                          | Required | Default
        x           | X-direction grid specification       | yes |
        y           | Y-direction grid specification       | yes |
        z           | Z-direction grid specification       | yes |
        scale       | Point scaling                        | no  | 1.0
        expansion   | Type of expansion (ratio/relative)   | no  | ratio
        boundary    | Boundary patches                     | yes |
        defaultPatch | Default patch specification         | no  |
    \endtable

    Grid coordinate controls
    \table
        Property| Description                               | Required | Default
        points  | Locations defining the mesh segment       | yes |
        nCells  | Divisions per mesh segment                | yes |
        ratios  | Expansion values per segment              | no  |
    \endtable

    A negative expansion value is trapped and treated as its reciprocal.
    by default, the expansion is as per blockMesh and represents the ratio
    of end-size / start-size for the section.
    Alternatively, the relative size can be given.

SourceFiles
    PDRblockI.H
    PDRblock.C
    PDRblockCreate.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_PDRblock_H
#define Foam_PDRblock_H

#include "ijkMesh.H"
#include "boundBox.H"
#include "pointField.H"
#include "faceList.H"
#include "Enum.H"
#include "vector2D.H"
#include "labelVector2D.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class IOobject;
class blockMesh;
class polyMesh;
class gradingDescriptors;

/*---------------------------------------------------------------------------*\
                           Class PDRblock Declaration
\*---------------------------------------------------------------------------*/

class PDRblock
:
    public ijkMesh
{
public:

    // Data Types

        //- The expansion type
        enum expansionType : uint8_t
        {
            EXPAND_UNIFORM,     //!< Uniform expansion (ie, no expansion)
            EXPAND_RATIO,       //!< End/start ratio
            EXPAND_RELATIVE     //!< Relative expansion ratio
        };

        //- Named enumerations for the expansion type
        const static Enum<expansionType> expansionNames_;


    // Public Classes

        //- Grid locations in an axis direction.
        //  The number of points is one larger than the number of elements
        //  it represents
        class location
        :
            public scalarList
        {
        public:

            //- Reset with min/max range and number of divisions
            void reset(const scalar low, const scalar upp, const label nCells);

            //- The location list is valid if it contains 2 or more points
            inline bool good() const noexcept;

            //- The number of cells in this direction.
            inline label nCells() const noexcept;

            //- The number of points in this direction.
            inline label nPoints() const noexcept;

            //- True if the location is within the range
            inline bool contains(const scalar p) const;

            //- The first() value is considered the min value.
            inline const scalar& min() const;

            //- The last() value is considered the max value.
            inline const scalar& max() const;

            //- Mid-point location, zero for an empty list.
            inline scalar centre() const;

            //- The difference between min/max values, zero for an empty list.
            inline scalar length() const;

            //- Check that element index is within valid range.
            inline void checkIndex(const label i) const;

            //- Cell size at element position.
            inline scalar width(const label i) const;

            //- Cell centre at element position.
            //  Treats -1 and nCells positions like a halo cell.
            inline scalar C(const label i) const;

            //- Return min/max edge lengths
            scalarMinMax edgeLimits() const;

            //- Return edge grading descriptors for the locations
            //  \see Foam::gradingDescriptor
            gradingDescriptors grading() const;

            //- Find the cell index enclosing this location
            //  \return -1 for out-of-bounds
            label findCell(const scalar p) const;

            //- Find the grid index, within the given tolerance
            //  Return -1 for out-of-bounds and -2 for a point that is
            //  within bounds, but not aligned with a grid point.
            label findIndex(const scalar p, const scalar tol) const;

            //- If out of range, return the respective min/max limits,
            //- otherwise return the value itself.
            //  If the range is invalid, always return the value.
            inline const scalar& clip(const scalar& val) const;
        };


        //- The begin/end nodes for each segment,
        //- with divisions and expansion for each segment
        //  Not normally used outside of PDRblock
        struct gridControl
        :
            public scalarList
        {
            //- The number of division per segment
            labelList divisions_;

            //- The expansion ratio per segment
            scalarList expansion_;

            //- Total number of cells in this direction
            label nCells() const;

            //- Return edge grading descriptors for the locations
            //  \see Foam::gradingDescriptor
            gradingDescriptors grading() const;

            //- Resize lists
            void resize(label len);

            //- Add point/divisions/expand to end of list (push_back)
            void append(const scalar p, label nDiv, scalar expRatio=1);

            //- Add point/divisions/expand to front of list (push_front)
            void prepend(const scalar p, label nDiv, scalar expRatio=1);

            //- Write as dictionary contents for specified vector direction
            void writeDict(Ostream& os, const direction cmpt) const;
        };


private:

    // Private Classes

        //- Extracted patch settings
        struct boundaryEntry
        {
            //- The patch name
            word name_;

            //- The patch type
            word type_;

            //- The patch size
            label size_;

            //- The associated block face ids [0,5]
            labelList faces_;
        };

        //- The begin/end nodes for each segment,
        //- with divisions and expansion for each segment
        struct outerControl
        {
            //- The control type
            enum controlType : uint8_t
            {
                OUTER_NONE = 0,     //!< No outer region
                OUTER_EXTEND,       //!< Extend inner region (orthogonal)
                OUTER_BOX,          //!< Cuboid
                OUTER_SPHERE        //!< Spherical
            };

            //- Named enumerations for the control type
            const static Enum<controlType> controlNames_;

            //- The control type
            controlType type_;

            //- The expansion type
            expansionType expandType_;

            //- True if on the ground
            bool onGround_;

            //- Relative size(s) for the outer region
            vector2D relSize_;

            //- Number of cells in outer region
            //  Generally only single component is used
            labelVector2D nCells_;

            //- Expansion ratio(s) for the outer region
            vector2D expansion_;


        // Constructors

            //- Default construct. NONE
            outerControl();


        // Member Functions

            //- Reset to default (NONE) values
            void clear();

            //- Is enabled (not NONE)
            bool active() const;

            //- Project on to sphere (is SPHERE)
            bool isSphere() const;

            //- Is the outer region on the ground?
            bool onGround() const;

            //- Define that the outer region is on the ground or not
            //  \return the old value
            bool onGround(const bool on);

            //- Read content from dictionary
            void read(const dictionary& dict);

            //- Report information about outer region
            void report(Ostream& os) const;
        };


    // Private Data

        //- Reference to mesh dictionary
        const dictionary& meshDict_;

        //- The grid controls in (i,j,k / x,y,z) directions.
        Vector<gridControl> control_;

        //- The grid points in all (i,j,k / x,y,z) directions,
        //- after applying the internal subdivisions.
        Vector<location> grid_;

        //- Control for the outer-region (if any)
        outerControl outer_;

        //- The mesh bounding box
        boundBox bounds_;

        //- The boundary patch information
        PtrList<boundaryEntry> patches_;

        //- The min/max edge lengths
        scalarMinMax edgeLimits_;

        //- Verbosity
        bool verbose_;


    // Private Member Functions

        //- Check that points increase monotonically
        static bool checkMonotonic
        (
            const direction cmpt,
            const UList<scalar>& pts
        );

        //- Add zero-sized patches with default naming/types
        void addDefaultPatches();

        //- Adjust sizing for updated grid points
        void adjustSizes();

        //- Read and define grid points in given direction
        void readGridControl
        (
            const direction cmpt,
            const dictionary& dict,
            const scalar scaleFactor = -1,
            expansionType expandType = expansionType::EXPAND_RATIO
        );

        //- Read "boundary" information
        void readBoundary(const dictionary& dict);

        //- Populate point field for the block
        void createPoints(pointField& pts) const;

        //- Add internal faces to lists.
        //  Lists must be properly sized!
        //  \return the number of faces added
        label addInternalFaces
        (
            faceList::iterator& faceIter,
            labelList::iterator& ownIter,
            labelList::iterator& neiIter
        ) const;

        //- Add boundary faces for the shape face to lists
        //  Lists must be properly sized!
        //  \return the number of faces added
        label addBoundaryFaces
        (
            const direction shapeFacei,
            faceList::iterator& faceIter,
            labelList::iterator& ownIter
        ) const;

        //- Obtain i,j,k index for cell enclosing this location
        //  \return false for out-of-bounds
        bool findCell(const point& pt, labelVector& pos) const;

        //- Obtain i,j,k grid index for point location
        //  \return false for out-of-bounds and off-grid
        bool gridIndex
        (
            const point& pt,
            labelVector& pos,
            const scalar tol
        ) const;

        //- The bounding box of the grid points
        static boundBox bounds
        (
            const scalarList& x, //!< X-points, monotonically increasing
            const scalarList& y, //!< Y-points, monotonically increasing
            const scalarList& z  //!< T-points, monotonically increasing
        );

        //- Equivalent edge grading descriptors in (x,y,z) directions.
        static Vector<gradingDescriptors> grading
        (
            const Vector<gridControl>& ctrl
        );

        //- Mesh sizes based on the controls
        static labelVector sizes
        (
            const Vector<gridControl>& ctrl
        );


    // Mesh Generation

        //- Create a blockMesh
        autoPtr<blockMesh> createBlockMesh(const IOobject& io) const;

        //- Create polyMesh via blockMesh
        autoPtr<polyMesh> meshBlockMesh(const IOobject& io) const;


public:

    // Static Member Functions

        //- Return a null PDRblock (reference to a nullObject).
        static const PDRblock& null() noexcept
        {
            return NullObjectRef<PDRblock>();
        }


    // Constructors

        //- Default construct, zero-size, inverted bounds etc
        PDRblock();

        //- Construct from components
        PDRblock
        (
            const UList<scalar>& xgrid,
            const UList<scalar>& ygrid,
            const UList<scalar>& zgrid
        );

        //- Construct from cube with specified griding
        PDRblock(const boundBox& box, const labelVector& nCells);

        //- Construct from dictionary
        explicit PDRblock(const dictionary& dict, bool verboseOutput=false);


    // Member Functions

        //- Read dictionary
        bool read(const dictionary& dict);

        //- Reset grid locations and mesh i-j-k sizing
        void reset
        (
            const UList<scalar>& xgrid,
            const UList<scalar>& ygrid,
            const UList<scalar>& zgrid
        );

        //- Reset cube and mesh i-j-k sizing
        void reset(const boundBox& box, const labelVector& nCells);


    // Access

        //- The grid point locations in the i,j,k (x,y,z) directions.
        const Vector<location>& grid() const noexcept { return grid_; }

        //- Equivalent edge grading descriptors in (x,y,z) directions.
        Vector<gradingDescriptors> grading() const;

        //- Equivalent edge grading descriptors in specified (x,y,z) direction.
        gradingDescriptors grading(const direction cmpt) const;


    // Mesh Information

        //- Mesh sizing as per ijkMesh
        using ijkMesh::sizes;

        //- The mesh bounding box
        const boundBox& bounds() const noexcept { return bounds_; }

        //- The min/max edge length
        const scalarMinMax& edgeLimits() const noexcept { return edgeLimits_; }

        //- Cell size in x-direction at i position.
        inline scalar dx(const label i) const;

        //- Cell size in x-direction at i position.
        inline scalar dx(const labelVector& ijk) const;

        //- Cell size in y-direction at j position.
        inline scalar dy(const label j) const;

        //- Cell size in y-direction at j position.
        inline scalar dy(const labelVector& ijk) const;

        //- Cell size in z-direction at k position.
        inline scalar dz(const label k) const;

        //- Cell size in z-direction at k position.
        inline scalar dz(const labelVector& ijk) const;

        //- Cell dimensions at i,j,k position.
        inline vector span(const label i, const label j, const label k) const;

        //- Cell dimensions at i,j,k position.
        inline vector span(const labelVector& ijk) const;

        //- Grid point at i,j,k position.
        inline point grid(const label i, const label j, const label k) const;

        //- Grid point at i,j,k position.
        inline point grid(const labelVector& ijk) const;

        //- Cell centre at i,j,k position.
        inline point C(const label i, const label j, const label k) const;

        //- Cell centre at i,j,k position.
        inline point C(const labelVector& ijk) const;

        //- Cell volume at i,j,k position.
        inline scalar V(const label i, const label j, const label k) const;

        //- Cell volume at i,j,k position.
        inline scalar V(const labelVector& ijk) const;

        //- Characteristic cell size at i,j,k position.
        //  This is the cubic root of the volume
        inline scalar width(const label i, const label j, const label k) const;

        //- Characteristic cell size at i,j,k position.
        //  This is the cubic root of the volume
        inline scalar width(const labelVector& ijk) const;


    // Searching

        //- Return i,j,k index for cell enclosing this location
        //  The value (-1,-1,-1) is returned for out-of-bounds (not found).
        labelVector findCell(const point& pt) const;

        //- Obtain i,j,k grid index for point location within specified
        //  relative tolerance of the min edge length
        //  The value (-1,-1,-1) is returned for out-of-bounds (not found).
        //  and off-grid
        labelVector gridIndex(const point& pt, const scalar relTol=0.01) const;


    // Mesh Generation

        //- Output content for an equivalent blockMeshDict
        //  Optionally generate header/footer content
        Ostream& blockMeshDict(Ostream& os, const bool withHeader=false) const;

        //- Content for an equivalent blockMeshDict
        dictionary blockMeshDict() const;

        //- Write an equivalent blockMeshDict
        void writeBlockMeshDict(const IOobject& io) const;

        //- Create polyMesh for grid definition and patch information
        autoPtr<polyMesh> mesh(const IOobject& io) const;

        //- Create polyMesh for inner-mesh only,
        //- ignore any outer block definitions
        autoPtr<polyMesh> innerMesh(const IOobject& io) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PDRblockI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
