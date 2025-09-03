/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2024 OpenCFD Ltd.
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
    Foam::coordSystem

Description
    Namespace for coordinate systems.

Class
    Foam::coordinateSystem

Description
    Base class for coordinate system specification,
    the default coordinate system type is
    \link coordSystem::cartesian cartesian \endlink.

    All systems are defined by an origin point and a coordinate rotation
    By default, the \link coordinateRotations::axes axes \endlink
    specification can be used directly as part of the
    coordinate system specification. For example,

    \verbatim
    coordinateSystem
    {
        origin  (0 0 0);
        e1      (0 1 0);
        e3      (1 0 0);
    }
    \endverbatim

    However, a more verbose format with rotation provided as a dictionary entry
    is possible:
    \verbatim
    coordinateSystem
    {
        type    cartesian;
        origin  (0 0 0);
        rotation
        {
            type    axes;
            e1      (0 1 0);
            e3      (1 0 0);
        }
    }
    \endverbatim

    It also also possible to use the compact (single-dictionary) form
    and specific a different type of rotation:
    \verbatim
    coordinateSystem
    {
        type        cartesian;
        origin      (0 0 0);
        rotation    euler;
        angles      (90 0 0);
    }
    \endverbatim

    This last form can be particularly readable for an identity rotation:
    coordinateSystem
    {
        type        cartesian;
        origin      (0 0 0);
        rotation    none;
    }
    \endverbatim

    Types of coordinateRotation:
      -# \link coordinateRotations::identity none \endlink
      -# \link coordinateRotations::axes axes \endlink
      -# \link coordinateRotations::axisAngle axisAngle \endlink
      -# \link coordinateRotations::euler euler \endlink
      -# \link coordinateRotations::starcd starcd \endlink

    Type of coordinateSystem:
      -# \link coordSystem::cartesian cartesian \endlink
      -# \link coordSystem::cylindrical cylindrical \endlink
      -# \link coordSystem::indirect indirect \endlink (references
      an entry in coordinateSystems).

SourceFiles
    coordinateSystem.C
    coordinateSystemNew.C
    coordinateSystemTemplates.C
    coordinateSystemTransform.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_coordinateSystem_H
#define Foam_coordinateSystem_H

#include "vector.H"
#include "point.H"
#include "tensor.H"
#include "vectorField.H"
#include "pointField.H"
#include "tensorField.H"
#include "pointIndList.H"
#include "coordinateRotation.H"
#include "autoPtr.H"
#include "tmp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class coordinateSystem;
class objectRegistry;

namespace coordSystem
{
class indirect;
}


/*---------------------------------------------------------------------------*\
                      Class coordinateSystem Declaration
\*---------------------------------------------------------------------------*/

class coordinateSystem
{
    // Private Member Functions

        //- Select construct specified coordinate system type
        //
        //  \param modelType Eg, cartesian, cylindrical, indirect
        //      An empty modelType will be treated as "cartesian".
        //  \param dict The specifications
        //  \param readOrigin The preferred handling for reading 'origin'
        //  \param obrPtr Optional reference to an objectRegistry
        //         (for indirect entries).
        static autoPtr<coordinateSystem> New
        (
            const word& modelType,
            const dictionary& dict,
            IOobjectOption::readOption readOrigin,
            const objectRegistry* obrPtr
        );

        //- Select constructed from dictionary
        //
        //  \param dict The top-level dictionary contents
        //  \param dictName The sub-dictionary name to use for the
        //      coordinate system specifications. An empty name invokes
        //      an implicit search for a "coordinateSystem" sub-dictionary
        //      which is convenient and provides compatibility with previous
        //      versions (1806 and earlier).
        //  \param readOrigin The preferred handling for reading 'origin'
        //  \param obrPtr The objectRegistry for lookup of indirect entries.
        static autoPtr<coordinateSystem> New
        (
            const dictionary& dict,
            const word& dictName,
            IOobjectOption::readOption readOrigin,
            const objectRegistry* obrPtr
        );

        //- Optional select construct from dictionary
        //
        //  \param dict The top-level dictionary to search
        //  \param dictName The sub-dictionary name to select.
        //     Return nullptr if it does not exist.
        //  \param obrPtr The objectRegistry for lookup of indirect entries.
        //
        //  Since the specifications are isolated within a sub-dictionary,
        //  the 'origin' entry is treated as optional
        static autoPtr<coordinateSystem> NewIfPresent
        (
            const dictionary& dict,
            const word& dictName,
            const objectRegistry* obrPtr
        );


protected:

    //- Friendship with indirect for dispatching to its underlying system
    friend coordSystem::indirect;


    // Protected Data

        //- User specification of the coordinate rotation
        //  May be invalid after a move assignment or transfer
        autoPtr<coordinateRotation> spec_;

        //- The coordinate system origin
        point origin_;

        //- The rotation tensor
        tensor rot_;

        //- The name of the coordinate system (optional)
        word name_;

        //- An optional note describing the coordinate system
        string note_;

        //- Dummy coordinate system for suppressed manipulation
        static coordinateSystem dummy_;


    // Protected Member Functions

        //- Implementation for R() methods
        template<class PointField>
        tmp<tensorField> rotationsImpl(const PointField& global) const;

        //- Implementation for transformPoint() methods
        template<class PointField>
        tmp<pointField> transformPointImpl(const PointField& localCart) const;

        //- Implementation for transformPosition() methods
        template<class PointField>
        tmp<pointField> invTransformPointImpl(const PointField& global) const;

        //- Apply single transform tensor for multiple inputs
        template<class RetType, class Type, class BinaryOp>
        static tmp<Field<RetType>> manyTimesImpl
        (
            const tensor& tt,
            const UList<Type>& input,
            const BinaryOp& bop
        );

        //- Use position-dependent transform tensors for multiple inputs
        template<class RetType, class PointField, class Type, class BinaryOp>
        tmp<Field<RetType>> oneToOneImpl
        (
            const PointField& global,
            const UList<Type>& input,
            const BinaryOp& bop
        ) const;

        //- Use position-dependent transform tensors for single input
        template<class RetType, class PointField, class Type, class BinaryOp>
        tmp<Field<RetType>> oneToManyImpl
        (
            const PointField& global,
            const Type& input,
            const BinaryOp& bop
        ) const;


        //- From local coordinate system to the global Cartesian system
        //- with optional translation for the origin
        virtual vector localToGlobal
        (
            const vector& local,
            bool translate
        ) const;

        //- From local coordinate system to the global Cartesian system
        //- with optional translation for the origin
        virtual tmp<vectorField> localToGlobal
        (
            const vectorField& local,
            bool translate
        ) const;

        //- From global Cartesian system to the local coordinate system
        //- with optional translation for the origin
        virtual vector globalToLocal
        (
            const vector& global,
            bool translate
        ) const;

        //- From global Cartesian system to the local coordinate system
        //- with optional translation for the origin
        virtual tmp<vectorField> globalToLocal
        (
            const vectorField& global,
            bool translate
        ) const;

        //- Assign from dictionary content with specified read handling
        //- of the 'origin' entry
        void assign
        (
            const dictionary& dict,
            IOobjectOption::readOption readOrigin = IOobjectOption::MUST_READ
        );


    // Constructors

        //- Construct null, without allocating a coordinateRotation
        //- specification.
        coordinateSystem(std::nullptr_t);


public:

    //- Runtime type information
    TypeName("coordinateSystem");

        //- Helper for construction of coordinateSystem PtrList
        //  The Istream contains a word followed by a dictionary.
        struct iNew
        {
            autoPtr<coordinateSystem> operator()(Istream& is) const
            {
                return coordinateSystem::New(is);
            }
        };


    // Constructors

        //- Default construct. This is an identity coordinate system
        coordinateSystem();

        //- Copy construct from rotation with origin=0
        explicit coordinateSystem(const coordinateRotation& crot);

        //- Move construct from rotation with origin=0
        explicit coordinateSystem(coordinateRotation&& crot);

        //- Copy construct
        coordinateSystem(const coordinateSystem& csys);

        //- Move construct
        coordinateSystem(coordinateSystem&& csys);

        //- Move construct from autoPtr
        explicit coordinateSystem(autoPtr<coordinateSystem>&& csys);

        //- Copy construct with a different name
        coordinateSystem
        (
            const word& name,
            const coordinateSystem& csys
        );

        //- Construct from origin and rotation
        coordinateSystem
        (
            const point& origin,
            const coordinateRotation& crot
        );

        //- Construct from origin and 2 axes
        coordinateSystem
        (
            const point& origin,
            const vector& axis,
            const vector& dirn
        );

        //- Construct from origin and rotation
        coordinateSystem
        (
            const word& name,
            const point& origin,
            const coordinateRotation& crot
        );

        //- Construct named from origin and 2 axes
        coordinateSystem
        (
            const word& name,
            const point& origin,
            const vector& axis,
            const vector& dirn
        );

        //- Construct from dictionary with optional
        //- read handling for the 'origin' entry (default: MUST_READ).
        //
        //  \note The readOrigin is downgraded to READ_IF_PRESENT
        //  if the dictionary itself is "coordinateSystem"
        explicit coordinateSystem
        (
            const dictionary& dict,
            IOobjectOption::readOption readOrigin = IOobjectOption::MUST_READ
        );

        //- Construct from dictionary with optional subDict lookup and optional
        //- read handling for the 'origin' entry (default: MUST_READ).
        //
        //  \param dictName If non-empty, mandatory sub-dictionary to use.
        //
        //  \note The readOrigin is downgraded to READ_IF_PRESENT
        //  if the dictionary itself is "coordinateSystem"
        //  or if a sub-dictionary is being used
        coordinateSystem
        (
            const dictionary& dict,
            const word& dictName,
            IOobjectOption::readOption readOrigin = IOobjectOption::MUST_READ
        );


    //- Return clone
    virtual autoPtr<coordinateSystem> clone() const
    {
        return autoPtr<coordinateSystem>::New(*this);
    }


    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        coordinateSystem,
        dictionary,
        (
            const dictionary& dict,
            IOobjectOption::readOption readOrigin
        ),
        (dict, readOrigin)
    );

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        coordinateSystem,
        registry,
        (
            const objectRegistry& obr,
            const dictionary& dict,
            IOobjectOption::readOption readOrigin
        ),
        (obr, dict, readOrigin)
    );


    // Factory Methods

        //- Clone a coordinate system
        template<class Derived>
        static autoPtr<coordinateSystem> Clone(const Derived& csys)
        {
            return autoPtr<coordinateSystem>(new Derived(csys));
        }

        //- Select construct the specified coordinate system type
        //- with reference to objectRegistry for indirect entries.
        //
        //  An empty modelType will be treated as "cartesian"
        static autoPtr<coordinateSystem> Clone
        (
            const word& modelType,
            const objectRegistry& obr,
            const dictionary& dict,
            IOobjectOption::readOption readOrigin = IOobjectOption::MUST_READ
        );

        //- Select construct the specified coordinate system type
        //- with reference to objectRegistry for indirect entries.
        //
        //  An empty modelType will be treated as "cartesian"
        static autoPtr<coordinateSystem> New
        (
            const word& modelType,
            const objectRegistry& obr,
            const dictionary& dict,
            IOobjectOption::readOption readOrigin = IOobjectOption::MUST_READ
        );

        //- Select construct the specified coordinate system type
        //
        //  An empty modelType will be treated as "cartesian"
        static autoPtr<coordinateSystem> New
        (
            const word& modelType,
            const dictionary& dict,
            IOobjectOption::readOption readOrigin = IOobjectOption::MUST_READ
        );

        //- Select construct from dictionary with reference to objectRegistry
        //- for indirect entries.
        //
        //  \param dictName If non-empty, the sub-dictionary name to use
        //      for the coordinate system description.
        //
        //  \note When the dictName is empty, it includes an implicit search
        //      for a "coordinateSystem" sub-dictionary for convenience and
        //      compatibility with previous versions (1806 and earlier).
        static autoPtr<coordinateSystem> New
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& dictName = word::null,
            IOobjectOption::readOption readOrigin = IOobjectOption::MUST_READ
        );

        //- Select constructed from dictionary
        //  \param dictName If non-empty, the sub-dictionary name to use
        //      for the coordinate system description.
        //
        //  \note When the dictName is empty, it includes an implicit search
        //      for a "coordinateSystem" sub-dictionary for convenience and
        //      compatibility with previous versions (1806 and earlier).
        static autoPtr<coordinateSystem> New
        (
            const dictionary& dict,
            const word& dictName = word::null,
            IOobjectOption::readOption readOrigin = IOobjectOption::MUST_READ
        );

        //- Select constructed from Istream
        //  Expects a name/dictionary as input
        static autoPtr<coordinateSystem> New
        (
            Istream& is,
            IOobjectOption::readOption readOrigin = IOobjectOption::MUST_READ
        );


        //- Optional select construct from dictionary, with registry reference
        //
        //  \param obr The objectRegistry (for lookup of indirect entries)
        //  \param dict The top-level dictionary to search
        //  \param dictName The sub-dictionary name to select the
        //     for coordinate system specification
        //     (default is 'coordinateSystem').
        //
        //  Since the specifications are isolated within a sub-dictionary,
        //  the 'origin' entry is treated as optional
        static autoPtr<coordinateSystem> NewIfPresent
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& dictName = coordinateSystem::typeName
        );

        //- Optional select construct from dictionary
        //
        //  \param dict The top-level dictionary to search
        //  \param dictName The sub-dictionary name to select the
        //     for coordinate system specification
        //     (default is 'coordinateSystem').
        //
        //  Since the specifications are isolated within a sub-dictionary,
        //  the 'origin' entry is treated as optional
        static autoPtr<coordinateSystem> NewIfPresent
        (
            const dictionary& dict,
            const word& dictName = coordinateSystem::typeName
        );


    //- Destructor
    virtual ~coordinateSystem() = default;


    // Member Functions

    // Characteristics

        //- Consider good if it has a specification
        virtual bool good() const { return bool(spec_); }

        //- True if the rotation tensor is uniform for all locations
        virtual bool uniform() const { return true; }

        //- Same as good() - 2023-07
        virtual bool valid() const { return this->good(); }


    // Access

        //- Return origin
        virtual const point& origin() const
        {
            return origin_;
        }

        //- The rotation specification
        virtual const coordinateRotation& rotation() const
        {
            return *spec_;
        }

        //- Return the name
        virtual const word& name() const
        {
            return name_;
        }

        //- Return the optional note
        virtual const string& note() const
        {
            return note_;
        }

        //- Return const reference to the rotation tensor
        virtual const tensor& R() const
        {
            return rot_;
        }

        //- The local Cartesian x-axis in global coordinates
        virtual const vector e1() const
        {
            return rot_.cx();
        }

        //- The local Cartesian y-axis in global coordinates
        virtual const vector e2() const
        {
            return rot_.cy();
        }

        //- The local Cartesian z-axis in global coordinates
        virtual const vector e3() const
        {
            return rot_.cz();
        }


    // Edit

        //- Rename
        virtual void rename(const word& newName)
        {
            name_ = newName;
        }

        //- Edit access to optional note
        virtual string& note()
        {
            return note_;
        }

        //- Edit access to origin
        virtual point& origin()
        {
            return origin_;
        }

        //- Reset origin and rotation to an identity coordinateSystem
        //  Also resets the note
        virtual void clear();

        //- Change the rotation
        virtual void rotation(autoPtr<coordinateRotation>&& crot);


    // Write

        //- Write
        virtual void write(Ostream& os) const;

        //- Write 'coordinateSystem' dictionary entry
        virtual void writeEntry(Ostream& os) const;

        //- Write dictionary entry
        virtual void writeEntry(const word& keyword, Ostream& os) const;


    // Member Operators

        //- Copy assignment
        void operator=(const coordinateSystem& csys);

        //- Move assignment
        void operator=(coordinateSystem&& csys);

        //- Copy assignment from autoPtr
        void operator=(const autoPtr<coordinateSystem>& csys);

        //- Move assignment from autoPtr
        void operator=(autoPtr<coordinateSystem>&& csys);


    // Rotation

        //- Position-dependent rotation tensor (when uniform = false)
        //- \return tensor
        virtual tensor R(const point& global) const;

        //- Position-dependent rotation tensors (when uniform = false)
        //- \return tensorField
        virtual tmp<tensorField> R(const UList<point>& global) const;

        //- Position-dependent rotation tensors (when uniform = false)
        //- \return tensorField
        virtual tmp<tensorField> R(const pointUIndList& global) const;


    // Position

        //- Transform point and add origin offset.
        //  Corresponds to a local-to-global transformation using Cartesian
        //  coordinates for both local and global.
        point transformPoint(const point& localCart) const;

        //- Transform points and add origin offset.
        tmp<pointField> transformPoint(const UList<point>& localCart) const;

        //- Transform points and add origin offset.
        tmp<pointField> transformPoint(const pointUIndList& localCart) const;


        //- Remove origin offset and inverse transform point.
        //  Corresponds to a global-to-local transformation using Cartesian
        //  coordinates for both local and global.
        point invTransformPoint(const point& global) const;

        //- Remove origin offset and inverse transform points.
        tmp<pointField> invTransformPoint(const UList<point>& global) const;

        //- Remove origin offset and inverse transform points.
        tmp<pointField> invTransformPoint(const pointUIndList& global) const;


    // Transformations with change of coordinate types

        //- From local coordinate position to global (cartesian) position
        point globalPosition(const point& local) const
        {
            return localToGlobal(local, true);
        }

        //- From local coordinate position to global (cartesian) position
        tmp<pointField> globalPosition(const pointField& local) const
        {
            return localToGlobal(local, true);
        }

        //- From global (cartesian) position to local coordinate position
        point localPosition(const point& global) const
        {
            return globalToLocal(global, true);
        }

        //- From global (cartesian) position to local coordinate position
        tmp<pointField> localPosition(const pointField& global) const
        {
            return globalToLocal(global, true);
        }



        //- From local to global (cartesian) vector components
        vector globalVector(const vector& local) const
        {
            return localToGlobal(local, false);
        }

        //- From local to global (cartesian) vector components
        tmp<vectorField> globalVector(const vectorField& local) const
        {
            return localToGlobal(local, false);
        }

        //- From global (cartesian) to local vector components
        vector localVector(const vector& global) const
        {
            return globalToLocal(global, false);
        }

        //- From global (cartesian) to local vector components
        tmp<vectorField> localVector(const vectorField& global) const
        {
            return globalToLocal(global, false);
        }


    // Transformations (input and output are Cartesian)

#undef  defineCoordinateSystemTransform
#define defineCoordinateSystemTransform(Op, RetType, Type)                    \
                                                                              \
    /*! \brief With constant rotation tensor */                               \
    virtual RetType Op(const Type& input) const;                              \
                                                                              \
    /*! \brief With constant rotation tensor */                               \
    virtual tmp<Field<RetType>> Op(const UList<Type>& input) const;           \
                                                                              \
    /*! \brief With rotation tensor at given global position */               \
    virtual RetType Op(const point& global, const Type& input) const;         \
                                                                              \
    /*! \brief With rotation tensors at given global positions */             \
    virtual tmp<Field<RetType>> Op                                            \
    (                                                                         \
        const UList<point>& global,                                           \
        const Type& input                                                     \
    ) const;                                                                  \
                                                                              \
    /*! \brief With rotation tensors at given global positions */             \
    virtual tmp<Field<RetType>> Op                                            \
    (                                                                         \
        const pointUIndList& global,                                          \
        const Type& input                                                     \
    ) const;                                                                  \
                                                                              \
    /*! \brief With rotation tensors at given global positions */             \
    virtual tmp<Field<RetType>> Op                                            \
    (                                                                         \
        const UList<point>& global,                                           \
        const UList<Type>& input                                              \
    ) const;                                                                  \
                                                                              \
    /*! \brief With rotation tensors at given global positions */             \
    virtual tmp<Field<RetType>> Op                                            \
    (                                                                         \
        const pointUIndList& global,                                          \
        const UList<Type>& input                                              \
    ) const;


    defineCoordinateSystemTransform(transformPrincipal, symmTensor, vector);

    defineCoordinateSystemTransform(transform, scalar, scalar);
    defineCoordinateSystemTransform(transform, vector, vector);
    defineCoordinateSystemTransform
    (
        transform,
        sphericalTensor,
        sphericalTensor
    );
    defineCoordinateSystemTransform(transform, symmTensor, symmTensor);
    defineCoordinateSystemTransform(transform, tensor, tensor);

    defineCoordinateSystemTransform(invTransform, scalar, scalar);
    defineCoordinateSystemTransform(invTransform, vector, vector);
    defineCoordinateSystemTransform
    (
        invTransform,
        sphericalTensor,
        sphericalTensor
    );
    defineCoordinateSystemTransform(invTransform, symmTensor, symmTensor);
    defineCoordinateSystemTransform(invTransform, tensor, tensor);

    #undef defineCoordinateSystemTransform
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Global Operators

//- Compare inequality
bool operator!=(const coordinateSystem& a, const coordinateSystem& b);

//- Output operator
Ostream& operator<<(Ostream& os, const coordinateSystem& csys);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "coordinateSystemTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
