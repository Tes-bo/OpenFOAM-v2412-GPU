/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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
    Foam::fv::snGradScheme

Group
    grpFvSnGradSchemes

Description
    Abstract base class for runtime selected \c snGrad surface
    normal gradient schemes.

    A surface normal gradient is evaluated at a cell face. It
    is the normal-to-face component of the gradient of
    values at the centres of two cells that the face connects.

    Unit-surface-normal vector decomposition is based on the
    so-called over-relaxed approach. Orthogonal components are
    treated implicitly and non-orthogonal components are treated
    explicitly with (or without) various limiters.

SourceFiles
    snGradScheme.C

\*---------------------------------------------------------------------------*/

#ifndef snGradScheme_H
#define snGradScheme_H

#include "tmp.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                        Class snGradScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class snGradScheme
:
    public refCount
{
    // Private Data

        //- Hold const reference to mesh
        const fvMesh& mesh_;


    // Private Member Functions

        //- No copy construct
        snGradScheme(const snGradScheme&) = delete;

        //- No copy assignment
        void operator=(const snGradScheme&) = delete;


public:

    //- Runtime type information
    virtual const word& type() const = 0;


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            tmp,
            snGradScheme,
            Mesh,
            (const fvMesh& mesh, Istream& schemeData),
            (mesh, schemeData)
        );


    // Constructors

        //- Construct from mesh
        snGradScheme(const fvMesh& mesh)
        :
            mesh_(mesh)
        {}


    // Selectors

        //- Return new tmp interpolation scheme
        static tmp<snGradScheme<Type>> New
        (
            const fvMesh& mesh,
            Istream& schemeData
        );


    //- Destructor
    virtual ~snGradScheme() = default;


    // Member Functions

        //- Return const reference to mesh
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        //- Return the snGrad of the given cell field
        //- by using the given deltaCoeffs
        static tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        snGrad
        (
            const GeometricField<Type, fvPatchField, volMesh>&,
            const tmp<surfaceScalarField>&,
            const word& snGradName = "snGrad"
        );

        //- Return the sndGrad of the given cell field
        static tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        sndGrad
        (
            const GeometricField<Type, fvPatchField, volMesh>&,
            const word& snGradName = "sndGrad"
        );

        //- Return the interpolation weighting factors for the given field
        virtual tmp<surfaceScalarField> deltaCoeffs
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const = 0;

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const noexcept
        {
            return false;
        }

        //- Return the explicit correction to the snGrad for the given field
        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        correction(const GeometricField<Type, fvPatchField, volMesh>&) const
        {
            return nullptr;
        }

        //- Return the snGrad of the given cell field
        //- with explicit correction
        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        snGrad(const GeometricField<Type, fvPatchField, volMesh>&) const;

        //- Return the snGrad of the given tmp cell field
        //- with explicit correction
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        snGrad
        (
            const tmp<GeometricField<Type, fvPatchField, volMesh>>&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables

#define makeSnGradTypeScheme(SS, Type)                                         \
    defineNamedTemplateTypeNameAndDebug(Foam::fv::SS<Foam::Type>, 0);          \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        namespace fv                                                           \
        {                                                                      \
            snGradScheme<Type>::addMeshConstructorToTable<SS<Type>>            \
                add##SS##Type##MeshConstructorToTable_;                        \
        }                                                                      \
    }

#define makeSnGradScheme(SS)                                                   \
                                                                               \
makeSnGradTypeScheme(SS, scalar)                                               \
makeSnGradTypeScheme(SS, vector)                                               \
makeSnGradTypeScheme(SS, sphericalTensor)                                      \
makeSnGradTypeScheme(SS, symmTensor)                                           \
makeSnGradTypeScheme(SS, tensor)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "snGradScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
