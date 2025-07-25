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
    Foam::fa::divScheme

Description
    Abstract base class for finite area calculus div schemes.

SourceFiles
    faDivScheme.C
    faDivSchemes.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_faDivScheme_H
#define Foam_faDivScheme_H

#include "areaFieldsFwd.H"
#include "edgeFieldsFwd.H"
#include "linearEdgeInterpolation.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class Type> class faMatrix;
class faMesh;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

/*---------------------------------------------------------------------------*\
                          Class divScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class divScheme
:
    public refCount
{
protected:

    // Protected Data

        //- Reference to mesh
        const faMesh& mesh_;

        //- Edge interpolation scheme
        tmp<edgeInterpolationScheme<Type>> tinterpScheme_;


public:

    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            tmp,
            divScheme,
            Istream,
            (const faMesh& mesh, Istream& schemeData),
            (mesh, schemeData)
        );


    // Generated Methods

        //- No copy construct
        divScheme(const divScheme&) = delete;

        //- No copy assignment
        void operator=(const divScheme&) = delete;


    // Constructors

        //- Construct from mesh
        divScheme(const faMesh& mesh)
        :
            mesh_(mesh),
            tinterpScheme_(new linearEdgeInterpolation<Type>(mesh))
        {}

        //- Construct from mesh and Istream
        divScheme(const faMesh& mesh, Istream& is)
        :
            mesh_(mesh)
        {
            if (is.eof())
            {
                tinterpScheme_.reset
                (
                    new linearEdgeInterpolation<Type>(mesh)
                );
            }
            else
            {
                tinterpScheme_.reset
                (
                    edgeInterpolationScheme<Type>::New(mesh, is)
                );
            }
        }


    // Selectors

        //- Return a pointer to a new divScheme created on freestore
        static tmp<divScheme<Type>> New
        (
            const faMesh& mesh,
            Istream& schemeData
        );


    //- Destructor
    virtual ~divScheme() = default;


    // Member Functions

        //- Return mesh reference
        const faMesh& mesh() const noexcept { return mesh_; }

        virtual tmp
        <
            GeometricField
            <typename innerProduct<vector, Type>::type, faPatchField, areaMesh>
        > facDiv
        (
            const GeometricField<Type, faPatchField, areaMesh>&
        ) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables

#define makeFaDivTypeScheme(SS, Type)                                          \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(Foam::fa::SS<Foam::Type>, 0);          \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        namespace fa                                                           \
        {                                                                      \
            divScheme<Type>::addIstreamConstructorToTable<SS<Type>>            \
                add##SS##Type##IstreamConstructorToTable_;                     \
        }                                                                      \
    }

#define makeFaDivScheme(SS)                                                    \
                                                                               \
makeFaDivTypeScheme(SS, vector)                                                \
makeFaDivTypeScheme(SS, tensor)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "faDivScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
