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
    Foam::fa::laplacianScheme

Description
    Abstract base class for finite area calculus laplacian schemes.

SourceFiles
    faLaplacianScheme.C
    faLaplacianSchemes.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_faLaplacianScheme_H
#define Foam_faLaplacianScheme_H

#include "areaFieldsFwd.H"
#include "edgeFieldsFwd.H"
#include "correctedLnGrad.H"
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
                       Class laplacianScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class laplacianScheme
:
    public refCount
{
protected:

    // Protected Data

        const faMesh& mesh_;

        tmp<edgeInterpolationScheme<scalar>> tinterpGammaScheme_;
        tmp<lnGradScheme<Type>> tlnGradScheme_;


public:

    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            tmp,
            laplacianScheme,
            Istream,
            (const faMesh& mesh, Istream& schemeData),
            (mesh, schemeData)
        );


    // Generated Methods

        //- No copy construct
        laplacianScheme(const laplacianScheme&) = delete;

        //- No copy assignment
        void operator=(const laplacianScheme&) = delete;


    // Constructors

        //- Construct from mesh
        laplacianScheme(const faMesh& mesh)
        :
            mesh_(mesh),
            tinterpGammaScheme_(new linearEdgeInterpolation<scalar>(mesh)),
            tlnGradScheme_(new correctedLnGrad<Type>(mesh))
        {}

        //- Construct from mesh and Istream
        laplacianScheme(const faMesh& mesh, Istream& is)
        :
            mesh_(mesh)
        {
            if (is.eof())
            {
                tinterpGammaScheme_.reset
                (
                    new linearEdgeInterpolation<scalar>(mesh)
                );

                tlnGradScheme_.reset
                (
                    new correctedLnGrad<Type>(mesh)
                );
            }
            else
            {
                tinterpGammaScheme_.reset
                (
                    edgeInterpolationScheme<scalar>::New(mesh, is)
                );

                tlnGradScheme_.reset
                (
                    lnGradScheme<Type>::New(mesh, is)
                );
            }
        }


    // Selectors

        //- Return a pointer to a new laplacianScheme created on freestore
        static tmp<laplacianScheme<Type>> New
        (
            const faMesh& mesh,
            Istream& schemeData
        );


    //- Destructor
    virtual ~laplacianScheme();


    // Member Functions

        //- Return mesh reference
        const faMesh& mesh() const
        {
            return mesh_;
        }

        virtual tmp<faMatrix<Type>> famLaplacian
        (
            const edgeScalarField&,
            const GeometricField<Type, faPatchField, areaMesh>&
        ) = 0;

        virtual tmp<faMatrix<Type>> famLaplacian
        (
            const areaScalarField&,
            const GeometricField<Type, faPatchField, areaMesh>&
        );

        virtual tmp<GeometricField<Type, faPatchField, areaMesh>> facLaplacian
        (
            const GeometricField<Type, faPatchField, areaMesh>&
        ) = 0;

        virtual tmp<GeometricField<Type, faPatchField, areaMesh>> facLaplacian
        (
            const edgeScalarField&,
            const GeometricField<Type, faPatchField, areaMesh>&
        ) = 0;

        virtual tmp<GeometricField<Type, faPatchField, areaMesh>> facLaplacian
        (
            const areaScalarField&,
            const GeometricField<Type, faPatchField, areaMesh>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables

#define makeFaLaplacianTypeScheme(SS, Type)                                    \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(Foam::fa::SS<Foam::Type>, 0);          \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
    namespace fa                                                               \
    {                                                                          \
        laplacianScheme<Type>::addIstreamConstructorToTable<SS<Type>>          \
            add##SS##Type##IstreamConstructorToTable_;                         \
    }                                                                          \
    }


#define makeFaLaplacianScheme(SS)                                              \
                                                                               \
makeFaLaplacianTypeScheme(SS, scalar)                                          \
makeFaLaplacianTypeScheme(SS, vector)                                          \
makeFaLaplacianTypeScheme(SS, tensor)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "faLaplacianScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
