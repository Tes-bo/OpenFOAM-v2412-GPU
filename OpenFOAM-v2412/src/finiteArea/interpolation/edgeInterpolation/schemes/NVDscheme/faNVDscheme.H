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
    Foam::faNVDscheme

Description
    Class to create the weighting-factors based on the NVD
    (Normalised Variable Diagram).

    The particular differencing scheme class is supplied as a template argument,
    the weight function of which is called by the weight function of this class
    for the internal edges as well as edges of coupled patches
    (e.g. processor-processor patches). The weight function is supplied the
    central-differencing weighting factor, the edge-flux, the cell and edge
    gradients (from which the normalised variable distribution may be created)
    and the cell centre distance.

    This code organisation is both neat and efficient, allowing for convenient
    implementation of new schemes to run on parallelised cases.

SourceFiles
    faNVDscheme.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_faNVDscheme_H
#define Foam_faNVDscheme_H

#include "edgeInterpolationScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class faNVDscheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class NVDweight>
class faNVDscheme
:
    public edgeInterpolationScheme<Type>,
    public NVDweight
{
protected:

    // Protected Data

        const edgeScalarField& edgeFlux_;


public:

    //- Define a typedef for the NVDweight
    typedef NVDweight Weight;


    //- Runtime type information
    TypeName("faNVDscheme");


    // Generated Methods

        //- No copy construct
        faNVDscheme(const faNVDscheme&) = delete;

        //- No copy assignment
        void operator=(const faNVDscheme&) = delete;


    // Constructors

        //- Construct from mesh and edgeFlux and blendingFactor
        faNVDscheme
        (
            const faMesh& mesh,
            const edgeScalarField& edgeFlux,
            const NVDweight& weight
        )
        :
            edgeInterpolationScheme<Type>(mesh),
            NVDweight(weight),
            edgeFlux_(edgeFlux)
        {}

        //- Construct from mesh and Istream.
        //  The name of the flux field is read from the Istream and looked-up
        //  from the database
        faNVDscheme
        (
            const faMesh& mesh,
            Istream& is
        )
        :
            edgeInterpolationScheme<Type>(mesh),
            NVDweight(is),
            edgeFlux_
            (
                mesh.thisDb().lookupObject<edgeScalarField>
                (
                    word(is)
                )
            )
        {}

        //- Construct from mesh, edgeFlux and Istream
        faNVDscheme
        (
            const faMesh& mesh,
            const edgeScalarField& edgeFlux,
            Istream& is
        )
        :
            edgeInterpolationScheme<Type>(mesh),
            NVDweight(is),
            edgeFlux_(edgeFlux)
        {}


    // Member Functions

        //- Return the interpolation weighting factors
        virtual tmp<edgeScalarField> weights
        (
            const GeometricField<Type, faPatchField, areaMesh>&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables

#define makeNVDedgeInterpolationTypeScheme(SS, WEIGHT, NAME, TYPE)             \
                                                                               \
typedef faNVDscheme<TYPE, WEIGHT> faNVDscheme##TYPE##WEIGHT_;                  \
defineTemplateTypeNameAndDebugWithName(faNVDscheme##TYPE##WEIGHT_, NAME, 0);   \
                                                                               \
edgeInterpolationScheme<TYPE>::addMeshConstructorToTable                       \
<faNVDscheme<TYPE, WEIGHT>>                                                   \
    add##SS##TYPE##MeshConstructorToTable_;                                    \
                                                                               \
edgeInterpolationScheme<TYPE>::addMeshFluxConstructorToTable                   \
<faNVDscheme<TYPE, WEIGHT>>                                                   \
    add##SS##TYPE##MeshFluxConstructorToTable_;


#define makeNVDedgeInterpolationScheme(SS, WEIGHT, NAME)                       \
                                                                               \
makeNVDedgeInterpolationTypeScheme(SS, WEIGHT, NAME, scalar)                   \
makeNVDedgeInterpolationTypeScheme(SS, WEIGHT, NAME, vector)                   \
makeNVDedgeInterpolationTypeScheme(SS, WEIGHT, NAME, tensor)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "faNVDscheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
