/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

Class
    Foam::PatchFunction1Types::Sampled

Description
    PatchFunction1 to sample an existing field.

    It is the exact functionality from mappedField boundary condition
    with the following differences:
    - the field name is specified. So any derived fields will have the
    same field name to sample.
    - if used with uniformFixedValue boundary condition there is the problem
    that that re-evaluates instead of reading/mapping.

Usage
    \table
        Property            | Description             | Required | Default value
        field               | Field name              | yes      |
        sampleMode          | how to obtain samples   | yes      |
        sampleRegion        | mesh to sample          | no       | ""
        samplePatch         | patch to sample         | no       | ""
        offsetMode          | how to offset samples   | no       | uniform
        offset              | uniform offset vector   | no       | (0 0 0)
        interpolationScheme | interpolation method    | yes      | cell
        setAverage          | optional average adjustment | no   | false
        average             | optional average value  | no       |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Field to sample
        field           U;

        // Geometric/mapping info (equivalent of 'mappedPatch' patch type)
        sampleMode      nearestCell;
        offset          (0 -0.001 0);

        // Field specific info (equivalent of 'mappedField' patch field type)
        interpolationScheme cell;
    }
    \endverbatim

SourceFiles
    Sampled.C

\*---------------------------------------------------------------------------*/

#ifndef PatchFunction1Types_Sampled_H
#define PatchFunction1Types_Sampled_H

#include "PatchFunction1.H"
#include "mappedPatchBase.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace PatchFunction1Types
{

/*---------------------------------------------------------------------------*\
                           Class Sampled Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Sampled
:
    public PatchFunction1<Type>,
    public mappedPatchBase
{
    // Private Member Functions

        //- Selective retrieval of "average" entry from the dictionary
        static Type getAverage(const dictionary& dict, bool mandatory);


protected:

    // Protected data

        //- Name of the field
        word fieldName_;

        //- If true adjust the mapped field to maintain average value average_
        const bool setAverage_;

        //- Average value the mapped field is adjusted to maintain if
        //- setAverage_ is set true
        const Type average_;

        //- Interpolation scheme to use for nearestcell mode
        word interpolationScheme_;


    // Private Member Functions

        //- Field to sample. Either on my or nbr mesh
        bool haveSampleField() const;


public:

    // Runtime type information
    TypeName("sampled");


    // Generated Methods

        //- No copy assignment
        void operator=(const Sampled<Type>&) = delete;


    // Constructors

        //- Construct from entry name and dictionary
        Sampled
        (
            const polyPatch& pp,
            const word& redirectType,
            const word& entryName,
            const dictionary& dict,
            const bool faceValues = true
        );

        //- Copy construct setting patch
        explicit Sampled
        (
            const Sampled<Type>& rhs,
            const polyPatch& pp
        );

        //- Copy construct
        explicit Sampled(const Sampled<Type>& rhs);

        //- Return a clone
        virtual tmp<PatchFunction1<Type>> clone() const
        {
            return PatchFunction1<Type>::Clone(*this);
        }

        //- Return a clone, setting the patch
        virtual tmp<PatchFunction1<Type>> clone(const polyPatch& pp) const
        {
            return PatchFunction1<Type>::Clone(*this, pp);
        }


    //- Destructor
    virtual ~Sampled() = default;


    // Member Functions

        //- Field to sample. Either on my or nbr mesh
        const GeometricField<Type, fvPatchField, volMesh>& sampleField() const;

        // Evaluation

            //- Return sampled value
            virtual tmp<Field<Type>> value(const scalar x) const;

            //- Is value constant (i.e. independent of x)
            virtual inline bool constant() const
            {
                return false;
            }

            //- Is value uniform (i.e. independent of coordinate)
            virtual inline bool uniform() const
            {
                return false;
            }

            //- Integrate between two values
            virtual tmp<Field<Type>> integrate
            (
                const scalar x1,
                const scalar x2
            ) const;


        // I-O

            //- Write in dictionary format
            virtual void writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace PatchFunction1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "Sampled.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
