/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2019 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
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
    Foam::distributionModels::general

Description
    Particle-size distribution model wherein random samples are
    drawn from a given arbitrary probability density function
    or cumulative distribution function. Input distributions are
    specified as pairs of (size, probability) (i.e. (point, value) ).

Usage
    Minimal example by using \c constant/\<CloudProperties\>:
    \verbatim
    subModels
    {
        injectionModels
        {
            <name>
            {
                ...

                sizeDistribution
                {
                    type        general;
                    generalDistribution
                    {
                        cumulative      false;

                        distribution
                        (
                            (<size1>    <probability1>)
                            (<size2>    <probability2>)
                            ...
                            (<sizeN>    <probabilityN>)
                        );
                    }
                }
            }
        }
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                       | Type | Reqd | Deflt
      type         | Type name: general                | word | yes  | -
      generalDistribution | Distribution settings      | dict | yes  | -
      distribution | \<size\>-\<probability\> pairs    | dict | yes  | -
      \<size\>     | Particle size                     | scalar | yes | -
      \<probability\> | Volume fraction/probability    | scalar | yes | -
      cumulative   | Flag to determine if input distribution is specified <!--
                   --> as cumulative or as density     | bool | no   | false
    \endtable

Note
  - An example for a pair within \c distribution subdictionary
  can be given as follows: "(1e-07 1.3)" means 1.3\% of particles
  are modelled with a particle diameter of "1e-7" [m], and so on.
  - Variation between input pairs is assumed to be linear.
  - Elements in the second column (i.e. probability) are normalised.
  - Elements in the second column for cumulative distribution
  functions must start from zero and must be non-decreasing (i.e. monotonic).
  - Elements in the first column (i.e. size) must be specified
  in an ascending order.
  - Input pairs cannot contain any negative element.

SourceFiles
    general.C

\*---------------------------------------------------------------------------*/

#ifndef distributionModels_general_H
#define distributionModels_general_H

#include "distributionModel.H"
#include "Vector.H"
#include "VectorSpace.H"
#include "Field.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace distributionModels
{
    // Forward Declarations
    class general;
}

// Forward Declarations
Istream& operator>>(Istream&, distributionModels::general&);
Ostream& operator<<(Ostream&, const distributionModels::general&);

namespace distributionModels
{

/*---------------------------------------------------------------------------*\
                          Class general Declaration
\*---------------------------------------------------------------------------*/

class general
:
    public distributionModel
{
    typedef VectorSpace<Vector<scalar>, scalar, 2> pair;

    // Private Data

        //- List of (x, y=f(x)) pairs
        List<pair> xy_;

        //- Number of entries in the xy_ list
        label nEntries_;

        //- Mean of the distribution
        scalar meanValue_;

        //- Values of cumulative distribution function
        List<scalar> integral_;

        //- Flag to determine if input is specified as cumulative or as density
        bool cumulative_;


    // Private Member Functions

        //- Initialise the distribution parameters
        void initialise();


public:

    //- Runtime type information
    TypeName("general");


    // Constructors

        //- Construct from components
        general(const dictionary& dict, Random& rndGen);

        //- Construct from components
        //  Allows negative entries
        general
        (
            const UList<scalar>& sampleData,
            const scalar binWidth,
            Random& rndGen
        );

        //- Copy construct
        general(const general& p);

        //- Construct and return a clone
        virtual autoPtr<distributionModel> clone() const
        {
            return autoPtr<distributionModel>(new general(*this));
        }

        //- No copy assignment
        void operator=(const general&) = delete;


    //- Destructor
    virtual ~general() = default;


    // Member Functions

        //- Bin boundaries
        virtual tmp<scalarField> x() const;

        //- Probabilities
        virtual tmp<scalarField> y() const;

        //- Sample the distribution
        virtual scalar sample() const;

        //- Return the arithmetic mean of the distribution data
        virtual scalar meanValue() const;

        //- Write data to stream
        virtual void writeData(Ostream& os) const;

        //- Read data from stream
        virtual void readData(Istream& os);

        //- Write data in dictionary format
        virtual dictionary writeDict(const word& dictName) const;

        //- Read data from dictionary
        virtual void readDict(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace distributionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
