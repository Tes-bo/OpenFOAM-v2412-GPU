/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2021 OpenCFD Ltd.
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
    Foam::distributionModels::binned

Description
    Particle-size distribution model wherein random samples are
    drawn from a given discrete set of (\c bin, \c probability) pairs,
    where in terms of its meaning, \c bins correspond to particle sizes
    and \c probabilities correspond to (relative) probability of occurrences.

    The second column (i.e. \c probability) are normalised by the sum of all
    its values, resulting in a normalised column in the domain of [0,1].
    To generate a sample, first a sample drawn from the uniform probability
    density function on the unit interval (i.e. \c u), and then, the \c bin
    corresponding to the first \c probability larger than \c u is fetched
    as the particle size to be further processed.

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
                    type        binned;
                    binnedDistribution
                    {
                        distribution
                        (
                            (<bin1>    <probability1>)
                            (<bin2>    <probability2>)
                            ...
                            (<binN>    <probabilityN>)
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
      type         | Type name: binned                 | word | yes  | -
      binnedDistribution | Distribution settings       | dict | yes  | -
      distribution | \<bin\>-\<probability\> pairs     | dict | yes  | -
      \<bin\>      | Particle size                     | scalar | yes | -
      \<probability\> | Probability of occurrence      | scalar | yes | -
    \endtable

SourceFiles
    binned.C

\*---------------------------------------------------------------------------*/

#ifndef distributionModels_binned_H
#define distributionModels_binned_H

#include "distributionModel.H"
#include "Field.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace distributionModels
{
    // Forward Declarations
    class binned;
}

// Forward Declarations
Istream& operator>>(Istream&, distributionModels::binned&);
Ostream& operator<<(Ostream&, const distributionModels::binned&);

namespace distributionModels
{

/*---------------------------------------------------------------------------*\
                           Class binned Declaration
\*---------------------------------------------------------------------------*/

class binned
:
    public distributionModel
{
    typedef VectorSpace<Vector<scalar>, scalar, 2> pair;

    // Private Data

        // List of (bin probability) pairs
        List<pair> xy_;

        //- Mean of the distribution
        scalar meanValue_;


    // Private Member Functions

        //- Initialise the distribution parameters
        void initialise();


public:

    //- Runtime type information
    TypeName("binned");

    static const char* header;


    // Constructors

        //- Construct from dictionary
        binned(const dictionary& dict, Random& rndGen);

        //- Construct from components
        //  Allows negative entries
        binned
        (
            const UList<scalar>& sampleData,
            const scalar binWidth,
            Random& rndGen
        );

        //- Copy construct
        binned(const binned& p);

        //- Construct and return a clone
        virtual autoPtr<distributionModel> clone() const
        {
            return autoPtr<distributionModel>(new binned(*this));
        }

        //- No copy assignment
        void operator=(const binned&) = delete;


    //- Destructor
    virtual ~binned() = default;


    // Member Functions

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
