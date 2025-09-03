/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2023 OpenCFD Ltd.
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
    Foam::scotchDecomp

Description
    Scotch domain decomposition.

    When run in parallel will collect the whole graph on to the master,
    decompose and send back. Use ptscotchDecomp for proper distributed
    decomposition.

    Coefficients dictionary: \a scotchCoeffs, \a coeffs.

    Quoting from the Scotch forum, on the 2008-08-22 10:09, Francois
    PELLEGRINI posted the following details:
    \verbatim
    RE: Graph mapping 'strategy' string

    Strategy handling in Scotch is a bit tricky. In order
    not to be confused, you must have a clear view of how they are built.
    Here are some rules:

    1- Strategies are made up of "methods" which are combined by means of
    "operators".

    2- A method is of the form "m{param=value,param=value,...}", where "m"
    is a single character (this is your first error: "f" is a method name,
    not a parameter name).

    3- There exist different sort of strategies : bipartitioning strategies,
    mapping strategies, ordering strategies, which cannot be mixed. For
    instance, you cannot build a bipartitioning strategy and feed it to a
    mapping method (this is your second error).

    To use the "mapCompute" routine, you must create a mapping strategy, not
    a bipartitioning one, and so use stratGraphMap() and not
    stratGraphBipart(). Your mapping strategy should however be based on the
    "recursive bipartitioning" method ("b"). For instance, a simple (and
    hence not very efficient) mapping strategy can be :

    "b{sep=f}"

    which computes mappings with the recursive bipartitioning method "b",
    this latter using the Fiduccia-Mattheyses method "f" to compute its
    separators.

    If you want an exact partition (see your previous post), try
    "b{sep=fx}".

    However, these strategies are not the most efficient, as they do not
    make use of the multi-level framework.

    To use the multi-level framework, try for instance:

    "b{sep=m{vert=100,low=h,asc=f}x}"

    The current default mapping strategy in Scotch can be seen by using the
    "-vs" option of program gmap. It is, to date:

    r
    {
        job=t,
        map=t,
        poli=S,
        sep=
        (
            m
            {
                asc=b
                {
                    bnd=
                    (
                        d{pass=40,dif=1,rem=1}
                     |
                    )
                    f{move=80,pass=-1,bal=0.002491},
                    org=f{move=80,pass=-1,bal=0.002491},
                    width=3
                },
                low=h{pass=10}
                f{move=80,pass=-1,bal=0.002491},
                type=h,
                vert=80,
                rat=0.8
            }
          | m
            {
                asc=b
                {
                    bnd=
                    (
                        d{pass=40,dif=1,rem=1}
                      |
                    )
                    f{move=80,pass=-1,bal=0.002491},
                    org=f{move=80,pass=-1,bal=0.002491},
                    width=3
                },
                low=h{pass=10}
                f{move=80,pass=-1,bal=0.002491},
                type=h,
                vert=80,
                rat=0.8
            }
        )
    }
    \endverbatim

    Given that this information was written in 2008, this example strategy will
    unlikely work as-is with the more recent Scotch versions. Therefore, the
    steps for getting the current default strategy from within Scotch, is to do
    the following steps:

    <ol>
    <li> Edit the file <tt>system/decomposeParDict</tt> and use the following
    settings:

    \verbatim
    method          scotch;

    scotchCoeffs
    {
        writeGraph true;
    }
    \endverbatim
    </li>

    <li> Run \c decomposePar. For example, it will write a file named
         <tt>region0.grf</tt>.
    </li>

    <li> Now, instead of using \c gmap, run \c gpart with the following
    command structure to get the default strategy:

    \verbatim
    gpart \<nProcs\> -vs \<grfFile\>
    \endverbatim

    where:

    <ul>
        <li> \<grfFile\> is the file that was obtained with the option
             <tt>writeGraph=true</tt>, namely <tt>region0.grf</tt>.
        </li>
        <li> \<nProcs\> is the \c numberOfSubdomains defined in the dictionary
             file.
        </li>
    </ul>
    </li>

    <li> At the end of the execution will be shown a long string, similar to
    the following example (complete line was cropped at <tt>[...]</tt>):

    \verbatim
    S    Strat=m{asc=b{width=3,bnd=d{pass=40,dif=1,rem=0}[...],type=h}
    \endverbatim
    </li>

    <li> Edit the file <tt>system/decomposeParDict</tt> once again and add
    the \c strategy entry as exemplified:

    \verbatim
    method          scotch;

    scotchCoeffs
    {
        //writeGraph true;
        strategy "m{asc=b{width=3,bnd=d{pass=40,dif=1,rem=0}[...],type=h}";
    }
    \endverbatim
    </li>

    <li> Finally, run \c decomposePar once again, to at least test if it
         works as intended.
    </li>

    </ol>


    Also support for multi-level decomposition by specifying inter-level
    communication weights:

    \verbatim
    numberOfSubdomains  2048;
    coeffs
    {
        // Divide into 64 clusters, each of 32 cores
        domains (64 32);
        // Inside a cluster the communication weight is 1% of that inbetween
        // clusters
        domainWeights (1 0.01);
    }
    \endverbatim

    Alternatively the first-level decomposition can be left out by assuming
    weights are 1 for that level:

    \verbatim
    numberOfSubdomains  2048;
    coeffs
    {
        // Divide into 2048/32 clusters
        domains (32);
        // Inside a cluster the communication weight is 1% of that inbetween
        // clusters
        domainWeights (0.01);
    }
    \endverbatim



Note
    \c gpart can be found in the current search path by adding the respective
    \c bin folder from the Scotch installation, namely by running the following
    commands:

    \verbatim
    source $(foamEtcFile config.sh/scotch)
    export PATH=$PATH:$SCOTCH_ARCH_PATH/bin
    \endverbatim

SourceFiles
    scotchDecomp.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_scotchDecomp_H
#define Foam_scotchDecomp_H

#include "metisLikeDecomp.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class scotchDecomp Declaration
\*---------------------------------------------------------------------------*/

class scotchDecomp
:
    public metisLikeDecomp
{
    // Private Data

        //- Output path and name for optional grf file.
        mutable fileName graphPath_;


protected:

    // Protected Member Functions

        //- Decompose non-parallel
        virtual label decomposeSerial
        (
            const labelList& adjncy,
            const labelList& xadj,
            const List<scalar>& cWeights,
            labelList& decomp
        ) const;

public:

    // Generated Methods

        //- No copy construct
        scotchDecomp(const scotchDecomp&) = delete;

        //- No copy assignment
        void operator=(const scotchDecomp&) = delete;


    //- Runtime type information
    TypeName("scotch");


    // Constructors

        //- Construct with number of domains (no coefficients or constraints)
        explicit scotchDecomp(const label numDomains);

        //- Construct given decomposition dictionary and optional region name
        explicit scotchDecomp
        (
            const dictionary& decompDict,
            const word& regionName = ""
        );


    //- Destructor
    virtual ~scotchDecomp() = default;


    // Member Functions

        //- Knows about coupled boundaries
        virtual bool parallelAware() const
        {
            return true;
        }

        //- Inherit all decompose methods
        using decompositionMethod::decompose;

        //- Return for every coordinate the wanted processor number.
        virtual labelList decompose
        (
            const polyMesh& mesh,
            const pointField& points = pointField::null(),
            const scalarField& pointWeights = scalarField::null()
        ) const;

        //- Return for every coordinate the wanted processor number.
        virtual labelList decompose
        (
            const polyMesh& mesh,
            const labelList& agglom,  //!< agglomeration: fine-to-coarse
            const pointField& agglomPoints,
            const scalarField& agglomWeights = scalarField::null()
        ) const;

        //- Return the wanted processor number for every cell.
        virtual labelList decompose
        (
            const CompactListList<label>& globalCellCells,
            const pointField& cc = pointField::null(),
            const scalarField& cWeights = scalarField::null()
        ) const;

        //- Return the wanted processor number for every cell.
        virtual labelList decompose
        (
            const labelListList& globalCellCells,
            const pointField& cc = pointField::null(),
            const scalarField& cWeights = scalarField::null()
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
