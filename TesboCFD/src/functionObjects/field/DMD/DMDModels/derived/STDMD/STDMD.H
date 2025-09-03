/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2021 OpenCFD Ltd.
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
    Foam::DMDModels::STDMD

Description
    Streaming Total Dynamic Mode Decomposition (i.e. \c STDMD)
    is a variant of dynamic mode decomposition.

    Among other Dynamic Mode Decomposition (DMD) variants, \c STDMD is
    presumed to provide the general DMD method capabilities alongside
    economised and feasible memory and CPU usage.

    The code implementation corresponds to Figs. 15-16 of the first citation
    below, more broadly to Section 2.4.

    References:
    \verbatim
        DMD and mode-sorting algorithms (tags:K, HRDC, KZ, HWR):
            Kiewat, M. (2019).
            Streaming modal decomposition approaches for vehicle aerodynamics.
            PhD thesis. Munich: Technical University of Munich.
            URL:mediatum.ub.tum.de/doc/1482652/1482652.pdf

            Hemati, M. S., Rowley, C. W.,
            Deem, E. A., & Cattafesta, L. N. (2017).
            De-biasing the dynamic mode decomposition
            for applied Koopman spectral analysis of noisy datasets.
            Theoretical and Computational Fluid Dynamics, 31(4), 349-368.
            DOI:10.1007/s00162-017-0432-2

            Kou, J., & Zhang, W. (2017).
            An improved criterion to select
            dominant modes from dynamic mode decomposition.
            European Journal of Mechanics-B/Fluids, 62, 109-129.
            DOI:10.1016/j.euromechflu.2016.11.015

            Hemati, M. S., Williams, M. O., & Rowley, C. W. (2014).
            Dynamic mode decomposition for large and streaming datasets.
            Physics of Fluids, 26(11), 111701.
            DOI:10.1063/1.4901016

        Parallel classical Gram-Schmidt process (tag:Ka):
            Katagiri, T. (2003).
            Performance evaluation of parallel
            Gram-Schmidt re-orthogonalization methods.
            In: Palma J. M. L. M., Sousa A. A., Dongarra J., Hernández V. (eds)
            High Performance Computing for Computational Science — VECPAR 2002.
            Lecture Notes in Computer Science, vol 2565, p. 302-314.
            Berlin, Heidelberg: Springer.
            DOI:10.1007/3-540-36569-9_19

        Parallel direct tall-skinny QR decomposition (tags:BGD, DGHL):
            Benson, A. R., Gleich, D. F., & Demmel, J. (2013).
            Direct QR factorizations for
            tall-and-skinny matrices in MapReduce architectures.
            2013 IEEE International Conference on Big Data.
            DOI:10.1109/bigdata.2013.6691583

            Demmel, J., Grigori, L., Hoemmen, M., & Langou, J. (2012).
            Communication-optimal parallel
            and sequential QR and LU factorizations.
            SIAM Journal on Scientific Computing, 34(1), A206-A239.
            DOI:10.1137/080731992
    \endverbatim

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    DMD1
    {
        // Mandatory entries
        DMDModel            STDMD;

        // Conditional mandatory entries

            // Option-1
            interval            5.5;

            // Option-2
            executeInterval     10;

        // Optional entries
        modeSorter          kiewat;
        nGramSchmidt        5;
        maxRank             50;
        nModes              50;
        fMin                0;
        fMax                1000000000;
        nAgglomerationProcs 20;

        // Optional entries (not recommended to change)
        minBasis            0.00000001;
        minEVal             0.00000001;
        sortLimiter         500.0;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Reqd  | Deflt
      DMDModel     | Type: STDMD                        | word | yes   | -
      interval     | STDMD time-step size [s]           | scalar | cndtnl <!--
               --> | executeInterval*(current time-step of the simulation)
      modeSorter   | Mode-sorting algorithm model       | word | no    <!--
               --> | firstSnapshot
      nModes       | Number of output modes in input frequency range   <!--
               -->                                        | label | no | GREAT
      maxRank  | Max columns in accumulated matrices      | label | no | GREAT
      nGramSchmidt | Number of Gram-Schmidt iterations    | label | no | 5
      fMin         | Min (non-negative) output frequency  | label | no | 0
      fMax         | Max output frequency                 | label | no | GREAT
      nAgglomerationProcs | Number of processors at each agglomeration <!--
               --> unit during the computation of reduced Koopman      <!--
               --> operator                               | label | no | 20
      minBasis | Orthogonal basis expansion threshold     | scalar| no | 1e-8
      minEVal  | Min eigenvalue for below eigenvalues are omitted      <!--
               -->                                        | scalar| no | 1e-8
      sortLimiter | Max allowable magnitude for <!--
               -->  mag(eigenvalue)*(number of DMD steps) to be used in <!--
               -->  modeSorter={kiewat,kouZhang} to avoid overflow errors <!--
               -->                                        | scalar | no | 500.0
    \endtable

    Options for the \c modeSorter entry:
    \verbatim
      kiewat        | Modified weighted-amplitude scaling method
      kouZhang      | Weighted-amplitude scaling method
      firstSnapshot | First-snapshot amplitude magnitude method
    \endverbatim

Note
  - To specify the STDMD time-step size (not necessarily equal to the
    time step of the simulation), entries of either \c interval or
    \c executeInterval must be available (highest to lowest precedence). While
    \c interval allows users to directly specify the STDMD time-step size
    in seconds, in absence of \c interval, for convenience,
    \c executeInterval allows users to compute the STDMD time-step internally
    by multiplying itself with the current time-step size of the simulation.
  - Limitation: Restart is currently not available since intermediate writing
    of STDMD matrices are not supported.
  - Limitation: Non-physical input (e.g. full-zero fields) can upset STDMD.
  - Warning: In \c modeSorter algorithms of \c kiewat and \c kouZhang, a
    function of \f$power(x,y)\f$ exists where \c x is the magnitude of an
    eigenvalue, and \c y is the number of STDMD snapshots. This function poses
    a risk of overflow errors since, for example, if \c x ends up above 1.5 with
    500 snapshots or more, this function automatically throws the floating point
    error meaning overflow. Therefore, the domain-range of this function is
    heuristically constrained by the optional entry \c sortLimiter which skips
    the evaluation of this function for a given
    mag(eigenvalue)*(no. of STDMD steps), i.e. x*y, whose magnitude is larger
    than \c sortLimiter.

See also
  - Foam::functionObjects::DMD
  - Foam::DMDModel

SourceFiles
    STDMD.C
    STDMDTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef DMDModels_STDMD_H
#define DMDModels_STDMD_H

#include "DMDModel.H"
#include "SquareMatrix.H"
#include "RectangularMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace DMDModels
{

/*---------------------------------------------------------------------------*\
                              Class STDMD Declaration
\*---------------------------------------------------------------------------*/

class STDMD
:
    public DMDModel
{
    typedef SquareMatrix<scalar> SMatrix;
    typedef RectangularMatrix<scalar> RMatrix;
    typedef RectangularMatrix<complex> RCMatrix;


    // Private Enumerations

        //- Options for the mode-sorting algorithm
        enum modeSorterType : char
        {
            FIRST_SNAPSHOT = 0, //!< "First-snapshot amplitude magnitude method"
            KIEWAT,            //!< "Modified weighted-amplitude scaling method"
            KOU_ZHANG          //!< "Weighted-amplitude scaling method"
        };

        //- Names for modeSorterType
        static const Enum<modeSorterType> modeSorterTypeNames;


    // Private Data

        //- Mode-sorting algorithm
        enum modeSorterType modeSorter_;

        //- Accumulated-in-time unitary similarity matrix produced by the
        //- orthonormal decomposition of the augmented snapshot matrix 'z'
        //  (K:Eq. 60)
        RMatrix Q_;

        //- Accumulated-in-time (squared) upper triangular matrix produced by
        //- the orthonormal decomposition of the augmented snapshot matrix 'z'
        //  (K:Eq. 64)
        SMatrix G_;

        //- Upper half of 'Q'
        RMatrix Qupper_;

        //- Lower half of 'Q'
        RMatrix Qlower_;

        //- Moore-Penrose pseudo-inverse of 'R' produced by
        //- the QR decomposition of the last time-step 'Q'
        RMatrix RxInv_;

        //- Eigenvalues of modes
        List<complex> evals_;

        //- Eigenvectors of modes
        RCMatrix evecs_;

        //- (Non-negative) frequencies of modes
        List<scalar> freqs_;

        //- Indices of 'freqs' where frequencies are
        //- non-negative and within ['fMin', 'fMax']
        DynamicList<label> freqsi_;

        //- Amplitudes of modes
        List<complex> amps_;

        //- (Descending) magnitudes of (complex) amplitudes of modes
        List<scalar> mags_;

        //- Indices of 'mags'
        List<label> magsi_;

        //- Names of operand patches
        const wordRes patches_;

        //- Name of operand field
        const word fieldName_;

        //- First-processed snapshot required by the mode-sorting
        //- algorithms at the final output computations (K:p. 43)
        word timeName0_;

        //- Min value to execute orthogonal basis expansion of 'Q' and 'G'
        scalar minBasis_;

        //- Min eigenvalue magnitude below where 'evals' are omitted
        scalar minEval_;

        //- STDMD time-step size that equals to
        //- (executeInterval of DMD)*(deltaT of simulation) [s]
        scalar dt_;

        //- Maximum allowable magnitude for mag(eigenvalue)*(number of
        //- STDMD steps) to be used in modeSorter={kiewat,kouZhang} to
        //- avoid overflow errors
        scalar sortLimiter_;

        //- Min frequency: Output only entries corresponding to freqs >= 'fMin'
        label fMin_;

        //- Max frequency: Output only entries corresponding to freqs <= 'fMax'
        label fMax_;

        //- Number of maximum iterations of the classical
        //- Gram-Schmidt procedure for the orthonormalisation
        label nGramSchmidt_;

        //- Maximum allowable matrix column for 'Q' and 'G'
        //  'Q' is assumed to always have full-rank, thus 'Q.n() = rank'
        label maxRank_;

        //- Current execution-step index of STDMD,
        //- not necessarily that of the simulation
        label step_;

        //- Number of output modes within input frequency range
        //- starting from the most energetic mode
        label nModes_;

        //- Number of processors at each agglomeration unit
        //- during the computation of reduced Koopman operator
        label nAgglomerationProcs_;

        //- (Internal) Flag to tag snapshots which are effectively empty
        bool empty_;


    // Private Member Functions

        // Evaluation

            //- Return (parallel) L2-norm of a given column vector
            scalar L2norm(const RMatrix& z) const;

            //- Execute (parallel) classical Gram-Schmidt
            //- process to orthonormalise 'ez' (Ka:Fig. 5)
            RMatrix orthonormalise(RMatrix ez) const;

            //- Expand orthonormal bases 'Q' and 'G' by stacking a column
            //- '(ez/ezNorm)' to 'Q', and a row (Zero) and column (Zero)
            //- to 'G' if '(ezNorm/zNorm > minBasis)'
            void expand(const RMatrix& ez, const scalar ezNorm);

            //- Update 'G' before the potential orthonormal basis compression
            void updateG(const RMatrix& z);

            //- Compress orthonormal basis for 'Q' and 'G' if '(Q.n()>maxRank)'
            void compress();

            //- Return reduced Koopman operator 'Atilde' (K:Eq. 78)
            //  Also fills 'RxInv'.
            //  The function was not divided into subsections to ensure
            //  minimal usage of memory, hence the complexity of the function
            SMatrix reducedKoopmanOperator();

            //- Compute eigenvalues and eigenvectors of 'Atilde'
            //  Remove any eigenvalues whose magnitude is smaller than
            //  'minEVal' while keeping the order of elements the same
            bool eigendecomposition(SMatrix& Atilde);

            //- Compute and filter frequencies and its indices (K:Eq. 81)
            //  Indices of 'freqs' where frequencies are
            //  non-negative and within ['fMin', 'fMax']
            void frequencies();

            //- Compute amplitudes
            void amplitudes();

            //- Compute magnitudes and ordered magnitude indices
            //  Includes advanced sorting methods using
            //  the chosen weighted amplitude scaling method
            void magnitudes();

            //- Eigenvalue weighted amplitude scaling (KZ:Eq. 33)
            //- Modified eigenvalue weighted amplitude scaling (K)
            scalar sorter
            (
                const List<scalar>& weight,
                const complex& amplitude,
                const complex& eigenvalue,
                const scalar modeNorm
            ) const;

            //- Compute and write mode dynamics
            virtual bool dynamics();

            //- Compute and write modes
            virtual bool modes();

            //- Select field type for modes
            template<class Type>
            bool modes();

            //- Compute modes based on input field type
            template<class GeoFieldType>
            bool calcModes();

            //- Compute a mode for a given scalar-based input field
            template<class GeoFieldType>
            typename std::enable_if
            <
                std::is_same<scalar, typename GeoFieldType::value_type>::value,
                void
            >::type calcMode
            (
                GeoFieldType& modeRe,
                GeoFieldType& modeIm,
                const RMatrix& primitiveMode,
                const label magi,
                const label rowi = 0
            );

            //- Compute a mode for a given non-scalar-based input field
            template<class GeoFieldType>
            typename std::enable_if
            <
                !std::is_same<scalar, typename GeoFieldType::value_type>::value,
                void
            >::type calcMode
            (
                GeoFieldType& modeRe,
                GeoFieldType& modeIm,
                const RMatrix& primitiveMode,
                const label magi,
                const label rowi = 0
            );


        // I-O

            //- Write objects of dynamics
            void writeToFile(const word& fileName) const;

            // Notifying the compiler that we want both writeToFile functions
            using Foam::functionObjects::writeFile::writeToFile;

            //- Write file header information
            void writeFileHeader(Ostream& os) const;

            //- Filter objects of dynamics according to 'freqsi' and 'magsi'
            void filter();

            //- Return a new List containing elems of List at 'indices'
            template<class Type>
            void filterIndexed
            (
                List<Type>& lst,
                const UList<label>& indices
            );

            //- Return a new Matrix containing columns of Matrix at 'indices'
            template<class MatrixType>
            void filterIndexed
            (
                MatrixType& lst,
                const UList<label>& indices
            );


public:

    //- Runtime type information
    TypeName("STDMD");


    // Constructors

        //- Construct from components
        STDMD
        (
            const fvMesh& mesh,
            const word& name,
            const dictionary& dict
        );

        //- No copy construct
        STDMD(const STDMD&) = delete;

        //- No copy assignment
        void operator=(const STDMD&) = delete;


    //- Destructor
    virtual ~STDMD() = default;


    // Member Functions

        // Evaluation

            //- Initialise 'Q' and 'G' (both require the first two snapshots)
            virtual bool initialise(const RMatrix& z);

            //- Incremental orthonormal basis update (K:Fig. 15)
            virtual bool update(const RMatrix& z);

            //- Compute and write modes and
            //- mode dynamics of model data members
            virtual bool fit();


        // IO

            //- Read STDMD settings
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace DMDModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "STDMDTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
