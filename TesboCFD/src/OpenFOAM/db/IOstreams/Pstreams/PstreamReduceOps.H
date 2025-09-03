/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2024 OpenCFD Ltd.
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

InNamespace
    Foam

Description
    Inter-processor communication reduction functions.

\*---------------------------------------------------------------------------*/

#ifndef Foam_PstreamReduceOps_H
#define Foam_PstreamReduceOps_H

#include "Pstream.H"
#include "FixedList.H"
#include "ops.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Reduce inplace (cf. MPI Allreduce)
//- using linear/tree communication schedule
template<class T, class BinaryOp>
void reduce
(
    T& value,
    const BinaryOp& bop,
    const int tag = UPstream::msgType(),
    const label comm = UPstream::worldComm
)
{
    if (UPstream::is_parallel(comm))
    {
        if (UPstream::warnComm >= 0 && comm != UPstream::warnComm)
        {
            Perr<< "** reducing:" << value << " with comm:" << comm << endl;
            error::printStack(Perr);
        }
        Pstream::gather(value, bop, tag, comm);
        Pstream::broadcast(value, comm);
    }
}


//- Reduce inplace (cf. MPI Allreduce)
//- multiple values (same size on all ranks!)
template<class T, class BinaryOp>
void reduce
(
    T values[],
    const int size,
    const BinaryOp&,
    const int tag,
    const label comm
)
{
    NotImplemented;
}

//- Non-blocking reduce inplace (cf. MPI Iallreduce)
//- single value. Sets request.
template<class T, class BinaryOp>
void reduce
(
    T& Value,
    const BinaryOp&,
    const int tag,
    const label comm,
    UPstream::Request& req
)
{
    NotImplemented;
}

//- Non-blocking reduce inplace (cf. MPI Iallreduce)
//- single value. Sets request.
template<class T, class BinaryOp>
void reduce
(
    T& Value,
    const BinaryOp&,
    const int tag,
    const label comm,
    label& request
)
{
    NotImplemented;
}

//- Non-blocking reduce inplace (cf. MPI Iallreduce)
//- of multiple values (same size on all ranks!). Sets request.
template<class T, class BinaryOp>
void reduce
(
    T values[],
    const int size,
    const BinaryOp&,
    const int tag,
    const label comm,
    UPstream::Request& req
)
{
    NotImplemented;
}

//- Non-blocking reduce inplace (cf. MPI Iallreduce)
//- of multiple values (same size on all ranks!). Sets request.
template<class T, class BinaryOp>
void reduce
(
    T values[],
    const int size,
    const BinaryOp&,
    const int tag,
    const label comm,
    label& request
)
{
    NotImplemented;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Special reductions for bool

//- Logical (and) inplace reduction. Uses UPstream::reduceAnd
void reduce
(
    bool& value,
    const andOp<bool>&,
    const int tag = UPstream::msgType(),  /*!< (ignored) */
    const label comm = UPstream::worldComm
);

//- Logical (or) inplace reduction. Uses UPstream::reduceOr
void reduce
(
    bool& value,
    const orOp<bool>&,
    const int tag = UPstream::msgType(),  /*!< (ignored) */
    const label comm = UPstream::worldComm
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Common reductions

#undef  Pstream_CommonReductions
#define Pstream_CommonReductions(Native)                                      \
                                                                              \
/*! \brief Reduce (min) multiple Native values (same size on all ranks!) */   \
void reduce                                                                   \
(                                                                             \
    Native values[],                                                          \
    const int size,                                                           \
    const minOp<Native>&,                                                     \
    const int tag = UPstream::msgType(),  /*!< (ignored) */                   \
    const label comm = UPstream::worldComm                                    \
);                                                                            \
                                                                              \
/*! \brief Reduce (min) single Native value */                                \
void reduce                                                                   \
(                                                                             \
    Native& value,                                                            \
    const minOp<Native>&,                                                     \
    const int tag = UPstream::msgType(),  /*!< (ignored) */                   \
    const label comm = UPstream::worldComm                                    \
);                                                                            \
                                                                              \
/*! \brief Reduce (min) multiple Native values */                             \
template<unsigned N>                                                          \
inline void reduce                                                            \
(                                                                             \
    FixedList<Native, N>& values,                                             \
    const minOp<Native>&,                                                     \
    const int tag = UPstream::msgType(),  /*!< (ignored) */                   \
    const label comm = UPstream::worldComm                                    \
)                                                                             \
{                                                                             \
    reduce(values.data(), int(values.size()), minOp<Native>(), tag, comm);    \
}                                                                             \
                                                                              \
/*! \brief Reduce (max) multiple Native values (same size on all ranks!) */   \
void reduce                                                                   \
(                                                                             \
    Native values[],                                                          \
    const int size,                                                           \
    const maxOp<Native>&,                                                     \
    const int tag,  /*!< (ignored) */                                         \
    const label comm                                                          \
);                                                                            \
                                                                              \
/*! \brief Reduce (max) single Native value */                                \
void reduce                                                                   \
(                                                                             \
    Native& value,                                                            \
    const maxOp<Native>&,                                                     \
    const int tag = UPstream::msgType(),  /*!< (ignored) */                   \
    const label comm = UPstream::worldComm                                    \
);                                                                            \
                                                                              \
/*! \brief Reduce (max) multiple Native values */                             \
template<unsigned N>                                                          \
inline void reduce                                                            \
(                                                                             \
    FixedList<Native, N>& values,                                             \
    const maxOp<Native>&,                                                     \
    const int tag = UPstream::msgType(),  /*!< (ignored) */                   \
    const label comm = UPstream::worldComm                                    \
)                                                                             \
{                                                                             \
    reduce(values.data(), int(values.size()), maxOp<Native>(), tag, comm);    \
}                                                                             \
                                                                              \
/*! \brief Reduce (sum) multiple Native values (same size on all ranks!) */   \
void reduce                                                                   \
(                                                                             \
    Native values[],                                                          \
    const int size,                                                           \
    const sumOp<Native>&,                                                     \
    const int tag,  /*!< (ignored) */                                         \
    const label comm                                                          \
);                                                                            \
                                                                              \
/*! \brief Reduce (sum) single Native value */                                \
void reduce                                                                   \
(                                                                             \
    Native& value,                                                            \
    const sumOp<Native>&,                                                     \
    const int tag = UPstream::msgType(),  /*!< (ignored) */                   \
    const label comm = UPstream::worldComm                                    \
);                                                                            \
                                                                              \
/*! \brief Reduce (sum) multiple Native values */                             \
template<unsigned N>                                                          \
inline void reduce                                                            \
(                                                                             \
    FixedList<Native, N>& values,                                             \
    const sumOp<Native>&,                                                     \
    const int tag = UPstream::msgType(),  /*!< (ignored) */                   \
    const label comm = UPstream::worldComm                                    \
)                                                                             \
{                                                                             \
    reduce(values.data(), int(values.size()), sumOp<Native>(), tag, comm);    \
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Floating-point reductions

#undef  Pstream_FloatReductions
#define Pstream_FloatReductions(Native)                                       \
                                                                              \
Pstream_CommonReductions(Native);                                             \
                                                                              \
/*! \brief Non-blocking reduce (sum) multiple Native values. Sets request */  \
void reduce                                                                   \
(                                                                             \
    Native values[],                                                          \
    const int size,                                                           \
    const sumOp<Native>&,                                                     \
    const int tag,  /*!< (ignored) */                                         \
    const label comm,                                                         \
    UPstream::Request& req  /*!< [out] request information */                 \
);                                                                            \
                                                                              \
/*! \brief Non-blocking reduce (sum) multiple Native values. Sets request */  \
/*! \deprecated prefer version with UPstream::Request */                      \
void reduce                                                                   \
(                                                                             \
    Native values[],                                                          \
    const int size,                                                           \
    const sumOp<Native>&,                                                     \
    const int tag,  /*!< (ignored) */                                         \
    const label comm,                                                         \
    label& requestID                                                          \
);                                                                            \
                                                                              \
/*! \brief Non-blocking reduce (sum) single Native value. Sets request */     \
void reduce                                                                   \
(                                                                             \
    Native& value,                                                            \
    const sumOp<Native>&,                                                     \
    const int tag,  /*!< (ignored) */                                         \
    const label comm,                                                         \
    UPstream::Request& req  /*!< [out] request information */                 \
);                                                                            \
                                                                              \
/*! \brief Non-blocking reduce (sum) single Native value. Sets request */     \
/*! \deprecated prefer version with UPstream::Request */                      \
void reduce                                                                   \
(                                                                             \
    Native& value,                                                            \
    const sumOp<Native>&,                                                     \
    const int tag,  /*!< (ignored) */                                         \
    const label comm,                                                         \
    label& requestID                                                          \
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Bitwise reductions

#undef  Pstream_BitwiseReductions
#define Pstream_BitwiseReductions(Native)                                     \
                                                                              \
/*! \brief Reduce (bit-or) multiple Native values (same size on all ranks!)*/ \
void reduce                                                                   \
(                                                                             \
    Native values[],                                                          \
    const int size,                                                           \
    const bitOrOp<Native>&,                                                   \
    const int tag = UPstream::msgType(),  /*!< (ignored) */                   \
    const label comm = UPstream::worldComm                                    \
);                                                                            \
                                                                              \
/*! \brief Reduce (bit-or) single Native value */                             \
void reduce                                                                   \
(                                                                             \
    Native& value,                                                            \
    const bitOrOp<Native>&,                                                   \
    const int tag = UPstream::msgType(),  /*!< (ignored) */                   \
    const label comm = UPstream::worldComm                                    \
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Pstream_CommonReductions(int32_t);
Pstream_CommonReductions(int64_t);
Pstream_CommonReductions(uint32_t);
Pstream_CommonReductions(uint64_t);

Pstream_FloatReductions(float);
Pstream_FloatReductions(double);

Pstream_BitwiseReductions(unsigned char);
Pstream_BitwiseReductions(unsigned int);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef Pstream_CommonReductions
#undef Pstream_FloatReductions
#undef Pstream_BitwiseReductions


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Reduce inplace (cf. MPI Allreduce)
//- the sum of both value and count (for averaging)
template<class T>
void sumReduce
(
    T& value,
    label& count,
    const int tag = UPstream::msgType(),
    const label comm = UPstream::worldComm
)
{
    if (UPstream::is_parallel(comm))
    {
        Foam::reduce(value, sumOp<T>(), tag, comm);
        Foam::reduce(count, sumOp<label>(), tag, comm);
    }
}


// Floating-point sum-reduce

#undef  Pstream_SumReduce
#define Pstream_SumReduce(Native)                                             \
                                                                              \
/*! \brief Sum of both Native value and count (for averaging) */              \
void sumReduce                                                                \
(                                                                             \
    Native& value,                                                            \
    label& count,                                                             \
    const int tag = UPstream::msgType(),  /*!< (ignored) */                   \
    const label comm = UPstream::worldComm                                    \
);


Pstream_SumReduce(float);
Pstream_SumReduce(double);

#undef Pstream_SumReduce


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Convenience wrappers - defined after all specialisations are known

//- Perform reduction on a copy, using specified binary operation
//  \return the resulting value
template<class T, class BinaryOp>
T returnReduce
(
    const T& value,
    const BinaryOp& bop,
    const int tag = UPstream::msgType(),
    const label comm = UPstream::worldComm
)
{
    T work(value);
    Foam::reduce(work, bop, tag, comm);
    return work;
}


//- Perform logical (and) MPI Allreduce on a copy. Uses UPstream::reduceAnd
//  \return the resulting value
inline bool returnReduceAnd
(
    const bool value,
    const label comm = UPstream::worldComm
)
{
    bool work(value);
    UPstream::reduceAnd(work, comm);
    return work;
}


//- Perform logical (or) MPI Allreduce on a copy. Uses UPstream::reduceOr
//  \return the resulting value
inline bool returnReduceOr
(
    const bool value,
    const label comm = UPstream::worldComm
)
{
    bool work(value);
    UPstream::reduceOr(work, comm);
    return work;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
