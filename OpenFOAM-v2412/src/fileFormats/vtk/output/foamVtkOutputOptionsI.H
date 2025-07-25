/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2022 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::vtk::outputOptions::outputOptions()
:
    fmtType_(formatType::INLINE_ASCII),
    precision_(IOstream::defaultPrecision())
{}


inline Foam::vtk::outputOptions::outputOptions
(
    enum vtk::formatType fmtType
)
:
    fmtType_(fmtType),
    precision_(IOstream::defaultPrecision())
{}


inline Foam::vtk::outputOptions::outputOptions
(
    enum vtk::formatType fmtType,
    unsigned prec
)
:
    fmtType_(fmtType),
    precision_(prec)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

inline Foam::autoPtr<Foam::vtk::formatter>
Foam::vtk::outputOptions::newFormatter(std::ostream& os) const
{
    return vtk::newFormatter(os, fmtType_, precision_);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

inline Foam::vtk::formatType Foam::vtk::outputOptions::fmt() const noexcept
{
    return fmtType_;
}


inline Foam::word Foam::vtk::outputOptions::ext(vtk::fileTag contentType) const
{
    return
    (
        legacy()
      ? vtk::legacy::fileExtension
      : vtk::fileExtension[contentType]
    );
}


inline bool Foam::vtk::outputOptions::legacy() const noexcept
{
    return
    (
        fmtType_ == formatType::LEGACY_ASCII
     || fmtType_ == formatType::LEGACY_BINARY
    );
}


inline bool Foam::vtk::outputOptions::xml() const noexcept
{
    return !legacy();
}


inline bool Foam::vtk::outputOptions::append() const noexcept
{
    return
    (
        fmtType_ == formatType::APPEND_BASE64
     || fmtType_ == formatType::APPEND_BINARY
    );
}


inline bool Foam::vtk::outputOptions::insitu() const noexcept
{
    return !append();
}


inline bool Foam::vtk::outputOptions::ascii() const noexcept
{
    return !(unsigned(fmtType_) & 0x0F);
}


inline unsigned Foam::vtk::outputOptions::precision() const noexcept
{
    return precision_;
}


// ************************************************************************* //
