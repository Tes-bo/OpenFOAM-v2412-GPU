/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::IOPosition<CloudType>::IOPosition
(
    const CloudType& c,
    cloud::geometryType geomType
)
:
    regIOobject
    (
        IOobject
        (
            cloud::geometryTypeNames[geomType],
            c.time().timeName(),
            c,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    geometryType_(geomType),
    cloud_(c)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::IOPosition<CloudType>::write
(
    const bool /* writeOnProc */
) const
{
    return regIOobject::write(cloud_.size() > 0);
}


template<class CloudType>
bool Foam::IOPosition<CloudType>::writeData(Ostream& os) const
{
    os  << cloud_.size() << nl << token::BEGIN_LIST << nl;

    switch (geometryType_)
    {
        case cloud::geometryType::COORDINATES:
        {
            for (const auto& p : cloud_)
            {
                p.writeCoordinates(os);
                os  << nl;
            }
            break;
        }
        case cloud::geometryType::POSITIONS:
        {
            for (const auto& p : cloud_)
            {
                p.writePosition(os);
                os  << nl;
            }
            break;
        }
    }

    os  << token::END_LIST << endl;

    return os.good();
}


template<class CloudType>
void Foam::IOPosition<CloudType>::readData(Istream& is, CloudType& c)
{
    const polyMesh& mesh = c.pMesh();

    token tok(is);

    const bool newFormat = (geometryType_ == cloud::geometryType::COORDINATES);

    if (tok.isLabel())
    {
        const label len = tok.labelToken();

        // Read beginning of contents
        is.readBeginList(FUNCTION_NAME);

        for (label i=0; i<len; ++i)
        {
            // Read position only
            c.append
            (
                new typename CloudType::particleType
                (
                    mesh,
                    is,
                    false,
                    newFormat
                )
            );
        }

        // Read end of contents
        is.readEndList(FUNCTION_NAME);
    }
    else if (tok.isPunctuation(token::BEGIN_LIST))
    {
        is >> tok;
        while (!tok.isPunctuation(token::END_LIST))
        {
            is.putBack(tok);

            // Read position only
            c.append
            (
                new typename CloudType::particleType(mesh, is, false, newFormat)
            );
            is >> tok;
        }
    }
    else
    {
        FatalIOErrorInFunction(is)
            << "incorrect first token, expected <int> or '(', found "
            << tok.info() << nl
            << exit(FatalIOError);
    }

    is.check(FUNCTION_NAME);
}


// ************************************************************************* //
