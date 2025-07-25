/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "molecule.H"
#include "IOstreams.H"
#include "moleculeCloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const std::size_t Foam::molecule::sizeofFields
(
    offsetof(molecule, siteForces_) - offsetof(molecule, Q_)
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::molecule::molecule
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields,
    bool newFormat
)
:
    particle(mesh, is, readFields, newFormat),
    Q_(Zero),
    v_(Zero),
    a_(Zero),
    pi_(Zero),
    tau_(Zero),
    specialPosition_(Zero),
    potentialEnergy_(0.0),
    rf_(Zero),
    special_(0),
    id_(0),
    siteForces_(),
    sitePositions_()
{
    if (readFields)
    {
        if (is.format() == IOstreamOption::ASCII)
        {
            is  >> Q_
                >> v_
                >> a_
                >> pi_
                >> tau_
                >> specialPosition_
                >> potentialEnergy_
                >> rf_
                >> special_
                >> id_;
        }
        else if (!is.checkLabelSize<>() || !is.checkScalarSize<>())
        {
            // Non-native label or scalar size

            is.beginRawRead();

            readRawScalar(is, Q_.data(), tensor::nComponents);
            readRawScalar(is, v_.data(), vector::nComponents);
            readRawScalar(is, a_.data(), vector::nComponents);
            readRawScalar(is, pi_.data(), vector::nComponents);
            readRawScalar(is, tau_.data(), vector::nComponents);
            readRawScalar(is, specialPosition_.data(), vector::nComponents);
            readRawScalar(is, &potentialEnergy_);
            readRawScalar(is, rf_.data(), tensor::nComponents);
            readRawLabel(is, &special_);
            readRawLabel(is, &id_);

            is.endRawRead();
        }
        else
        {
            is.read(reinterpret_cast<char*>(&Q_), sizeofFields);
        }

        is  >> siteForces_ >> sitePositions_;
    }

    is.check(FUNCTION_NAME);
}


void Foam::molecule::readFields(Cloud<molecule>& mC)
{
    const bool readOnProc = mC.size();

    particle::readFields(mC);

    IOField<tensor> Q(mC.newIOobject("Q", IOobject::MUST_READ), readOnProc);
    mC.checkFieldIOobject(mC, Q);

    IOField<vector> v(mC.newIOobject("v", IOobject::MUST_READ), readOnProc);
    mC.checkFieldIOobject(mC, v);

    IOField<vector> a(mC.newIOobject("a", IOobject::MUST_READ), readOnProc);
    mC.checkFieldIOobject(mC, a);

    IOField<vector> pi(mC.newIOobject("pi", IOobject::MUST_READ), readOnProc);
    mC.checkFieldIOobject(mC, pi);

    IOField<vector> tau
    (
        mC.newIOobject("tau", IOobject::MUST_READ),
        readOnProc
    );
    mC.checkFieldIOobject(mC, tau);

    IOField<vector> specialPosition
    (
        mC.newIOobject("specialPosition", IOobject::MUST_READ),
        readOnProc
    );
    mC.checkFieldIOobject(mC, specialPosition);

    IOField<label> special
    (
        mC.newIOobject("special", IOobject::MUST_READ),
        readOnProc
    );
    mC.checkFieldIOobject(mC, special);

    IOField<label> id(mC.newIOobject("id", IOobject::MUST_READ), readOnProc);
    mC.checkFieldIOobject(mC, id);

    label i = 0;
    for (molecule& mol : mC)
    {
        mol.Q_ = Q[i];
        mol.v_ = v[i];
        mol.a_ = a[i];
        mol.pi_ = pi[i];
        mol.tau_ = tau[i];
        mol.specialPosition_ = specialPosition[i];
        mol.special_ = special[i];
        mol.id_ = id[i];

        ++i;
    }
}


void Foam::molecule::writeFields(const Cloud<molecule>& mC)
{
    particle::writeFields(mC);

    const label np = mC.size();
    const bool writeOnProc = mC.size();

    IOField<tensor> Q(mC.newIOobject("Q", IOobject::NO_READ), np);
    IOField<vector> v(mC.newIOobject("v", IOobject::NO_READ), np);
    IOField<vector> a(mC.newIOobject("a", IOobject::NO_READ), np);
    IOField<vector> pi(mC.newIOobject("pi", IOobject::NO_READ), np);
    IOField<vector> tau(mC.newIOobject("tau", IOobject::NO_READ), np);
    IOField<vector> specialPosition
    (
        mC.newIOobject("specialPosition", IOobject::NO_READ),
        np
    );
    IOField<label> special(mC.newIOobject("special", IOobject::NO_READ), np);
    IOField<label> id(mC.newIOobject("id", IOobject::NO_READ), np);

    // Post processing fields

    IOField<vector> piGlobal
    (
        mC.newIOobject("piGlobal", IOobject::NO_READ),
        np
    );

    IOField<vector> tauGlobal
    (
        mC.newIOobject("tauGlobal", IOobject::NO_READ),
        np
    );

    IOField<vector> orientation1
    (
        mC.newIOobject("orientation1", IOobject::NO_READ),
        np
    );

    IOField<vector> orientation2
    (
        mC.newIOobject("orientation2", IOobject::NO_READ),
        np
    );

    IOField<vector> orientation3
    (
        mC.newIOobject("orientation3", IOobject::NO_READ),
        np
    );

    label i = 0;
    for (const molecule& mol : mC)
    {
        Q[i] = mol.Q_;
        v[i] = mol.v_;
        a[i] = mol.a_;
        pi[i] = mol.pi_;
        tau[i] = mol.tau_;
        specialPosition[i] = mol.specialPosition_;
        special[i] = mol.special_;
        id[i] = mol.id_;

        piGlobal[i] = mol.Q_ & mol.pi_;
        tauGlobal[i] = mol.Q_ & mol.tau_;

        orientation1[i] = mol.Q_ & vector(1,0,0);
        orientation2[i] = mol.Q_ & vector(0,1,0);
        orientation3[i] = mol.Q_ & vector(0,0,1);

        ++i;
    }

    Q.write(writeOnProc);
    v.write(writeOnProc);
    a.write(writeOnProc);
    pi.write(writeOnProc);
    tau.write(writeOnProc);
    specialPosition.write(writeOnProc);
    special.write(writeOnProc);
    id.write(writeOnProc);

    piGlobal.write(writeOnProc);
    tauGlobal.write(writeOnProc);

    orientation1.write(writeOnProc);
    orientation2.write(writeOnProc);
    orientation3.write(writeOnProc);

    Info<< "writeFields " << mC.name() << endl;

    if (isA<moleculeCloud>(mC))
    {
        const moleculeCloud& m = dynamic_cast<const moleculeCloud&>(mC);

        m.writeXYZ
        (
            m.mesh().time().timePath()/cloud::prefix/"moleculeCloud.xmol"
        );
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const molecule& mol)
{
    if (os.format() == IOstreamOption::ASCII)
    {
        os  << token::SPACE << static_cast<const particle&>(mol)
            << token::SPACE << mol.Q_
            << token::SPACE << mol.v_
            << token::SPACE << mol.a_
            << token::SPACE << mol.pi_
            << token::SPACE << mol.tau_
            << token::SPACE << mol.specialPosition_
            << token::SPACE << mol.potentialEnergy_
            << token::SPACE << mol.rf_
            << token::SPACE << mol.special_
            << token::SPACE << mol.id_
            << token::SPACE << mol.siteForces_
            << token::SPACE << mol.sitePositions_;
    }
    else
    {
        os  << static_cast<const particle&>(mol);
        os.write
        (
            reinterpret_cast<const char*>(&mol.Q_),
            molecule::sizeofFields
        );
        os  << mol.siteForces_ << mol.sitePositions_;
    }

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
