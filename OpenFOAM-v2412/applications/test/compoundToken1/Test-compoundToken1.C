/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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

Description
    Test token construct assign etc.
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IOobject.H"
#include "IOstreams.H"
#include "IFstream.H"
#include "cpuTime.H"
#include "labelList.H"
#include "DynamicList.H"

namespace Foam
{

template<class OS>
OS& printTypeCode(OS& os, char typeCode)
{
    os << int(static_cast<unsigned char>(typeCode));
    return os;
}


/*---------------------------------------------------------------------------*\
                          Class IFstream Declaration
\*---------------------------------------------------------------------------*/

bool test_pending = false;

class IFstreamDelayed
:
    public IFstream
{
    virtual bool readCompoundToken(token& tok, const word& type) override
    {
        auto& is = *this;

        bool delay = true;

        // Low-level: get next valid character (after comments)
        // and branch based on it being a '{' or not

        char c = 0;
        if (is.read(c))
        {
            // Delay further reading?
            delay = (c == token::BEGIN_BLOCK);
            is.putback(c);

            if (c)
            {
                cerr<< "nextChar:" << c << " : delay read: " << delay << nl;
            }
        }

        // Caller already checked token::compound::isCompound(...)
        // but use readCompoundToken anyhow for convenience

        if (tok.readCompoundToken(type, is, !delay))
        {
            cerr<< "readCompound(" << type << ")\n";
            cerr<< "typeCode: ";
            printTypeCode(cerr, tok.compoundToken().typeCode()) << nl;

            if (test_pending && delay)
            {
                InfoErr<< "pending read "
                    << tok.compoundToken().type() << endl;

                tok.refCompoundToken().pending(true);
            }

            return true;
        }

        return false;
    }


public:

    // Constructors
    using IFstream::IFstream;

    //- Destructor
    ~IFstreamDelayed() = default;

    // Testing deprecation warnings
    FOAM_DEPRECATED_STRICT(2023-08, "direct calling")
    Istream& operator()() const
    {
        return const_cast<IFstreamDelayed&>(*this);
    }
};


} // End namespace Foam


using namespace Foam;

void populateCompound(token::compound& ct, const dictionary& dict)
{
    Info<< "populateCompound: " << nl;

    // This is where runTime dispatch, eg based on transport type
    // could be used...

    switch (ct.typeCode())
    {
        #undef  fillComponents
        #define fillComponents(Type, Variable, Value)          \
        {                                                      \
            ct.pending(false);                                 \
            ct.resize(10);                                     \
            UList<Type> Variable                               \
            (                                                  \
                reinterpret_cast<Type*>(ct.data_bytes()),      \
                label(ct.size_bytes() / sizeof(Type))          \
            );                                                 \
            Variable = Value;                                  \
        }

        case token::tokenType::PUNCTUATION :
        {
            fillComponents(char, cmpts, '@');
        }
        break;

        case token::tokenType::BOOL :
        {
            fillComponents(bool, cmpts, false);
        }
        break;

        case token::tokenType::LABEL :
        {
            fillComponents(label, cmpts, 123);
        }
        break;

        case token::tokenType::FLOAT :
        {
            fillComponents(float, cmpts, 2.7);
        }
        break;

        case token::tokenType::DOUBLE :
        {
            fillComponents(double, cmpts, 3.1415);
        }
        break;

        default:
            break;

        #undef fillComponents
    }

    if (!ct.pending())
    {
        Info<< "assigned values:" << endl;
    }
}


void rewriteCompounds(ITstream& is)
{
    Info<< "rewrite: " << flatOutput(is) << endl;

    for (label toki = 0; toki < is.size(); ++toki)
    {
        if (is[toki].isCompound() && is[toki].compoundToken().pending())
        {
            Info<< "replace : " << is[toki].info() << endl;

            if (is.peekToken(toki+1).isPunctuation(token::BEGIN_BLOCK))
            {
                labelRange slice
                (
                    is.find(token::BEGIN_BLOCK, token::END_BLOCK, toki+1)
                );

                if (slice.good() && (slice.start() == toki+1))
                {
                    Info<< "Compound at:" << toki
                        << " dict:" << slice << endl;

                    ITstream substream(is.extract(slice));

                    dictionary dict(substream);

                    populateCompound(is[toki].refCompoundToken(), dict);
                }
            }
        }
    }
}


void rewriteDict(dictionary& dict)
{
    for (entry& e : dict)
    {
        if (e.isDict())
        {
            rewriteDict(e.dict());
        }
        else if (e.isStream())
        {
            rewriteCompounds(e.stream());
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noParallel();
    argList::addBoolOption("std", "standard reading (no delayed compounds)");
    argList::addBoolOption("pending", "read with pending");

    argList args(argc, argv, false, true);

    Info<< "typeCodes:" << nl;
    Info<< "    bool=";
    printTypeCode(Info, token::tokenType::BOOL) << nl;
    Info<< "    label=";
    printTypeCode(Info, token::tokenType::LABEL) << nl;
    Info<< "    float=";
    printTypeCode(Info, token::tokenType::FLOAT) << nl;
    Info<< "    double=";
    printTypeCode(Info, token::tokenType::DOUBLE) << nl;
    Info<< nl;

    if (args.found("pending"))
    {
        test_pending = true;
    }

    if (args.found("std"))
    {
        for (label argi = 1; argi < args.size(); ++argi)
        {
            Info<< "Read: " << args[argi] << endl;
            IFstream is(args[argi]);

            dictionary dict(is);

            Info<< "read: " << dict << nl;
        }
    }
    else
    {
        for (label argi = 1; argi < args.size(); ++argi)
        {
            Info<< "Read delay: " << args[argi] << endl;

            IFstreamDelayed is(args[argi]);

            // Trigger strict warning?
            Info<< "stream: " << is().name() << nl;

            dictionary dict(is);
            Info<< "read: " << dict << nl;

            rewriteDict(dict);

            Info<< "modified: " << dict << nl;
        }
    }

    return 0;
}


// ************************************************************************* //
