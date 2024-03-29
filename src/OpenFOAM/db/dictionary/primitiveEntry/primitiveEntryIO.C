/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Released 2004-2011 OpenCFD Ltd.
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Modified code Copyright (C) 2017-2018 OpenCFD Ltd.
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

#include "primitiveEntry.H"
#include "functionEntry.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{
    // This is akin to a SafeIOWarning, which does not yet exist
    inline void safeIOWarning
    (
        const Foam::IOstream& is,
        const std::string& msg
    )
    {
        std::cerr
            << "--> FOAM Warning :\n"
            << "    Reading \"" << is.name() << "\" at line "
            << is.lineNumber() << '\n'
            << "    " << msg << std::endl;
    }

} // End anonymous namespace


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::primitiveEntry::acceptToken
(
    const token& tok,
    const dictionary& dict,
    Istream& is
)
{
    bool accept = tok.good();

    if (tok.isWord())
    {
        const word& key = tok.wordToken();

        accept =
        (
            disableFunctionEntries
         || key.size() == 1
         || (
                !(key[0] == '$' && expandVariable(key.substr(1), dict))
             && !(key[0] == '#' && expandFunction(key.substr(1), dict, is))
            )
        );
    }
    else if (tok.isVariable())
    {
        const string& key = tok.stringToken();

        accept =
        (
            disableFunctionEntries
         || key.size() <= 3
         || !(
                key[0] == '$'
             && key[1] == token::BEGIN_BLOCK
             && expandVariable(key.substr(1), dict)
            )
        );
    }

    return accept;
}


bool Foam::primitiveEntry::expandFunction
(
    const word& functionName,
    const dictionary& dict,
    Istream& is
)
{
    return functionEntry::execute(functionName, dict, *this, is);
}


bool Foam::primitiveEntry::read(const dictionary& dict, Istream& is)
{
    is.fatalCheck(FUNCTION_NAME);

    // Track balanced bracket/brace pairs, with max stack depth of 60.
    // Use a bitmask to track the opening char: 0 = '()', 1 = '{}'
    //
    // Notes
    // - the bitmask is set *before* increasing the depth since the left
    //   shift implicitly carries a 1-offset with it.
    //   Eg, (1u << 0) already corresponds to depth=1 (the first bit)
    //
    // - similarly, the bitmask is tested *after* decreasing depth

    uint64_t balanced = 0u;
    label depth = 0;
    token tok;

    while
    (
        !is.read(tok).bad() && tok.good()
     && !(tok == token::END_STATEMENT && depth == 0)
    )
    {
        if (tok.isPunctuation())
        {
            const char c = tok.pToken();
            switch (c)
            {
                case token::BEGIN_LIST:
                {
                    if (depth >= 0 && depth < 61)
                    {
                        balanced &= ~(1u << depth); // clear bit
                    }
                    ++depth;
                }
                break;

                case token::BEGIN_BLOCK:
                {
                    if (depth >= 0 && depth < 61)
                    {
                        balanced |= (1u << depth); // set bit
                    }
                    ++depth;
                }
                break;

                case token::END_LIST:
                {
                    --depth;
                    if (depth < 0)
                    {
                        safeIOWarning
                        (
                            is,
                            "Too many closing ')' ... was a ';' forgotten?"
                        );
                    }
                    else if (depth < 61 && ((balanced >> depth) & 1u))
                    {
                        // Bit was set, but expected it to be unset.
                        safeIOWarning(is, "Imbalanced '{' with ')'");
                    }
                }
                break;

                case token::END_BLOCK:
                {
                    --depth;
                    if (depth < 0)
                    {
                        safeIOWarning
                        (
                            is,
                            "Too many closing '}' ... was a ';' forgotten?"
                        );
                    }
                    else if (depth < 61 && !((balanced >> depth) & 1u))
                    {
                        // Bit was unset, but expected it to be set.
                        safeIOWarning(is, "Imbalanced '(' with '}'");
                    }
                }
                break;
            }
        }

        if (acceptToken(tok, dict, is))
        {
            newElmt(tokenIndex()++) = std::move(tok);
        }

        // With/without move: clear any old content and force to have a
        // known good token so that we can rely on it for the return value.

        tok = token::punctuationToken::NULL_TOKEN;
    }

    if (depth)
    {
        safeIOWarning(is, "Imbalanced brackets");
    }

    is.fatalCheck(FUNCTION_NAME);
    return tok.good();
}


void Foam::primitiveEntry::readEntry(const dictionary& dict, Istream& is)
{
    const label keywordLineNumber = is.lineNumber();
    tokenIndex() = 0;

    if (read(dict, is))
    {
        setSize(tokenIndex());
        tokenIndex() = 0;
    }
    else
    {
        std::ostringstream os;
        os  << "ill defined primitiveEntry starting at keyword '"
            << keyword() << '\''
            << " on line " << keywordLineNumber
            << " and ending at line " << is.lineNumber();

        SafeFatalIOErrorInFunction
        (
            is,
            os.str()
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::primitiveEntry::primitiveEntry
(
    const keyType& key,
    const dictionary& dict,
    Istream& is
)
:
    entry(key),
    ITstream
    (
        is.name() + '.' + key,
        tokenList(10),
        is.format(),
        is.version()
    )
{
    readEntry(dict, is);
}


Foam::primitiveEntry::primitiveEntry(const keyType& key, Istream& is)
:
    primitiveEntry(key, dictionary::null, is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::primitiveEntry::write(Ostream& os, const bool contentsOnly) const
{
    if (!contentsOnly)
    {
        os.writeKeyword(keyword());
    }

    bool addSpace = false;  // Separate from previous tokens with a space
    for (const token& tok : *this)
    {
        if (addSpace) os << token::SPACE;

        // Try to output token directly, with special handling in Ostreams.

        if (!os.write(tok))
        {
            os  << tok;   // Revert to normal '<<' output operator
        }

        addSpace = true;  // Separate from following tokens
    }

    if (!contentsOnly)
    {
        os  << token::END_STATEMENT << endl;
    }
}


void Foam::primitiveEntry::write(Ostream& os) const
{
    this->write(os, false);
}


// * * * * * * * * * * * * * Ostream operator  * * * * * * * * * * * * * * * //

template<>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const InfoProxy<primitiveEntry>& ip
)
{
    const primitiveEntry& e = ip.t_;

    e.print(os);

    const label nPrintTokens = 10;

    os  << "    primitiveEntry '" << e.keyword() << "' comprises ";

    for (label i=0; i<min(e.size(), nPrintTokens); ++i)
    {
        os  << nl << "        " << e[i].info();
    }

    if (e.size() > nPrintTokens)
    {
        os  << " ...";
    }

    os  << endl;

    return os;
}


// ************************************************************************* //
