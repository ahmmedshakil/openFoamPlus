/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Released 2004-2011 OpenCFD Ltd.
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Modified code Copyright (C) 2017-2019 OpenCFD Ltd.
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
    Foam::keyType

Description
    A class for handling keywords in dictionaries.

    A keyType is the keyword of a dictionary.
    It differs from word in that it also accepts patterns (regular expressions).
    It is very similar to wordRe, but doesn't store a regular expression.

SourceFiles
    keyType.C

\*---------------------------------------------------------------------------*/

#ifndef keyType_H
#define keyType_H

#include "word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class Istream;
class Ostream;

/*---------------------------------------------------------------------------*\
                           Class keyType Declaration
\*---------------------------------------------------------------------------*/

class keyType
:
    public word
{
    // Private Data

        //- Treat keyType as a pattern (regular expression)
        bool isPattern_;


    // Private Member Functions

        //- No assignment where we cannot determine string/word type
        void operator=(const std::string&) = delete;


public:

    // Static Data Members

        //- An empty keyType
        static const keyType null;


    // Public Data Types

        //- Enumeration for search/match modes as bitmask
        //  eg, (keyType::REGEX | keyType::RECURSIVE)
        enum option
        {
            LITERAL = 0,        //!< String literal
            REGEX   = 1,        //!< Regular expression
            RECURSIVE = 0x10,   //!< Recursive search (eg, in dictionary)
            LITERAL_RECURSIVE = (LITERAL | RECURSIVE),
            REGEX_RECURSIVE   = (REGEX | RECURSIVE)
        };


    // Constructors

        //- Construct null
        inline keyType();

        //- Copy construct, retaining type (literal or regex)
        inline keyType(const keyType& s);

        //- Copy construct from word, treat as literal.
        inline keyType(const word& s);

        //- Copy construct from string, treat as regular expression.
        inline keyType(const string& s);

        //- Construct as copy of character array, treat as literal.
        inline keyType(const char* s);

        //- Copy construct from std::string with specified treatment
        inline keyType(const std::string& s, bool isPattern);

        //- Move construct, retaining type (literal or regex)
        inline keyType(keyType&& s);

        //- Move construct from word, treat as literal.
        inline keyType(word&& s);

        //- Move construct from string, treat as regular expression.
        inline keyType(string&& s);

        //- Move construct from std::string with specified treatment
        inline keyType(std::string&& s, bool isPattern);

        //- Construct from Istream
        //  Treat as regular expression if surrounded by quotation marks.
        explicit keyType(Istream& is);


    // Member Functions

        //- Is this character valid for a keyType?
        //  This is largely identical with what word accepts, but also
        //  permit brace-brackets, which are valid for some regexs.
        inline static bool valid(char c);


    // Access

        //- The keyType is treated as literal, not as pattern.
        inline bool isLiteral() const;

        //- The keyType is treated as a pattern, not as literal string.
        inline bool isPattern() const;


    // Infrastructure

        //- Mark as regular expression
        inline bool compile();

        //- Mark as literal, instead of a regular expression.
        //  Optionally strip invalid word characters.
        inline void uncompile(bool doStrip = false);


    // Editing

        //- Clear string and set as literal
        inline void clear();

        //- Swap contents. Self-swapping is a no-op.
        inline void swap(keyType& s);


    // Matching/Searching

        //- Smart match as regular expression or as a string.
        //  Optionally force a literal match only
        bool match(const std::string& text, bool literal=false) const;


    // Member Operators

        //- Perform smart match on text, as per match()
        //  Allows use as a predicate.
        inline bool operator()(const std::string& text) const;


        //- Copy assignment, retaining type (literal or regex)
        //  Self-assignment is a no-op.
        inline void operator=(const keyType& s);

        //- Move assignment, retaining type (literal or regex)
        //  Self-assignment is a no-op.
        inline void operator=(keyType&& s);

        //- Assign as word, treat as literal
        inline void operator=(const word& s);

        //- Assign from Foam::string, treat as regular expression
        inline void operator=(const string& s);

        //- Assign as word, treat as literal
        inline void operator=(const char* s);
};


// IOstream Operators

//- Read operator
Istream& operator>>(Istream& is, keyType& val);

//- Write operator
Ostream& operator<<(Ostream& os, const keyType& val);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "keyTypeI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
