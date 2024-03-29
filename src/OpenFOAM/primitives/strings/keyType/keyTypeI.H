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

\*---------------------------------------------------------------------------*/

#include <algorithm>

// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

inline bool Foam::keyType::valid(char c)
{
    return
    (
        !isspace(c)
     && c != '"'   // string quote
     && c != '\''  // string quote
     && c != '/'   // path separator
     && c != ';'   // end statement
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::keyType::keyType()
:
    word(),
    isPattern_(false)
{}


inline Foam::keyType::keyType(const keyType& s)
:
    word(s, false),
    isPattern_(s.isPattern())
{}


inline Foam::keyType::keyType(const word& s)
:
    word(s, false),
    isPattern_(false)
{}


inline Foam::keyType::keyType(const string& s)
:
    word(s, false),
    isPattern_(true)
{}


inline Foam::keyType::keyType(const char* s)
:
    word(s, false),
    isPattern_(false)
{}


inline Foam::keyType::keyType(const std::string& s, bool isPattern)
:
    word(s, false),
    isPattern_(isPattern)
{}


inline Foam::keyType::keyType(keyType&& s)
:
    word(std::move(static_cast<word&>(s)), false),
    isPattern_(s.isPattern())
{
    s.isPattern_ = false;
}


inline Foam::keyType::keyType(word&& s)
:
    word(std::move(s), false),
    isPattern_(false)
{}


inline Foam::keyType::keyType(string&& s)
:
    word(std::move(s), false),
    isPattern_(true)
{}


inline Foam::keyType::keyType(std::string&& s, bool isPattern)
:
    word(std::move(s), false),
    isPattern_(isPattern)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::keyType::isLiteral() const
{
    return !isPattern_;
}


inline bool Foam::keyType::isPattern() const
{
    return isPattern_;
}


inline bool Foam::keyType::compile()
{
    isPattern_ = true;
    return true;
}


inline void Foam::keyType::uncompile(bool doStrip)
{
    // Only strip when debug is active (potentially costly operation)
    if (isPattern_ && doStrip && word::debug)
    {
        string::stripInvalid<word>(*this);
    }

    isPattern_ = false;
}


inline void Foam::keyType::clear()
{
    word::clear();
    isPattern_ = false;
}


inline void Foam::keyType::swap(keyType& s)
{
    // Self-swapping is a no-op
    if (this == &s)
    {
        return;
    }

    word::swap(static_cast<word&>(s));
    std::swap(isPattern_, s.isPattern_);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::keyType::operator()(const std::string& text) const
{
    return match(text);  // Use smart match
}


inline void Foam::keyType::operator=(const keyType& s)
{
    // Self-assignment is a no-op
    if (this == &s)
    {
        return;
    }

    assign(s); // Bypasses char checking
    isPattern_ = s.isPattern_;
}


inline void Foam::keyType::operator=(keyType&& s)
{
    // Self-assignment is a no-op
    if (this == &s)
    {
        return;
    }

    clear();
    swap(s);
}


inline void Foam::keyType::operator=(const word& s)
{
    assign(s); // Bypasses char checking
    isPattern_ = false;
}


inline void Foam::keyType::operator=(const string& s)
{
    assign(s); // Bypasses char checking
    isPattern_ = true;
}


inline void Foam::keyType::operator=(const char* s)
{
    assign(s); // Bypasses char checking
    isPattern_ = false;
}


// ************************************************************************* //
