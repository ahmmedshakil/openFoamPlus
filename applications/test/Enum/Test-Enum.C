/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenCFD Ltd.
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
    Testing of Enum lookups.

\*---------------------------------------------------------------------------*/

#include "Enum.H"
#include "dictionary.H"
#include "FlatOutput.H"
#include "IOstreams.H"  // For 'Sin'

#include <array>

using namespace Foam;

struct testing
{
    enum class option { A, B, C, D };

    static const Foam::Enum<option> option1Names;
    static const Foam::Enum<option> option2Names;
};

// All names
const Foam::Enum<testing::option> testing::option1Names
({
    { testing::option::A, "a" },
    { testing::option::B, "b" },
    { testing::option::C, "c" },
    { testing::option::D, "d" },
});

// Subset of names
const Foam::Enum<testing::option> testing::option2Names
({
    { testing::option::C, "c" },
    { testing::option::D, "d" },
});


// Can use for integers as well, but not scalar etc.
const Foam::Enum<int> otherNames1
({
    { 0, "a" },
    { 2, "b" },
    { 3, "c" },
    { 3, "d" },
});


// Can use for integers as well, but not scalar etc.
const Foam::Enum<int> otherNames2
({
    { 0, "a" },
    { 2, "b" },
    { 3, "c" },
    { 3, "asdasd" },
});


std::array<const char*, 2> myarray{ "false", "true" };


// Verify compile-time warnings
// #include "NamedEnum.H"
// const Foam::NamedEnum<testing::option, 2> bad_legacy;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    Info<<"Enum 1" << nl
        << "    names:  " << testing::option1Names << nl
        << "    values: " << flatOutput(testing::option1Names.values())
        << nl << nl;

    Info<<"Enum 2" << nl
        << "    names:  " << testing::option2Names << nl
        << "    values: " << flatOutput(testing::option2Names.values())
        << nl << nl;

    Info<<"Other Enum" << nl
        << "    names:  " << otherNames2 << nl
        << "    values: " << flatOutput(otherNames2.values())
        << nl << nl;


    dictionary testDict;
    testDict.add("lookup1", "c");
    testDict.add("lookup2", "rubbish");

    Info<< nl
        << int(testing::option1Names["a"]) << nl
        << testing::option1Names[testing::option::A] << nl;

    Info<< "--- test dictionary lookup ---" << endl;
    {
        Info<< "dict: " << testDict << endl;

        Info<< "lookupOrDefault(notFound) = "
            <<  int
                (
                    testing::option1Names.lookupOrDefault
                    (
                        "notFound",
                        testDict,
                        testing::option::A
                    )
                )
            << nl;

        Info<< "lookupOrDefault(lookup1) = "
            <<  int
                (
                    testing::option1Names.lookupOrDefault
                    (
                        "lookup1",
                        testDict,
                        testing::option::A
                    )
                )
            << nl;

        Info<< "lookupOrDefault(lookup1) = "
            <<  int
                (
                    testing::option2Names.lookupOrDefault
                    (
                        "lookup1",
                        testDict,
                        testing::option::A
                    )
                )
            << nl;
    }

    Info<< "--- test read ---" << endl;

    testing::option dummy(testing::option1Names.read(Sin));
    Info<< testing::option1Names[dummy] << endl;

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
