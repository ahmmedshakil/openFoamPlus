/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Released 2004-2011 OpenCFD Ltd.
    Copyright (C) 2011 OpenFOAM Foundation
    Modified code Copyright (C) 2019 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "stringListOps.H"
#include "ListOps.H"
#include "FlatOutput.H"
#include "IOstreams.H"
#include "StringStream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    stringList strLst
    {
        "hello",
        "heello",
        "heeello",
        "bye",
        "bbye",
        "bbbye",
        "okey",
        "okkey",
        "okkkey",
    };

    wordRes reLst(IStringStream("( okey \"[hy]e+.*\" )")());

    Info<< "stringList " << strLst << nl;

    labelList matches = findStrings(regExp(".*ee.*"), strLst);

    Info<< "matches found for regexp .*ee.* :" << nl << matches << nl;

    forAll(matches, i)
    {
        Info<< " -> " << strLst[matches[i]] << nl;
    }

    Info<< "Match found using ListOps = "
        << ListOps::found(strLst, regExp(".*ee.*")) << nl;

    Info<< "First index = "
        << ListOps::find(strLst, regExp(".*ee.*")) << nl;

    Info<< endl;

    matches = findStrings(reLst, strLst);

    Info<< "matching " << flatOutput(reLst) << " => "
        << reLst.matching(strLst) << nl;
    Info<< "matches found for " << flatOutput(reLst) << " => "
        << matches << nl;
    forAll(matches, i)
    {
        Info<< " -> " << strLst[matches[i]] << nl;
    }
    Info<< endl;

    stringList subLst = subsetStrings(regExp(".*ee.*"), strLst);
    Info<< "subset stringList: " << subLst << nl;

    subLst = subsetStrings(reLst, strLst);
    Info<< "subset stringList: " << subLst << nl;

    inplaceSubsetStrings(reLst, strLst);
    Info<< "subsetted stringList: " << strLst << nl;

    inplaceSubsetStrings(regExp(".*l.*"), strLst);
    Info<< "subsetted stringList: " << strLst << nl;

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
