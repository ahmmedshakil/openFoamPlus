/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Released 2004-2011 OpenCFD Ltd.
    Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::valveBank

Description
    A list of valves.

\*---------------------------------------------------------------------------*/

#ifndef valveBank_H
#define valveBank_H

#include "PtrList.H"
#include "engineValve.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class valveBank Declaration
\*---------------------------------------------------------------------------*/

class valveBank
:
    public PtrList<engineValve>
{
private:

    // Private Member Functions

        //- No copy construct
        valveBank(const valveBank&) = delete;

        //- No copy assignment
        void operator=(const valveBank&) = delete;


public:

    // Constructors

        //- Construct from Istream
        valveBank
        (
            const polyMesh& mesh,
            Istream& is
        )
        {
            PtrList<entry> valveEntries(is);

            setSize(valveEntries.size());

            forAll(valveEntries, valveI)
            {
                valveI,
                set
                (
                    new engineValve
                    (
                        valveEntries[valveI].keyword(),
                        mesh,
                        valveEntries[valveI].dict()
                    )
                );
            }
        }


    // Destructor - default


    // Member Functions
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
