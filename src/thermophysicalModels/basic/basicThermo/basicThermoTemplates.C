/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2017 OpenFOAM Foundation
    Modified code Copyright (C) 2017 OpenCFD Ltd.
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

#include "basicThermo.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class Thermo, class Table>
typename Table::iterator Foam::basicThermo::lookupThermo
(
    const dictionary& thermoTypeDict,
    Table* tablePtr,
    std::initializer_list<const char*> cmptNames,
    const word& thermoTypeName
)
{
    // Lookup the thermo package
    auto cstrIter = tablePtr->find(thermoTypeName);

    // Print error message if package not found in the table
    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown " << Thermo::typeName << " type " << nl
            << "thermoType" << thermoTypeDict << nl << nl
            << "Valid " << Thermo::typeName << " types are:"
            << nl << nl;

        // Get the list of all the suitable thermo packages available
        wordList validThermoTypeNames
        (
            tablePtr->sortedToc()
        );

        // Build a table of the thermo packages constituent parts
        // Note: row-0 contains the names of constituent parts
        List<wordList> validThermoTypeNameCmpts
        (
            validThermoTypeNames.size() + 1
        );

        const int nCmpt = cmptNames.size();
        validThermoTypeNameCmpts[0].setSize(nCmpt);

        label j = 0;
        for (const char* cmptName : cmptNames)
        {
            validThermoTypeNameCmpts[0][j] = cmptName;
            ++j;
        }

        // Split the thermo package names into their constituent parts
        // Removing incompatible entries from the list
        j = 0;
        forAll(validThermoTypeNames, i)
        {
            wordList names
            (
                Thermo::splitThermoName(validThermoTypeNames[i], nCmpt)
            );

            if (names.size())
            {
                validThermoTypeNameCmpts[j++] = names;
            }
        }
        validThermoTypeNameCmpts.setSize(j);

        // Print the table of available packages
        // in terms of their constituent parts
        printTable(validThermoTypeNameCmpts, FatalError);

        FatalError<< exit(FatalError);
    }

    return cstrIter;
}


template<class Thermo, class Table>
typename Table::iterator Foam::basicThermo::lookupThermo
(
    const dictionary& thermoDict,
    Table* tablePtr
)
{
    if (thermoDict.isDict("thermoType"))
    {
        const dictionary& thermoTypeDict(thermoDict.subDict("thermoType"));

        Info<< "Selecting thermodynamics package " << thermoTypeDict << endl;

        if (thermoTypeDict.found("properties"))
        {
            std::initializer_list<const char*> cmptNames
            {
                "type",
                "mixture",
                "properties",
                "energy"
            };

            // Construct the name of the thermo package from the components
            const word thermoTypeName
            (
                thermoTypeDict.get<word>("type") + '<'
              + thermoTypeDict.get<word>("mixture") + '<'
              + thermoTypeDict.get<word>("properties") + ','
              + thermoTypeDict.get<word>("energy") + ">>"
            );

            return lookupThermo<Thermo, Table>
            (
                thermoTypeDict,
                tablePtr,
                cmptNames,
                thermoTypeName
            );
        }
        else
        {
            std::initializer_list<const char*> cmptNames
            {
                "type",
                "mixture",
                "transport",
                "thermo",
                "equationOfState",
                "specie",
                "energy"
            };

            // Construct the name of the thermo package from the components
            const word thermoTypeName
            (
                thermoTypeDict.get<word>("type") + '<'
              + thermoTypeDict.get<word>("mixture") + '<'
              + thermoTypeDict.get<word>("transport") + '<'
              + thermoTypeDict.get<word>("thermo") + '<'
              + thermoTypeDict.get<word>("equationOfState") + '<'
              + thermoTypeDict.get<word>("specie") + ">>,"
              + thermoTypeDict.get<word>("energy") + ">>>"
            );

            return lookupThermo<Thermo, Table>
            (
                thermoTypeDict,
                tablePtr,
                cmptNames,
                thermoTypeName
            );
        }
    }
    else
    {
        const word thermoTypeName(thermoDict.get<word>("thermoType"));

        Info<< "Selecting thermodynamics package " << thermoTypeName << endl;

        auto cstrIter = tablePtr->find(thermoTypeName);

        if (!cstrIter.found())
        {
            FatalErrorInFunction
                << "Unknown " << Thermo::typeName << " type "
                << thermoTypeName << nl << nl
                << "Valid " << Thermo::typeName << " types are:" << nl
                << tablePtr->sortedToc() << nl
                << exit(FatalError);
        }

        return cstrIter;
    }
}


template<class Thermo>
Foam::autoPtr<Thermo> Foam::basicThermo::New
(
    const fvMesh& mesh,
    const word& phaseName
)
{
    IOdictionary thermoDict
    (
        IOobject
        (
            phasePropertyName(dictName, phaseName),
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        )
    );

    auto cstrIter =
        lookupThermo<Thermo, typename Thermo::fvMeshConstructorTable>
        (
            thermoDict,
            Thermo::fvMeshConstructorTablePtr_
        );

    return autoPtr<Thermo>(cstrIter()(mesh, phaseName));
}


template<class Thermo>
Foam::autoPtr<Thermo> Foam::basicThermo::New
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
{
    auto cstrIter =
        lookupThermo<Thermo, typename Thermo::dictionaryConstructorTable>
        (
            dict,
            Thermo::dictionaryConstructorTablePtr_
        );

    return autoPtr<Thermo>(cstrIter()(mesh, dict, phaseName));
}


template<class Thermo>
Foam::autoPtr<Thermo> Foam::basicThermo::New
(
    const fvMesh& mesh,
    const word& phaseName,
    const word& dictName
)
{
    IOdictionary thermoDict
    (
        IOobject
        (
            dictName,
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        )
    );

    auto cstrIter =
        lookupThermo<Thermo, typename Thermo::fvMeshDictPhaseConstructorTable>
        (
            thermoDict,
            Thermo::fvMeshDictPhaseConstructorTablePtr_
        );

    return autoPtr<Thermo>(cstrIter()(mesh, phaseName, dictName));
}



// ************************************************************************* //
