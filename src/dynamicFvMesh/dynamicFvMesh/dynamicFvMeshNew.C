/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Released 2004-2011 OpenCFD Ltd.
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "staticFvMesh.H"
#include "simplifiedDynamicFvMesh.H"
#include "argList.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::dynamicFvMesh> Foam::dynamicFvMesh::New(const IOobject& io)
{
    // Note: - do not register the dictionary since dynamicFvMeshes themselves
    // do this.
    // - defaultRegion (region0) gets loaded from constant, other ones
    //   get loaded from constant/<regionname>. Normally we'd use
    //   polyMesh::dbDir() but we haven't got a polyMesh yet ...
    IOobject dictHeader
    (
        "dynamicMeshDict",
        io.time().constant(),
        (io.name() == polyMesh::defaultRegion ? "" : io.name()),
        io.db(),
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE,
        false
    );

    if (dictHeader.typeHeaderOk<IOdictionary>(true))
    {
        IOdictionary dict(dictHeader);

        const word modelType(dict.get<word>("dynamicFvMesh"));

        Info<< "Selecting dynamicFvMesh " << modelType << endl;

        const_cast<Time&>(io.time()).libs().open
        (
            dict,
            "dynamicFvMeshLibs",
            IOobjectConstructorTablePtr_
        );

        if (!IOobjectConstructorTablePtr_)
        {
            FatalErrorInFunction
                << "dynamicFvMesh table is empty"
                << exit(FatalError);
        }

        auto cstrIter = IOobjectConstructorTablePtr_->cfind(modelType);

        if (!cstrIter.found())
        {
            FatalIOErrorInFunction(dict)
                << "Unknown dynamicFvMesh type "
                << modelType << nl << nl
                << "Valid dynamicFvMesh types are :" << endl
                << IOobjectConstructorTablePtr_->sortedToc()
                << exit(FatalIOError);
        }

        return autoPtr<dynamicFvMesh>(cstrIter()(io));
    }


    return autoPtr<dynamicFvMesh>(new staticFvMesh(io));
}


Foam::autoPtr<Foam::dynamicFvMesh> Foam::dynamicFvMesh::New
(
    const argList& args,
    const Time& runTime
)
{
    if (args.found("dry-run") || args.found("dry-run-write"))
    {
        Info
            << "Operating in 'dry-run' mode: case will run for 1 time step.  "
            << "All checks assumed OK on a clean exit" << endl;

        FieldBase::allowConstructFromLargerSize = true;

        // Stop after 1 iteration of the simplified mesh

        if (args.found("dry-run-write"))
        {
            // Using saWriteNow triggers function objects execute(), write()
            runTime.stopAt(Foam::Time::saWriteNow);
        }
        else
        {
            // Using saNoWriteNow triggers function objects execute(),
            // but not write()
            runTime.stopAt(Foam::Time::saNoWriteNow);
        }

        functionObject::outputPrefix = "postProcessing-dry-run";

        return
            simplifiedMeshes::simplifiedDynamicFvMeshBase::New
            (
                IOobject
                (
                    dynamicFvMesh::defaultRegion,
                    runTime.timeName(),
                    runTime,
                    IOobject::MUST_READ
                )
            );
    }
    else
    {
        return
            New
            (
                IOobject
                (
                    dynamicFvMesh::defaultRegion,
                    runTime.timeName(),
                    runTime,
                    IOobject::MUST_READ
                )
            );
    }
}


// ************************************************************************* //
