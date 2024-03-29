/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Modified code Copyright (C) 2013-2019 FOSS GP
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

\*---------------------------------------------------------------------------*/

#include "primalSolver.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(primalSolver, 0);
    defineRunTimeSelectionTable(primalSolver, primalSolver);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::primalSolver::primalSolver
(
    fvMesh& mesh,
    const word& managerType,
    const dictionary& dict
)
:
    solver(mesh, managerType, dict)
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::primalSolver> Foam::primalSolver::New
(
    fvMesh& mesh,
    const word& managerType,
    const dictionary& dict
)
{
    const word primalSolverType(dict.get<word>("type"));

    auto cstrIter = primalSolverConstructorTablePtr_->cfind(primalSolverType);

    if (!cstrIter.found())
    {
        FatalIOErrorInFunction(dict)
            << "Unknown primalSolver type " << primalSolverType
            << nl << nl
            << "Valid primalSolver types are :" << nl
            << primalSolverConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return autoPtr<primalSolver>(cstrIter()(mesh, managerType, dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::primalSolver::readDict(const dictionary& dict)
{
    if (solver::readDict(dict))
    {
        return true;
    }

    return false;
}


void Foam::primalSolver::correctBoundaryConditions()
{
    // Do nothing
}


// ************************************************************************* //
