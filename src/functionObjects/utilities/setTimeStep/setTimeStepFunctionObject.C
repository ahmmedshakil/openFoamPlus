/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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

#include "setTimeStepFunctionObject.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(setTimeStepFunctionObject, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        setTimeStepFunctionObject,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::setTimeStepFunctionObject::setTimeStepFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObject(name),
    time_(runTime)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::setTimeStepFunctionObject::~setTimeStepFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::Time&
Foam::functionObjects::setTimeStepFunctionObject::time() const
{
    return time_;
}


bool Foam::functionObjects::setTimeStepFunctionObject::adjustTimeStep()
{
    // Wanted timestep
    scalar newDeltaT = timeStepPtr_().value(time_.timeOutputValue());

    static label index = -1;

    if (time_.timeIndex() != index)
    {
        // Store current time so we don't get infinite recursion (since
        // setDeltaT calls adjustTimeStep() again)
        index = time_.timeIndex();

        // Set time, allow deltaT to be adjusted for writeInterval purposes
        const_cast<Time&>(time_).setDeltaT(newDeltaT, true);
    }

    return true;
}


bool Foam::functionObjects::setTimeStepFunctionObject::read
(
    const dictionary& dict
)
{
    functionObject::read(dict);

    timeStepPtr_ = Function1<scalar>::New("deltaT", dict);

    // Check that adjustTimeStep is active
    const dictionary& controlDict = time_.controlDict();

    Switch adjust;
    if
    (
       !controlDict.readIfPresent<Switch>("adjustTimeStep", adjust)
    || !adjust
    )
    {
        FatalIOErrorInFunction(dict)
            << "Need to set 'adjustTimeStep' true to allow timestep control"
            << exit(FatalIOError);
    }

    return true;
}


bool Foam::functionObjects::setTimeStepFunctionObject::execute()
{
    return true;
}


bool Foam::functionObjects::setTimeStepFunctionObject::write()
{
    return true;
}


// ************************************************************************* //