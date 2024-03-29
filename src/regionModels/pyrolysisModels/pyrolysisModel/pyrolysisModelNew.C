/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Released 2009-2011 OpenCFD Ltd.
    Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "pyrolysisModel.H"
#include "fvMesh.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace pyrolysisModels
{

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<pyrolysisModel> pyrolysisModel::New
(
    const fvMesh& mesh,
    const word& regionType
)
{
    // get model name, but do not register the dictionary
    const word modelType
    (
        IOdictionary
        (
            IOobject
            (
                regionType + "Properties",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        ).get<word>("pyrolysisModel")
    );

    Info<< "Selecting pyrolysisModel " << modelType << endl;

    auto cstrIter = meshConstructorTablePtr_->cfind(modelType);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown pyrolysisModel type "
            << modelType << nl << nl
            << "Valid pyrolysisModel types :" << nl
            << meshConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<pyrolysisModel>(cstrIter()(modelType, mesh, regionType));
}


autoPtr<pyrolysisModel> pyrolysisModel::New
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& regionType
)
{

    const word modelType(dict.get<word>("pyrolysisModel"));

    Info<< "Selecting pyrolysisModel " << modelType << endl;

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(modelType);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown pyrolysisModel type "
            << modelType << nl << nl
            << "Valid pyrolysisModel types :" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<pyrolysisModel>
    (
        cstrIter()
        (
            modelType,
            mesh,
            dict,
            regionType
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
