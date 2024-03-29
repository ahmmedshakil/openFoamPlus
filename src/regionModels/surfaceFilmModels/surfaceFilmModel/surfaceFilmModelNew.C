/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Released 2009-2011 OpenCFD Ltd.
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

#include "surfaceFilmModel.H"
#include "noFilm.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<surfaceFilmModel> surfaceFilmModel::New
(
    const fvMesh& mesh,
    const dimensionedVector& g,
    const word& regionType
)
{
    word modelType;

    {
        IOobject surfaceFilmPropertiesDictHeader
        (
            regionType + "Properties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        );

        if (surfaceFilmPropertiesDictHeader.typeHeaderOk<IOdictionary>())
        {
            IOdictionary surfaceFilmPropertiesDict
            (
                surfaceFilmPropertiesDictHeader
            );

            surfaceFilmPropertiesDict.readEntry("surfaceFilmModel", modelType);
        }
        else
        {
            modelType = surfaceFilmModels::noFilm::typeName;
        }
    }

    Info<< "Selecting surfaceFilmModel " << modelType << endl;

    auto cstrIter = meshConstructorTablePtr_->cfind(modelType);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown surfaceFilmModel type "
            << modelType << nl << nl
            << "Valid surfaceFilmModel types :" << nl
            << meshConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<surfaceFilmModel>
    (
        cstrIter()
        (
            modelType,
            mesh,
            g,
            regionType
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
