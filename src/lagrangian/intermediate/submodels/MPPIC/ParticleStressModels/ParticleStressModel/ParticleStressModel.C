/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
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

#include "ParticleStressModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ParticleStressModel, 0);
    defineRunTimeSelectionTable(ParticleStressModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ParticleStressModel::ParticleStressModel
(
    const dictionary& dict
)
:
    alphaPacked_(dict.get<scalar>("alphaPacked"))
{}


Foam::ParticleStressModel::ParticleStressModel
(
    const ParticleStressModel& cm
)
:
    alphaPacked_(cm.alphaPacked_)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::ParticleStressModel> Foam::ParticleStressModel::New
(
    const dictionary& dict
)
{
    const word modelType(dict.get<word>("type"));

    Info<< "Selecting particle stress model " << modelType << endl;

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(modelType);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown particle stress model type "
            << modelType << nl << nl
            << "Valid particle stress model types :" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << abort(FatalError);
    }

    return autoPtr<ParticleStressModel>(cstrIter()(dict));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ParticleStressModel::~ParticleStressModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::ParticleStressModel::alphaPacked() const
{
    return alphaPacked_;
}


Foam::tmp<Foam::FieldField<Foam::Field, Foam::scalar>>
Foam::ParticleStressModel::tau
(
    const FieldField<Field, scalar>& alpha,
    const FieldField<Field, scalar>& rho,
    const FieldField<Field, scalar>& uRms
) const
{
    tmp<FieldField<Field, scalar>> value
    (
        new FieldField<Field, scalar>(alpha.size())
    );

    forAll(alpha, i)
    {
        value->set(i, tau(alpha[i], rho[i], uRms[i]));
    }

    return value;
}


// ************************************************************************* //
