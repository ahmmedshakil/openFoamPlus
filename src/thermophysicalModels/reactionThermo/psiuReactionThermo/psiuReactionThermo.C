/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "psiuReactionThermo.H"
#include "fvMesh.H"
#include "zeroGradientFvPatchFields.H"
#include "fixedUnburntEnthalpyFvPatchScalarField.H"
#include "gradientUnburntEnthalpyFvPatchScalarField.H"
#include "mixedUnburntEnthalpyFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

defineTypeNameAndDebug(psiuReactionThermo, 0);
defineRunTimeSelectionTable(psiuReactionThermo, fvMesh);
defineRunTimeSelectionTable(psiuReactionThermo, fvMeshDictPhase);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

wordList psiuReactionThermo::heuBoundaryTypes()
{
    const volScalarField::Boundary& tbf =
        this->Tu().boundaryField();

    wordList hbt = tbf.types();

    forAll(tbf, patchi)
    {
        if (isA<fixedValueFvPatchScalarField>(tbf[patchi]))
        {
            hbt[patchi] = fixedUnburntEnthalpyFvPatchScalarField::typeName;
        }
        else if
        (
            isA<zeroGradientFvPatchScalarField>(tbf[patchi])
         || isA<fixedGradientFvPatchScalarField>(tbf[patchi])
        )
        {
            hbt[patchi] = gradientUnburntEnthalpyFvPatchScalarField::typeName;
        }
        else if (isA<mixedFvPatchScalarField>(tbf[patchi]))
        {
            hbt[patchi] = mixedUnburntEnthalpyFvPatchScalarField::typeName;
        }
    }

    return hbt;
}

void psiuReactionThermo::heuBoundaryCorrection(volScalarField& heu)
{
    volScalarField::Boundary& hbf = heu.boundaryFieldRef();

    forAll(hbf, patchi)
    {
        if
        (
            isA<gradientUnburntEnthalpyFvPatchScalarField>(hbf[patchi])
        )
        {
            refCast<gradientUnburntEnthalpyFvPatchScalarField>(hbf[patchi])
                .gradient() = hbf[patchi].fvPatchField::snGrad();
        }
        else if
        (
            isA<mixedUnburntEnthalpyFvPatchScalarField>(hbf[patchi])
        )
        {
            refCast<mixedUnburntEnthalpyFvPatchScalarField>(hbf[patchi])
                .refGrad() = hbf[patchi].fvPatchField::snGrad();
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

psiuReactionThermo::psiuReactionThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    psiReactionThermo(mesh, phaseName)
{}


psiuReactionThermo::psiuReactionThermo
(
    const fvMesh& mesh,
    const word& phaseName,
    const word& dictName
)
:
    psiReactionThermo(mesh, phaseName, dictName)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::psiuReactionThermo> Foam::psiuReactionThermo::New
(
    const fvMesh& mesh,
    const word& phaseName
)
{
    return basicThermo::New<psiuReactionThermo>(mesh, phaseName);
}


Foam::autoPtr<Foam::psiuReactionThermo> Foam::psiuReactionThermo::New
(
    const fvMesh& mesh,
    const word& phaseName,
    const word& dictName
)
{
    return basicThermo::New<psiuReactionThermo>(mesh, phaseName, dictName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

psiuReactionThermo::~psiuReactionThermo()
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
