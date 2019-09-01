/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2013-2016 OpenFOAM Foundation
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

#include "blendingFactor.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(blendingFactor, 0);
    addToRunTimeSelectionTable(functionObject, blendingFactor, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::blendingFactor::writeFileHeader(Ostream& os) const
{
    writeHeader(os, "Blending factor");
    writeCommented(os, "Time");
    writeTabbed(os, "Scheme1");
    writeTabbed(os, "Scheme2");
    writeTabbed(os, "Blended");
    os  << nl;
}


bool Foam::functionObjects::blendingFactor::calc()
{
    bool processed = false;

    processed = processed || calcScheme<scalar>();
    processed = processed || calcScheme<vector>();

    return processed;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::blendingFactor::blendingFactor
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict),
    writeFile(obr_, name, typeName, dict),
    phiName_("phi"),
    tolerance_(0.001)
{
    read(dict);
    writeFileHeader(file());
    setResultName(typeName, "");

    auto indicatorPtr = tmp<volScalarField>::New
    (
        IOobject
        (
            resultName_,
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero),
        zeroGradientFvPatchScalarField::typeName
    );

    store(resultName_, indicatorPtr);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::blendingFactor::read(const dictionary& dict)
{
    if (fieldExpression::read(dict) && writeFile::read(dict))
    {
        phiName_ = dict.getOrDefault<word>("phi", "phi");

        auto allowedRange = [&](const scalar tol){ return tol < 0 || 1 < tol; };

        tolerance_ =
            dict.getCheckOrDefault
            (
                "tolerance",
                0.001,
                allowedRange
            );

        return true;
    }

    return false;
}


bool Foam::functionObjects::blendingFactor::write()
{
    if (fieldExpression::write()) // Add if the scheme is blended
    {
        const volScalarField& indicator =
            lookupObject<volScalarField>(resultName_);

        // Generate scheme statistics
        label nCellsScheme1 = Zero;
        label nCellsScheme2 = Zero;
        label nCellsBlended = Zero;

        for (const auto& i : indicator)
        {
            if (i < tolerance_)
            {
                ++nCellsScheme1;
            }
            else if ((1 - tolerance_) < i)
            {
                ++nCellsScheme2;
            }
            else
            {
                ++nCellsBlended;
            }
        }

        reduce(nCellsScheme1, sumOp<label>());
        reduce(nCellsScheme2, sumOp<label>());
        reduce(nCellsBlended, sumOp<label>());

        Log << "    scheme 1 cells :  " << nCellsScheme1 << nl
            << "    scheme 2 cells :  " << nCellsScheme2 << nl
            << "    blended cells  :  " << nCellsBlended << nl
            << endl;

        writeTime(file());

        file()
            << tab << nCellsScheme1
            << tab << nCellsScheme2
            << tab << nCellsBlended
            << nl;
    }

    return true;
}


// ************************************************************************* //
