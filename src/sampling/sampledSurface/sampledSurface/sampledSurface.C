/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Released 2004-2011 OpenCFD Ltd.
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Modified code Copyright (C) 2018-2019 OpenCFD Ltd.
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

#include "sampledSurface.H"
#include "polyMesh.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sampledSurface, 0);
    defineRunTimeSelectionTable(sampledSurface, word);
}


const Foam::wordList Foam::sampledSurface::surfaceFieldTypes
({
    "surfaceScalarField",
    "surfaceVectorField",
    "surfaceSphericalTensorField",
    "surfaceSymmTensorField",
    "surfaceTensorField"
});



// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::sampledSurface::clearGeom() const
{
    area_ = -1;
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::sampledSurface> Foam::sampledSurface::New
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
{
    const word sampleType(dict.get<word>("type"));

    DebugInfo
        << "Selecting sampledType " << sampleType << endl;

    auto cstrIter = wordConstructorTablePtr_->cfind(sampleType);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown sample type "
            << sampleType << nl << nl
            << "Valid sample types :" << endl
            << wordConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<sampledSurface>(cstrIter()(name, mesh, dict));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledSurface::sampledSurface(const word& name, std::nullptr_t)
:
    name_(name),
    mesh_(NullObjectRef<polyMesh>()),
    enabled_(true),
    interpolate_(false),
    area_(-1)
{}


Foam::sampledSurface::sampledSurface
(
    const word& name,
    const polyMesh& mesh,
    const bool interpolate
)
:
    name_(name),
    mesh_(mesh),
    enabled_(true),
    interpolate_(interpolate),
    area_(-1)
{}


Foam::sampledSurface::sampledSurface
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    name_(dict.lookupOrDefault<word>("name", name)),
    mesh_(mesh),
    enabled_(dict.lookupOrDefault("enabled", true)),
    interpolate_(dict.lookupOrDefault("interpolate", false)),
    area_(-1)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampledSurface::~sampledSurface()
{
    clearGeom();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::sampledSurface::area() const
{
    if (area_ < 0)
    {
        area_ = gSum(magSf());
    }

    return area_;
}


bool Foam::sampledSurface::withSurfaceFields() const
{
    return false;
}


Foam::tmp<Foam::scalarField> Foam::sampledSurface::sample
(
    const surfaceScalarField& sField
) const
{
    NotImplemented;
    return nullptr;
}


Foam::tmp<Foam::vectorField> Foam::sampledSurface::sample
(
    const surfaceVectorField& sField
) const
{
    NotImplemented;
    return nullptr;
}


Foam::tmp<Foam::sphericalTensorField> Foam::sampledSurface::sample
(
    const surfaceSphericalTensorField& sField
) const
{
    NotImplemented;
    return nullptr;
}


Foam::tmp<Foam::symmTensorField> Foam::sampledSurface::sample
(
    const surfaceSymmTensorField& sField
) const
{
    NotImplemented;
    return nullptr;
}


Foam::tmp<Foam::tensorField> Foam::sampledSurface::sample
(
    const surfaceTensorField& sField
) const
{
    NotImplemented;
    return nullptr;
}


void Foam::sampledSurface::print(Ostream& os) const
{
    os << type();
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const sampledSurface& s)
{
    s.print(os);
    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
