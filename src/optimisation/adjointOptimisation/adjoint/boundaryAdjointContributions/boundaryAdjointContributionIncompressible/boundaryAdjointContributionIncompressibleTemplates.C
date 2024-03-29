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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class returnType, class sourceType, class castType>
Foam::tmp<Foam::Field<returnType>>
Foam::boundaryAdjointContributionIncompressible::sumContributions
(
    PtrList<sourceType>& sourceList,
    const fvPatchField<returnType>&(castType::*boundaryFunction)(const label)
)
{
    // Objective function contribution
    tmp<Field<returnType>> tdJtotdvar
        (new Field<returnType>(patch_.size(), pTraits<returnType>::zero));
    Field<returnType>& dJtotdvar = tdJtotdvar.ref();

    // Get weighthed contribution
    for (sourceType& funcI : sourceList)
    {
        castType& cfuncI = refCast<castType>(funcI);
        const fvPatchField<returnType>& dJdvar =
            (cfuncI.*boundaryFunction)(patch_.index());
        dJtotdvar += cfuncI.weight()*dJdvar;
    }

    return tdJtotdvar;
}


// ************************************************************************* //
