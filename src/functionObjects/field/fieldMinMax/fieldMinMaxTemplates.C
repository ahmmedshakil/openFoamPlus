/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Released 2008-2011 OpenCFD Ltd.
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Modified code Copyright (C) 2015-2017 OpenCFD Ltd.
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

#include "fieldMinMax.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::functionObjects::fieldMinMax::output
(
    const word& fieldName,
    const word& outputName,
    const label minCell,
    const label maxCell,
    const vector& minC,
    const vector& maxC,
    const label minProci,
    const label maxProci,
    const Type& minValue,
    const Type& maxValue
)
{
    OFstream& file = this->file();

    if (location_)
    {
        writeTime(file());

        writeTabbed(file, fieldName);

        file<< token::TAB << minValue
            << token::TAB << minC;

        if (Pstream::parRun())
        {
            file<< token::TAB << minProci;
        }

        file<< token::TAB << maxValue
            << token::TAB << maxC;

        if (Pstream::parRun())
        {
            file<< token::TAB << maxProci;
        }

        file<< endl;

        Log << "    min(" << outputName << ") = " << minValue
            << " in cell " << minCell
            << " at location " << minC;

        if (Pstream::parRun())
        {
            Log << " on processor " << minProci;
        }

        Log << nl << "    max(" << outputName << ") = " << maxValue
            << " in cell " << maxCell
            << " at location " << maxC;

        if (Pstream::parRun())
        {
            Log << " on processor " << maxProci;
        }
    }
    else
    {
        file<< token::TAB << minValue << token::TAB << maxValue;

        Log << "    min/max(" << outputName << ") = "
            << minValue << ' ' << maxValue;
    }

    Log << endl;

    // Write state/results information
    word nameStr('(' + outputName + ')');
    this->setResult("min" + nameStr, minValue);
    this->setResult("min" + nameStr + "_cell", minCell);
    this->setResult("min" + nameStr + "_position", minC);
    this->setResult("min" + nameStr + "_processor", minProci);
    this->setResult("max" + nameStr, maxValue);
    this->setResult("max" + nameStr + "_cell", maxCell);
    this->setResult("max" + nameStr + "_position", maxC);
    this->setResult("max" + nameStr + "_processor", maxProci);
}


template<class Type>
void Foam::functionObjects::fieldMinMax::calcMinMaxFieldType
(
    const GeometricField<Type, fvPatchField, volMesh>& field,
    const word& outputFieldName
)
{
    const label proci = Pstream::myProcNo();

    // Find min internal field value info
    List<Type> minVs(Pstream::nProcs());
    labelList minCells(Pstream::nProcs());
    List<vector> minCs(Pstream::nProcs());

    label minProci = findMin(field);
    if (minProci != -1)
    {
        minVs[proci] = field[minProci];
        minCells[proci] = minProci;
        minCs[proci] = mesh_.C()[minProci];
    }
    else
    {
        minVs[proci] = pTraits<Type>::max;
        minCells[proci] = -1;
        minCs[proci] = vector::max;
    }

    // Find max internal field value info
    List<Type> maxVs(Pstream::nProcs());
    labelList maxCells(Pstream::nProcs());
    List<vector> maxCs(Pstream::nProcs());

    label maxProci = findMax(field);
    if (maxProci != -1)
    {
        maxVs[proci] = field[maxProci];
        maxCells[proci] = maxProci;
        maxCs[proci] = mesh_.C()[maxProci];
    }
    else
    {
        maxVs[proci] = pTraits<Type>::min;
        maxCells[proci] = -1;
        maxCs[proci] = vector::max;
    }

    // Find min and max boundary field info
    const auto& fieldBoundary = field.boundaryField();
    const auto& CfBoundary = mesh_.C().boundaryField();

    forAll(fieldBoundary, patchi)
    {
        const Field<Type>& fp = fieldBoundary[patchi];
        if (fp.size())
        {
            const vectorField& Cfp = CfBoundary[patchi];

            const labelList& faceCells =
                fieldBoundary[patchi].patch().faceCells();

            label minPi = findMin(fp);
            if (fp[minPi] < minVs[proci])
            {
                minVs[proci] = fp[minPi];
                minCells[proci] = faceCells[minPi];
                minCs[proci] = Cfp[minPi];
            }

            label maxPi = findMax(fp);
            if (fp[maxPi] > maxVs[proci])
            {
                maxVs[proci] = fp[maxPi];
                maxCells[proci] = faceCells[maxPi];
                maxCs[proci] = Cfp[maxPi];
            }
        }
    }

    // Collect info from all processors and output
    Pstream::gatherList(minVs);
    Pstream::scatterList(minVs);
    Pstream::gatherList(minCells);
    Pstream::scatterList(minCells);
    Pstream::gatherList(minCs);
    Pstream::scatterList(minCs);

    Pstream::gatherList(maxVs);
    Pstream::scatterList(maxVs);
    Pstream::gatherList(maxCells);
    Pstream::scatterList(maxCells);
    Pstream::gatherList(maxCs);
    Pstream::scatterList(maxCs);

    label mini = findMin(minVs);
    const Type& minValue = minVs[mini];
    const label minCell = minCells[mini];
    const vector& minC = minCs[mini];

    label maxi = findMax(maxVs);
    const Type& maxValue = maxVs[maxi];
    const label maxCell = maxCells[maxi];
    const vector& maxC = maxCs[maxi];

    output
    (
        field.name(),
        outputFieldName,
        minCell,
        maxCell,
        minC,
        maxC,
        mini,
        maxi,
        minValue,
        maxValue
    );
}


template<class Type>
void Foam::functionObjects::fieldMinMax::calcMinMaxFields
(
    const word& fieldName,
    const modeType& mode
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    if (obr_.foundObject<fieldType>(fieldName))
    {
        const fieldType& field = lookupObject<fieldType>(fieldName);

        switch (mode)
        {
            case mdMag:
            {
                calcMinMaxFieldType<scalar>
                (
                    mag(field),
                    word("mag(" + fieldName + ")")
                );
                break;
            }
            case mdCmpt:
            {
                calcMinMaxFieldType(field, fieldName);
                break;
            }
            default:
            {
                FatalErrorInFunction
                    << "Unknown min/max mode: " << modeTypeNames_[mode_]
                    << exit(FatalError);
            }
        }
    }
}


// ************************************************************************* //
