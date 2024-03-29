/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Released 2010-2011 OpenCFD Ltd.
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "solidReaction.H"
#include "DynamicList.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::solidReaction<ReactionThermo>::solidReaction
(
    const Reaction<ReactionThermo>& reaction,
    const speciesTable& pyrolisisGases,
    const List<specieCoeffs>& glhs,
    const List<specieCoeffs>& grhs
)
:
    Reaction<ReactionThermo>(reaction),
    pyrolisisGases_(pyrolisisGases),
    glhs_(glhs),
    grhs_(grhs)
{}


template<class ReactionThermo>
Foam::solidReaction<ReactionThermo>::solidReaction
(
    const solidReaction<ReactionThermo>& r,
    const speciesTable& pyrolisisGases
)
:
    Reaction<ReactionThermo>(r),
    pyrolisisGases_(pyrolisisGases),
    glhs_(r.glhs_),
    grhs_(r.grhs_)
{}


template<class ReactionThermo>
Foam::solidReaction<ReactionThermo>::solidReaction
(
    const speciesTable& species,
    const HashPtrTable<ReactionThermo>& thermoDatabase,
    const dictionary& dict
)
:
    Reaction<ReactionThermo>(species, thermoDatabase, dict, false),
    pyrolisisGases_(dict.parent().parent().lookup("gaseousSpecies")),
    glhs_(),
    grhs_()
{
    this->setLRhs
    (
        IStringStream(dict.getString("reaction"))(),
        pyrolisisGases_,
        glhs_,
        grhs_
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ReactionThermo>
const Foam::List<typename Foam::solidReaction<ReactionThermo>::specieCoeffs>&
Foam::solidReaction<ReactionThermo>::glhs() const
{
    return glhs_;
}


template<class ReactionThermo>
const Foam::List<typename Foam::Reaction<ReactionThermo>::specieCoeffs>&
Foam::solidReaction<ReactionThermo>::grhs() const
{
    return grhs_;
}


template<class ReactionThermo>
const Foam::speciesTable& Foam::solidReaction<ReactionThermo>::
gasSpecies() const
{
    return pyrolisisGases_;
}


template<class ReactionThermo>
void Foam::solidReaction<ReactionThermo>::write(Ostream& os) const
{
    OStringStream reaction;
    os.writeEntry("reaction", solidReactionStr(reaction));
}


template<class ReactionThermo>
Foam::string Foam::solidReaction<ReactionThermo>::solidReactionStr
(
    OStringStream& reaction
) const
{
    this->reactionStrLeft(reaction);
    if (glhs().size() > 0)
    {
        reaction << " + ";
        solidReactionStrLeft(reaction);
    }
    reaction << " = ";
    this->reactionStrRight(reaction);
    if (grhs().size() > 0)
    {
        reaction << " + ";
        solidReactionStrRight(reaction);
    }
    return reaction.str();

}


template<class ReactionThermo>
void Foam::solidReaction<ReactionThermo>::solidReactionStrLeft
(
    OStringStream& reaction
) const
{
    for (label i = 0; i < glhs().size(); ++i)
    {
        if (i > 0)
        {
            reaction << " + ";
        }
        if (mag(glhs()[i].stoichCoeff - 1) > SMALL)
        {
            reaction << glhs()[i].stoichCoeff;
        }
        reaction << gasSpecies()[glhs()[i].index];
        if (mag(glhs()[i].exponent - glhs()[i].stoichCoeff) > SMALL)
        {
            reaction << "^" << glhs()[i].exponent;
        }
    }
}


template<class ReactionThermo>
void Foam::solidReaction<ReactionThermo>::solidReactionStrRight
(
    OStringStream& reaction
) const
{

    for (label i = 0; i < grhs().size(); ++i)
    {
        if (i > 0)
        {
            reaction << " + ";
        }
        if (mag(grhs()[i].stoichCoeff - 1) > SMALL)
        {
            reaction << grhs()[i].stoichCoeff;
        }
        reaction << gasSpecies()[grhs()[i].index];
        if (mag(grhs()[i].exponent - grhs()[i].stoichCoeff) > SMALL)
        {
            reaction << "^" << grhs()[i].exponent;
        }
    }
}

// ************************************************************************* //
