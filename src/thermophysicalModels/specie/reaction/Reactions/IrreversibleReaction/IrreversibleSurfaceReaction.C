/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "IrreversibleSurfaceReaction.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo, class ReactionRate>
Foam::IrreversibleSurfaceReaction<ReactionThermo, ReactionRate>::
IrreversibleSurfaceReaction
(
    const SurfaceReaction<ReactionThermo>& reaction,
    const ReactionRate& k
)
:
    SurfaceReaction<ReactionThermo>(reaction),
    k_(k)
{}


template<class ReactionThermo, class ReactionRate>
Foam::IrreversibleSurfaceReaction<ReactionThermo, ReactionRate>::
IrreversibleSurfaceReaction
(
    const speciesTable& species,
    const HashPtrTable<ReactionThermo>& thermoDatabase,
    const dictionary& dict
)
:
    SurfaceReaction<ReactionThermo>(species, thermoDatabase, dict),
    k_(species, dict)
{}


template<class ReactionThermo, class ReactionRate>
Foam::IrreversibleSurfaceReaction<ReactionThermo, ReactionRate>::
IrreversibleSurfaceReaction
(
    const speciesTable& species,
    const HashPtrTable<ReactionThermo>& thermoDatabase,
    const objectRegistry& ob,
    const dictionary& dict
)
:
    SurfaceReaction<ReactionThermo>(species, thermoDatabase, dict),
    k_(species, ob, dict)
{}


template<class ReactionThermo, class ReactionRate>
Foam::IrreversibleSurfaceReaction<ReactionThermo, ReactionRate>::
IrreversibleSurfaceReaction
(
    const IrreversibleSurfaceReaction<ReactionThermo,ReactionRate>& irr,
    const speciesTable& species
)
:
    SurfaceReaction<ReactionThermo>(irr, species),
    k_(irr.k_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ReactionThermo, class ReactionRate>
Foam::scalar Foam::IrreversibleSurfaceReaction<ReactionThermo, ReactionRate>::kf
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return k_(p, T, c, li);
}


template<class ReactionThermo, class ReactionRate>
Foam::scalar Foam::IrreversibleSurfaceReaction<ReactionThermo, ReactionRate>::kr
(
    const scalar kfwd,
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return 0;
}


template<class ReactionThermo, class ReactionRate>
Foam::scalar Foam::IrreversibleSurfaceReaction<ReactionThermo, ReactionRate>::kr
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return 0;
}


template<class ReactionThermo, class ReactionRate>
Foam::scalar Foam::IrreversibleSurfaceReaction<ReactionThermo, ReactionRate>::dkfdT
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return k_.ddT(p, T, c, li);
}


template<class ReactionThermo, class ReactionRate>
Foam::scalar Foam::IrreversibleSurfaceReaction<ReactionThermo, ReactionRate>::dkrdT
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    const scalar dkfdT,
    const scalar kr
) const
{
    return 0;
}


template<class ReactionThermo, class ReactionRate>
const Foam::List<Foam::Tuple2<Foam::label, Foam::scalarField>>&
Foam::IrreversibleSurfaceReaction<ReactionThermo, ReactionRate>::covdep() const
{
    return k_.covdep();
}


template<class ReactionThermo, class ReactionRate>
const Foam::List<Foam::Tuple2<Foam::label, Foam::scalar>>&
Foam::IrreversibleSurfaceReaction<ReactionThermo, ReactionRate>::beta() const
{
    return k_.beta();
}


template<class ReactionThermo, class ReactionRate>
void Foam::IrreversibleSurfaceReaction<ReactionThermo, ReactionRate>::dcidc
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    scalarField& dcidc
) const
{
    k_.dcidc(p, T, c, li, dcidc);
}


template<class ReactionThermo, class ReactionRate>
Foam::scalar Foam::IrreversibleSurfaceReaction<ReactionThermo, ReactionRate>::dcidT
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return k_.dcidT(p, T, c, li);
}


template<class ReactionThermo, class ReactionRate>
void Foam::IrreversibleSurfaceReaction<ReactionThermo, ReactionRate>::write
(
    Ostream& os
) const
{
    SurfaceReaction<ReactionThermo>::write(os);
    k_.write(os);
}


// ************************************************************************* //
