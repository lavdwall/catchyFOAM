/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "SurfaceReactionProxy.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::SurfaceReactionProxy<ReactionThermo>::SurfaceReactionProxy
(
    const Reaction<ReactionThermo>& reaction,
    const speciesTable& gases,
    const List<specieCoeffs>& glhs,
    const List<specieCoeffs>& grhs
)
:
    SurfaceReaction<ReactionThermo>
    (
        reaction, gases, glhs, grhs
    )
{}


template<class ReactionThermo>
Foam::SurfaceReactionProxy<ReactionThermo>::SurfaceReactionProxy
(
    const SurfaceReactionProxy<ReactionThermo>& r,
    const speciesTable& gases
)
:
    SurfaceReaction<ReactionThermo>
    (
        r, gases
    )
{}


template<class ReactionThermo>
Foam::SurfaceReactionProxy<ReactionThermo>::SurfaceReactionProxy
(
    const speciesTable& species,
    const HashPtrTable<ReactionThermo>& thermoDatabase,
    const dictionary& dict
)
:
    SurfaceReaction<ReactionThermo>
    (
        species, thermoDatabase, dict
    )
{}


template<class ReactionThermo>
Foam::autoPtr<Foam::Reaction<ReactionThermo>>
Foam::SurfaceReactionProxy<ReactionThermo>::clone() const
{
    NotImplemented;
    return autoPtr<Reaction<ReactionThermo>>();
}


template<class ReactionThermo>
Foam::autoPtr<Foam::Reaction<ReactionThermo>>
Foam::SurfaceReactionProxy<ReactionThermo>::clone
(
    const speciesTable& species
) const
{
    NotImplemented;
    return autoPtr<Reaction<ReactionThermo>>();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



template<class ReactionThermo>
Foam::scalar Foam::SurfaceReactionProxy<ReactionThermo>::kf
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    NotImplemented;
    return 0;
}


template<class ReactionThermo>
Foam::scalar Foam::SurfaceReactionProxy<ReactionThermo>::kr
(
    const scalar kfwd,
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    NotImplemented;
    return 0;
}


template<class ReactionThermo>
Foam::scalar Foam::SurfaceReactionProxy<ReactionThermo>::kr
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    NotImplemented;
    return 0;
}


template<class ReactionThermo>
Foam::scalar Foam::SurfaceReactionProxy<ReactionThermo>::dkfdT
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    NotImplemented;
    return 0;
}


template<class ReactionThermo>
Foam::scalar Foam::SurfaceReactionProxy<ReactionThermo>::dkrdT
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    const scalar dkfdT,
    const scalar kr
) const
{
    NotImplemented;
    return 0;
}


template<class ReactionThermo>
const Foam::List<Foam::Tuple2<Foam::label, Foam::scalar>>&
Foam::SurfaceReactionProxy<ReactionThermo>::beta() const
{
    NotImplemented;
    return NullObjectRef<List<Tuple2<label, scalar>>>();
}


template<class ReactionThermo>
const Foam::List<Foam::Tuple2<Foam::label, Foam::scalarField>>&
Foam::SurfaceReactionProxy<ReactionThermo>::covdep() const
{
    NotImplemented;
    return NullObjectRef<List<Tuple2<label, scalarField>>>();
}


template<class ReactionThermo>
void Foam::SurfaceReactionProxy<ReactionThermo>::dcidc
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    scalarField& dcidc
) const
{
    NotImplemented;
}


template<class ReactionThermo>
Foam::scalar Foam::SurfaceReactionProxy<ReactionThermo>::dcidT
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    NotImplemented;
    return 0;
}


// ************************************************************************* //
