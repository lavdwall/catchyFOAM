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

InClass
    Foam::makeSurfaceReaction

Description
    Macros for instantiating reactions on given thermo packages

\*---------------------------------------------------------------------------*/

#ifndef makeSurfaceReaction_H
#define makeSurfaceReaction_H

#include "SurfaceReaction.H"

#include "IrreversibleSurfaceReaction.H"

#include "thermo.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define defineSurfaceReaction(nullArg, ThermoPhysics)                          \
    typedef SurfaceReaction<ThermoPhysics> SurfaceReaction##ThermoPhysics;     \
                                                                               \
    defineTemplateTypeNameAndDebug(SurfaceReaction##ThermoPhysics, 0);         \
                                                                               \
    defineTemplateRunTimeSelectionTable(SurfaceReaction##ThermoPhysics, dictionary) \
    defineTemplateRunTimeSelectionTable                                        \
    (                                                                          \
        SurfaceReaction##ThermoPhysics,                                        \
        objectRegistry                                                         \
    )


#define makeSurfaceReaction(ReactionType, ReactionRate, ThermoPhysics)         \
                                                                               \
    typedef ReactionType<ThermoPhysics, ReactionRate>                          \
        ReactionType##ThermoPhysics##ReactionRate;                             \
    typedef SurfaceReaction<ThermoPhysics>                                     \
        ReactionType##ThermoPhysics##ReactionRate##Base;                       \
                                                                               \
    template<>                                                                 \
    const word ReactionType##ThermoPhysics##ReactionRate::typeName             \
    (                                                                          \
        ReactionType::typeName_()                                              \
      + ReactionRate::type().capitalise()                                      \
    );                                                                         \
                                                                               \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        ReactionType##ThermoPhysics##ReactionRate##Base,                       \
        ReactionType##ThermoPhysics##ReactionRate,                             \
        dictionary                                                             \
    )

#define makeISurfaceReactions(ReactionRate, ThermoPhysics)                     \
                                                                               \
    makeSurfaceReaction(IrreversibleSurfaceReaction, ReactionRate, ThermoPhysics)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
