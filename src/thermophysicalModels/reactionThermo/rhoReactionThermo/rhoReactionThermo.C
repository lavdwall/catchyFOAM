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

#include "rhoReactionThermo.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(rhoReactionThermo, 0);
    defineRunTimeSelectionTable(rhoReactionThermo, fvMesh);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rhoReactionThermo::rhoReactionThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    rhoThermo(mesh, phaseName),
    reactingCells_(mesh.nCells(), true)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::rhoReactionThermo> Foam::rhoReactionThermo::New
(
    const fvMesh& mesh,
    const word& phaseName
)
{
    return basicThermo::New<rhoReactionThermo>(mesh, phaseName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::rhoReactionThermo::~rhoReactionThermo()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::rhoReactionThermo::setCellReacting
(
    const label celli,
    const bool active
)
{
    reactingCells_[celli] = active;
}

// ************************************************************************* //
