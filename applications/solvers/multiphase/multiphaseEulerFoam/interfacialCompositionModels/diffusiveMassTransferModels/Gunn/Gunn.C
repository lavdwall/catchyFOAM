/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2019 OpenFOAM Foundation
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

#include "Gunn.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diffusiveMassTransferModels
{
    defineTypeNameAndDebug(Gunn, 0);
    addToRunTimeSelectionTable
    (
        diffusiveMassTransferModel,
        Gunn,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diffusiveMassTransferModels::Gunn::Gunn
(
    const dictionary& dict,
    const phasePair& pair
)
:
    diffusiveMassTransferModel(dict, pair),
    Le_("Le", dimless, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diffusiveMassTransferModels::Gunn::~Gunn()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::diffusiveMassTransferModels::Gunn::K() const
{
    volScalarField Nu
    (
        (scalar(7.0) - 10.0*pair_.continuous() + 5.0*sqr(pair_.continuous()))
            *(scalar(1.0) + 0.7*pow(pair_.Re(), 0.2)*cbrt(pair_.Pr()))
      + (scalar(1.33) - 2.4*pair_.continuous() + 1.2*sqr(pair_.continuous()))
            *pow(pair_.Re(), 0.7)*cbrt(pair_.Pr())
    );

    volScalarField Sh(Nu*cbrt(Le_));

    return 6.0*pair_.dispersed()*Sh/sqr(pair_.dispersed().d()); 
}


// ************************************************************************* //
