/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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

#ifndef forGasSolids_H
#define forGasSolids_H

#include "specie.H"

#include "rhoConst.H"
#include "perfectGas.H"

#include "eConstThermo.H"
#include "hConstThermo.H"
#include "janafThermo.H"

#include "sensibleEnthalpy.H"
#include "sensibleInternalEnergy.H"

#include "constTransport.H"
#include "sutherlandTransport.H"
#include "logPolynomialTransport.H"

#include "thermo.H"

#include "forGSThermo.H"
#include "forThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define forGasSolidEquations(GMu, He, GCp, SCp, Macro, Args...)                \
    forGSThermo(constTransport, He, SCp, rhoConst, specie, GMu, He, GCp, perfectGas, specie, Macro, Args)

#define forGasSolidEnergiesAndThermos(Mu, Macro, Args...)                      \
    forGasSolidEquations(Mu, sensibleEnthalpy, hConstThermo, hConstThermo, Macro, Args);\
    forGasSolidEquations(Mu, sensibleEnthalpy, janafThermo, hConstThermo, Macro, Args);\
    forGasSolidEquations(Mu, sensibleInternalEnergy, eConstThermo, eConstThermo, Macro, Args);\
    forGasSolidEquations(Mu, sensibleInternalEnergy, janafThermo, eConstThermo, Macro, Args)

#define forGasSolidTransports(Macro, Args...)                                  \
    forGasSolidEnergiesAndThermos(sutherlandTransport, Macro, Args);           \
    forGasSolidEnergiesAndThermos(logPolynomialTransport, Macro, Args)

#define forGasSolids(Macro, Args...)                                           \
    forGasSolidTransports(Macro, Args)

#define forReactingSolids(Macro, Args...)                                      \
    forThermo(constTransport, sensibleEnthalpy, hConstThermo, rhoConst, specie, Macro, Args);\
    forThermo(constTransport, sensibleInternalEnergy, eConstThermo, rhoConst, specie, Macro, Args)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
