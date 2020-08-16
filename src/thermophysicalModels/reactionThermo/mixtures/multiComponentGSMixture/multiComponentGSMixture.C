/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "multiComponentGSMixture.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class SolidThermoType, class GasThermoType>
Foam::PtrList<GasThermoType>
Foam::multiComponentGSMixture<SolidThermoType, GasThermoType>::readSpeciesData
(
    const dictionary& thermoDict
) const
{
    PtrList<GasThermoType> speciesData(species_.size());

    forAll(species_, i)
    {
        speciesData.set
        (
            i,
            new GasThermoType(thermoDict.subDict(species_[i]))
        );
    }

    return speciesData;
}


template<class SolidThermoType, class GasThermoType>
Foam::PtrList<SolidThermoType>
Foam::multiComponentGSMixture<SolidThermoType, GasThermoType>::readSolidSpeciesData
(
    const dictionary& thermoDict
) const
{
    PtrList<SolidThermoType> solidSpeciesData(solidSpecies_.size());

    forAll(solidSpecies_, i)
    {
        solidSpeciesData.set
        (
            i,
            new SolidThermoType(thermoDict.subDict(solidSpecies_[i]))
        );
    }

    return solidSpeciesData;
}


template<class SolidThermoType, class GasThermoType>
typename Foam::multiComponentGSMixture<SolidThermoType, GasThermoType>::speciesCompositionTable
Foam::multiComponentGSMixture<SolidThermoType, GasThermoType>::readSpeciesComposition
(
    const dictionary& thermoDict,
    const speciesTable& species
) const
{
    speciesCompositionTable speciesComposition_;

    // Loop through all species in thermoDict to retrieve
    // the species composition
    forAll(species, si)
    {
        if (thermoDict.subDict(species[si]).isDict("elements"))
        {
            dictionary currentElements
            (
                thermoDict.subDict(species[si]).subDict("elements")
            );

            wordList currentElementsName(currentElements.toc());
            List<specieElement> currentComposition(currentElementsName.size());

            forAll(currentElementsName, eni)
            {
                currentComposition[eni].name() = currentElementsName[eni];

                currentComposition[eni].nAtoms() =
                    currentElements.lookupOrDefault
                    (
                        currentElementsName[eni],
                        0
                    );
            }

            // Add current specie composition to the hash table
            speciesCompositionTable::iterator specieCompositionIter
            (
                speciesComposition_.find(species[si])
            );

            if (specieCompositionIter != speciesComposition_.end())
            {
                speciesComposition_.erase(specieCompositionIter);
            }

            speciesComposition_.insert(species[si], currentComposition);
        }
    }

    return speciesComposition_;
}


template<class SolidThermoType, class GasThermoType>
void Foam::multiComponentGSMixture<SolidThermoType, GasThermoType>::correctMassFractions()
{
    // Multiplication by 1.0 changes Yt patches to "calculated"
    volScalarField Yt("Yt", 1.0*Y_[0]);

    for (label n=1; n<Y_.size(); n++)
    {
        Yt += Y_[n];
    }

    if (mag(max(Yt).value()) < rootVSmall)
    {
        FatalErrorInFunction
            << "Sum of mass fractions is zero for species " << this->species()
            << exit(FatalError);
    }

    forAll(Y_, n)
    {
        Y_[n] /= Yt;
    }

    // Multiplication by 1.0 changes Yst patches to "calculated"
    volScalarField Yst("Yst", 1.0*Ys_[0]);

    for (label n=1; n<Ys_.size(); n++)
    {
        Yst += Ys_[n];
    }

    if (mag(max(Yst).value()) < rootVSmall)
    {
        FatalErrorInFunction
            << "Sum of coverages is zero for solid species " << this->solidSpecies()
            << exit(FatalError);
    }

    forAll(Ys_, n)
    {
        Ys_[n] /= Yst;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class SolidThermoType, class GasThermoType>
Foam::multiComponentGSMixture<SolidThermoType, GasThermoType>::multiComponentGSMixture
(
    const dictionary& thermoDict,
    const wordList& solidSpecieNames,
    const HashPtrTable<SolidThermoType>& solidThermoData,
    const wordList& specieNames,
    const HashPtrTable<GasThermoType>& thermoData,
    const fvMesh& mesh,
    const word& phaseName
)
:
    basicSpecieMixture
    (
        thermoDict,
        specieNames,
        solidSpecieNames,
        mesh,
        phaseName
    ),
    solidSpeciesData_(readSolidSpeciesData(thermoDict.subDict("solid"))),
    speciesData_(readSpeciesData(thermoDict.subDict("gas"))),
    speciesComposition_(readSpeciesComposition(thermoDict.subDict("gas"), species())),
    mixture_("mixture", solidSpeciesData_[0]),
    mixtureVol_("volMixture", solidSpeciesData_[0]),
    gasMixture_("gasMixture", speciesData_[0]),
    gasMixtureVol_("gasVolMixture", speciesData_[0])
{
    correctMassFractions();
}


template<class SolidThermoType, class GasThermoType>
Foam::multiComponentGSMixture<SolidThermoType, GasThermoType>::multiComponentGSMixture
(
    const dictionary& thermoDict,
    const fvMesh& mesh,
    const word& phaseName
)
:
    basicSpecieMixture
    (
        thermoDict,
        thermoDict.subDict("gas").lookup("species"),
        thermoDict.subDict("solid").lookup("species"),
        mesh,
        phaseName
    ),
    solidSpeciesData_(readSolidSpeciesData(thermoDict.subDict("solid"))),
    speciesData_(readSpeciesData(thermoDict.subDict("gas"))),
    speciesComposition_(readSpeciesComposition(thermoDict.subDict("gas"), species())),
    mixture_("mixture", solidSpeciesData_[0]),
    mixtureVol_("volMixture", solidSpeciesData_[0]),
    gasMixture_("gasMixture", speciesData_[0]),
    gasMixtureVol_("gasVolMixture", speciesData_[0])
{
    correctMassFractions();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class SolidThermoType, class GasThermoType>
const SolidThermoType& Foam::multiComponentGSMixture<SolidThermoType, GasThermoType>::cellMixture
(
    const label celli
) const
{
    return solidSpeciesData_[0];
}


template<class SolidThermoType, class GasThermoType>
const SolidThermoType& Foam::multiComponentGSMixture<SolidThermoType, GasThermoType>::patchFaceMixture
(
    const label patchi,
    const label facei
) const
{
    return solidSpeciesData_[0];
}


template<class SolidThermoType, class GasThermoType>
const SolidThermoType& Foam::multiComponentGSMixture<SolidThermoType, GasThermoType>::cellVolMixture
(
    const scalar p,
    const scalar T,
    const label celli
) const
{
    return solidSpeciesData_[0];
}


template<class SolidThermoType, class GasThermoType>
const SolidThermoType& Foam::multiComponentGSMixture<SolidThermoType, GasThermoType>::
patchFaceVolMixture
(
    const scalar p,
    const scalar T,
    const label patchi,
    const label facei
) const
{
    return solidSpeciesData_[0];
}


template<class SolidThermoType, class GasThermoType>
const GasThermoType& Foam::multiComponentGSMixture<SolidThermoType, GasThermoType>::gasCellMixture
(
    const label celli
) const
{
    gasMixture_ = Y_[0][celli]*speciesData_[0];

    for (label n=1; n<Y_.size(); n++)
    {
        gasMixture_ += Y_[n][celli]*speciesData_[n];
    }

    return gasMixture_;
}


template<class SolidThermoType, class GasThermoType>
const GasThermoType& Foam::multiComponentGSMixture<SolidThermoType, GasThermoType>::gasPatchFaceMixture
(
    const label patchi,
    const label facei
) const
{
    gasMixture_ = Y_[0].boundaryField()[patchi][facei]*speciesData_[0];

    for (label n=1; n<Y_.size(); n++)
    {
        gasMixture_ += Y_[n].boundaryField()[patchi][facei]*speciesData_[n];
    }

    return gasMixture_;
}


template<class SolidThermoType, class GasThermoType>
const GasThermoType& Foam::multiComponentGSMixture<SolidThermoType, GasThermoType>::gasCellVolMixture
(
    const scalar p,
    const scalar T,
    const label celli
) const
{
    scalar rhoInv = 0.0;
    forAll(speciesData_, i)
    {
        rhoInv += Y_[i][celli]/speciesData_[i].rho(p, T);
    }

    gasMixtureVol_ =
        Y_[0][celli]/speciesData_[0].rho(p, T)/rhoInv*speciesData_[0];

    for (label n=1; n<Y_.size(); n++)
    {
        gasMixtureVol_ +=
            Y_[n][celli]/speciesData_[n].rho(p, T)/rhoInv*speciesData_[n];
    }

    return gasMixtureVol_;
}


template<class SolidThermoType, class GasThermoType>
const GasThermoType& Foam::multiComponentGSMixture<SolidThermoType, GasThermoType>::
gasPatchFaceVolMixture
(
    const scalar p,
    const scalar T,
    const label patchi,
    const label facei
) const
{
    scalar rhoInv = 0.0;
    forAll(speciesData_, i)
    {
        rhoInv +=
            Y_[i].boundaryField()[patchi][facei]/speciesData_[i].rho(p, T);
    }

    gasMixtureVol_ =
        Y_[0].boundaryField()[patchi][facei]/speciesData_[0].rho(p, T)/rhoInv
      * speciesData_[0];

    for (label n=1; n<Y_.size(); n++)
    {
        gasMixtureVol_ +=
            Y_[n].boundaryField()[patchi][facei]/speciesData_[n].rho(p,T)
          / rhoInv*speciesData_[n];
    }

    return gasMixtureVol_;
}


template<class SolidThermoType, class GasThermoType>
void Foam::multiComponentGSMixture<SolidThermoType, GasThermoType>::read
(
    const dictionary& thermoDict
)
{
    forAll(species_, i)
    {
        speciesData_[i] = GasThermoType(thermoDict.subDict(species_[i]));
    }

    solidSpeciesData_[0] = SolidThermoType(thermoDict.subDict("solidMixture"));
}


// ************************************************************************* //
