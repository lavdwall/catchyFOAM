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

Class
    Foam::multiComponentGSMixture

Description
    Foam::multiComponentGSMixture

SourceFiles
    multiComponentGSMixture.C

\*---------------------------------------------------------------------------*/

#ifndef multiComponentGSMixture_H
#define multiComponentGSMixture_H

#include "basicSpecieMixture.H"
#include "HashPtrTable.H"
#include "specieElement.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class multiComponentGSMixture Declaration
\*---------------------------------------------------------------------------*/

template<class SolidThermoType, class GasThermoType>
class multiComponentGSMixture
:
    public basicSpecieMixture
{
    // Private Data

        //- Solid species data
        PtrList<SolidThermoType> solidSpeciesData_;

        //- Gas species data
        PtrList<GasThermoType> speciesData_;

        typedef HashTable<List<specieElement>> speciesCompositionTable;

        //- Table of species composition
        speciesCompositionTable speciesComposition_;

        //- Temporary storage for the cell/face mixture thermo data
        mutable SolidThermoType mixture_;

        //- Temporary storage for the volume weighted
        //  cell/face mixture thermo data
        mutable SolidThermoType mixtureVol_;

        //- Temporary storage for the cell/face gas mixture thermo data
        mutable GasThermoType gasMixture_;

        //- Temporary storage for the volume weighted
        //  cell/face gas mixture thermo data
        mutable GasThermoType gasMixtureVol_;


    // Private Member Functions

        //- Read the species data from the given dictionary and return
        PtrList<GasThermoType> readSpeciesData(const dictionary& thermoDict) const;

        //- Read the solid species data from the given dictionary and return
        PtrList<SolidThermoType> readSolidSpeciesData(const dictionary& thermoDict) const;

        //- Read the species composition from the given dictionary and return
        speciesCompositionTable readSpeciesComposition
        (
            const dictionary& thermoDict,
            const speciesTable& speciesTable_
        ) const;

        //- Correct the mass fractions to sum to 1
        void correctMassFractions();


public:

    //- The type of thermodynamics this mixture is instantiated for
    typedef SolidThermoType thermoType;
    typedef GasThermoType gthermoType;


    // Constructors

        //- Construct from dictionary, specie names, thermo database,
        //  mesh and phase name
        multiComponentGSMixture
        (
            const dictionary&,
            const wordList& solidSpecieNames,
            const HashPtrTable<SolidThermoType>& solidThermoData,
            const wordList& specieNames,
            const HashPtrTable<GasThermoType>& thermoData,
            const fvMesh&,
            const word&
        );

        //- Construct from dictionary, mesh and phase name
        multiComponentGSMixture(const dictionary&, const fvMesh&, const word&);

        //- Disallow default bitwise copy construction
        multiComponentGSMixture
        (
            const multiComponentGSMixture<SolidThermoType, GasThermoType>&
        ) = delete;


    //- Destructor
    virtual ~multiComponentGSMixture()
    {}


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "multiComponentGSMixture<" + SolidThermoType::typeName()
                                        + ',' + GasThermoType::typeName() + '>';
        }

        //- Return the raw specie thermodynamic data
        const PtrList<GasThermoType>& specieThermos() const
        {
            return speciesData_;
        }

        //- Return the raw solid specie thermodynamic data
        const PtrList<SolidThermoType>& solidSpecieThermos() const
        {
            return solidSpeciesData_;
        }

        //- Return the table of species composition
        const speciesCompositionTable& specieComposition() const
        {
            return speciesComposition_;
        }

        const SolidThermoType& cellMixture(const label celli) const;

        const SolidThermoType& patchFaceMixture
        (
            const label patchi,
            const label facei
        ) const;

        const SolidThermoType& cellVolMixture
        (
            const scalar p,
            const scalar T,
            const label celli
        ) const;

        const SolidThermoType& patchFaceVolMixture
        (
            const scalar p,
            const scalar T,
            const label patchi,
            const label facei
        ) const;

        const GasThermoType& gasCellMixture(const label celli) const;

        const GasThermoType& gasPatchFaceMixture
        (
            const label patchi,
            const label facei
        ) const;

        const GasThermoType& gasCellVolMixture
        (
            const scalar p,
            const scalar T,
            const label celli
        ) const;

        const GasThermoType& gasPatchFaceVolMixture
        (
            const scalar p,
            const scalar T,
            const label patchi,
            const label facei
        ) const;

        //- Read dictionary
        void read(const dictionary&);

        //- Return thermo based on index
        inline const GasThermoType& specieThermo(const label speciei) const
        {
            return speciesData_[speciei];
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "multiComponentGSMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
