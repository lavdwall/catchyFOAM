/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------


\*---------------------------------------------------------------------------*/

#include "catalyticWallSpeciesFluxFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "basicGSChemistryModel.H"
#include "rhoReactionThermophysicalTransportModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::catalyticWallSpeciesFluxFvPatchScalarField::
catalyticWallSpeciesFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    catalystName_("catalyst")
{}


Foam::catalyticWallSpeciesFluxFvPatchScalarField::
catalyticWallSpeciesFluxFvPatchScalarField
(
    const catalyticWallSpeciesFluxFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    catalystName_(ptf.catalystName_)
{}


Foam::catalyticWallSpeciesFluxFvPatchScalarField::
catalyticWallSpeciesFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF, dict),
    catalystName_(dict.lookupOrDefault<word>("catalystName", "catalyst"))
{}


Foam::catalyticWallSpeciesFluxFvPatchScalarField::
catalyticWallSpeciesFluxFvPatchScalarField
(
    const catalyticWallSpeciesFluxFvPatchScalarField& tppsf
)
:
    fixedGradientFvPatchScalarField(tppsf),
    catalystName_(tppsf.catalystName_)
{}


Foam::catalyticWallSpeciesFluxFvPatchScalarField::
catalyticWallSpeciesFluxFvPatchScalarField
(
    const catalyticWallSpeciesFluxFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(tppsf, iF),
    catalystName_(tppsf.catalystName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::catalyticWallSpeciesFluxFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const label patchi = patch().index();
    const word speciei = IOobject::member(internalField().name());

    const rhoReactionThermo& thermo = db().lookupObject<rhoReactionThermo>
    (
        IOobject::groupName(basicThermo::dictName, IOobject::group(internalField().name()))
    );
    const basicGSChemistryModel& reaction = db().lookupObject<basicGSChemistryModel>
    (
        IOobject::groupName("chemistryProperties",catalystName_)
    );
    const rhoReactionThermophysicalTransportModel& thermophysicalTransport = db().lookupObject<rhoReactionThermophysicalTransportModel>
    (
        IOobject::groupName("thermophysicalTransport", IOobject::group(internalField().name()))
    );

    label i = thermo.composition().species()[speciei];
    const scalarField Rgw = reaction.catalyticWallSpeciesFlux(speciei, patchi); // [kg/(m2.s)]
    const scalarField rhoDw = thermophysicalTransport.alphaEff(patchi);

    // rhoDw*gradient = Rgw
    // [kg/m3]*[m2/s]*[1/m] = [kg/(m2.s)]
    gradient() = Rgw/rhoDw; // [1/m]

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::catalyticWallSpeciesFluxFvPatchScalarField::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);
    writeEntry(os, "value", *this);

    os.writeKeyword("catalystName") << catalystName_
        << token::END_STATEMENT << nl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        catalyticWallSpeciesFluxFvPatchScalarField
    );
}

// ************************************************************************* //
