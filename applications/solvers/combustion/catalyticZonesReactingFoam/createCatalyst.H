Info << "Reading catalyst zone properties\n" << endl;
IOdictionary catDict
(
    IOobject
    (
        "catalyticZonesProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);
word catalystName(catDict.lookup("catalystName"));
List<word> catZones(catDict.lookup("catalyticCellZones"));

// Write required thermo fields for catalyst 'phase'
{
    forAll(Y,i)
    {
        volScalarField Ysi
        (
            IOobject
            (
                IOobject::groupName(IOobject::member(Y[i].name()), catalystName),
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            Y[i],
            calculatedFvPatchScalarField::typeName
        );
        Ysi.write();
    }

    volScalarField Tsi
    (
        IOobject
        (
            IOobject::groupName(IOobject::member(thermo.T().name()), catalystName),
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        thermo.T(),
        calculatedFvPatchScalarField::typeName
    );
    Tsi.write();
}

Info<< "Creating catalyst thermophysical properties\n" << endl;
autoPtr<rhoReactionThermo> catThermo(rhoReactionThermo::New(mesh, catalystName));
rhoReactionThermo& cThermo = catThermo();

autoPtr<basicGSChemistryModel> catReaction
(
    basicGSChemistryModel::New(cThermo)
);

forAll(mesh.C(), celli)
{
    cThermo.setCellReacting(celli, false);
    thermo.setCellReacting(celli, true);
}
forAll(catZones, czi)
{
    const cellZone& cz = mesh.cellZones()[catZones[czi]];
    forAll(cz, i)
    {
        cThermo.setCellReacting(cz[i], true);
        thermo.setCellReacting(cz[i], false);
    }
}

catReaction->disableSurfaceSpeciesEqn();
