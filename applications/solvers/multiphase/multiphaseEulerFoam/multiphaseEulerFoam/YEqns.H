for (int Ycorr=0; Ycorr<nSpecieCorrectors; Ycorr++)
{
    fluid.correctReactions();

    autoPtr<phaseSystem::specieTransferTable>
        specieTransferPtr(fluid.specieTransfer());

    phaseSystem::specieTransferTable&
        specieTransfer(specieTransferPtr());

    forAll(fluid.multiComponentPhases(), multiComponentPhasei)
    {
        phaseModel& phase = fluid.multiComponentPhases()[multiComponentPhasei];

        UPtrList<volScalarField>& Y = phase.YActiveRef();
        const volScalarField& alpha = phase;
        const volScalarField& rho = phase.rho();

        forAll(Y, i)
        {
            fvScalarMatrix YiEqn
            (
                phase.YiEqn(Y[i])
             ==
               *specieTransfer[Y[i].name()]
              + fvOptions(alpha, rho, Y[i])
            );

            YiEqn.relax();
            YiEqn.solve("Yi");
        }

        UPtrList<volScalarField>& Ys = phase.YsActiveRef();

        forAll(Ys, i)
        {
            fvScalarMatrix YsiEqn
            (
                phase.YsiEqn(Ys[i])
            );

            YsiEqn.relax();
            YsiEqn.solve("Ysi");
        }
    }

    fluid.correctSpecies();
}
