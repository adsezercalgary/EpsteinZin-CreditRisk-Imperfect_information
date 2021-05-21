# EpsteinZin-CreditRisk-Imperfect_information
This script generates the numerical results in the paper:  ``Credit risk pricing in a consumption based equilibrium framework with incomplete accounting information'',
by M. Ogunsolu, J. Ma, J. Qiu, A.D. Sezer.
The main script is main.m.

To execute the script make sure to have the following files in the same folder as main.m:
ccoefficent.m
BK.m
rhofunction1bc.m
rhofunction1ic.m
rhofunction1pde.m
trial.mat
trialX.mat

trial.mat and trialX.mat contain the scenarios for the consumption, volatility and the firm value processes.  If the parameters of consumption and volatility cahenge in the main.m, the data in trial.mat must be regenerated with the new parameters.  For this, use the script scenarios.m.  Make sure also to update the data in trialX.mat, if the consumption and volatility scenarios change or the firm parameters change.  This can be done using the script Xscenario.m.   

