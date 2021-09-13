#include "..\..\include\phiNNLibSettings.hpp"

using namespace std;

int main()
{
    phiNNParameters pNN;

    pNN.phiInOutFileName = "nnInputOutputFile.txt";

    pNN.phiInOutReadFileFromText();

    pNN.phiInOutWriteDataSetMatrices();

    /////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////
    // NN lib is started!

    printf("\nPress enter to start phiLM!\n");
    getchar();

    // randomize the whole process!!
    srand((unsigned int)time(NULL));

    pNN.phiNNPhiInitializeNeuralNetwork(3, 4, 1);

    pNN.phiNNPhiTrainingParametersSettings(1e-9, 100000, 1e200, 1e-12);

    pNN.phiNNPhiTrainingNeuralNetworkWithLM();

    printf("%lf ", pNN.phiNNErrorNowValue);

    pNN.phiNNWriteCoefficientMatrices();

    printf("Press enter to exit!\n");
    getchar();

    // NN lib is started!
    /////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////
    // exit functions

    pNN.phiInOutPhiExitInOut();

    pNN.phiNNPhiExitNNLib();

    // exit functions
    ////////////////////////////////////////////////////////

    return 0;
}