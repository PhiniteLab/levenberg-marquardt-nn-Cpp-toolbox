#ifndef __PHI_NN_LIB_H__
#define __PHI_NN_LIB_H__

#ifdef __cplusplus
extern "C"
{
#endif

    ////////////////////////////////////////////////////////////
    // including some libraries for using input/output functions

#include "stdio.h"
#include "string.h"
#include "stdlib.h"

    // including some libraries for using input/output functions
    ////////////////////////////////////////////////////////////

#ifdef __cplusplus
}
#endif

////////////////////////////////////////////////////////////
// including some libraries for using input/output functions

#include "phiMathParameters.hpp"

// including some libraries for using input/output functions
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
// CLASS DEFINITION

class phiNNParameters : public phiMathParameters
{
private:
protected:
public:
    //////////////////////////////////////////////////////
    //////////////////////////////////////////////////////
    // GENERAL VARIABLES

    ////////////////////////////////////////////////////////////
    // input output data structure to store every variables

    /////////////////////////////////////////////////
    // standard NN parameters

    int phiNNTrainingNumber; /*
                             * length of dataset in training process
                             */

    int phiNNNumberOfInputLayerNode; /*
                                     * number of input layer node
                                     */

    int phiNNNumberOfHiddenLayerNode; /*
                                     * number of hidden layer node
                                     */

    int phiNNNumberOfOutputLayerNode; /*
                                     * number of output layer node
                                     */

    int phiNNI, phiNNH, phiNNK; // for shorthining code

    double phiNNDevInvMatrices; /*
                             * Jacobian check condition
                             */

    double phiNNMu; /*
                   * learning rate for Levenberg Marquardt algorithm
                   */

    int phiNNNnTrainingCondition; /*
                                  * control variable to check NN training loop
                                  * it is used to terminate loop if the conditions are
                                  * satisfied!
                                  */

    FILE *phiNNNnOutputFileId; /*
                               * storage variable to keep the whole output of LM algorithm
                               */

    // standard NN parameters
    /////////////////////////////////////////////////

    /////////////////////////////////////////////////
    // neural network activation functions

    double **phiNNZActFunc; /*
                           * hidden activation function values
                           */

    double **phiNNYActFunc; /*
                           * output activation function values
                           */

    // neural network activation functions
    /////////////////////////////////////////////////

    /////////////////////////////////////////////////
    // neural network coefficients

    double **phiNNPhiWmatrices; /*
                               * input to hidden vector coefficients
                               */

    double **phiNNPhiVmatrices; /*
                               * hidden to output vector coefficients
                               */

    // neural network coefficients
    /////////////////////////////////////////////////

    /////////////////////////////////////////////////
    // vector parameters

    double **phiNNErrorNow; /*
                           * error vector for the present values
                           */

    double **phiNNErrorNowJac; /*
                           * error Jacobian vector for the present values
                           */

    double **phiNNErrorPre; /*
                           * error vector for the past values
                           */

    double phiNNErrorNowValue; /*
                              * the sum value of errorNow parameters
                              */

    double phiNNErrorPreValue; /*
                              * the sum value of errorPre parameters
                              */

    // vector parameters
    /////////////////////////////////////////////////

    /////////////////////////////////////////////////
    // NN training termination parameters

    int phiNNIterationMax; /*
                           * maximum iteration number of loop
                           */

    double phiNNEpsilonError; /*
                             * minimum value of error function
                             */

    double phiNNMuMin; /*
                      * minimum value of mu parameter
                      */

    double phiNNMaxDetInvMatrices; /*
                                  * maximum value of Jacobian determinant value
                                  */

    // NN training termination parameters
    /////////////////////////////////////////////////

    /////////////////////////////////////////////////
    // NN training Matrices

    double **phiNNJacobianTerm; /*
         * storing the total Jacobian matrices to update the coefficients
         */

    // NN training Matrices
    /////////////////////////////////////////////////

    // input output data structure to store every variables
    ////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////
    // CONSTRUCTOR

    phiNNParameters()
    {
        /////////////////////////////////////////////////
        // standard NN parameters

        phiNNTrainingNumber = 0;

        phiNNNumberOfInputLayerNode = 0;

        phiNNNumberOfHiddenLayerNode = 0;

        phiNNNumberOfOutputLayerNode = 0;

        phiNNI = 0;
        phiNNH = 0;
        phiNNK = 0;

        phiNNDevInvMatrices = 0;

        phiNNMu = 0;

        phiNNNnTrainingCondition = 0;

        phiNNNnOutputFileId = NULL;

        // standard NN parameters
        /////////////////////////////////////////////////

        /////////////////////////////////////////////////
        // neural network activation functions

        phiNNZActFunc = NULL;

        phiNNYActFunc = NULL;

        // neural network activation functions
        /////////////////////////////////////////////////

        /////////////////////////////////////////////////
        // neural network coefficients

        phiNNPhiWmatrices = NULL;

        phiNNPhiVmatrices = NULL;

        // neural network coefficients
        /////////////////////////////////////////////////

        /////////////////////////////////////////////////
        // vector parameters

        phiNNErrorNow = NULL;

        phiNNErrorNowJac = NULL;

        phiNNErrorPre = NULL;

        phiNNErrorNowValue = 0;

        phiNNErrorPreValue = 0;

        // vector parameters
        /////////////////////////////////////////////////

        /////////////////////////////////////////////////
        // NN training termination parameters

        phiNNIterationMax = 0;

        phiNNEpsilonError = 0;

        phiNNMuMin = 0;

        phiNNMaxDetInvMatrices = 0;

        // NN training termination parameters
        /////////////////////////////////////////////////

        /////////////////////////////////////////////////
        // NN training Matrices

        phiNNJacobianTerm = NULL;

        // NN training Matrices
        /////////////////////////////////////////////////
    }

    ~phiNNParameters()
    {
    }

    // CONSTRUCTOR
    ////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////

    // GENERAL VARIABLES
    //////////////////////////////////////////////////////
    //////////////////////////////////////////////////////

    //////////////////////////////////////////////////////
    //////////////////////////////////////////////////////
    // GENERAL FUNCTIONS

    void phiNNPhiInitializeNeuralNetwork(int inputLayerNode,
                                         int hiddenLayerNode,
                                         int outputLayerNode)
    {
        /*
         * initializing neural network for Levenberg marquadt structure
         * output -> returns nothing
         * input  -> value of input layer node
         *        -> value of hidden layer node
         *        -> value of output layer node
         */

        // set training number value
        this->phiNNTrainingNumber = this->phiInOutLengthNNDataset;

        // set the ptrNN parameters to use in other sections
        this->phiNNNumberOfHiddenLayerNode = hiddenLayerNode;
        this->phiNNNumberOfInputLayerNode = this->phiInOutNumberNNInputDataset;
        this->phiNNNumberOfOutputLayerNode = this->phiInOutNumberNNOutputDataset;

        // for summary operation NN library

        this->phiNNI = this->phiNNNumberOfInputLayerNode;
        this->phiNNH = this->phiNNNumberOfHiddenLayerNode;
        this->phiNNK = this->phiNNNumberOfOutputLayerNode;

        // initialize NN coefficients
        this->phiNNPhiWmatrices = this->phiNNPhiRandInitialization(-10, 10, this->phiNNH, this->phiNNI);
        this->phiNNPhiVmatrices = this->phiNNPhiRandInitialization(-10, 10, this->phiNNK, this->phiNNH + 1);
        printf("NN coefficients are created with random values\n");

        //writeCoefficientMatrices(ptrNN);

        // initialize error matrices
        this->phiNNErrorNow = this->phiNNPhiRandInitialization(0, 1, this->phiNNTrainingNumber, this->phiNNK);
        this->phiNNErrorPre = this->phiNNPhiRandInitialization(0, 1, this->phiNNTrainingNumber, this->phiNNK);
        this->phiNNErrorNowJac = this->phiNNPhiRandInitialization(0, 1, this->phiNNTrainingNumber, this->phiNNK);
        printf("Error vectors are created with random values\n");

        this->phiNNErrorNowValue = this->phiMathPhiMeanMatrices(this->phiNNErrorNow, this->phiNNTrainingNumber, this->phiNNK);
        this->phiNNErrorPreValue = this->phiMathPhiMeanMatrices(this->phiNNErrorPre, this->phiNNTrainingNumber, this->phiNNK);

        // creating NN stucture with proper activation functions
        this->phiNNCreatingActivationFunctions();

        // training condition is set to the start position!
        this->phiNNNnTrainingCondition = 1;

        this->phiNNMu = 0.08;

        // Jacobian update function is created
        this->phiNNJacobianTerm = this->phiInOutCreatingEmptyMatrices(this->phiNNTrainingNumber,
                                                                      (this->phiNNH * this->phiNNI + this->phiNNK * (this->phiNNH + 1)));

        printf("Activation functions are created with random values\n\n\n");

        printf("CONGRATULATIONS, NEURAL NETWORK STRUCTURE HAS BEEN ESTABLISHED!\n\n");

        this->phiNNWriteCoefficientMatrices();
    }

    void phiNNCreatingActivationFunctions()
    {
        /*
         * initializing NN activation functions
         * output -> returns nothing
         * input  -> taking nothing
         */

        this->phiNNZActFunc = this->phiNNPhiRandInitialization(0, 0, this->phiNNH + 1, this->phiNNTrainingNumber);
        this->phiNNYActFunc = this->phiNNPhiRandInitialization(0, 0, this->phiNNK, this->phiNNTrainingNumber);
    }

    double **phiNNPhiRandInitialization(double min, double max, int firstNodeNumber, int secondNodeNumber)
    {
        /*
         * initializing random NN coefficients
         * output -> returns address of W,V
         * input  -> value of minimum value
         *        -> value of maximum value
         *        -> value of first node number
         *        -> value of second node number
         */

        double **pd = this->phiInOutCreatingEmptyMatrices(firstNodeNumber, secondNodeNumber);

        for (int i = 0; i < firstNodeNumber; i++)
        {
            for (int j = 0; j < secondNodeNumber; j++)
            {
                pd[i][j] = (int)this->phiMathPhiRand(min, max);
                //pd[i][j] = phiRand(min, max);
            }
        }

        return pd;
    }

    void phiNNWriteCoefficientMatrices()
    {
        /*
         * write the coefficient matrices W,V
         * output -> returns nothings
         * input  -> taking nothing
         */

        printf("\n\n");
        printf("W input to hidden NN coefficients are written!\n");
        printf("\n");

        for (int i = 0; i < this->phiNNH; i++)
        {
            for (int j = 0; j < this->phiNNI; j++)
            {
                printf("%lf ", this->phiNNPhiWmatrices[i][j]);
            }
            printf("\n");
        }

        printf("\n");
        printf("V hidden to output NN coefficients are written!\n");
        printf("\n");

        for (int i = 0; i < this->phiNNK; i++)
        {
            for (int j = 0; j < this->phiNNH + 1; j++)
            {
                printf("%lf ", this->phiNNPhiVmatrices[i][j]);
            }
            printf("\n");
        }
        printf("\n\n");
    }

    void phiNNPhiTrainingParametersSettings(double epsilon,
                                            int iterationMax,
                                            double maxDetInvMatrices,
                                            double muMin)
    {
        /*
         * parameters settings of training Neural Network
         * output -> returns nothings
         * input  -> max error value
         *           max iteration step
         *           max Jacobian check
         *           max learning rate
         */

        this->phiNNEpsilonError = epsilon;
        this->phiNNIterationMax = iterationMax;
        this->phiNNMaxDetInvMatrices = maxDetInvMatrices;
        this->phiNNMuMin = muMin;

        printf("Neural Network Training Parameters:\n");
        printf("Minimum Error : %lf, Maximum Iteration : %d, Jacobian Check : %g, Minimum Mu : %lf",
               epsilon,
               iterationMax,
               maxDetInvMatrices,
               muMin);

        printf("\n\nThe whole parameters settings is done!\n\n");
    }

    void phiNNPhiTrainingNeuralNetworkWithLM()
    {
        /*
         * training Neural Network with Levenberg-Marquardt optimization technique
         * output -> returns nothings
         * input  -> taking nothing
         * */
        // storing iteration number
        int iteration = 0;

        // transpose storing Jacobian terms / indices are changed!
        double **JacobianTermTrans = this->phiInOutCreatingEmptyMatrices((this->phiNNH * this->phiNNI + this->phiNNK * (this->phiNNH + 1)),
                                                                         this->phiNNTrainingNumber);

        // multiplication of JacobianTerms
        double **JacobianMultiplicationResult = this->phiInOutCreatingEmptyMatrices(
            (this->phiNNH * this->phiNNI + this->phiNNK * (this->phiNNH + 1)),
            (this->phiNNH * this->phiNNI + this->phiNNK * (this->phiNNH + 1)));

        double **eyeMatrices = this->phiInOutCreatingEmptyMatrices(
            (this->phiNNH * this->phiNNI + this->phiNNK * (this->phiNNH + 1)),
            (this->phiNNH * this->phiNNI + this->phiNNK * (this->phiNNH + 1)));

        double **muEyeMatrices = this->phiInOutCreatingEmptyMatrices(
            (this->phiNNH * this->phiNNI + this->phiNNK * (this->phiNNH + 1)),
            (this->phiNNH * this->phiNNI + this->phiNNK * (this->phiNNH + 1)));

        double **JacobianFullTerm = this->phiInOutCreatingEmptyMatrices(
            (this->phiNNH * this->phiNNI + this->phiNNK * (this->phiNNH + 1)),
            (this->phiNNH * this->phiNNI + this->phiNNK * (this->phiNNH + 1)));

        double **InvJacobianFullTerm = this->phiInOutCreatingEmptyMatrices(
            (this->phiNNH * this->phiNNI + this->phiNNK * (this->phiNNH + 1)),
            (this->phiNNH * this->phiNNI + this->phiNNK * (this->phiNNH + 1)));

        double **coeffIntUpdate = this->phiInOutCreatingEmptyMatrices(this->phiNNH * this->phiNNI + this->phiNNK * (this->phiNNH + 1),
                                                                      this->phiNNTrainingNumber);

        double **coeffUpdate = this->phiInOutCreatingEmptyMatrices(
            (this->phiNNH * this->phiNNI + this->phiNNK * (this->phiNNH + 1)),
            1);

        while (this->phiNNNnTrainingCondition != 0)
        {
            // increasing iteration number
            iteration++;
            this->phiMathPhiMatrixAssignment(this->phiNNErrorPre, this->phiNNErrorNow, this->phiNNTrainingNumber, this->phiNNK);

            for (int i = 0; i < this->phiNNTrainingNumber; i++)
            {
                this->phiNNPhiActivationCalculation(i);

                this->phiNNPhiOutputFunctionCalculation(i);

                this->phiNNPhiLmAlgorithm(i);

                this->phiNNCostFunctionCalculation(i);

                this->phiNNCostFunctionCalculationJacobian(i);
            }

            // determinant calculation for devIntMatrices
            this->phiMathPhiTranspose(this->phiNNJacobianTerm, this->phiNNTrainingNumber,
                                      (this->phiNNH * this->phiNNI + this->phiNNK * (this->phiNNH + 1)),
                                      JacobianTermTrans);

            //printf("Phitranspose is done!\n\n");
            JacobianMultiplicationResult = this->phiMathPhiVectorMatrixMultiplication(JacobianTermTrans,
                                                                                      this->phiNNJacobianTerm,
                                                                                      (this->phiNNH * this->phiNNI + this->phiNNK * (this->phiNNH + 1)),
                                                                                      this->phiNNTrainingNumber,
                                                                                      this->phiNNTrainingNumber,
                                                                                      (this->phiNNH * this->phiNNI + this->phiNNK * (this->phiNNH + 1)));
            //printf("Jacobian multiplication is done!\n\n");

            eyeMatrices = this->phiMathPhiEyeCreationMatrices(this->phiNNH * this->phiNNI + this->phiNNK * (this->phiNNH + 1));

            //printf("Eye matrices is done!\n\n");

            muEyeMatrices = this->phiMathPhiSkalarMatrixMultiplication(this->phiNNMu, eyeMatrices,
                                                                       this->phiNNH * this->phiNNI + this->phiNNK * (this->phiNNH + 1),
                                                                       this->phiNNH * this->phiNNI + this->phiNNK * (this->phiNNH + 1));

            //printf("Mu Eye matrices is done!\n\n");

            JacobianFullTerm = this->phiMathPhiMatrixSummation(JacobianMultiplicationResult,
                                                               muEyeMatrices,
                                                               this->phiNNH * this->phiNNI + this->phiNNK * (this->phiNNH + 1),
                                                               this->phiNNH * this->phiNNI + this->phiNNK * (this->phiNNH + 1));

            //printf("Jacobian Full Term is done!\n\n");

            //printf("Inv matrices\n\n");
            InvJacobianFullTerm = this->phiMathPhiInverseMatrices(JacobianFullTerm,
                                                                  this->phiNNH * this->phiNNI + this->phiNNK * (this->phiNNH + 1),
                                                                  this->phiNNH * this->phiNNI + this->phiNNK * (this->phiNNH + 1));

            //printf("Internal update coeff is being calculated..\n");

            coeffIntUpdate = this->phiMathPhiVectorMatrixMultiplication(InvJacobianFullTerm,
                                                                        JacobianTermTrans,
                                                                        this->phiNNH * this->phiNNI + this->phiNNK * (this->phiNNH + 1),
                                                                        this->phiNNH * this->phiNNI + this->phiNNK * (this->phiNNH + 1),
                                                                        (this->phiNNH * this->phiNNI + this->phiNNK * (this->phiNNH + 1)),
                                                                        this->phiNNTrainingNumber);

            coeffUpdate = this->phiMathPhiVectorMatrixMultiplication(coeffIntUpdate, this->phiNNErrorNowJac,
                                                                     this->phiNNH * this->phiNNI + this->phiNNK * (this->phiNNH + 1),
                                                                     this->phiNNTrainingNumber,
                                                                     this->phiNNTrainingNumber, 1);

            // internal terms
            this->phiInOutPhiFree(coeffIntUpdate, this->phiNNH * this->phiNNI + this->phiNNK * (this->phiNNH + 1), this->phiNNTrainingNumber);
            //printf("Total update coeff is being calculated..\n");

            this->phiNNDevInvMatrices = this->phiMathPhiDeterminationCalculation(InvJacobianFullTerm,
                                                                                 this->phiNNH * this->phiNNI + this->phiNNK * (this->phiNNH + 1),
                                                                                 this->phiNNH * this->phiNNI + this->phiNNK * (this->phiNNH + 1));
            /*
        for (int i = 0; i < ptrNN->trainingNumber; i++)
        {
            for (int j = 0; j < ptrNN->K; j++)
            {
                printf("%lf ", ptrNN->errorNow[i][j]);
            }
            printf("\n");
        }*/

            /////////////////////////////////////////////////////////////
            // LM update to coefficients
            int counter = 0;

            for (int k = 0; k < this->phiNNK; k++)
            {
                for (int h = 0; h < this->phiNNH + 1; h++)
                {
                    this->phiNNPhiVmatrices[k][h] -= coeffUpdate[counter][0];
                    counter++;
                }
            }

            for (int h = 0; h < this->phiNNH; h++)
            {
                for (int in = 0; in < this->phiNNI; in++)
                {
                    this->phiNNPhiWmatrices[h][in] -= coeffUpdate[counter][0];

                    counter++;
                }
            }

            // LM update to coefficients
            /////////////////////////////////////////////////////////////

            this->phiNNErrorNowValue = this->phiMathPhiMeanMatrices(this->phiNNErrorNow, this->phiNNTrainingNumber, this->phiNNK);
            this->phiNNErrorPreValue = this->phiMathPhiMeanMatrices(this->phiNNErrorPre, this->phiNNTrainingNumber, this->phiNNK);

            if ((this->phiNNErrorPreValue - this->phiNNErrorNowValue) > 0)
            {
                double internalAssesmentMu = (this->phiNNErrorPreValue - this->phiNNErrorNowValue);

                if (fabs(internalAssesmentMu) > (1e-4 / (this->phiNNTrainingNumber)))
                {
                    this->phiNNMu = this->phiNNMu - this->phiNNMu * 0.01;
                }
                else
                {
                    this->phiNNMu = this->phiNNMu + this->phiNNMu * 0.01;
                }
            }

            if ((iteration % 1) == 0)
            {
                printf("\nIteration : %d, error: %e, Jac: %e, mu: %e\n",
                       iteration,
                       this->phiNNErrorNowValue,
                       this->phiNNDevInvMatrices,
                       this->phiNNMu);
                //writeCoefficientMatrices(ptrNN);
                /*
            for (int i = 0; i < 1; i++)
            {
                for (int j = 0; j < ptrNN->H * ptrNN->I + ptrNN->K * (ptrNN->H + 1); j++)
                {
                    printf("%e ", coeffUpdate[j][i]);
                }
                printf("\n");
            }*/
            }

            if ((this->phiNNErrorPreValue - this->phiNNErrorNowValue) >= 0)
            {
                double internalAssesmentMu = (this->phiNNErrorPreValue - this->phiNNErrorNowValue);

                if (fabs(internalAssesmentMu) > ((1e-1) / (this->phiNNTrainingNumber)))
                {
                    this->phiNNMu = this->phiNNMu + this->phiNNMu * 0.1;
                }
                else
                {
                    this->phiNNMu = this->phiNNMu - this->phiNNMu * 0.1;
                }
            }

            this->phiNNNnTrainingCondition = (iteration < this->phiNNIterationMax) &&
                                             (this->phiNNErrorNowValue > this->phiNNEpsilonError);

            double coeffControl = this->phiMathPhiSumMatrices(coeffUpdate, this->phiNNH * this->phiNNI + this->phiNNK * (this->phiNNH + 1), 1);

            if ((this->phiNNErrorPreValue - this->phiNNErrorNowValue) < (this->phiNNEpsilonError * this->phiNNEpsilonError) &&
                    (iteration > 1000) ||
                (_isnan(coeffControl) != 0) || (this->phiNNErrorNowValue > 500 * this->phiNNH))
            {
                this->phiInOutPhiFree(this->phiNNPhiWmatrices, this->phiNNH, this->phiNNI);
                this->phiInOutPhiFree(this->phiNNPhiVmatrices, this->phiNNK, this->phiNNH + 1);
                this->phiInOutPhiFree(this->phiNNErrorNow, this->phiNNTrainingNumber, this->phiNNK);
                this->phiInOutPhiFree(this->phiNNErrorPre, this->phiNNTrainingNumber, this->phiNNK);
                this->phiInOutPhiFree(this->phiNNErrorNowJac, this->phiNNTrainingNumber, this->phiNNK);

                this->phiInOutPhiFree(this->phiNNZActFunc, this->phiNNH + 1, this->phiNNTrainingNumber);
                this->phiInOutPhiFree(this->phiNNYActFunc, this->phiNNK, this->phiNNTrainingNumber);

                iteration = 0;

                this->phiNNPhiInitializeNeuralNetwork(this->phiNNI, this->phiNNH, this->phiNNK);
                printf("\n\nReinitialize!!!\n\n");
            }
        }
    }

    ////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////
    // training internal functions

    void phiNNPhiActivationCalculation(int index)
    {
        /*
         * calculation of hidden layer outputs
         * output -> returns nothings
         * input  -> index of activation function
         */

        // creating a loop to search the whole hidden node
        for (int h = 0; h < this->phiNNH + 1; h++)
        {
            // assign 1 if the hiddes node equals to bias term
            if (h == this->phiNNH)
            {
                this->phiNNZActFunc[h][index] = 1;
            }
            else
            {
                // calculatin input times coefficients
                double internalSum = 0.0;

                for (int i = 0; i < this->phiNNI; i++)
                {
                    // the whole inputs are summed with the terms of W matrices * inputMatrices
                    internalSum = internalSum + this->phiNNPhiWmatrices[h][i] * this->phiInOutInputDataArray[i][index];
                }
                // sigmoid activation function is used in here!
                this->phiNNZActFunc[h][index] = 1.0 / (1.0 + exp(-internalSum));
                /*
            if (internalSum > 20)
            {

                ptrNN->zActFunc[h][index] = 1.0;
            }
            else
            {
                if (internalSum < -20)
                {

                    ptrNN->zActFunc[h][index] = 0.0;
                }
                else
                {
                    ptrNN->zActFunc[h][index] = 1.0 / (1.0 + exp(-internalSum));
                }
            }*/
            }
        }
    }

    void phiNNPhiOutputFunctionCalculation(int index)
    {
        /*
         * calculation of output layer outputs
         * output -> returns nothings
         * input  -> index of output function
         */

        // search for the whole outputs
        for (int k = 0; k < this->phiNNK; k++)
        {
            double internalSum = 0.0;

            for (int h = 0; h < this->phiNNH + 1; h++)
            {
                // multiplying the whole terms with the hidden output values
                // linear neuron is used in here!
                internalSum = internalSum + this->phiNNPhiVmatrices[k][h] * this->phiNNZActFunc[h][index];
            }
            this->phiNNYActFunc[k][index] = internalSum;
            /*
        // assigned internalSum to the output function
        if (internalSum > 15)
        {
            ptrNN->yActFunc[k][index] = 15;
        }
        else
        {
            if (internalSum < -15)
            {
                ptrNN->yActFunc[k][index] = -15;
            }
            else
            {
                ptrNN->yActFunc[k][index] = internalSum;
            }
        }*/
        }
    }

    void phiNNPhiLmAlgorithm(int index)
    {
        /*
         * calculation total deflection in NN coefficient with Levenberg Marquardt
         * output -> returns nothings
         * input  -> index of Lm algorithm
         */

        int jacobianCounter = 0;

        // calculation of increment values in V matrices
        for (int k = 0; k < this->phiNNK; k++)
        {
            for (int h = 0; h < this->phiNNH + 1; h++)
            {
                // Jacobian terms is assigned to dV value since it is updated in the next session!
                this->phiNNJacobianTerm[index][jacobianCounter] = this->phiNNZActFunc[h][index] * (-1.0);

                jacobianCounter++;
            }
        }

        for (int k = 0; k < this->phiNNK; k++)
        {
            for (int h = 0; h < this->phiNNH; h++)
            {
                for (int in = 0; in < this->phiNNI; in++)
                {

                    this->phiNNJacobianTerm[index][jacobianCounter] = this->phiNNPhiVmatrices[k][h] *
                                                                      this->phiNNZActFunc[h][index] *
                                                                      (1.0 - this->phiNNZActFunc[h][index]) *
                                                                      this->phiInOutInputDataArray[in][index] *
                                                                      (-1.0);

                    jacobianCounter++;
                }
            }
        }
    }

    void phiNNCostFunctionCalculation(int index)
    {
        /*
         * calculation of total error
         * output -> returns nothings
         * input  -> index of cost function
         */

        // calculate the total error value on training/output dataset
        for (int k = 0; k < this->phiNNK; k++)
        {
            // 1/2*(y_tr - y)^2
            this->phiNNErrorNow[index][k] = fabs(1.0 / 2 *
                                                 (this->phiInOutOutputDataArray[k][index] - this->phiNNYActFunc[k][index]) *
                                                 (this->phiInOutOutputDataArray[k][index] - this->phiNNYActFunc[k][index]));
        }
    }

    void phiNNCostFunctionCalculationJacobian(int index)
    {
        /*
         * calculation of single error
         * output -> returns nothings
         * input  -> index of training process
         */

        // calculate the total error value on training/output dataset
        for (int k = 0; k < this->phiNNK; k++)
        {
            // 1/2*(y_tr - y)^2
            this->phiNNErrorNowJac[index][k] = (this->phiInOutOutputDataArray[k][index] - this->phiNNYActFunc[k][index]);
        }
    }

    // training internal functions
    ////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////

    void phiNNPhiExitNNLib()
    {
        /*
         * exiting the whole variables and functions
         * output -> returns nothing
         * input -> taking nothing
         */

        this->phiInOutPhiFree(this->phiNNPhiWmatrices, this->phiNNH, this->phiNNI);
        this->phiInOutPhiFree(this->phiNNPhiVmatrices, this->phiNNK, this->phiNNH + 1);
        this->phiInOutPhiFree(this->phiNNErrorNow, this->phiNNTrainingNumber, 1);
        this->phiInOutPhiFree(this->phiNNErrorPre, this->phiNNTrainingNumber, 1);
        this->phiInOutPhiFree(this->phiNNErrorNowJac, this->phiNNTrainingNumber, 1);
        this->phiInOutPhiFree(this->phiNNZActFunc, this->phiNNH + 1, this->phiNNTrainingNumber);
        this->phiInOutPhiFree(this->phiNNYActFunc, this->phiNNK, this->phiNNTrainingNumber);
        this->phiInOutPhiFree(this->phiNNJacobianTerm, this->phiNNTrainingNumber, (this->phiNNH * this->phiNNI + this->phiNNK * (this->phiNNH + 1)));
    }

    // type decleration for usage of NN Library
    ////////////////////////////////////////////////////////////

    // GENERAL FUNCTIONS
    //////////////////////////////////////////////////////
    //////////////////////////////////////////////////////
};

typedef phiNNParameters *phiNNParametersPtr;

// CLASS DEFINITION
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

#endif