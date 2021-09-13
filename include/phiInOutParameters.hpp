#ifndef __PHI_IN_OUT_PARAMETERS_HPP__
#define __PHI_IN_OUT_PARAMETERS_HPP__

#ifdef __cplusplus
extern "C"
{
#endif

    ////////////////////////////////////////////////////////////
    // including some libraries for using input/output functions

// standard library usage
#include "stdio.h"
#include "string.h"
#include "stdlib.h"
#include "time.h"
#include "math.h"

    // including some libraries for using input/output functions
    ////////////////////////////////////////////////////////////

#ifdef __cplusplus
}
#endif

////////////////////////////////////////////////////////////
// including some libraries for using input/output functions

// including some libraries for using input/output functions
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// define error types

#define ALLOCATION_ERROR 1
#define INCONSISTENT_ROW_COLUMN 2
#define FILE_OPEN_ERROR 3
#define SAMPLING_RATE_ERROR 4

// define error types
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
// CLASS DEFINITION

class phiInOutParameters
{
private:
protected:
public:
    ////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////
    // input file information

    FILE *phiInOutInputFile; /*
                          * creating an inputfile to read the file
                          */
    char *phiInOutFileName;  /*
                          * storing the name of file
                          */

    int phiInOutLengthNNDataset; /*
                           * row value of NN dataset
                           */

    int phiInOutNumberNNInputDataset;  /*
                                   * column value of input NN dataset
                                   */
    int phiInOutNumberNNOutputDataset; /*
                                   * column value of output NN dataset
                                   */

    // input file information
    ////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////
    // storing the whole elements in one pointer

    double **phiInOutInputDataArray; /*
                                 * storing the whole input elements in one dynamic array
                                 */

    double **phiInOutOutputDataArray; /*
                                 * storing the whole output elements in one dynamic array
                                 */

    // storing the whole elements in one pointer
    ////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////
    // CONSTRUCTOR

    phiInOutParameters()
    {
        this->phiInOutInputFile = NULL;

        this->phiInOutFileName = NULL;

        this->phiInOutLengthNNDataset = 0;

        this->phiInOutNumberNNInputDataset = 0;
        this->phiInOutNumberNNOutputDataset = 0;

        this->phiInOutInputDataArray = NULL;

        this->phiInOutOutputDataArray = NULL;
    }

    phiInOutParameters(char *fileName)
    {
        this->phiInOutInputFile = NULL;

        this->phiInOutFileName = fileName;

        this->phiInOutLengthNNDataset = 0;

        this->phiInOutNumberNNInputDataset = 0;
        this->phiInOutNumberNNOutputDataset = 0;

        this->phiInOutInputDataArray = NULL;

        this->phiInOutOutputDataArray = NULL;
    }

    ~phiInOutParameters()
    {
    }

    // CONSTRUCTOR
    ////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////
    // type decleration for usage of input output section library

    /*
    * This structure represents the basic variables of input/output
    sessions in the other files such as excel, text, etc.
    */

    void phiInOutPrintPhi()
    {
        /*
        * simple test code to laod the whole library to the project
        * output -> returns nothing
        * input -> taking nothing
        */

        printf("\n\nPhi Libraries are loaded to the project!\n");
        printf("Process is started...\n\n");
        Sleep(5);
    }

    void phiInOutReadFileFromText()
    {
        /*
        * reading the whole data for a given text file
        * output -> returns nothing
        * input -> taking nothing 
        */

        this->phiInOutInputFile = fopen(this->phiInOutFileName, "r");

        if (this->phiInOutInputFile == NULL)
        {
            printf("Text file is not opened!\n");
            exit(EXIT_FAILURE);
        }

        double emptyReading1, emptyReading2, emptyReading3;
        double emptyReading;
        double fillReading;

        char endOfFileControl = 0;

        int varTextInfo = 0;
        int rowCounter = 0;

        while (endOfFileControl != EOF)
        {
            if (varTextInfo == 0)
            {
                varTextInfo = 1;
                fscanf(this->phiInOutInputFile, "%lf %lf %lf", &emptyReading1,
                       &emptyReading2,
                       &emptyReading3);

                printf("Row %d Input Column %d Output Column %d\n", (int)emptyReading1,
                       (int)emptyReading2,
                       (int)emptyReading3);

                this->phiInOutLengthNNDataset = (int)emptyReading1;
                this->phiInOutNumberNNInputDataset = (int)emptyReading2;
                this->phiInOutNumberNNOutputDataset = (int)emptyReading3;

                for (int j = 0; j < (this->phiInOutNumberNNOutputDataset + this->phiInOutNumberNNInputDataset - 3); j++)
                {
                    fscanf(this->phiInOutInputFile, "%lf", &emptyReading);
                }

                this->phiInOutFillInputAndOutputMatrices();

                endOfFileControl = getc(this->phiInOutInputFile);
            }
            else
            {
                if (rowCounter < this->phiInOutLengthNNDataset)
                {
                    for (int j = 0; j < this->phiInOutNumberNNInputDataset; j++)
                    {
                        fscanf(this->phiInOutInputFile, "%lf ", &fillReading);
                        this->phiInOutInputDataArray[j][rowCounter] = fillReading;
                    }

                    for (int j = 0; j < this->phiInOutNumberNNOutputDataset; j++)
                    {
                        fscanf(this->phiInOutInputFile, "%lf", &fillReading);
                        this->phiInOutOutputDataArray[j][rowCounter] = fillReading;
                    }
                }

                rowCounter++;
                endOfFileControl = getc(this->phiInOutInputFile);
            }
        }

        fclose(this->phiInOutInputFile);
    }

    void phiInOutReadFileFromCSV()
    {
        /*
         * reading the whole data for a given CSV file
         * output -> returns nothing
         * input -> address of InOutparameters structure 
         */

        this->phiInOutInputFile = fopen(this->phiInOutFileName, "r");

        if (this->phiInOutInputFile == NULL)
        {
            printf("CSV file is not opened!\n");
            exit(EXIT_FAILURE);
        }

        double emptyReading1, emptyReading2, emptyReading3;
        double emptyReading;
        double fillReading;

        char endOfFileControl = 0;

        int varTextInfo = 0;
        int rowCounter = 0;

        while (endOfFileControl != EOF)
        {
            if (varTextInfo == 0)
            {
                varTextInfo = 1;
                fscanf(this->phiInOutInputFile, "%lf,%lf,%lf,", &emptyReading1,
                       &emptyReading2,
                       &emptyReading3);

                printf("Row %d Input Column %d Output Column %d\n", (int)emptyReading1,
                       (int)emptyReading2,
                       (int)emptyReading3);

                this->phiInOutLengthNNDataset = (int)emptyReading1;
                this->phiInOutNumberNNInputDataset = (int)emptyReading2;
                this->phiInOutNumberNNOutputDataset = (int)emptyReading3;

                for (int j = 0; j < (this->phiInOutNumberNNOutputDataset + this->phiInOutNumberNNInputDataset - 3); j++)
                {
                    fscanf(this->phiInOutInputFile, "%lf,", &emptyReading);
                }

                this->phiInOutFillInputAndOutputMatrices();

                endOfFileControl = getc(this->phiInOutInputFile);
            }
            else
            {
                if (rowCounter < this->phiInOutLengthNNDataset)
                {
                    for (int j = 0; j < this->phiInOutNumberNNInputDataset; j++)
                    {
                        fscanf(this->phiInOutInputFile, "%lf,", &fillReading);
                        this->phiInOutInputDataArray[j][rowCounter] = fillReading;
                    }

                    for (int j = 0; j < this->phiInOutNumberNNOutputDataset; j++)
                    {
                        fscanf(this->phiInOutInputFile, "%lf,", &fillReading);
                        this->phiInOutOutputDataArray[j][rowCounter] = fillReading;
                    }
                }

                rowCounter++;
                endOfFileControl = getc(this->phiInOutInputFile);
            }
        }

        fclose(this->phiInOutInputFile);
    }

    /////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////

    double **phiInOutCreatingEmptyInputMatrices()
    {
        /*
         * creating dynamically empty input matrices for data storage
         * output -> address of matrices (float **)
         * input  -> address of phinioutparameters
         */

        double **pd = (double **)malloc(this->phiInOutNumberNNInputDataset * sizeof(double *));

        if (pd == NULL)
        {
            this->phiInOutPhiErrorHandler(ALLOCATION_ERROR);
            exit(EXIT_FAILURE);
        }

        for (int i = 0; i < this->phiInOutNumberNNInputDataset; i++)
            pd[i] = (double *)malloc(this->phiInOutLengthNNDataset * sizeof(double));

        return pd;
    }

    double **phiInOutCreatingEmptyOutputMatrices()
    {
        /*
         * creating dynamically empty output matrices for data storage
         * output -> address of matrices (float **)
         * input  -> taking nothing
         */

        double **pd = (double **)malloc(this->phiInOutNumberNNOutputDataset * sizeof(double *));

        if (pd == NULL)
        {
            this->phiInOutPhiErrorHandler(ALLOCATION_ERROR);
            exit(EXIT_FAILURE);
        }

        for (int i = 0; i < this->phiInOutNumberNNOutputDataset; i++)
            pd[i] = (double *)malloc(this->phiInOutLengthNNDataset * sizeof(double));

        return pd;
    }

    double **phiInOutCreatingEmptyMatrices(int rows, int cols)
    {
        /*
         * creating dynamically empty matrices for general usage
         * output -> address of empty matrices
         * input  -> rows and columns values of matrices
         */

        double **pd = (double **)malloc(rows * sizeof(double *));

        if (pd == NULL)
        {
            this->phiInOutPhiErrorHandler(ALLOCATION_ERROR);
            exit(EXIT_FAILURE);
        }

        for (int i = 0; i < rows; i++)
            pd[i] = (double *)malloc(cols * sizeof(double));

        // initialize the matrices
        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < cols; j++)
            {
                pd[i][j] = 0.0;
            }
        }

        return pd;
    }

    int **phiInOutCreatingEmptyMatricesIntegralType(int rows,
                                                    int cols)
    {
        /*
         * creating dynamically empty matrices for general usage
         * output -> address of empty matrices
         * input  -> rows and columns values of matrices
         */

        int **pd = (int **)malloc(rows * sizeof(int *));

        if (pd == NULL)
        {
            this->phiInOutPhiErrorHandler(ALLOCATION_ERROR);
            exit(EXIT_FAILURE);
        }

        for (int i = 0; i < rows; i++)
            pd[i] = (int *)malloc(cols * sizeof(int));

        // initialize the matrices
        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < cols; j++)
            {
                pd[i][j] = 0;
            }
        }

        return pd;
    }

    void phiInOutFillInputAndOutputMatrices()
    {
        /*
         * filling the whole data to the matrices
         * output -> returns nothing
         * input  -> taking nothing
         */
        this->phiInOutInputDataArray = this->phiInOutCreatingEmptyInputMatrices();
        this->phiInOutOutputDataArray = this->phiInOutCreatingEmptyOutputMatrices();
    }

    void phiInOutWriteDataSetMatrices()
    {
        /*
         * writing the whole data to the matrices
         * output -> returns nothing
         * input  -> taking nothing
         */

        printf("Printing dataset...\n");

        printf("In1 In2 In3 Out1\n");

        for (int i = 0; i < this->phiInOutLengthNNDataset; i++)
        {
            for (int j = 0; j < this->phiInOutNumberNNInputDataset; j++)
            {
                printf("%lf ", this->phiInOutInputDataArray[j][i]);
            }

            for (int j = 0; j < this->phiInOutNumberNNOutputDataset; j++)
            {
                printf("%lf ", this->phiInOutOutputDataArray[j][i]);
            }

            printf("\n");
        }
    }

    /////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////

    void phiInOutPhiFree(double **pd,
                         int row,
                         int col)
    {
        /*
         * free the whole memory allocation
         * output -> return nothing
         * input  -> address of memory
         *           row value of matrices
         *           column value of matrices
         */

        for (int i = 0; i < row; i++)
            free(pd[i]);

        free(pd);
    }

    void phiInOutPhiFreeIntegralType(int **pd,
                                     int row,
                                     int col)
    {
        /*
         * free the whole memory allocation
         * output -> return nothing
         * input  -> address of memory
         *           row value of matrices
         *           column value of matrices
         */

        for (int i = 0; i < row; i++)
            free(pd[i]);

        free(pd);
    }

    void phiInOutPhiExitInOut()
    {
        /*
         * exiting the whole variables and functions
         * output -> returns nothing
         * input -> address of ptrInOut
         */
        this->phiInOutPhiFree(this->phiInOutInputDataArray, this->phiInOutNumberNNInputDataset, this->phiInOutLengthNNDataset);
        this->phiInOutPhiFree(this->phiInOutOutputDataArray, this->phiInOutNumberNNOutputDataset, this->phiInOutLengthNNDataset);
    }

    void phiInOutPhiErrorHandler(int errorType)
    {
        /*
         * notifying the error result
         * output -> return nothing
         * input  -> type of error
         */

        switch (errorType)
        {
        case FILE_OPEN_ERROR:
            printf("System Dynamic Parameter files cannot be created!\n");
            break;
        case INCONSISTENT_ROW_COLUMN:
            printf("The rows and columns are not consistent!\n");
            break;
        case ALLOCATION_ERROR:
            printf("Memory allocation cannot be done!\n");
            break;
        case SAMPLING_RATE_ERROR:
            printf("Sampling period cannot be assigned to either negative or zero value!\n");
            break;

        default:
            break;
        }
    }

    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////

    void phiInOutReadingDemoCSV()
    {
        /*
         * reading the data from csv file
         */

        this->phiInOutFileName = "nnInputOutputFile.csv";

        this->phiInOutReadFileFromCSV();

        this->phiInOutWriteDataSetMatrices();

        this->phiInOutPhiExitInOut();
    }

    void phiInOutReadingDemoText()
    {
        /*
         * reading the data from text file
         */

        this->phiInOutFileName = "nnInputOutputFile.txt";

        this->phiInOutReadFileFromText();

        this->phiInOutWriteDataSetMatrices();

        this->phiInOutPhiExitInOut();
    }

    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////

    // type decleration for usage of NN Library
    ////////////////////////////////////////////////////////////
};

typedef phiInOutParameters *phiInOutParameterPtr;

// CLASS DEFINITION
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

#endif