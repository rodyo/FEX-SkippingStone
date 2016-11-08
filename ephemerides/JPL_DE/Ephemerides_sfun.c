#define S_FUNCTION_NAME  Ephemerides_sfun
#define S_FUNCTION_LEVEL 2

#include "Ephemerides_wrapper.h"
#include "simstruc.h"
#include "Utils.h"
#include "singletonSfunction.h"


//#define MDL_SET_WORK_WIDTHS
#define MDL_START
#define MDL_INITIAL_SIZES
#define MDL_OUTPUTS
#define MDL_INITIALIZE_SAMPLE_TIMES
#define MDL_SET_INPUT_PORT_FRAME_DATA
#define MDL_SET_INPUT_PORT_DATA_TYPE
#define MDL_SET_OUTPUT_PORT_DATA_TYPE
#define MDL_SET_DEFAULT_PORT_DATA_TYPES


#ifdef MDL_SET_WORK_WIDTHS
// TODO: (Rody Oldenhuis) crashes in R2010a
static
void
mdlSetWorkWidths(SimStruct *S)
{
    ssSetNumRunTimeParams(S,0);
}
#endif

#ifdef MDL_INITIAL_SIZES
static
void
mdlInitializeSizes(SimStruct *S)
{
    unsigned int
            i = 0u,

            num_params = 0u,

            num_inputs   = 1u,
            input_dims[] = {1u},

            num_outputs   = 1u,
            output_dims[] = {6u};

    /* Check consistency of signal size and chosen bodies */
    Ephem_check_params(output_dims[0]);

    /* Ensure that there can be only 1 instance of this S-function */
    ssVerifySingleInstance(S);

    /* Parameters */
    ssSetNumSFcnParams(S, num_params);
    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S))
        return;

    /* Inputs */
    if (!ssSetNumInputPorts(S, num_inputs))
        return;

    for (i=0u; i<num_inputs; ++i) {
        ssSetInputPortWidth(S, 0, input_dims[i]);
        ssSetInputPortDataType(S, 0, SS_DOUBLE);
        ssSetInputPortComplexSignal(S, 0, COMPLEX_NO);
        ssSetInputPortDirectFeedThrough(S, 0, 1);
        ssSetInputPortRequiredContiguous(S, 0, 1);
    }

    /* Outputs */
    if (!ssSetNumOutputPorts(S, num_outputs))
        return;

    for (i=0u; i<num_outputs; ++i) {
        ssSetOutputPortWidth(S, 0, output_dims[i]);
        ssSetOutputPortDataType(S, 0, SS_DOUBLE);
        ssSetOutputPortComplexSignal(S, 0, COMPLEX_NO);
    }

    /* S-Function options */
    ssSetNumContStates(S, 0);  ssSetNumSampleTimes(S, 1);
    ssSetNumDiscStates(S, 0);  ssSetNumNonsampledZCs(S, 0);

    ssSetNumRWork(S, 0);  ssSetNumPWork(S, 0);
    ssSetNumIWork(S, 0);  ssSetNumModes(S, 0);

    ssSetModelReferenceNormalModeSupport(S,
            MDL_START_AND_MDL_PROCESS_PARAMS_OK);

    ssSetOptions(S,
                 SS_OPTION_REQ_INPUT_SAMPLE_TIME_MATCH |
                 SS_OPTION_CALL_TERMINATE_ON_EXIT   |
                 SS_OPTION_USE_TLC_WITH_ACCELERATOR |
                 SS_OPTION_CAN_BE_CALLED_CONDITIONALLY |
                 SS_OPTION_SFUNCTION_INLINED_FOR_RTW );
}
#endif

#ifdef MDL_SET_INPUT_PORT_FRAME_DATA
static
void
mdlSetInputPortFrameData(SimStruct  *S,
                         int_T      port,
                         Frame_T    frameData)
{
    ssSetInputPortFrameData(S, port, frameData);
}
#endif


#ifdef MDL_INITIALIZE_SAMPLE_TIMES
static
void
mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, INHERITED_SAMPLE_TIME);
    ssSetOffsetTime(S, 0, 0.0);
    ssSetModelReferenceSampleTimeDefaultInheritance(S);
}
#endif

#ifdef MDL_SET_INPUT_PORT_DATA_TYPE
static
void
mdlSetInputPortDataType(SimStruct *S,
                        int port,
                        DTypeId dType)
{
    UNREFERENCED_FORMAL_PARAMETER(port);
    ssSetInputPortDataType(S, 0, dType);
}
#endif

#ifdef MDL_SET_OUTPUT_PORT_DATA_TYPE
static
void
mdlSetOutputPortDataType(SimStruct *S,
                        int port,
                        DTypeId dType)
{
    UNREFERENCED_FORMAL_PARAMETER(port);
    ssSetOutputPortDataType(S, 0, dType);
}
#endif

#ifdef MDL_SET_DEFAULT_PORT_DATA_TYPES
static
void
mdlSetDefaultPortDataTypes(SimStruct *S)
{
    ssSetInputPortDataType (S, 0, SS_DOUBLE);
    ssSetOutputPortDataType(S, 0, SS_DOUBLE);
}
#endif




#ifdef MDL_START
static
void
mdlStart(SimStruct *S)
{
    UNREFERENCED_FORMAL_PARAMETER(S);
    Ephem_initialize(QUOTE(S_FUNCTION_NAME));
}
#endif /*  MDL_START */


#ifdef MDL_OUTPUTS
static
void
mdlOutputs(SimStruct *S,
           int_T tid)
{
    const double  JD  = *((const double*) ssGetInputPortSignal     (S,0));
    double *positions =         (double*) ssGetOutputPortRealSignal(S,0);

    UNREFERENCED_FORMAL_PARAMETER(tid);

    Ephem_step(JD, positions);
}
#endif


static
void
mdlTerminate(SimStruct *S)
{
    ssKillSingleton(S);
    Ephem_cleanup();
}


#if RUNNING_IN_SIMULINK
    #include "simulink.c" /* MEX-file interface mechanism */
#else
    #include "cg_sfun.h"  /* Code generation registration function */
#endif

