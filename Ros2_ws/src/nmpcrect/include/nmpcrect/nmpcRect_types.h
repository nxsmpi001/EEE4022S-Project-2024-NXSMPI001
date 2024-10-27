//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: nmpcRect_types.h
//
// Code generated for Simulink model 'nmpcRect'.
//
// Model version                  : 1.104
// Simulink Coder version         : 24.2 (R2024b) 21-Jun-2024
// C/C++ source code generated on : Sun Oct 27 09:47:01 2024
//
// Target selection: ert.tlc
// Embedded hardware selection: Intel->x86-64 (Linux 64)
// Code generation objectives: Unspecified
// Validation result: Not run
//
#ifndef nmpcRect_types_h_
#define nmpcRect_types_h_
#include "rtwtypes.h"
#include "coder_bounded_array.h"
#ifndef DEFINED_TYPEDEF_FOR_SL_Bus_geometry_msgs_Vector3_
#define DEFINED_TYPEDEF_FOR_SL_Bus_geometry_msgs_Vector3_

// MsgType=geometry_msgs/Vector3
struct SL_Bus_geometry_msgs_Vector3
{
  real_T x;
  real_T y;
  real_T z;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_SL_Bus_geometry_msgs_Twist_
#define DEFINED_TYPEDEF_FOR_SL_Bus_geometry_msgs_Twist_

// MsgType=geometry_msgs/Twist
struct SL_Bus_geometry_msgs_Twist
{
  // MsgType=geometry_msgs/Vector3
  SL_Bus_geometry_msgs_Vector3 linear;

  // MsgType=geometry_msgs/Vector3
  SL_Bus_geometry_msgs_Vector3 angular;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_SL_Bus_std_msgs_Int32_
#define DEFINED_TYPEDEF_FOR_SL_Bus_std_msgs_Int32_

// MsgType=std_msgs/Int32
struct SL_Bus_std_msgs_Int32
{
  int32_T data;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_SL_Bus_geometry_msgs_Point_
#define DEFINED_TYPEDEF_FOR_SL_Bus_geometry_msgs_Point_

// MsgType=geometry_msgs/Point
struct SL_Bus_geometry_msgs_Point
{
  real_T x;
  real_T y;
  real_T z;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_SL_Bus_catarob_interfaces_StateEstimate_
#define DEFINED_TYPEDEF_FOR_SL_Bus_catarob_interfaces_StateEstimate_

// MsgType=catarob_interfaces/StateEstimate
struct SL_Bus_catarob_interfaces_StateEstimate
{
  real_T x;
  real_T y;
  real_T psi;
  real_T u;
  real_T v;
  real_T r;
  real_T covariance[36];
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_SL_Bus_builtin_interfaces_Time_
#define DEFINED_TYPEDEF_FOR_SL_Bus_builtin_interfaces_Time_

// MsgType=builtin_interfaces/Time
struct SL_Bus_builtin_interfaces_Time
{
  int32_T sec;
  uint32_T nanosec;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_SL_Bus_ROSVariableLengthArrayInfo_
#define DEFINED_TYPEDEF_FOR_SL_Bus_ROSVariableLengthArrayInfo_

struct SL_Bus_ROSVariableLengthArrayInfo
{
  uint32_T CurrentLength;
  uint32_T ReceivedLength;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_SL_Bus_std_msgs_Header_
#define DEFINED_TYPEDEF_FOR_SL_Bus_std_msgs_Header_

// MsgType=std_msgs/Header
struct SL_Bus_std_msgs_Header
{
  // MsgType=builtin_interfaces/Time
  SL_Bus_builtin_interfaces_Time stamp;

  // PrimitiveROSType=string:IsVarLen=1:VarLenCategory=data:VarLenElem=frame_id_SL_Info:TruncateAction=warn 
  uint8_T frame_id[128];

  // IsVarLen=1:VarLenCategory=length:VarLenElem=frame_id
  SL_Bus_ROSVariableLengthArrayInfo frame_id_SL_Info;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_SL_Bus_sensor_msgs_NavSatStatus_
#define DEFINED_TYPEDEF_FOR_SL_Bus_sensor_msgs_NavSatStatus_

// MsgType=sensor_msgs/NavSatStatus
struct SL_Bus_sensor_msgs_NavSatStatus
{
  int8_T status;
  uint16_T service;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_SL_Bus_sensor_msgs_NavSatFix_
#define DEFINED_TYPEDEF_FOR_SL_Bus_sensor_msgs_NavSatFix_

// MsgType=sensor_msgs/NavSatFix
struct SL_Bus_sensor_msgs_NavSatFix
{
  // MsgType=std_msgs/Header
  SL_Bus_std_msgs_Header header;

  // MsgType=sensor_msgs/NavSatStatus
  SL_Bus_sensor_msgs_NavSatStatus status;
  real_T latitude;
  real_T longitude;
  real_T altitude;
  real_T position_covariance[9];
  uint8_T position_covariance_type;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_SL_Bus_std_msgs_Float64_
#define DEFINED_TYPEDEF_FOR_SL_Bus_std_msgs_Float64_

// MsgType=std_msgs/Float64
struct SL_Bus_std_msgs_Float64
{
  real_T data;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_SL_Bus_geometry_msgs_Quaternion_
#define DEFINED_TYPEDEF_FOR_SL_Bus_geometry_msgs_Quaternion_

// MsgType=geometry_msgs/Quaternion
struct SL_Bus_geometry_msgs_Quaternion
{
  real_T x;
  real_T y;
  real_T z;
  real_T w;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_SL_Bus_sensor_msgs_Imu_
#define DEFINED_TYPEDEF_FOR_SL_Bus_sensor_msgs_Imu_

// MsgType=sensor_msgs/Imu
struct SL_Bus_sensor_msgs_Imu
{
  // MsgType=std_msgs/Header
  SL_Bus_std_msgs_Header header;

  // MsgType=geometry_msgs/Quaternion
  SL_Bus_geometry_msgs_Quaternion orientation;
  real_T orientation_covariance[9];

  // MsgType=geometry_msgs/Vector3
  SL_Bus_geometry_msgs_Vector3 angular_velocity;
  real_T angular_velocity_covariance[9];

  // MsgType=geometry_msgs/Vector3
  SL_Bus_geometry_msgs_Vector3 linear_acceleration;
  real_T linear_acceleration_covariance[9];
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_ParamBus_
#define DEFINED_TYPEDEF_FOR_ParamBus_

struct ParamBus
{
  real_T Constant;
  real_T Constant1;
  real_T Constant8;
};

#endif

// Custom Type definition for MATLABSystem: '<S6>/SourceBlock'
#include "rmw/qos_profiles.h"
#ifndef struct_sJ4ih70VmKcvCeguWN0mNVF
#define struct_sJ4ih70VmKcvCeguWN0mNVF

struct sJ4ih70VmKcvCeguWN0mNVF
{
  real_T sec;
  real_T nsec;
};

#endif                                 // struct_sJ4ih70VmKcvCeguWN0mNVF

#ifndef struct_ros_slros2_internal_block_Pub_T
#define struct_ros_slros2_internal_block_Pub_T

struct ros_slros2_internal_block_Pub_T
{
  boolean_T matlabCodegenIsDeleted;
  int32_T isInitialized;
  boolean_T isSetupComplete;
  boolean_T QOSAvoidROSNamespaceConventions;
};

#endif                                // struct_ros_slros2_internal_block_Pub_T

#ifndef struct_e_robotics_slcore_internal_bl_T
#define struct_e_robotics_slcore_internal_bl_T

struct e_robotics_slcore_internal_bl_T
{
  int32_T __dummy;
};

#endif                                // struct_e_robotics_slcore_internal_bl_T

#ifndef struct_ros_slros2_internal_block_Get_T
#define struct_ros_slros2_internal_block_Get_T

struct ros_slros2_internal_block_Get_T
{
  boolean_T matlabCodegenIsDeleted;
  int32_T isInitialized;
  boolean_T isSetupComplete;
  e_robotics_slcore_internal_bl_T SampleTimeHandler;
};

#endif                                // struct_ros_slros2_internal_block_Get_T

#ifndef struct_matlabshared_tracking_interna_T
#define struct_matlabshared_tracking_interna_T

struct matlabshared_tracking_interna_T
{
  int32_T isInitialized;
};

#endif                                // struct_matlabshared_tracking_interna_T

#ifndef struct_matlabshared_tracking_inter_a_T
#define struct_matlabshared_tracking_inter_a_T

struct matlabshared_tracking_inter_a_T
{
  int32_T isInitialized;
};

#endif                                // struct_matlabshared_tracking_inter_a_T

// Custom Type definition for MATLAB Function: '<S24>/NLMPC'
#ifndef struct_sG8JZ69axY52WWR6RKyApQC_nmpcR_T
#define struct_sG8JZ69axY52WWR6RKyApQC_nmpcR_T

struct sG8JZ69axY52WWR6RKyApQC_nmpcR_T
{
  real_T penaltyParam;
  real_T threshold;
  int32_T nPenaltyDecreases;
  real_T linearizedConstrViol;
  real_T initFval;
  real_T initConstrViolationEq;
  real_T initConstrViolationIneq;
  real_T phi;
  real_T phiPrimePlus;
  real_T phiFullStep;
  real_T feasRelativeFactor;
  real_T nlpPrimalFeasError;
  real_T nlpDualFeasError;
  real_T nlpComplError;
  real_T firstOrderOpt;
  boolean_T hasObjective;
};

#endif                                // struct_sG8JZ69axY52WWR6RKyApQC_nmpcR_T

#ifndef struct_s7RdrPWkr8UPAUyTdDJkLaG_nmpcR_T
#define struct_s7RdrPWkr8UPAUyTdDJkLaG_nmpcR_T

struct s7RdrPWkr8UPAUyTdDJkLaG_nmpcR_T
{
  boolean_T gradOK;
  boolean_T fevalOK;
  boolean_T done;
  boolean_T stepAccepted;
  boolean_T failedLineSearch;
  int32_T stepType;
};

#endif                                // struct_s7RdrPWkr8UPAUyTdDJkLaG_nmpcR_T

#ifndef struct_ros_slros2_internal_block_Sub_T
#define struct_ros_slros2_internal_block_Sub_T

struct ros_slros2_internal_block_Sub_T
{
  boolean_T matlabCodegenIsDeleted;
  int32_T isInitialized;
  boolean_T isSetupComplete;
  boolean_T QOSAvoidROSNamespaceConventions;
};

#endif                                // struct_ros_slros2_internal_block_Sub_T

// Custom Type definition for MATLAB Function: '<S24>/NLMPC'
#ifndef struct_sttYSJM5GCi2c1Eu0R50efC_nmpcR_T
#define struct_sttYSJM5GCi2c1Eu0R50efC_nmpcR_T

struct sttYSJM5GCi2c1Eu0R50efC_nmpcR_T
{
  real_T iterations;
  real_T funcCount;
  char_T algorithm[3];
  real_T constrviolation;
  real_T stepsize;
  real_T lssteplength;
  real_T firstorderopt;
};

#endif                                // struct_sttYSJM5GCi2c1Eu0R50efC_nmpcR_T

#ifndef struct_somzaGboVhDG7PNQS6E98jD_nmpcR_T
#define struct_somzaGboVhDG7PNQS6E98jD_nmpcR_T

struct somzaGboVhDG7PNQS6E98jD_nmpcR_T
{
  char_T SolverName[7];
  int32_T MaxIterations;
  real_T StepTolerance;
  real_T OptimalityTolerance;
  real_T ConstraintTolerance;
  real_T ObjectiveLimit;
  real_T PricingTolerance;
  real_T ConstrRelTolFactor;
  real_T ProbRelTolFactor;
  boolean_T RemainFeasible;
  boolean_T IterDisplayQP;
};

#endif                                // struct_somzaGboVhDG7PNQS6E98jD_nmpcR_T

#ifndef struct_s_YHCj6MDTSR65so1fxJsjSD_nmpc_T
#define struct_s_YHCj6MDTSR65so1fxJsjSD_nmpc_T

struct s_YHCj6MDTSR65so1fxJsjSD_nmpc_T
{
  real_T constants[67];
};

#endif                                // struct_s_YHCj6MDTSR65so1fxJsjSD_nmpc_T

// Custom Type definition for MATLAB Function: '<S24>/NLMPC'
#ifndef struct_s_2COE1uYisQtyPYvPjrXP9G_nmpc_T
#define struct_s_2COE1uYisQtyPYvPjrXP9G_nmpc_T

struct s_2COE1uYisQtyPYvPjrXP9G_nmpc_T
{
  int32_T nVarMax;
  int32_T mNonlinIneq;
  int32_T mNonlinEq;
  int32_T mIneq;
  int32_T mEq;
  int32_T iNonIneq0;
  int32_T iNonEq0;
  real_T sqpFval;
  real_T sqpFval_old;
  real_T xstarsqp[125];
  real_T xstarsqp_old[125];
  coder::bounded_array<real_T, 320U, 1U> cIneq;
  coder::bounded_array<real_T, 320U, 1U> cIneq_old;
  real_T cEq[120];
  real_T cEq_old[120];
  coder::bounded_array<real_T, 686U, 1U> grad;
  coder::bounded_array<real_T, 686U, 1U> grad_old;
  int32_T FunctionEvaluations;
  int32_T sqpIterations;
  int32_T sqpExitFlag;
  coder::bounded_array<real_T, 1251U, 1U> lambdasqp;
  coder::bounded_array<real_T, 1251U, 1U> lambdaStopTest;
  coder::bounded_array<real_T, 1251U, 1U> lambdaStopTestPrev;
  real_T steplength;
  coder::bounded_array<real_T, 1251U, 1U> delta_x;
  coder::bounded_array<real_T, 686U, 1U> socDirection;
  coder::bounded_array<int32_T, 1251U, 1U> workingset_old;
  coder::bounded_array<real_T, 109760U, 2U> JacCineqTrans_old;
  coder::bounded_array<real_T, 82320U, 2U> JacCeqTrans_old;
  coder::bounded_array<real_T, 686U, 1U> gradLag;
  coder::bounded_array<real_T, 686U, 1U> delta_gradLag;
  coder::bounded_array<real_T, 1251U, 1U> xstar;
  real_T fstar;
  real_T firstorderopt;
  coder::bounded_array<real_T, 1251U, 1U> lambda;
  int32_T state;
  real_T maxConstr;
  int32_T iterations;
  coder::bounded_array<real_T, 1251U, 1U> searchDir;
};

#endif                                // struct_s_2COE1uYisQtyPYvPjrXP9G_nmpc_T

// Custom Type definition for MATLAB Function: '<S24>/NLMPC'
#ifndef struct_s_kmYqIq13KlaOrGCTq3ShMG_nmpc_T
#define struct_s_kmYqIq13KlaOrGCTq3ShMG_nmpc_T

struct s_kmYqIq13KlaOrGCTq3ShMG_nmpc_T
{
  coder::bounded_array<real_T, 858186U, 2U> workspace_float;
  coder::bounded_array<int32_T, 1251U, 1U> workspace_int;
  coder::bounded_array<int32_T, 1251U, 1U> workspace_sort;
};

#endif                                // struct_s_kmYqIq13KlaOrGCTq3ShMG_nmpc_T

// Custom Type definition for MATLAB Function: '<S24>/NLMPC'
#ifndef struct_s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T
#define struct_s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T

struct s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T
{
  int32_T mConstr;
  int32_T mConstrOrig;
  int32_T mConstrMax;
  int32_T nVar;
  int32_T nVarOrig;
  int32_T nVarMax;
  int32_T ldA;
  coder::bounded_array<real_T, 219520U, 1U> Aineq;
  coder::bounded_array<real_T, 320U, 1U> bineq;
  coder::bounded_array<real_T, 82320U, 1U> Aeq;
  real_T beq[120];
  coder::bounded_array<real_T, 686U, 1U> lb;
  coder::bounded_array<real_T, 686U, 1U> ub;
  coder::bounded_array<int32_T, 686U, 1U> indexLB;
  coder::bounded_array<int32_T, 686U, 1U> indexUB;
  coder::bounded_array<int32_T, 686U, 1U> indexFixed;
  int32_T mEqRemoved;
  int32_T indexEqRemoved[120];
  coder::bounded_array<real_T, 858186U, 1U> ATwset;
  coder::bounded_array<real_T, 1251U, 1U> bwset;
  int32_T nActiveConstr;
  coder::bounded_array<real_T, 1251U, 1U> maxConstrWorkspace;
  int32_T sizes[5];
  int32_T sizesNormal[5];
  int32_T sizesPhaseOne[5];
  int32_T sizesRegularized[5];
  int32_T sizesRegPhaseOne[5];
  int32_T isActiveIdx[6];
  int32_T isActiveIdxNormal[6];
  int32_T isActiveIdxPhaseOne[6];
  int32_T isActiveIdxRegularized[6];
  int32_T isActiveIdxRegPhaseOne[6];
  coder::bounded_array<boolean_T, 1251U, 1U> isActiveConstr;
  coder::bounded_array<int32_T, 1251U, 1U> Wid;
  coder::bounded_array<int32_T, 1251U, 1U> Wlocalidx;
  int32_T nWConstr[5];
  int32_T probType;
  real_T SLACK0;
};

#endif                                // struct_s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T

// Custom Type definition for MATLAB Function: '<S24>/NLMPC'
#ifndef struct_s_xtSBzQGTZuMYOTjcuMqLQH_nmpc_T
#define struct_s_xtSBzQGTZuMYOTjcuMqLQH_nmpc_T

struct s_xtSBzQGTZuMYOTjcuMqLQH_nmpc_T
{
  coder::bounded_array<real_T, 686U, 1U> grad;
  coder::bounded_array<real_T, 685U, 1U> Hx;
  boolean_T hasLinear;
  int32_T nvar;
  int32_T maxVar;
  real_T beta;
  real_T rho;
  int32_T objtype;
  int32_T prev_objtype;
  int32_T prev_nvar;
  boolean_T prev_hasLinear;
  real_T gammaScalar;
};

#endif                                // struct_s_xtSBzQGTZuMYOTjcuMqLQH_nmpc_T

// Custom Type definition for MATLAB Function: '<S24>/NLMPC'
#ifndef struct_s_0RmwrXfzGd5lqbHvgKQe2_nmpcR_T
#define struct_s_0RmwrXfzGd5lqbHvgKQe2_nmpcR_T

struct s_0RmwrXfzGd5lqbHvgKQe2_nmpcR_T
{
  int32_T ldq;
  coder::bounded_array<real_T, 1565001U, 2U> QR;
  coder::bounded_array<real_T, 1565001U, 2U> Q;
  coder::bounded_array<int32_T, 1251U, 1U> jpvt;
  int32_T mrows;
  int32_T ncols;
  coder::bounded_array<real_T, 1251U, 1U> tau;
  int32_T minRowCol;
  boolean_T usedPivoting;
};

#endif                                // struct_s_0RmwrXfzGd5lqbHvgKQe2_nmpcR_T

#ifndef struct_s_mDApTYzDBpxuvxemclsuEF_nmpc_T
#define struct_s_mDApTYzDBpxuvxemclsuEF_nmpc_T

struct s_mDApTYzDBpxuvxemclsuEF_nmpc_T
{
  coder::bounded_array<real_T, 1565001U, 2U> FMat;
  int32_T ldm;
  int32_T ndims;
  int32_T info;
  real_T scaleFactor;
  boolean_T ConvexCheck;
  real_T regTol_;
  real_T workspace_;
  real_T workspace2_;
};

#endif                                // struct_s_mDApTYzDBpxuvxemclsuEF_nmpc_T

#ifndef struct_s_jex761Cl1dvQqVqRqjms8C_nmpc_T
#define struct_s_jex761Cl1dvQqVqRqjms8C_nmpc_T

struct s_jex761Cl1dvQqVqRqjms8C_nmpc_T
{
  real_T x[6];
  real_T lastMV[2];
  real_T ref[80];
  real_T OutputWeights[80];
  real_T MVWeights[40];
  real_T MVRateWeights[40];
  real_T ECRWeight;
  real_T OutputMin[80];
  real_T OutputMax[80];
  real_T StateMin[120];
  real_T StateMax[120];
  real_T MVMin[40];
  real_T MVMax[40];
  real_T MVRateMin[40];
  real_T MVRateMax[40];
  real_T MVScaledTarget[40];
  real_T Parameters[3];
};

#endif                                // struct_s_jex761Cl1dvQqVqRqjms8C_nmpc_T

#ifndef struct_sAc4bxvmmjmxjQV9i3feLrE_nmpcR_T
#define struct_sAc4bxvmmjmxjQV9i3feLrE_nmpcR_T

struct sAc4bxvmmjmxjQV9i3feLrE_nmpcR_T
{
  real_T Ts;
  real_T CurrentStates[6];
  real_T LastMV[2];
  real_T References[80];
  real_T MVTarget[40];
  real_T PredictionHorizon;
  real_T NumOfStates;
  real_T NumOfOutputs;
  real_T NumOfInputs;
  real_T MVIndex[2];
  real_T UDIndex[2];
  real_T InputPassivityIndex;
  real_T OutputPassivityIndex;
  boolean_T PassivityUsePredictedX;
};

#endif                                // struct_sAc4bxvmmjmxjQV9i3feLrE_nmpcR_T

#ifndef struct_s_I4XPpWw7d7shktLagLlNtD_nmpc_T
#define struct_s_I4XPpWw7d7shktLagLlNtD_nmpc_T

struct s_I4XPpWw7d7shktLagLlNtD_nmpc_T
{
  s_jex761Cl1dvQqVqRqjms8C_nmpc_T runtimedata;
  sAc4bxvmmjmxjQV9i3feLrE_nmpcR_T userdata;
};

#endif                                // struct_s_I4XPpWw7d7shktLagLlNtD_nmpc_T

#ifndef struct_anonymous_function_nmpcRect_T
#define struct_anonymous_function_nmpcRect_T

struct anonymous_function_nmpcRect_T
{
  s_I4XPpWw7d7shktLagLlNtD_nmpc_T workspace;
};

#endif                                 // struct_anonymous_function_nmpcRect_T

#ifndef struct_coder_internal_stickyStruct_1_T
#define struct_coder_internal_stickyStruct_1_T

struct coder_internal_stickyStruct_1_T
{
  anonymous_function_nmpcRect_T b_value;
};

#endif                                // struct_coder_internal_stickyStruct_1_T

#ifndef struct_coder_internal_stickyStruct_j_T
#define struct_coder_internal_stickyStruct_j_T

struct coder_internal_stickyStruct_j_T
{
  anonymous_function_nmpcRect_T b_value;
  coder_internal_stickyStruct_1_T next;
};

#endif                                // struct_coder_internal_stickyStruct_j_T

#ifndef struct_coder_internal_stickyStruct_f_T
#define struct_coder_internal_stickyStruct_f_T

struct coder_internal_stickyStruct_f_T
{
  coder_internal_stickyStruct_j_T next;
};

#endif                                // struct_coder_internal_stickyStruct_f_T

#ifndef struct_coder_internal_stickyStruct_m_T
#define struct_coder_internal_stickyStruct_m_T

struct coder_internal_stickyStruct_m_T
{
  int32_T b_value;
  coder_internal_stickyStruct_f_T next;
};

#endif                                // struct_coder_internal_stickyStruct_m_T

#ifndef struct_coder_internal_stickyStruct_d_T
#define struct_coder_internal_stickyStruct_d_T

struct coder_internal_stickyStruct_d_T
{
  coder_internal_stickyStruct_m_T next;
};

#endif                                // struct_coder_internal_stickyStruct_d_T

#ifndef struct_coder_internal_stickyStruct_a_T
#define struct_coder_internal_stickyStruct_a_T

struct coder_internal_stickyStruct_a_T
{
  coder_internal_stickyStruct_d_T next;
};

#endif                                // struct_coder_internal_stickyStruct_a_T

#ifndef struct_coder_internal_stickyStruct_i_T
#define struct_coder_internal_stickyStruct_i_T

struct coder_internal_stickyStruct_i_T
{
  coder_internal_stickyStruct_a_T next;
};

#endif                                // struct_coder_internal_stickyStruct_i_T

#ifndef struct_coder_internal_stickyStruc_mn_T
#define struct_coder_internal_stickyStruc_mn_T

struct coder_internal_stickyStruc_mn_T
{
  coder_internal_stickyStruct_i_T next;
};

#endif                                // struct_coder_internal_stickyStruc_mn_T

#ifndef struct_coder_internal_stickyStruct_2_T
#define struct_coder_internal_stickyStruct_2_T

struct coder_internal_stickyStruct_2_T
{
  coder_internal_stickyStruc_mn_T next;
};

#endif                                // struct_coder_internal_stickyStruct_2_T

#ifndef struct_s_TO7fEqhVtfGgzjgQpVOcHG_nmpc_T
#define struct_s_TO7fEqhVtfGgzjgQpVOcHG_nmpc_T

struct s_TO7fEqhVtfGgzjgQpVOcHG_nmpc_T
{
  anonymous_function_nmpcRect_T objfun;
  anonymous_function_nmpcRect_T nonlin;
  real_T f_1;
  coder::bounded_array<real_T, 160U, 1U> cIneq_1;
  real_T cEq_1[120];
  real_T f_2;
  coder::bounded_array<real_T, 160U, 1U> cIneq_2;
  real_T cEq_2[120];
  int32_T nVar;
  int32_T mIneq;
  int32_T mEq;
  int32_T numEvals;
  boolean_T SpecifyObjectiveGradient;
  boolean_T SpecifyConstraintGradient;
  boolean_T isEmptyNonlcon;
  boolean_T hasLB[125];
  boolean_T hasUB[125];
  boolean_T hasBounds;
  int32_T FiniteDifferenceType;
};

#endif                                // struct_s_TO7fEqhVtfGgzjgQpVOcHG_nmpc_T

// Forward declaration for rtModel
typedef struct tag_RTM_nmpcRect_T RT_MODEL_nmpcRect_T;

#endif                                 // nmpcRect_types_h_

//
// File trailer for generated code.
//
// [EOF]
//
