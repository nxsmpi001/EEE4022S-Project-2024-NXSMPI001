//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: nmpcSine.h
//
// Code generated for Simulink model 'nmpcSine'.
//
// Model version                  : 1.102
// Simulink Coder version         : 24.2 (R2024b) 21-Jun-2024
// C/C++ source code generated on : Sun Oct 27 09:42:42 2024
//
// Target selection: ert.tlc
// Embedded hardware selection: Intel->x86-64 (Linux 64)
// Code generation objectives: Unspecified
// Validation result: Not run
//
#ifndef nmpcSine_h_
#define nmpcSine_h_
#include <stdio.h>
#include <string.h>
#include "rtwtypes.h"
#include "slros2_initialize.h"
#include "nmpcSine_types.h"

extern "C"
{

#include "rt_nonfinite.h"

}

extern "C"
{

#include "rtGetInf.h"

}

extern "C"
{

#include "rtGetNaN.h"

}

#include <stddef.h>

// Block signals (default storage)
struct B_nmpcSine_T {
  s_2COE1uYisQtyPYvPjrXP9G_nmpc_T TrialState;
  s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T WorkingSet;
  s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T b_obj;
  s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T c_WorkingSet;
  s_0RmwrXfzGd5lqbHvgKQe2_nmpcS_T QRManager;
  s_kmYqIq13KlaOrGCTq3ShMG_nmpc_T memspace;
  s_xtSBzQGTZuMYOTjcuMqLQH_nmpc_T QPObjective;
  s_mDApTYzDBpxuvxemclsuEF_nmpc_T CholManager;
  real_T y_data[1565001];
  real_T B_data[858186];
  real_T A_data[25760];
  real_T varargin_1_data[25600];
  real_T JacCineqTrans_data[20000];
  real_T b_varargin_1_data[20000];
  real_T a__4_data[20000];
  real_T Jx[19200];
  real_T Jx_data[19200];
  real_T varargin_1_data_m[19200];
  real_T unusedExpr[15625];
  real_T JacCeqTrans[15000];
  real_T JacEqTrans_tmp[15000];
  real_T Jx_c[14400];
  real_T Auf_data[6400];
  real_T tmp_data[6400];
  real_T Jmv[4800];
  coder_internal_stickyStruct_2_T FcnEvaluator;
  real_T y_data_k[1251];
  real_T y_data_c[1251];
  real_T work_data[1251];
  real_T work_data_b[1251];
  real_T vn1_data[1251];
  real_T vn2_data[1251];
  real_T work_data_p[1251];
  real_T b_data[1251];
  real_T y_data_cv[1251];
  real_T y_data_f[1251];
  s_jex761Cl1dvQqVqRqjms8C_nmpc_T expl_temp;
  s_jex761Cl1dvQqVqRqjms8C_nmpc_T expl_temp_g;
  int8_T Au[6400];
  real_T tmp_data_g[640];
  real_T varargin_2_data[640];
  real_T tmp_data_m[640];
  real_T Jmv_n[480];
  real_T B_data_p[160];
  real_T Cineq_data[160];
  real_T a__3_data[160];
  real_T varargin_1_data_l[160];
  real_T varargin_1_data_j[160];
  real_T c[160];
  real_T b_c[160];
  real_T Je_data[160];
  real_T b_Bu[160];
  sAc4bxvmmjmxjQV9i3feLrE_nmpcS_T expl_temp_d;
  real_T x_Delay[126];                 // '<S25>/x_Delay'
  real_T X[126];
  real_T b_X[126];
  real_T X_g[126];
  real_T b_X_l[126];
  real_T X_d[126];
  real_T z0[125];
  real_T zUB[125];
  real_T z[125];
  real_T dv[125];
  real_T Gfxu_data[125];
  real_T Selector[120];
  real_T Ceq[120];
  real_T b_x[120];
  real_T reshapes_f1[120];
  real_T z_d[120];
  real_T Gfxu_data_l[120];
  real_T Gfxu_data_o[120];
  real_T Selector_b[114];
  real_T U0[84];
  real_T U[84];
  real_T b_U[84];
  real_T U_n[84];
  real_T b_U_b[84];
  real_T U_l[84];
  int32_T ineqRange_data[160];
  real_T a__1[72];
  real_T y[72];
  real_T b_J[60];
  real_T c_J[60];
  real_T arg478[60];
  SL_Bus_sensor_msgs_Imu In1;          // '<S55>/In1'
  SL_Bus_sensor_msgs_Imu rtb_SourceBlock_o2_h;
  real_T a__1_b[48];
  real_T y_d[48];
  SL_Bus_catarob_interfaces_StateEstimate BusAssignment;// '<S17>/Bus Assignment' 
  real_T a__1_e[42];
  real_T y_b[42];
  real_T Umv[42];
  real_T c_A[42];
  real_T Selector1[40];
  real_T dv1[40];
  real_T dv2[40];
  real_T dv3[40];
  real_T y_j[40];
  real_T Gfuu_data[40];
  real_T Gfuu_data_f[40];
  real_T err_data[40];
  real_T err[40];
  real_T Selector1_a[38];
  real_T A[36];
  real_T y_ju[36];
  real_T arg278[36];
  real_T arg287[36];
  real_T arg296[36];
  real_T arg301[36];
  real_T arg305[36];
  real_T arg309[36];
  real_T arg313[36];
  real_T arg317[36];
  real_T arg321[36];
  real_T arg324[36];
  real_T arg328[36];
  real_T A_j[36];
  real_T S[36];
  SL_Bus_sensor_msgs_NavSatFix In1_j;  // '<S53>/In1'
  real_T VectorConcatenate1[80];       // '<S3>/Vector Concatenate1'
  SL_Bus_sensor_msgs_NavSatFix rtb_SourceBlock_o2_e_o;
  real_T arg283[24];
  real_T arg292[24];
  real_T arg337[24];
  real_T arg340[24];
  real_T arg349[24];
  real_T arg357[24];
  real_T arg396[24];
  real_T arg401[24];
  real_T arg425[24];
  real_T arg430[24];
  boolean_T icf[160];
  uint8_T tmp_data_n[160];
  real_T err_x_data[20];
  real_T err_y_data[20];
  boolean_T icf_i[160];
  int8_T Je[160];
  uint8_T tmp_data_o[160];
  uint8_T ii_data[160];
  boolean_T x[160];
  real_T x_data[19];
  real_T y_n[19];
  real_T a__1_m[16];
  real_T K[16];
  sG8JZ69axY52WWR6RKyApQC_nmpcS_T MeritFunction;
  real_T b_dHdx[12];
  real_T K_c[12];
  real_T C[12];
  real_T C_m[12];
  real_T b_X_m[10];
  boolean_T bv[80];
  boolean_T bv1[80];
  somzaGboVhDG7PNQS6E98jD_nmpcS_T expl_temp_j;
  somzaGboVhDG7PNQS6E98jD_nmpcS_T expl_temp_h;
  real_T c_A_c[7];
  real_T M[7];
  sttYSJM5GCi2c1Eu0R50efC_nmpcS_T Out;
  SL_Bus_geometry_msgs_Twist BusAssignment_c;// '<S1>/Bus Assignment'
  real_T imvec[6];
  real_T Pxy[6];
  real_T z_c[6];
  real_T ic[6];
  real_T a__4[6];
  real_T b_X_p[6];
  real_T arg281[6];
  real_T arg288_tmp[6];
  real_T arg296_tmp[6];
  real_T arg296_tmp_p[6];
  real_T ic_a[6];
  real_T b_X_e[6];
  real_T b_X_a[6];
  real_T K_a[6];
  real_T b_tau[6];
  real_T work[6];
  real_T imvec_i[6];
  real_T k1[6];
  real_T k2[6];
  real_T k3[6];
  real_T x_l[6];
  int32_T icf_tmp[8];
  int32_T icf_tmp_o[8];
  real_T Sy[4];
  real_T R[4];
  real_T b[4];
  SL_Bus_geometry_msgs_Point BusAssignment_h;// '<S16>/Bus Assignment'
  real_T params[3];                    // '<Root>/Control Logic'
  real_T TmpSignalConversionAtMATLAB[3];
  real_T mv[2];                        // '<S24>/NLMPC'
  real_T gps0[3];                      // '<Root>/Control Logic'
  real_T ref[2];                       // '<S3>/MATLAB Function1'
  real_T dv4[2];
  real_T dv5[2];
  real_T dv6[2];
  real_T dv7[2];
  real_T dv8[2];
  real_T dv9[2];
  real_T dv10[2];
  int32_T Jx_size[3];
  s7RdrPWkr8UPAUyTdDJkLaG_nmpcS_T Flags;
  boolean_T icf_o[8];
  boolean_T icf_ip[8];
  real_T s1;
  real_T c1;
  real_T s2;
  real_T c2;
  real_T w1;
  real_T w2;
  real_T b_g;
  real_T c_c;
  real_T e;
  real_T epsilon;
  real_T l_data;
  real_T o;
  real_T psi;                          // '<S3>/HeadingToPsi'
  real_T y_data_idx_0;
  real_T resToWrap_data_idx_0;
  real_T b_c_o;
  real_T scale;
  real_T absxk;
  real_T t;
  real_T TrialState_lambdasqp;
  real_T e_l;
  real_T t2;
  real_T t3;
  real_T b_t2;
  real_T b_t3;
  real_T b_X_mv;
  real_T b_X_mj;
  real_T b_X_c;
  real_T b_U_f;
  real_T b_U_p;
  real_T b_U_e;
  real_T b_U_o;
  real_T b_X_tmp;
  real_T b_X_tmp_h;
  real_T b_X_tmp_l;
  real_T b_X_tmp_h2;
  real_T b_X_tmp_m;
  real_T b_X_tmp_mc;
  real_T b_X_tmp_h3;
  real_T arg45;
  real_T arg96;
  real_T arg121;
  real_T arg132;
  real_T arg140;
  real_T arg150;
  real_T arg163;
  real_T arg176;
  real_T arg212;
  real_T arg222;
  real_T arg232;
  real_T arg244;
  real_T arg258;
  real_T arg272;
  real_T arg10_tmp;
  real_T arg6_tmp;
  real_T extraParams;
  real_T extraParams_c;
  real_T extraParams_k;
  real_T arg288_tmp_tmp;
  real_T phi_alpha;
  real_T e_p;
  real_T t2_p;
  real_T t3_p;
  real_T b_t2_a;
  real_T b_t3_j;
  real_T b_X_ek;
  real_T b_X_o;
  real_T b_U_bb;
  real_T b_U_a;
  real_T b_U_g;
  real_T b_U_ex;
  real_T b_yk_idx_0;
  real_T b_yk_idx_1;
  real_T b_yk_idx_2;
  real_T b_yk_idx_3;
  real_T b_ic_idx_0;
  real_T b_X_tmp_f;
  real_T b_X_tmp_h22;
  real_T b_X_tmp_e;
  real_T z_ch;
  real_T z_a;
  real_T z_da;
  real_T z_af;
  real_T e_pb;
  real_T R_m;
  real_T obj_next_next_next_next_next_ne;
  real_T obj_next_next_next_next_next__o;
  real_T data_References;
  real_T data_References_n;
  real_T varargin_1;
  real_T J_tmp;
  real_T yk_idx_0;
  real_T yk_idx_1;
  real_T yk_idx_2;
  real_T yk_idx_3;
  real_T ic_idx_0;
  real_T ic_idx_1;
  real_T ic_idx_2;
  real_T ic_idx_3;
  real_T runtimedata_OutputMin;
  real_T runtimedata_OutputMin_l;
  real_T runtimedata_OutputMin_p;
  real_T runtimedata_OutputMin_pt;
  real_T runtimedata_OutputMax;
  real_T runtimedata_OutputMax_f;
  real_T runtimedata_OutputMax_i;
  real_T runtimedata_OutputMax_o;
  real_T runtimedata_MVRateMin;
  real_T runtimedata_MVRateMin_k;
  real_T runtimedata_MVRateMax;
  real_T runtimedata_MVRateMax_i;
  real_T runtimedata_MVMin;
  real_T runtimedata_MVMin_o;
  real_T runtimedata_MVMax;
  real_T runtimedata_MVMax_m;
  real_T s;
  real_T beta1;
  real_T c_cu;
  real_T scale_f;
  real_T absxk_h;
  real_T t_m;
  real_T nrmGradInf;
  real_T nrmDirInf;
  real_T u1;
  real_T b_c_a;
  real_T maxConstr_new;
  real_T normDelta;
  real_T solution_lambda;
  real_T tempMaxConstr;
  real_T constrViolation_basicX;
  real_T temp;
  real_T temp2;
  real_T smax;
  real_T s_k;
  real_T xnorm;
  real_T a;
  real_T scale_p;
  real_T absxk_b;
  real_T t_c;
  real_T data_References_nb;
  real_T data_References_i;
  real_T U_m;
  real_T epsilon_j;
  real_T beta1_e;
  real_T scale_m;
  real_T absxk_m;
  real_T t_j;
  real_T beta;
  real_T qpfvalQuadExcess;
  real_T smax_f;
  real_T s_a;
  real_T qpfvalQuadExcess_tmp;
  real_T normDelta_g;
  real_T solution_lambda_n;
  real_T smax_d;
  real_T s_n;
  real_T smax_c;
  real_T s_f;
  real_T y_p;
  real_T penaltyParamTrial;
  real_T constrViolationEq;
  real_T constrViolationIneq;
  real_T TrialState_cIneq;
  real_T denomTol;
  real_T phaseOneCorrectionX;
  real_T phaseOneCorrectionP;
  real_T pk_corrected;
  real_T ratio;
  real_T b_c_p;
  real_T c_n;
  real_T tempMaxConstr_k;
  real_T c_n3;
  real_T t2_o;
  real_T t3_g;
  real_T k1_tmp;
  real_T k1_tmp_c;
  real_T k1_tmp_cj;
  real_T k1_tmp_m;
  real_T k1_tmp_j;
  real_T k1_tmp_k;
  real_T k1_tmp_mx;
  real_T constrViolation;
  real_T tol;
  real_T maxDiag;
  real_T u1_p;
  real_T b_atmp;
  real_T tau;
  real_T a_d;
  real_T b_g4;
  real_T c_c_c;
  real_T b_s;
  real_T b_temp;
  real_T roe;
  real_T absa;
  real_T absb;
  real_T scale_c;
  real_T ads;
  real_T bds;
  real_T u1_i;
  real_T c_d;
  real_T u1_g;
  real_T c_l;
  real_T alpha1;
  real_T y_f;
  real_T temp_d;
  real_T ssq;
  real_T c_j;
  real_T u1_i3;
  real_T obj_maxConstrWorkspace;
  real_T c_h;
  real_T temp_n;
  SL_Bus_std_msgs_Float64 In1_c;       // '<S54>/In1'
  SL_Bus_std_msgs_Float64 SourceBlock_o2_b;// '<S5>/SourceBlock'
  int32_T A_size[2];
  int32_T Cineq_size[2];
  int32_T JacCineqTrans_size[2];
  int32_T a__3_size[2];
  int32_T a__4_size[2];
  int32_T varargin_1_size[2];
  int32_T b_varargin_1_size[2];
  int32_T Jx_f[2];
  int32_T tmp_size[2];
  int32_T b_i[2];
  int32_T tmp_size_f[2];
  int32_T b_j;
  int32_T coffset;
  int32_T aoffset;
  int32_T b_k;
  int32_T controlOn;                   // '<Root>/Control Logic'
  int32_T i;
  int32_T mIneq;
  int32_T mConstrMax;
  int32_T maxDims;
  int32_T mLB;
  int32_T mUB;
  int32_T mFixed;
  int32_T iw0;
  int32_T iEq0;
  int32_T b_iy;
  int32_T i_k;
  int32_T mLinIneq_tmp;
  int32_T loop_ub;
  int32_T mFixed_d;
  int32_T mIneq_i;
  int32_T mLB_g;
  int32_T mUB_n;
  int32_T mConstr;
  int32_T mLinIneq;
  int32_T k;
  int32_T ix;
  int32_T h;
  int32_T o_l;
  int32_T d_ix;
  int32_T loop_ub_c;
  int32_T u1_n;
  int32_T vectorUB;
  int32_T loop_ub_p;
  int32_T y_size_idx_0;
  int32_T nVar_tmp_tmp;
  int32_T d_ix_tmp;
  int32_T h_tmp;
  int32_T row;
  int32_T col;
  int32_T col_end;
  int32_T idx_mat;
  int32_T i_d;
  int32_T loop_ub_o;
  int32_T loop_ub_j;
  int32_T i_c;
  int32_T k_h;
  int32_T b_U_tmp;
  int32_T Jx_tmp;
  int32_T c_tmp;
  int32_T i1;
  int32_T arg305_tmp;
  int32_T i2;
  int32_T k_d;
  int32_T mLinIneq_c;
  int32_T idx;
  int32_T k_p;
  int32_T scalarLB;
  int32_T vectorUB_p;
  int32_T scalarLB_tmp;
  int32_T vectorUB_tmp;
  int32_T n;
  int32_T yk;
  int32_T k_a;
  int32_T vectorUB_o;
  int32_T ineqRange_size_idx_1;
  int32_T i3;
  int32_T Umv_tmp;
  int32_T idx_current;
  int32_T Gfxu_data_tmp;
  int32_T k_j;
  int32_T k_pi;
  int32_T c_k;
  int32_T i_o;
  int32_T i4;
  int32_T vectorUB_l;
  int32_T tmp_size_idx_0;
  int32_T Jx_data_tmp;
  int32_T Jx_data_tmp_k;
  int32_T i5;
  int32_T i6;
  int32_T loop_ub_jk;
  int32_T loop_ub_f;
  int32_T loop_ub_cm;
  int32_T mc;
  int32_T coffset_f;
  int32_T boffset;
  int32_T aoffset_n;
  int32_T bkj;
  int32_T j;
  int32_T i_i;
  int32_T b_i_l;
  int32_T scalarLB_i;
  int32_T vectorUB_k;
  int32_T i7;
  int32_T b_i_f;
  int32_T i_a;
  int32_T ic_idx_0_d;
  int32_T ic_idx_1_e;
  int32_T b_i_e;
  int32_T aoffset_b;
  int32_T ii;
  int32_T knt;
  int32_T lastc;
  int32_T e_a;
  int32_T jA;
  int32_T nVar;
  int32_T k_i;
  int32_T idxStartIneq;
  int32_T idxEndIneq;
  int32_T vectorUB_f;
  int32_T nVar_j;
  int32_T mConstrMax_o;
  int32_T idxIneqOffset;
  int32_T idx_Aineq;
  int32_T idx_lower;
  int32_T idx_Partition;
  int32_T nWIneq_old;
  int32_T nWLower_old;
  int32_T nWUpper_old;
  int32_T iy;
  int32_T l;
  int32_T idxStartIneq_f;
  int32_T idxStartIneq_tmp;
  int32_T idx_lower_tmp;
  int32_T nVar_o;
  int32_T c_ln;
  int32_T ixlast;
  int32_T scalarLB_l;
  int32_T vectorUB_g;
  int32_T PROBTYPE_ORIG;
  int32_T mConstr_d;
  int32_T idxStartIneq_d;
  int32_T idxEndIneq_j;
  int32_T nVar_tmp;
  int32_T idxStartIneq_tmp_f;
  int32_T activeSetChangeID;
  int32_T nVar_js;
  int32_T globalActiveConstrIdx;
  int32_T idxMinLambda;
  int32_T k_ho;
  int32_T iQR0;
  int32_T g;
  int32_T nVar_c;
  int32_T mWConstr;
  int32_T nVar_n;
  int32_T rankQR;
  int32_T ldq;
  int32_T jBcol;
  int32_T iAcol;
  int32_T h_k;
  int32_T br;
  int32_T n_a;
  int32_T b_br;
  int32_T iQR0_f;
  int32_T b_iQR0;
  int32_T i_k_j;
  int32_T idx_k;
  int32_T ix0;
  int32_T iy0;
  int32_T b_b;
  int32_T k_hm;
  int32_T minmn;
  int32_T nfxd;
  int32_T b_j_e;
  int32_T ma_tmp;
  int32_T ma;
  int32_T minmn_h;
  int32_T ii_k;
  int32_T nmi;
  int32_T mmi;
  int32_T pvt;
  int32_T itemp;
  int32_T j_j;
  int32_T idxmax;
  int32_T knt_o;
  int32_T d;
  int32_T scalarLB_c;
  int32_T vectorUB_h;
  int32_T vectorUB_tmp_i;
  int32_T kend;
  int32_T k_pl;
  int32_T k_f;
  int32_T b_j_ew;
  int32_T nVarOrig;
  int32_T idx_max;
  int32_T mIneq_n;
  int32_T mLBOrig;
  int32_T mFiniteLBOrig;
  int32_T mIneq_tmp;
  int32_T activeSetChangeID_h;
  int32_T nVar_h;
  int32_T globalActiveConstrIdx_f;
  int32_T idxMinLambda_i;
  int32_T k_f4;
  int32_T iQR0_c;
  int32_T g_n;
  int32_T nVar_he;
  int32_T b_idx;
  int32_T nVars;
  int32_T LDimSizeP1;
  int32_T A_maxDiag_idx;
  int32_T k_k;
  int32_T order;
  int32_T c_ix;
  int32_T mNull_tmp;
  int32_T A_maxDiag_idx_tmp;
  int32_T LDimSizeP1_h;
  int32_T A_maxDiag_idx_b;
  int32_T k_o;
  int32_T order_n;
  int32_T iy0_m;
  int32_T LDimSizeP1_k;
  int32_T subRows;
  int32_T LD_diagOffset;
  int32_T subBlockSize;
  int32_T k_jk;
  int32_T ix_h;
  int32_T c_f;
  int32_T d_d;
  int32_T ia;
  int32_T e_li;
  int32_T lastColC;
  int32_T br_k;
  int32_T g_i;
  int32_T h_h;
  int32_T b_m;
  int32_T k_g;
  int32_T k_l;
  int32_T iyend;
  int32_T b_iy_m;
  int32_T f;
  int32_T g_nt;
  int32_T ia_g;
  int32_T e_tmp;
  int32_T b_d;
  int32_T b_iy_mq;
  int32_T scalarLB_f;
  int32_T vectorUB_gd;
  int32_T y_tmp;
  int32_T nVar_jc;
  int32_T ia_c;
  int32_T i_e;
  int32_T itau;
  int32_T c_m;
  int32_T b_k_o;
  int32_T scalarLB_a;
  int32_T vectorUB_j;
  int32_T lastv;
  int32_T lastc_g;
  int32_T coltop;
  int32_T ia_j;
  int32_T iy_e;
  int32_T b_iy_j;
  int32_T b_jb;
  int32_T d_g;
  int32_T ia_o;
  int32_T b_h;
  int32_T idxStartIneq_c;
  int32_T idxEndIneq_a;
  int32_T idxStartIneq_tmp_l;
  int32_T nFixedConstr;
  int32_T nDepIneq;
  int32_T idxDiag;
  int32_T b_idx_j;
  int32_T ix0_i;
  int32_T iy0_mi;
  int32_T d_f;
  int32_T idxDiag_tmp;
  int32_T lda;
  int32_T ii_o;
  int32_T mmi_i;
  int32_T i_e0;
  int32_T loop_ub_j0;
  int32_T i_o4;
  int32_T k_fr;
  int32_T idxRotGCol;
  int32_T QRk0;
  int32_T ix_m;
  int32_T b_n;
  int32_T b_ix;
  int32_T d_temp_tmp;
  int32_T QRk0_tmp;
  int32_T mIneq_a;
  int32_T b_mIneq;
  int32_T f_h;
  int32_T b_o;
  int32_T b_iy_h;
  int32_T e_j;
  int32_T scalarLB_g;
  int32_T vectorUB_jz;
  int32_T mIneq_l;
  int32_T b_mIneq_k;
  int32_T f_d;
  int32_T b_np;
  int32_T b_iy_j5;
  int32_T e_a3;
  int32_T scalarLB_h;
  int32_T vectorUB_i;
  int32_T LDimSizeP1_d;
  int32_T LD_diagOffset_b;
  int32_T subMatrixDim;
  int32_T k_hj;
  int32_T b_k_p;
  int32_T jA_n;
  int32_T j_jz;
  int32_T b_ot;
  int32_T ijA;
  int32_T alpha1_tmp;
  int32_T idx_b;
  int32_T ix0_j;
  int32_T iy0_e;
  int32_T b_k_i;
  int32_T idxA1j;
  int32_T idxAjj;
  int32_T nmj;
  int32_T b_j_n;
  int32_T ia0;
  int32_T iy_i;
  int32_T b_p;
  int32_T d_o;
  int32_T ia_m;
  int32_T b_k_ot;
  int32_T b_gz;
  int32_T idx_e;
  int32_T mIneq_iz;
  int32_T b_gb;
  int32_T k_ge;
  int32_T b_gbr;
  int32_T b_iy_g;
  int32_T e_c;
  int32_T scalarLB_k;
  int32_T vectorUB_d;
  int32_T jA_k;
  int32_T j_p;
  int32_T b_p5;
  int32_T ijA_m;
  int32_T lastColC_k;
  int32_T br_a;
  int32_T ar;
  int32_T cr;
  int32_T b_f;
  int32_T ic_c;
  int32_T c_jk;
  int32_T d_k;
  int32_T scalarLB_hj;
  int32_T vectorUB_d1;
  int32_T b_j1;
  int32_T i_n;
  int32_T ixlast_j;
  int32_T k_lc;
  int32_T iy_p;
  int32_T i_tmp;
  int32_T ix_p;
  int32_T b_iy_l;
  int32_T b_l;
  int32_T c_hb;
  int32_T ia_cg;
  int32_T i8;
  uint32_T len;
  SL_Bus_std_msgs_Int32 BusAssignment_cw;// '<S15>/Bus Assignment'
  int32_T B_size[1];
  int32_T tmp_size_o[1];
  int32_T tmp_size_c[1];
  int32_T Je_size[1];
  boolean_T b_x_b[4];
  boolean_T x_e[4];
  uint8_T sizes[2];
  uint8_T sizes_g[2];
  uint8_T Jx_e[2];
  uint8_T varargin_2[2];
  uint8_T Je_n[2];
  int8_T sizes_idx_0;
  uint8_T u;
  boolean_T SourceBlock_o1;            // '<S5>/SourceBlock'
  boolean_T SourceBlock_o1_o;          // '<S4>/SourceBlock'
  boolean_T MATLABSystem_o3_i;         // '<S18>/MATLAB System'
  boolean_T b_varargout_1;
  boolean_T unnamed_idx_0;
  boolean_T allFinite;
  boolean_T sizes_idx_1_tmp;
  boolean_T tooSmallX;
  boolean_T y_fn;
  boolean_T allFinite_n;
  boolean_T y_e;
  boolean_T checkBoundViolation;
  boolean_T subProblemChanged;
  boolean_T updateFval;
  boolean_T nonDegenerateWset;
  boolean_T tf;
  boolean_T b_tf;
  boolean_T subProblemChanged_b;
  boolean_T updateFval_a;
  boolean_T isEqAndIneqFeasible;
  boolean_T nonDegenerateWset_i;
  boolean_T okWorkingSet;
};

// Block states (default storage) for system '<Root>'
struct DW_nmpcSine_T {
  s_YHCj6MDTSR65so1fxJsjSD_nmpc_T ADdata;// '<S24>/NLMPC'
  ros_slros2_internal_block_Get_T obj; // '<Root>/Get Parameter'
  ros_slros2_internal_block_Pub_T obj_o;// '<S52>/SinkBlock'
  ros_slros2_internal_block_Pub_T obj_i;// '<S50>/SinkBlock'
  ros_slros2_internal_block_Pub_T obj_k;// '<S48>/SinkBlock'
  ros_slros2_internal_block_Pub_T obj_c;// '<S9>/SinkBlock'
  ros_slros2_internal_block_Sub_T obj_p;// '<S6>/SourceBlock'
  ros_slros2_internal_block_Sub_T obj_l;// '<S5>/SourceBlock'
  ros_slros2_internal_block_Sub_T obj_km;// '<S4>/SourceBlock'
  real_T UnitDelay_DSTATE[2];          // '<S3>/Unit Delay'
  real_T mv_Delay_DSTATE[42];          // '<S25>/mv_Delay'
  real_T x_Delay_DSTATE[126];          // '<S25>/x_Delay'
  real_T slack_delay_DSTATE;           // '<S25>/slack_delay'
  real_T P[36];                        // '<S10>/DataStoreMemory - P'
  real_T x[6];                         // '<S10>/DataStoreMemory - x'
  real_T k;                            // '<S3>/MATLAB Function1'
  real_T en1;                          // '<Root>/Control Logic'
  real_T en2;                          // '<Root>/Control Logic'
  real_T en3;                          // '<Root>/Control Logic'
  uint16_T temporalCounter_i1;         // '<Root>/Control Logic'
  uint8_T is_active_c4_nmpcSine;       // '<Root>/Control Logic'
  uint8_T is_c4_nmpcSine;              // '<Root>/Control Logic'
  boolean_T icLoad;                    // '<S25>/mv_Delay'
  boolean_T icLoad_j;                  // '<S25>/x_Delay'
  boolean_T icLoad_m;                  // '<S25>/slack_delay'
  boolean_T ADdata_not_empty;          // '<S24>/NLMPC'
};

// Constant parameters (default storage)
struct ConstP_nmpcSine_T {
  // Expression: min(3,PredictionHorizon+1):(PredictionHorizon+1)
  //  Referenced by: '<S25>/Constant'

  real_T Constant_Value_a[19];

  // Expression: 2:max(2,PredictionHorizon)
  //  Referenced by: '<S25>/Constant1'

  real_T Constant1_Value[19];

  // Expression: p.R{1}
  //  Referenced by: '<S10>/R1'

  real_T R1_Value[4];

  // Expression: [x_ref2, y_ref2]
  //  Referenced by: '<S3>/Constant5'

  real_T Constant5_Value[3002];

  // Expression: p.Q
  //  Referenced by: '<S10>/Q'

  real_T Q_Value[36];

  // Expression: p.InitialCovariance
  //  Referenced by: '<S10>/DataStoreMemory - P'

  real_T DataStoreMemoryP_InitialValue[36];
};

// Real-time Model Data Structure
struct tag_RTM_nmpcSine_T {
  const char_T * volatile errorStatus;

  //
  //  Timing:
  //  The following substructure contains information regarding
  //  the timing information for the model.

  struct {
    struct {
      uint8_T TID[3];
    } TaskCounters;
  } Timing;

  const char_T* getErrorStatus() const;
  void setErrorStatus(const char_T* const volatile aErrorStatus);
};

// Constant parameters (default storage)
extern const ConstP_nmpcSine_T nmpcSine_ConstP;

// Class declaration for model nmpcSine
class nmpcSine
{
  // public data and function members
 public:
  // Real-Time Model get method
  RT_MODEL_nmpcSine_T * getRTM();

  // model initialize function
  void initialize();

  // model step function
  void step();

  // model terminate function
  void terminate();

  // Constructor
  nmpcSine();

  // Destructor
  ~nmpcSine();

  // private data and function members
 private:
  // Block signals
  B_nmpcSine_T nmpcSine_B;

  // Block states
  DW_nmpcSine_T nmpcSine_DW;

  // private member function(s) for subsystem '<Root>'
  void nmpcSine_sind(real_T *x);
  void nmpcSine_cosd(real_T *x);
  real_T nmpcSine_xnrm2(int32_T n, const real_T x[16], int32_T ix0);
  real_T nmpcSine_rt_hypotd_snf(real_T u0, real_T u1);
  void nmpcSine_qr(const real_T A[16], real_T Q[16], real_T R[4]);
  void nmpcSine_trisolve(const real_T A[4], real_T B[12]);
  void nmpcSine_trisolve_m(const real_T A[4], real_T B[12]);
  real_T nmpcSine_xnrm2_m(int32_T n, const real_T x[48], int32_T ix0);
  void nmpcSine_qr_m(const real_T A[48], real_T Q[48], real_T R[36]);
  real_T nmpcSine_mod(real_T x, real_T y);
  void nmpcSine_trisolve_m3(real_T A, real_T B[6]);
  real_T nmpcSine_xnrm2_m3(int32_T n, const real_T x[42], int32_T ix0);
  void EKFCorrector_correctStateAndSqr(const real_T x[6], const real_T S[36],
    real_T residue, const real_T Pxy[6], real_T Sy, const real_T H[6], real_T
    Rsqrt, real_T b_x[6], real_T b_S[36]);
  real_T nmpcSine_xnrm2_m32(int32_T n, const real_T x[7], int32_T ix0);
  void EKFCorrectorAdditive_getMeasure(real_T Rs, const real_T x[6], const
    real_T S[36], real_T *zEstimated, real_T Pxy[6], real_T *Sy, real_T dHdx[6],
    real_T *Rsqrt);
  real_T nmpcSine_xnrm2_m32u(int32_T n, const real_T x[42], int32_T ix0);
  void nmpcSine_qr_m3(const real_T A[42], real_T Q[42], real_T R[36]);
  void nmpcSine_getXUe(const real_T z[125], const real_T x[6], real_T X[126],
                       real_T U[84], real_T *e);
  real_T nmpcSine_costFcn(const real_T X[126], const real_T U[84], const real_T
    data_References[80], real_T Q, real_T R, real_T Qt);
  void nmpcSine_mtimes(const real_T A_data[], const int32_T A_size[2], real_T
                       C_data[], int32_T C_size[2]);
  void nmpcSine_getUBounds(const real_T runtimedata_lastMV[2], const real_T
    runtimedata_MVMin[40], const real_T runtimedata_MVMax[40], const real_T
    runtimedata_MVRateMin[40], const real_T runtimedata_MVRateMax[40], real_T
    A_data[], int32_T A_size[2], real_T Bu_data[], int32_T Bu_size[1]);
  void nmpcSine_eye(real_T b_I[36]);
  void nm_stateTransitionFcnJacobianAD(const real_T inputVariables[10], const
    real_T extraParams[67], real_T obj[6], real_T grad[60]);
  void nmpcSine_stateEvolution(const real_T X[126], const real_T U[84], real_T
    c[120], real_T J[15000]);
  void nmpcSine_all(const boolean_T x[80], boolean_T y[4]);
  boolean_T nmpcSine_any(const boolean_T x[8]);
  void nmpcSine_reformJacobian(const real_T Jx_data[], const int32_T Jx_size[3],
    const real_T Jmv_data[], const real_T Je_data[], const int32_T Je_size[1],
    real_T Jc_data[], int32_T Jc_size[2]);
  void nmpcSine_outputBounds(const real_T runtimedata_OutputMin[80], const
    real_T runtimedata_OutputMax[80], const real_T X[126], real_T e, real_T
    c_data[], int32_T c_size[2], real_T Jc_data[], int32_T Jc_size[2]);
  void nmpcSine_c4_mpclib_anonFcn2(const real_T runtimedata_x[6], const real_T
    runtimedata_OutputMin[80], const real_T runtimedata_OutputMax[80], const
    real_T z[125], real_T varargout_1_data[], int32_T varargout_1_size[2],
    real_T varargout_2[120], real_T varargout_3_data[], int32_T
    varargout_3_size[2], real_T varargout_4[15000]);
  void nmpcSine_factoryConstruct(int32_T nVarMax, int32_T mConstrMax, int32_T
    mIneq, int32_T mNonlinIneq, s_2COE1uYisQtyPYvPjrXP9G_nmpc_T *obj);
  void nmpcSine_factoryConstruct_h53b(int32_T MaxVars, int32_T obj_grad_size[1],
    int32_T obj_Hx_size[1], boolean_T *obj_hasLinear, int32_T *obj_nvar, int32_T
    *obj_maxVar, real_T *obj_beta, real_T *obj_rho, int32_T *obj_objtype,
    int32_T *obj_prev_objtype, int32_T *obj_prev_nvar, boolean_T
    *obj_prev_hasLinear, real_T *obj_gammaScalar);
  void nmpcSine_factoryConstruct_h53bm(int32_T mIneqMax, int32_T nVarMax,
    int32_T mConstrMax, s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T *obj);
  real_T nmpcSine_costFcn_d(const real_T X[126], const real_T U[84], const
    real_T data_References[80], real_T Q, real_T R, real_T Qt);
  void computeObjectiveAndUserGradient(const s_I4XPpWw7d7shktLagLlNtD_nmpc_T
    *obj_next_next_next_next_next_ne, const real_T x[125], real_T
    grad_workspace_data[], real_T *fval, int32_T *status);
  int32_T nmpcSine_checkVectorNonFinite(int32_T N, const real_T vec_data[],
    int32_T iv0);
  int32_T nmpcSine_checkVectorNonFinite_n(const real_T vec[120]);
  int32_T computeConstraintsAndUserJacobi(int32_T
    obj_next_next_next_next_next_b_, const s_jex761Cl1dvQqVqRqjms8C_nmpc_T
    *obj_next_next_next_next_next_ne, const real_T x[125], real_T
    Cineq_workspace_data[], int32_T ineq0, real_T Ceq_workspace[120], real_T
    JacIneqTrans_workspace_data[], int32_T iJI_col, int32_T ldJI, real_T
    JacEqTrans_workspace_data[], int32_T ldJE);
  void evalObjAndConstrAndDerivatives(int32_T obj_next_next_next_next_next_b_,
    const s_jex761Cl1dvQqVqRqjms8C_nmpc_T *obj_next_next_next_next_next_ne,
    const s_I4XPpWw7d7shktLagLlNtD_nmpc_T *obj_next_next_next_next_next__0,
    const real_T x[125], real_T grad_workspace_data[], real_T
    Cineq_workspace_data[], int32_T ineq0, real_T Ceq_workspace[120], real_T
    JacIneqTrans_workspace_data[], int32_T iJI_col, int32_T ldJI, real_T
    JacEqTrans_workspace_data[], int32_T ldJE, real_T *fval, int32_T *status);
  void nmpcSin_modifyOverheadPhaseOne_(s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T *obj);
  void nmpcSine_setProblemType(s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T *obj, int32_T
    PROBLEM_TYPE);
  void nmpcSine_initActiveSet(s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T *obj);
  void nmpcSine_factoryConstruct_h5(int32_T maxRows, int32_T maxCols, int32_T
    *obj_ldq, int32_T obj_QR_size[2], real_T obj_Q_data[], int32_T obj_Q_size[2],
    int32_T obj_jpvt_data[], int32_T obj_jpvt_size[1], int32_T *obj_mrows,
    int32_T *obj_ncols, int32_T obj_tau_size[1], int32_T *obj_minRowCol,
    boolean_T *obj_usedPivoting);
  void nmpcSine_factoryConstruct_h53(int32_T MaxDims, int32_T obj_FMat_size[2],
    int32_T *obj_ldm, int32_T *obj_ndims, int32_T *obj_info, real_T
    *obj_scaleFactor, boolean_T *obj_ConvexCheck, real_T *obj_regTol_, real_T
    *obj_workspace_, real_T *obj_workspace2_);
  void nmpcSine_computeGradLag(real_T workspace_data[], int32_T ldA, int32_T
    nVar, const real_T grad_data[], int32_T mIneq, const real_T AineqTrans_data[],
    const real_T AeqTrans_data[], const int32_T finiteFixed_data[], int32_T
    mFixed, const int32_T finiteLB_data[], int32_T mLB, const int32_T
    finiteUB_data[], int32_T mUB, const real_T lambda_data[]);
  real_T nmpcSine_computePrimalFeasError(const real_T x[125], int32_T mLinIneq,
    int32_T mNonlinIneq, const real_T cIneq_data[], const real_T cEq[120], const
    int32_T finiteLB_data[], int32_T mLB, const real_T lb[125], const int32_T
    finiteUB_data[], int32_T mUB, const real_T ub[125]);
  void nmpcSine_computeDualFeasError(int32_T nVar, const real_T gradLag_data[],
    boolean_T *gradOK, real_T *val);
  void nmpcSine_test_exit(sG8JZ69axY52WWR6RKyApQC_nmpcS_T *MeritFunction, const
    s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T *WorkingSet, s_2COE1uYisQtyPYvPjrXP9G_nmpc_T
    *TrialState, const real_T lb[125], const real_T ub[125], boolean_T
    *Flags_gradOK, boolean_T *Flags_fevalOK, boolean_T *Flags_done, boolean_T
    *Flags_stepAccepted, boolean_T *Flags_failedLineSearch, int32_T
    *Flags_stepType);
  void nmpcSine_saveJacobian(s_2COE1uYisQtyPYvPjrXP9G_nmpc_T *obj, int32_T nVar,
    int32_T mIneq, const real_T JacCineqTrans_data[], int32_T ineqCol0, const
    real_T JacCeqTrans_data[], int32_T ldJ);
  real_T nmpcSine_computeComplError(const int32_T fscales_lineq_constraint_size
    [1], const int32_T fscales_cineq_constraint_size[1], const real_T xCurrent
    [125], int32_T mIneq, const real_T cIneq_data[], const int32_T
    finiteLB_data[], int32_T mLB, const real_T lb[125], const int32_T
    finiteUB_data[], int32_T mUB, const real_T ub[125], const real_T
    lambda_data[], int32_T iL0);
  void nmpcSine_computeGradLag_e(real_T workspace_data[], int32_T ldA, int32_T
    nVar, const real_T grad_data[], int32_T mIneq, const real_T AineqTrans_data[],
    const real_T AeqTrans_data[], const int32_T finiteFixed_data[], int32_T
    mFixed, const int32_T finiteLB_data[], int32_T mLB, const int32_T
    finiteUB_data[], int32_T mUB, const real_T lambda_data[]);
  void nmpcSine_computeDualFeasError_k(int32_T nVar, const real_T gradLag_data[],
    boolean_T *gradOK, real_T *val);
  void nmpcSi_updateWorkingSetForNewQP(const real_T xk[125],
    s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T *WorkingSet, int32_T mIneq, int32_T
    mNonlinIneq, const real_T cIneq_data[], const real_T cEq[120], int32_T mLB,
    const real_T lb[125], int32_T mUB, const real_T ub[125], int32_T mFixed);
  void nmpcSine_xswap(int32_T n, real_T x_data[], int32_T ix0, int32_T iy0);
  real_T nmpcSine_xnrm2_n(int32_T n, const real_T x_data[], int32_T ix0);
  real_T nmpcSine_xzlarfg(int32_T n, real_T *alpha1, real_T x_data[], int32_T
    ix0);
  void nmpcSine_xgemv(int32_T m, int32_T n, const real_T A_data[], int32_T ia0,
                      int32_T lda, const real_T x_data[], int32_T ix0, real_T
                      y_data[]);
  void nmpcSine_xgerc(int32_T m, int32_T n, real_T alpha1, int32_T ix0, const
                      real_T y_data[], real_T A_data[], int32_T ia0, int32_T lda);
  void nmpcSine_xzlarf(int32_T m, int32_T n, int32_T iv0, real_T tau, real_T
                       C_data[], int32_T ic0, int32_T ldc, real_T work_data[]);
  void nmpcSine_qrf(real_T A_data[], const int32_T A_size[2], int32_T m, int32_T
                    n, int32_T nfxd, real_T tau_data[]);
  void nmpcSine_qrpf(real_T A_data[], const int32_T A_size[2], int32_T m,
                     int32_T n, int32_T nfxd, real_T tau_data[], int32_T
                     jpvt_data[]);
  void nmpcSine_xgeqp3(real_T A_data[], const int32_T A_size[2], int32_T m,
                       int32_T n, int32_T jpvt_data[], real_T tau_data[],
                       int32_T tau_size[1]);
  void nmpcSine_factorQRE(s_0RmwrXfzGd5lqbHvgKQe2_nmpcS_T *obj, const real_T
    A_data[], int32_T mrows, int32_T ncols, int32_T ldA);
  void nmpcSine_xorgqr(int32_T m, int32_T n, int32_T k, real_T A_data[], const
                       int32_T A_size[2], int32_T lda, const real_T tau_data[]);
  void nmpcSine_sortLambdaQP(real_T lambda_data[], int32_T
    WorkingSet_nActiveConstr, const int32_T WorkingSet_sizes[5], const int32_T
    WorkingSet_isActiveIdx[6], const int32_T WorkingSet_Wid_data[], const
    int32_T WorkingSet_Wlocalidx_data[], real_T workspace_data[]);
  void nmpcSine_test_exit_n(s7RdrPWkr8UPAUyTdDJkLaG_nmpcS_T *Flags,
    s_kmYqIq13KlaOrGCTq3ShMG_nmpc_T *memspace, sG8JZ69axY52WWR6RKyApQC_nmpcS_T
    *MeritFunction, const int32_T fscales_lineq_constraint_size[1], const
    int32_T fscales_cineq_constraint_size[1], s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T
    *WorkingSet, s_2COE1uYisQtyPYvPjrXP9G_nmpc_T *TrialState,
    s_0RmwrXfzGd5lqbHvgKQe2_nmpcS_T *QRManager, const real_T lb[125], const
    real_T ub[125]);
  boolean_T nmpcSine_BFGSUpdate(int32_T nvar, real_T Bk[15625], const real_T
    sk_data[], real_T yk_data[], real_T workspace_data[]);
  void nmpcSine_factorQRE_o(s_0RmwrXfzGd5lqbHvgKQe2_nmpcS_T *obj, int32_T mrows,
    int32_T ncols);
  void nmpcSine_countsort(int32_T x_data[], int32_T xLen, int32_T
    workspace_data[], int32_T xMin, int32_T xMax);
  void nmpcSine_removeConstr(s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T *obj, int32_T
    idx_global);
  int32_T nmpcSine_RemoveDependentEq_(s_kmYqIq13KlaOrGCTq3ShMG_nmpc_T *memspace,
    s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T *workingset, s_0RmwrXfzGd5lqbHvgKQe2_nmpcS_T
    *qrmanager);
  void nmpcSine_RemoveDependentIneq_(s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T *workingset,
    s_0RmwrXfzGd5lqbHvgKQe2_nmpcS_T *qrmanager, s_kmYqIq13KlaOrGCTq3ShMG_nmpc_T *
    memspace);
  int32_T nmpcSine_rank(const real_T qrmanager_QR_data[], const int32_T
                        qrmanager_QR_size[2], int32_T qrmanager_mrows, int32_T
                        qrmanager_ncols);
  void nmpcSine_xgemv_d(int32_T m, int32_T n, const real_T A_data[], int32_T lda,
                        const real_T x_data[], real_T y_data[]);
  real_T nmpcSi_maxConstraintViolation_i(s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T *obj,
    const real_T x_data[]);
  void nmpcSine_xgemv_df(int32_T m, int32_T n, const real_T A_data[], int32_T
    lda, const real_T x_data[], int32_T ix0, real_T y_data[]);
  real_T nmpcS_maxConstraintViolation_in(s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T *obj,
    const real_T x_data[], int32_T ix0);
  boolean_T nmpcSin_feasibleX0ForWorkingSet(real_T workspace_data[], const
    int32_T workspace_size[2], real_T xCurrent_data[],
    s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T *workingset, s_0RmwrXfzGd5lqbHvgKQe2_nmpcS_T
    *qrmanager);
  void nmpcSine_RemoveDependentIneq__b(s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T
    *workingset, s_0RmwrXfzGd5lqbHvgKQe2_nmpcS_T *qrmanager,
    s_kmYqIq13KlaOrGCTq3ShMG_nmpc_T *memspace);
  void nmpcSine_xgemv_dff(int32_T m, int32_T n, const real_T A_data[], int32_T
    lda, const real_T x_data[], real_T y_data[]);
  real_T maxConstraintViolation_AMats_no(s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T *obj,
    const real_T x_data[]);
  real_T maxConstraintViolation_AMats_re(s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T *obj,
    const real_T x_data[]);
  real_T nmpc_maxConstraintViolation_ini(s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T *obj,
    const real_T x_data[]);
  void nmpcSine_PresolveWorkingSet(s_2COE1uYisQtyPYvPjrXP9G_nmpc_T *solution,
    s_kmYqIq13KlaOrGCTq3ShMG_nmpc_T *memspace, s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T
    *workingset, s_0RmwrXfzGd5lqbHvgKQe2_nmpcS_T *qrmanager);
  void nmpcSine_xgemv_dffo(int32_T m, int32_T n, const real_T A[15625], int32_T
    lda, const real_T x_data[], real_T y_data[]);
  void nmpcSine_computeGrad_StoreHx(s_xtSBzQGTZuMYOTjcuMqLQH_nmpc_T *obj, const
    real_T H[15625], const real_T f_data[], const real_T x_data[]);
  real_T nmpcSine_computeFval_ReuseHx(const s_xtSBzQGTZuMYOTjcuMqLQH_nmpc_T *obj,
    real_T workspace_data[], const real_T f_data[], const real_T x_data[]);
  void nmpcSine_xgeqrf(real_T A_data[], const int32_T A_size[2], int32_T m,
                       int32_T n, real_T tau_data[], int32_T tau_size[1]);
  void nmpcSine_factorQR(s_0RmwrXfzGd5lqbHvgKQe2_nmpcS_T *obj, const real_T
    A_data[], int32_T mrows, int32_T ncols, int32_T ldA);
  void nmpcSine_xrotg(real_T *a, real_T *b, real_T *c, real_T *s);
  void nmpcSine_squareQ_appendCol(s_0RmwrXfzGd5lqbHvgKQe2_nmpcS_T *obj, const
    real_T vec_data[], int32_T iv0);
  void nmpcSine_deleteColMoveEnd(s_0RmwrXfzGd5lqbHvgKQe2_nmpcS_T *obj, int32_T
    idx);
  boolean_T nmpcSine_strcmp(const char_T a[7]);
  void nmpcSine_xgemm(int32_T m, int32_T n, int32_T k, const real_T A[15625],
                      int32_T lda, const real_T B_data[], int32_T ib0, int32_T
                      ldb, real_T C_data[], int32_T ldc);
  void nmpcSine_xgemm_k(int32_T m, int32_T n, int32_T k, const real_T A_data[],
                        int32_T ia0, int32_T lda, const real_T B_data[], int32_T
                        ldb, real_T C_data[], int32_T ldc);
  void nmpcSine_fullColLDL2_(s_mDApTYzDBpxuvxemclsuEF_nmpc_T *obj, int32_T
    LD_offset, int32_T NColsRemain);
  void nmpcSine_partialColLDL3_(s_mDApTYzDBpxuvxemclsuEF_nmpc_T *obj, int32_T
    LD_offset, int32_T NColsRemain);
  int32_T nmpcSine_xpotrf(int32_T n, real_T A_data[], int32_T lda);
  void nmpcSine_xgemv_dffoy(int32_T m, int32_T n, const real_T A_data[], int32_T
    ia0, int32_T lda, const real_T x_data[], real_T y_data[]);
  void nmpcSine_factor_g(s_mDApTYzDBpxuvxemclsuEF_nmpc_T *obj, const real_T A
    [15625], int32_T ndims, int32_T ldA);
  void nmpcSine_factor(s_mDApTYzDBpxuvxemclsuEF_nmpc_T *obj, const real_T A
                       [15625], int32_T ndims, int32_T ldA);
  void nmpcSine_solve_i(const s_mDApTYzDBpxuvxemclsuEF_nmpc_T *obj, real_T
                        rhs_data[]);
  void nmpcSine_solve(const s_mDApTYzDBpxuvxemclsuEF_nmpc_T *obj, real_T
                      rhs_data[]);
  void nmpcSine_compute_deltax(const real_T H[15625],
    s_2COE1uYisQtyPYvPjrXP9G_nmpc_T *solution, s_kmYqIq13KlaOrGCTq3ShMG_nmpc_T
    *memspace, const s_0RmwrXfzGd5lqbHvgKQe2_nmpcS_T *qrmanager,
    s_mDApTYzDBpxuvxemclsuEF_nmpc_T *cholmanager, const
    s_xtSBzQGTZuMYOTjcuMqLQH_nmpc_T *objective, boolean_T alwaysPositiveDef);
  real_T nmpcSine_xnrm2_nl(int32_T n, const real_T x_data[]);
  void nmpcSine_xgemv_dffoyl(int32_T m, int32_T n, const real_T A_data[],
    int32_T lda, const real_T x_data[], real_T y_data[]);
  void nmpcSine_feasibleratiotest(const real_T solution_xstar_data[], const
    real_T solution_searchDir_data[], real_T workspace_data[], const int32_T
    workspace_size[2], int32_T workingset_nVar, int32_T workingset_ldA, const
    real_T workingset_Aineq_data[], const real_T workingset_bineq_data[], const
    real_T workingset_lb_data[], const real_T workingset_ub_data[], const
    int32_T workingset_indexLB_data[], const int32_T workingset_indexUB_data[],
    const int32_T workingset_sizes[5], const int32_T workingset_isActiveIdx[6],
    const boolean_T workingset_isActiveConstr_data[], const int32_T
    workingset_nWConstr[5], boolean_T isPhaseOne, real_T *alpha, boolean_T
    *newBlocking, int32_T *constrType, int32_T *constrIdx);
  void nmpcSi_checkUnboundedOrIllPosed(s_2COE1uYisQtyPYvPjrXP9G_nmpc_T *solution,
    const s_xtSBzQGTZuMYOTjcuMqLQH_nmpc_T *objective);
  void nmpc_addBoundToActiveSetMatrix_(s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T *obj,
    int32_T TYPE, int32_T idx_local);
  void nmpcSine_addAineqConstr(s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T *obj, int32_T
    idx_local);
  void nmpcSine_compute_lambda(real_T workspace_data[],
    s_2COE1uYisQtyPYvPjrXP9G_nmpc_T *solution, const
    s_xtSBzQGTZuMYOTjcuMqLQH_nmpc_T *objective, const
    s_0RmwrXfzGd5lqbHvgKQe2_nmpcS_T *qrmanager);
  void nm_checkStoppingAndUpdateFval_f(int32_T *activeSetChangeID,
    s_2COE1uYisQtyPYvPjrXP9G_nmpc_T *solution, s_kmYqIq13KlaOrGCTq3ShMG_nmpc_T
    *memspace, const s_xtSBzQGTZuMYOTjcuMqLQH_nmpc_T *objective,
    s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T *workingset, s_0RmwrXfzGd5lqbHvgKQe2_nmpcS_T
    *qrmanager, int32_T runTimeOptions_MaxIterations, boolean_T *updateFval);
  void nmpcSine_iterate_d(const real_T H[15625], const real_T f_data[],
    s_2COE1uYisQtyPYvPjrXP9G_nmpc_T *solution, s_kmYqIq13KlaOrGCTq3ShMG_nmpc_T
    *memspace, s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T *workingset,
    s_0RmwrXfzGd5lqbHvgKQe2_nmpcS_T *qrmanager, s_mDApTYzDBpxuvxemclsuEF_nmpc_T *
    cholmanager, s_xtSBzQGTZuMYOTjcuMqLQH_nmpc_T *objective, const char_T
    options_SolverName[7], int32_T runTimeOptions_MaxIterations);
  void nmpc_checkStoppingAndUpdateFval(int32_T *activeSetChangeID, const real_T
    f_data[], s_2COE1uYisQtyPYvPjrXP9G_nmpc_T *solution,
    s_kmYqIq13KlaOrGCTq3ShMG_nmpc_T *memspace, const
    s_xtSBzQGTZuMYOTjcuMqLQH_nmpc_T *objective, s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T *
    workingset, s_0RmwrXfzGd5lqbHvgKQe2_nmpcS_T *qrmanager, int32_T
    runTimeOptions_MaxIterations, const boolean_T *updateFval);
  void nmpcSine_iterate(const real_T H[15625], const real_T f_data[],
                        s_2COE1uYisQtyPYvPjrXP9G_nmpc_T *solution,
                        s_kmYqIq13KlaOrGCTq3ShMG_nmpc_T *memspace,
                        s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T *workingset,
                        s_0RmwrXfzGd5lqbHvgKQe2_nmpcS_T *qrmanager,
                        s_mDApTYzDBpxuvxemclsuEF_nmpc_T *cholmanager,
                        s_xtSBzQGTZuMYOTjcuMqLQH_nmpc_T *objective, const char_T
                        options_SolverName[7], int32_T
                        runTimeOptions_MaxIterations);
  void nmpcSine_phaseone(const real_T H[15625], const real_T f_data[],
    s_2COE1uYisQtyPYvPjrXP9G_nmpc_T *solution, s_kmYqIq13KlaOrGCTq3ShMG_nmpc_T
    *memspace, s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T *workingset,
    s_0RmwrXfzGd5lqbHvgKQe2_nmpcS_T *qrmanager, s_mDApTYzDBpxuvxemclsuEF_nmpc_T *
    cholmanager, s_xtSBzQGTZuMYOTjcuMqLQH_nmpc_T *objective, const char_T
    options_SolverName[7], const somzaGboVhDG7PNQS6E98jD_nmpcS_T *runTimeOptions);
  void nmpcSine_linearForm_(boolean_T obj_hasLinear, int32_T obj_nvar, real_T
    workspace_data[], const real_T H[15625], const real_T f_data[], const real_T
    x_data[]);
  void nmpcSine_driver_o(const real_T H[15625], const real_T f_data[],
    s_2COE1uYisQtyPYvPjrXP9G_nmpc_T *solution, s_kmYqIq13KlaOrGCTq3ShMG_nmpc_T
    *memspace, s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T *workingset,
    s_0RmwrXfzGd5lqbHvgKQe2_nmpcS_T *qrmanager, s_mDApTYzDBpxuvxemclsuEF_nmpc_T *
    cholmanager, s_xtSBzQGTZuMYOTjcuMqLQH_nmpc_T *objective, const
    somzaGboVhDG7PNQS6E98jD_nmpcS_T *options, const
    somzaGboVhDG7PNQS6E98jD_nmpcS_T *runTimeOptions);
  void nmpcSine_addAeqConstr(s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T *obj, int32_T
    idx_local);
  boolean_T nmpcSine_soc(const real_T Hessian[15625], const real_T grad_data[],
    s_2COE1uYisQtyPYvPjrXP9G_nmpc_T *TrialState, s_kmYqIq13KlaOrGCTq3ShMG_nmpc_T
    *memspace, s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T *WorkingSet,
    s_0RmwrXfzGd5lqbHvgKQe2_nmpcS_T *QRManager, s_mDApTYzDBpxuvxemclsuEF_nmpc_T *
    CholManager, s_xtSBzQGTZuMYOTjcuMqLQH_nmpc_T *QPObjective, const
    somzaGboVhDG7PNQS6E98jD_nmpcS_T *qpoptions);
  real_T nmpcSine_maxConstraintViolation(const s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T
    *obj, const real_T x_data[]);
  void nmpcSine_normal(const real_T Hessian[15625], const real_T grad_data[],
                       s_2COE1uYisQtyPYvPjrXP9G_nmpc_T *TrialState,
                       sG8JZ69axY52WWR6RKyApQC_nmpcS_T *MeritFunction,
                       s_kmYqIq13KlaOrGCTq3ShMG_nmpc_T *memspace,
                       s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T *WorkingSet,
                       s_0RmwrXfzGd5lqbHvgKQe2_nmpcS_T *QRManager,
                       s_mDApTYzDBpxuvxemclsuEF_nmpc_T *CholManager,
                       s_xtSBzQGTZuMYOTjcuMqLQH_nmpc_T *QPObjective, const
                       somzaGboVhDG7PNQS6E98jD_nmpcS_T *qpoptions,
                       s7RdrPWkr8UPAUyTdDJkLaG_nmpcS_T *stepFlags);
  void nmpcSine_relaxed(const real_T Hessian[15625], const real_T grad_data[],
                        s_2COE1uYisQtyPYvPjrXP9G_nmpc_T *TrialState,
                        sG8JZ69axY52WWR6RKyApQC_nmpcS_T *MeritFunction,
                        s_kmYqIq13KlaOrGCTq3ShMG_nmpc_T *memspace,
                        s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T *WorkingSet,
                        s_0RmwrXfzGd5lqbHvgKQe2_nmpcS_T *QRManager,
                        s_mDApTYzDBpxuvxemclsuEF_nmpc_T *CholManager,
                        s_xtSBzQGTZuMYOTjcuMqLQH_nmpc_T *QPObjective,
                        somzaGboVhDG7PNQS6E98jD_nmpcS_T *qpoptions);
  void nmpcSine_step_k(s7RdrPWkr8UPAUyTdDJkLaG_nmpcS_T *stepFlags, real_T
                       Hessian[15625], const real_T lb[125], const real_T ub[125],
                       s_2COE1uYisQtyPYvPjrXP9G_nmpc_T *TrialState,
                       sG8JZ69axY52WWR6RKyApQC_nmpcS_T *MeritFunction,
                       s_kmYqIq13KlaOrGCTq3ShMG_nmpc_T *memspace,
                       s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T *WorkingSet,
                       s_0RmwrXfzGd5lqbHvgKQe2_nmpcS_T *QRManager,
                       s_mDApTYzDBpxuvxemclsuEF_nmpc_T *CholManager,
                       s_xtSBzQGTZuMYOTjcuMqLQH_nmpc_T *QPObjective,
                       somzaGboVhDG7PNQS6E98jD_nmpcS_T *qpoptions);
  void nmpcSine_evalObjAndConstr(int32_T obj_next_next_next_next_next_b_, const
    s_jex761Cl1dvQqVqRqjms8C_nmpc_T *obj_next_next_next_next_next_ne, const
    s_I4XPpWw7d7shktLagLlNtD_nmpc_T *obj_next_next_next_next_next__0, const
    real_T x[125], real_T Cineq_workspace_data[], int32_T ineq0, real_T
    Ceq_workspace[120], real_T *fval, int32_T *status);
  void nmpcSine_computeLinearResiduals(const real_T x[125], int32_T nVar, real_T
    workspaceIneq_data[], const int32_T workspaceIneq_size[1], int32_T mLinIneq,
    const real_T AineqT_data[], const real_T bineq_data[], int32_T ldAi);
  real_T nmpcSine_computeMeritFcn(real_T obj_penaltyParam, real_T fval, const
    real_T Cineq_workspace_data[], int32_T mIneq, const real_T Ceq_workspace[120],
    boolean_T evalWellDefined);
  void nmpcSine_linesearch(boolean_T *evalWellDefined, const real_T bineq_data[],
    int32_T WorkingSet_nVar, int32_T WorkingSet_ldA, const real_T
    WorkingSet_Aineq_data[], s_2COE1uYisQtyPYvPjrXP9G_nmpc_T *TrialState, real_T
    MeritFunction_penaltyParam, real_T MeritFunction_phi, real_T
    MeritFunction_phiPrimePlus, real_T MeritFunction_phiFullStep, int32_T
    FcnEvaluator_next_next_next_nex, const s_jex761Cl1dvQqVqRqjms8C_nmpc_T
    *FcnEvaluator_next_next_next_n_0, const s_I4XPpWw7d7shktLagLlNtD_nmpc_T
    *FcnEvaluator_next_next_next_n_1, boolean_T socTaken, real_T *alpha, int32_T
    *exitflag);
  void nmpcSine_driver(const real_T bineq_data[], const real_T lb[125], const
                       real_T ub[125], s_2COE1uYisQtyPYvPjrXP9G_nmpc_T
                       *TrialState, sG8JZ69axY52WWR6RKyApQC_nmpcS_T
                       *MeritFunction, const coder_internal_stickyStruct_2_T
                       *FcnEvaluator, s_kmYqIq13KlaOrGCTq3ShMG_nmpc_T *memspace,
                       s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T *WorkingSet,
                       s_0RmwrXfzGd5lqbHvgKQe2_nmpcS_T *QRManager,
                       s_mDApTYzDBpxuvxemclsuEF_nmpc_T *CholManager,
                       s_xtSBzQGTZuMYOTjcuMqLQH_nmpc_T *QPObjective, const
                       int32_T fscales_lineq_constraint_size[1], const int32_T
                       fscales_cineq_constraint_size[1], real_T Hessian[15625]);
  void nmpcSine_fmincon(const s_jex761Cl1dvQqVqRqjms8C_nmpc_T
                        *fun_workspace_runtimedata, const
                        sAc4bxvmmjmxjQV9i3feLrE_nmpcS_T *fun_workspace_userdata,
                        const real_T x0[125], const real_T Aineq_data[], const
                        real_T bineq_data[], const int32_T bineq_size[1], const
                        real_T lb[125], const real_T ub[125], const
                        s_jex761Cl1dvQqVqRqjms8C_nmpc_T
                        *nonlcon_workspace_runtimedata, real_T x[125], real_T
                        *fval, real_T *exitflag, sttYSJM5GCi2c1Eu0R50efC_nmpcS_T
                        *output);
  void nmpcSine_stateTransitionFcnDT(const real_T xk[6], const real_T u[3],
    real_T xk1[6]);
  real_T nmpcSine_xnrm2_m32uc(int32_T n, const real_T x[72], int32_T ix0);
  void nmpcSine_qr_m32(const real_T A[72], real_T Q[72], real_T R[36]);
  void nmpcSine_Publisher_setupImpl_m(const ros_slros2_internal_block_Pub_T *obj);
  void nmpcSine_Publisher_setupImpl_m3(const ros_slros2_internal_block_Pub_T
    *obj);
  void nmpcSin_Publisher_setupImpl_m32(const ros_slros2_internal_block_Pub_T
    *obj);
  void nmpcSine_Subscriber_setupImpl(const ros_slros2_internal_block_Sub_T *obj);
  void nmpcSine_Subscriber_setupImpl_m(const ros_slros2_internal_block_Sub_T
    *obj);
  void nmpcSin_Subscriber_setupImpl_m3(const ros_slros2_internal_block_Sub_T
    *obj);
  void nmpcSine_Publisher_setupImpl(const ros_slros2_internal_block_Pub_T *obj);

  // Real-Time Model
  RT_MODEL_nmpcSine_T nmpcSine_M;
};

extern volatile boolean_T stopRequested;
extern volatile boolean_T runModel;

//-
//  These blocks were eliminated from the model due to optimizations:
//
//  Block '<S10>/checkMeasurementFcn1Signals' : Unused code path elimination
//  Block '<S10>/checkMeasurementFcn2Signals' : Unused code path elimination
//  Block '<S10>/checkMeasurementFcn3Signals' : Unused code path elimination
//  Block '<S10>/checkStateTransitionFcnSignals' : Unused code path elimination
//  Block '<S26>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S27>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S28>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S29>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S30>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S31>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S32>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S33>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S34>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S35>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S36>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S37>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S38>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S39>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S40>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S41>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S42>/Vector Dimension Check' : Unused code path elimination
//  Block '<S43>/Vector Dimension Check' : Unused code path elimination
//  Block '<S44>/Vector Dimension Check' : Unused code path elimination
//  Block '<S45>/Vector Dimension Check' : Unused code path elimination
//  Block '<S14>/mv.init_zero' : Unused code path elimination
//  Block '<S14>/x.init_zero' : Unused code path elimination
//  Block '<S10>/DataTypeConversion_Enable1' : Eliminate redundant data type conversion
//  Block '<S10>/DataTypeConversion_Enable2' : Eliminate redundant data type conversion
//  Block '<S10>/DataTypeConversion_Enable3' : Eliminate redundant data type conversion
//  Block '<S10>/DataTypeConversion_Q' : Eliminate redundant data type conversion
//  Block '<S10>/DataTypeConversion_R1' : Eliminate redundant data type conversion
//  Block '<S10>/DataTypeConversion_R2' : Eliminate redundant data type conversion
//  Block '<S10>/DataTypeConversion_R3' : Eliminate redundant data type conversion
//  Block '<S10>/DataTypeConversion_uMeas1' : Eliminate redundant data type conversion
//  Block '<S10>/DataTypeConversion_uMeas2' : Eliminate redundant data type conversion
//  Block '<S10>/DataTypeConversion_uMeas3' : Eliminate redundant data type conversion
//  Block '<S10>/DataTypeConversion_uState' : Eliminate redundant data type conversion
//  Block '<S10>/DataTypeConversion_y1' : Eliminate redundant data type conversion
//  Block '<S10>/DataTypeConversion_y2' : Eliminate redundant data type conversion
//  Block '<S10>/DataTypeConversion_y3' : Eliminate redundant data type conversion
//  Block '<S24>/Reshape' : Reshape block reduction
//  Block '<S24>/Reshape1' : Reshape block reduction
//  Block '<S24>/mo or x Conversion' : Eliminate redundant data type conversion
//  Block '<S24>/mo or x Conversion1' : Eliminate redundant data type conversion
//  Block '<S24>/mo or x Conversion10' : Eliminate redundant data type conversion
//  Block '<S24>/mo or x Conversion11' : Eliminate redundant data type conversion
//  Block '<S24>/mo or x Conversion12' : Eliminate redundant data type conversion
//  Block '<S24>/mo or x Conversion13' : Eliminate redundant data type conversion
//  Block '<S24>/mo or x Conversion14' : Eliminate redundant data type conversion
//  Block '<S24>/mo or x Conversion15' : Eliminate redundant data type conversion
//  Block '<S24>/mo or x Conversion16' : Eliminate redundant data type conversion
//  Block '<S24>/mo or x Conversion17' : Eliminate redundant data type conversion
//  Block '<S24>/mo or x Conversion18' : Eliminate redundant data type conversion
//  Block '<S24>/mo or x Conversion19' : Eliminate redundant data type conversion
//  Block '<S24>/mo or x Conversion2' : Eliminate redundant data type conversion
//  Block '<S24>/mo or x Conversion3' : Eliminate redundant data type conversion
//  Block '<S24>/mo or x Conversion4' : Eliminate redundant data type conversion
//  Block '<S24>/mo or x Conversion5' : Eliminate redundant data type conversion
//  Block '<S24>/mo or x Conversion6' : Eliminate redundant data type conversion
//  Block '<S24>/mo or x Conversion7' : Eliminate redundant data type conversion
//  Block '<S24>/mo or x Conversion8' : Eliminate redundant data type conversion
//  Block '<S24>/mo or x Conversion9' : Eliminate redundant data type conversion
//  Block '<S25>/reshape_mv' : Reshape block reduction
//  Block '<S25>/reshape_x' : Reshape block reduction
//  Block '<S17>/Reshape' : Reshape block reduction
//  Block '<S3>/Zero-Order Hold1' : Eliminated since input and output rates are identical
//  Block '<S3>/Zero-Order Hold2' : Eliminated since input and output rates are identical


//-
//  The generated code includes comments that allow you to trace directly
//  back to the appropriate location in the model.  The basic format
//  is <system>/block_name, where system is the system number (uniquely
//  assigned by Simulink) and block_name is the name of the block.
//
//  Use the MATLAB hilite_system command to trace the generated code back
//  to the model.  For example,
//
//  hilite_system('<S3>')    - opens system 3
//  hilite_system('<S3>/Kp') - opens and selects block Kp which resides in S3
//
//  Here is the system hierarchy for this model
//
//  '<Root>' : 'nmpcSine'
//  '<S1>'   : 'nmpcSine/Command Velocity Publisher'
//  '<S2>'   : 'nmpcSine/Control Logic'
//  '<S3>'   : 'nmpcSine/NMPC Controller'
//  '<S4>'   : 'nmpcSine/Subscribe'
//  '<S5>'   : 'nmpcSine/Subscribe1'
//  '<S6>'   : 'nmpcSine/Subscribe2'
//  '<S7>'   : 'nmpcSine/Command Velocity Publisher/Blank Message'
//  '<S8>'   : 'nmpcSine/Command Velocity Publisher/ConvertToCmdVel'
//  '<S9>'   : 'nmpcSine/Command Velocity Publisher/Publish'
//  '<S10>'  : 'nmpcSine/NMPC Controller/Extended Kalman Filter'
//  '<S11>'  : 'nmpcSine/NMPC Controller/HeadingToPsi'
//  '<S12>'  : 'nmpcSine/NMPC Controller/LatLonToXY'
//  '<S13>'  : 'nmpcSine/NMPC Controller/MATLAB Function1'
//  '<S14>'  : 'nmpcSine/NMPC Controller/Nonlinear MPC Controller'
//  '<S15>'  : 'nmpcSine/NMPC Controller/Optimization Status Publisher'
//  '<S16>'  : 'nmpcSine/NMPC Controller/Ref path publisher'
//  '<S17>'  : 'nmpcSine/NMPC Controller/State Publisher'
//  '<S18>'  : 'nmpcSine/NMPC Controller/Extended Kalman Filter/Correct1'
//  '<S19>'  : 'nmpcSine/NMPC Controller/Extended Kalman Filter/Correct2'
//  '<S20>'  : 'nmpcSine/NMPC Controller/Extended Kalman Filter/Correct3'
//  '<S21>'  : 'nmpcSine/NMPC Controller/Extended Kalman Filter/Output'
//  '<S22>'  : 'nmpcSine/NMPC Controller/Extended Kalman Filter/Predict'
//  '<S23>'  : 'nmpcSine/NMPC Controller/Extended Kalman Filter/Output/MATLAB Function'
//  '<S24>'  : 'nmpcSine/NMPC Controller/Nonlinear MPC Controller/MPC'
//  '<S25>'  : 'nmpcSine/NMPC Controller/Nonlinear MPC Controller/xmvs_router'
//  '<S26>'  : 'nmpcSine/NMPC Controller/Nonlinear MPC Controller/MPC/MPC Preview Signal Check'
//  '<S27>'  : 'nmpcSine/NMPC Controller/Nonlinear MPC Controller/MPC/MPC Preview Signal Check1'
//  '<S28>'  : 'nmpcSine/NMPC Controller/Nonlinear MPC Controller/MPC/MPC Preview Signal Check10'
//  '<S29>'  : 'nmpcSine/NMPC Controller/Nonlinear MPC Controller/MPC/MPC Preview Signal Check11'
//  '<S30>'  : 'nmpcSine/NMPC Controller/Nonlinear MPC Controller/MPC/MPC Preview Signal Check12'
//  '<S31>'  : 'nmpcSine/NMPC Controller/Nonlinear MPC Controller/MPC/MPC Preview Signal Check13'
//  '<S32>'  : 'nmpcSine/NMPC Controller/Nonlinear MPC Controller/MPC/MPC Preview Signal Check14'
//  '<S33>'  : 'nmpcSine/NMPC Controller/Nonlinear MPC Controller/MPC/MPC Preview Signal Check15'
//  '<S34>'  : 'nmpcSine/NMPC Controller/Nonlinear MPC Controller/MPC/MPC Preview Signal Check16'
//  '<S35>'  : 'nmpcSine/NMPC Controller/Nonlinear MPC Controller/MPC/MPC Preview Signal Check3'
//  '<S36>'  : 'nmpcSine/NMPC Controller/Nonlinear MPC Controller/MPC/MPC Preview Signal Check4'
//  '<S37>'  : 'nmpcSine/NMPC Controller/Nonlinear MPC Controller/MPC/MPC Preview Signal Check5'
//  '<S38>'  : 'nmpcSine/NMPC Controller/Nonlinear MPC Controller/MPC/MPC Preview Signal Check6'
//  '<S39>'  : 'nmpcSine/NMPC Controller/Nonlinear MPC Controller/MPC/MPC Preview Signal Check7'
//  '<S40>'  : 'nmpcSine/NMPC Controller/Nonlinear MPC Controller/MPC/MPC Preview Signal Check8'
//  '<S41>'  : 'nmpcSine/NMPC Controller/Nonlinear MPC Controller/MPC/MPC Preview Signal Check9'
//  '<S42>'  : 'nmpcSine/NMPC Controller/Nonlinear MPC Controller/MPC/MPC Scalar Signal Check1'
//  '<S43>'  : 'nmpcSine/NMPC Controller/Nonlinear MPC Controller/MPC/MPC Scalar Signal Check2'
//  '<S44>'  : 'nmpcSine/NMPC Controller/Nonlinear MPC Controller/MPC/MPC Vector Signal Check1'
//  '<S45>'  : 'nmpcSine/NMPC Controller/Nonlinear MPC Controller/MPC/MPC Vector Signal Check11'
//  '<S46>'  : 'nmpcSine/NMPC Controller/Nonlinear MPC Controller/MPC/NLMPC'
//  '<S47>'  : 'nmpcSine/NMPC Controller/Optimization Status Publisher/Blank Message'
//  '<S48>'  : 'nmpcSine/NMPC Controller/Optimization Status Publisher/Publish'
//  '<S49>'  : 'nmpcSine/NMPC Controller/Ref path publisher/Blank Message'
//  '<S50>'  : 'nmpcSine/NMPC Controller/Ref path publisher/Publish'
//  '<S51>'  : 'nmpcSine/NMPC Controller/State Publisher/Blank Message'
//  '<S52>'  : 'nmpcSine/NMPC Controller/State Publisher/Publish'
//  '<S53>'  : 'nmpcSine/Subscribe/Enabled Subsystem'
//  '<S54>'  : 'nmpcSine/Subscribe1/Enabled Subsystem'
//  '<S55>'  : 'nmpcSine/Subscribe2/Enabled Subsystem'

#endif                                 // nmpcSine_h_

//
// File trailer for generated code.
//
// [EOF]
//
