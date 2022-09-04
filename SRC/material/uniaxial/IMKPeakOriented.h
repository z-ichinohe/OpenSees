/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

//Modified Ibarra-Medina-Krawinkler with Peak-Oriented Hysteretic Response

//**********************************************************************                                                                     
// Code Developed by: Ahmed Elkady and Hammad ElJisr
// Last Updated: July 2020
//**********************************************************************

#ifndef IMKPeakOriented_h
#define IMKPeakOriented_h

#include <UniaxialMaterial.h>

class IMKPeakOriented : public UniaxialMaterial
{
public:
    IMKPeakOriented(int tag, double Ke,
        double Uy_pos, double Umax_pos, double Uu_pos, double Fy_pos, double FmaxFy_pos, double ResF_pos,
        double Uy_neg, double Umax_neg, double Uu_neg, double Fy_neg, double FmaxFy_neg, double ResF_neg,
        double LAMBDA_S, double LAMBDA_C, double LAMBDA_A, double LAMBDA_K, double c_S, double c_C, double c_A, double c_K, double D_pos, double D_neg);
    IMKPeakOriented();
    ~IMKPeakOriented();
    const char *getClassType(void) const { return "IMKPeakOriented"; };
    int setTrialStrain(double strain, double strainRate = 0.0);
    double getStrain(void);
    double getStress(void);
    double getTangent(void);
    double getInitialTangent(void);
    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);
    UniaxialMaterial *getCopy(void);
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    void Print(OPS_Stream &s, int flag = 0);


protected:

private:
    //my functions

    //Fixed input material parameters
    double Ke;
    double Up_pos;
    double Upc_pos;
    double Uu_pos;
    double Fy_pos;
    double FcapFy_pos;
    double ResF_pos;
    double Up_neg;
    double Upc_neg;
    double Uu_neg;
    double Fy_neg;
    double FcapFy_neg;
    double ResF_neg;
    double LAMBDA_S;
    double LAMBDA_C;
    double LAMBDA_A;
    double LAMBDA_K;
    double c_S;
    double c_C;
    double c_A;
    double c_K;
    double D_pos;
    double D_neg;
    
    //State variables 
    double U,               cU;

    //History variables 

    double Ui,              cUi;
    double Fi,              cFi;
    double Ui_1,            cUi_1;
    double Fi_1,            cFi_1;
    double du_i_1,          cdu_i_1;

    double posUy,           cPosUy;
    double posUcap,         cPosUcap;
    double posFy,           cPosFy;
    double posFcap,         cPosFcap;
    double posUglobal,      cPosUglobal;
    double posFglobal,      cPosFglobal;
    double posUres,         cPosUres;
    double posFres,         cPosFres;
    double posKp,           cPosKp;
    double posKpc,          cPosKpc;


    double negUy,           cNegUy;
    double negUcap,         cNegUcap;
    double negFy,           cNegFy;
    double negFcap,         cNegFcap;
    double negUglobal,      cNegUglobal;
    double negFglobal,      cNegFglobal;

    double negUres,         cNegUres;
    double negFres,         cNegFres;
    double negKp,           cNegKp;
    double negKpc,          cNegKpc;

    double K_unload,        cK_unload;

    double engAcml,         cEngAcml;
    double engDspt,         cEngDspt;
    
    double u0, cu0;

    double dU;
    double dF;

    bool   FailS;
    bool   FailC;
    bool   FailA;
    bool   FailK;
    bool   FailPp;
    bool   FailPn;
    bool   FailRp;
    bool   FailRn;
    bool   FailDp;
    bool   FailDn;

    double Ei;
    double dEi;
    double Epj;
    double EpjK;
    double EiK;

    double c_cS;
    double c_cC;
    double c_cA;
    double c_cK;

    double EtS;
    double EtC;
    double EtA;
    double EtK;

    double betaS;
    double betaC;
    double betaA;
    double betaK;

    double sPCsp,           sPCsn;
    double sPCpcp,          sPCpcn;

    double TangentK,        cTangentK, ki;

    double Uy_pos,          Uy_neg;
    double Ucap_pos,        Ucap_neg;
    double Fcap_pos,        Fcap_neg;
    double Kpc_pos,         Kpc_neg;
    double Kp_pos,          Kp_neg;

    double posUlocal,     cPosUlocal;
    double posFlocal,     cPosFlocal;
    double negUlocal,     cNegUlocal;
    double negFlocal,     cNegFlocal;

    bool   Failure_Flag,   cFailure_Flag;
    bool   Excursion_Flag, cExcursion_Flag;
    int    Branch,         cBranch;
    bool   Yield_Flag,     cYield_Flag;
    bool   Reversal_Flag,  cReversal_Flag;
    int    exBranch,       cExBranch;

    double K_reload,       cK_reload;

    double K_local;
    double K_global;

};

#endif
