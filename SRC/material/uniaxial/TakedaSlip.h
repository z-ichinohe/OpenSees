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
** redistribution, and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

// Takeda-Slip Model

//**********************************************************************
// Code Developed by: Kazuki Ichinohe
// Last Updated: September 2023
//**********************************************************************

#ifndef TakedaSlip_h
#define TakedaSlip_h

#include <UniaxialMaterial.h>
#include <array>
#include <tuple>

class TakedaSlip : public UniaxialMaterial
{
public:
    TakedaSlip(int tag,
        double Uc, double Uy,
        double Kc, double Ky, double Kp,
        double unload_from_global_factor, double unload_from_local_factor);
    TakedaSlip();
    ~TakedaSlip();
    const char *getClassType(void) const { return "TakedaSlip"; };
    int setTrialStrain(double strain, double strainRate = 0.0);
    double  getStrain(void);
    double  getStress(void);
    double  getTangent(void);
    double  getInitialTangent(void);
    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);
    UniaxialMaterial *getCopy(void);
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    void Print(OPS_Stream &s, int flag = 0);
protected:

private:
    int branch, cbranch;

    double k_tangent, ck_tangent;
    double k_unload, ck_unload;

    double f_new, cf_new;
    double d_new, cd_new;

    std::array<double, 3> f_global, cf_global;
    std::array<double, 3> d_global, cd_global;

    double pos_f_global, cpos_f_global;
    double pos_d_global, cpos_d_global;
    double neg_f_global, cneg_f_global;
    double neg_d_global, cneg_d_global;
    double pos_f_crack;
    double pos_d_crack;
    double neg_f_crack;
    double neg_d_crack;
    double pos_f_yield;
    double pos_d_yield;
    double neg_f_yield;
    double neg_d_yield;

    double f_local, cf_local;
    double d_local, cd_local;

    double f_pinch, cf_pinch;
    double d_pinch, cd_pinch;

    double d_zero, cd_zero;

    double f_crack;
    double f_yield;
    const double d_crack;
    const double d_yield;
    const double k_crack;
    const double k_yield;
    const double k_plastic;
    const double unload_from_global_factor;
    const double unload_from_local_factor;
};

#endif