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
    std::tuple<int, float, float>Tslip_120(double d_yield, double d_new, int sign, double f_crack, double d_crack, double k_yield, double k_plastic, double f_yield);
    int branch, cbranch;
    double d_old, cd_old;
    double d_new, cd_new;
    double f_old, cf_old;
    double f_new, cf_new;
    double k_tangent, ck_tangent;
    std::array<double, 3> k_unload, ck_unload;
    double k_local, ck_local;
    double f_reload, cf_reload;
    double f_local, cf_local;
    double d_reload, cd_reload;
    double d_local, cd_local;
    double d_zero, cd_zero;
    double d_zero_from_local, cd_zero_from_local;
    std::array<double, 3> f_global, cf_global;
    std::array<double, 3> d_global, cd_global;
    double d_pinch, cd_pinch;
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