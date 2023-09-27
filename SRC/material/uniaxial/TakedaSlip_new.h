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

class TakedaSlip : public UniaxialMaterial
{
public:
    TakedaSlip(int tag,
        double Uc, double Uy,
        double Kc, double Ky, double Kp,
        double b0, double b1);
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
// 7 Fixed input material parameters
    const double Uc;
    const double Uy;
    const double Kc;
    const double Ky;
    const double Kp;
    const double b0;
    const double b1;
// 4 Initial Variables
    const double posUc, negUc;
    const double posUy, negUy;
    const double posFc, negFc;
    const double posFy, negFy;
// History Variables
// 10 U and F
    double posUlocal, cPosUlocal;
    double posFlocal, cPosFlocal;
    double posUglobal, cPosUglobal;
    double posFglobal, cPosFglobal;
    double negUlocal, cNegUlocal;
    double negFlocal, cNegFlocal;
    double negUglobal, cNegUglobal;
    double negFglobal, cNegFglobal;
    double Fpinch, cFpinch;
    double Upinch, cUpinch;
// 3 State Variables
    double U, cU;
    double Ui, cUi;
    double Fi, cFi;
// 2 Stiffness
    double Kreload, cKreload, KgetTangent;
    double Kunload, cKunload;
// 1 Flag
    int    Branch,         cBranch;
};

#endif