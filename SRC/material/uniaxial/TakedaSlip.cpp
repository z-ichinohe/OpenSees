/*********************************************************************
** OpenSees - Open System for Earthquake Engineering Simulation **
** Pacific Earthquake Engineering Research Center **
****
****
** (C) Copyright 1999, The Regents of the University of California **
** All Rights Reserved. **
****
** Commercial use of this program without express permission of the **
** University of California, Berkeley, is strictly prohibited.  See **
** file 'COPYRIGHT'  in main directory for information on usage and **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES. **
****
** Developed by: **
** Frank McKenna (fmckenna@ce.berkeley.edu) **
** Gregory L. Fenves (fenves@ce.berkeley.edu) **
** Filip C. Filippou (filippou@ce.berkeley.edu) **
****
*********************************************************************/

#include <math.h>
#include <TakedaSlip.h>
#include <elementAPI.h>
#include <Vector.h>
#include <Channel.h>
#include <OPS_Globals.h>
#include <algorithm>

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
using namespace std;

static int numTakedaSlipMaterials = 0;

void *
OPS_TakedaSlip()
{
    if (numTakedaSlipMaterials == 0) {
        numTakedaSlipMaterials++;
        opserr << "Takeda with Slipping Response - Code by Kazuki Ichinohe (Sep23)\n";
    }

 // Pointer to a uniaxial material that will be returned
    UniaxialMaterial * theMaterial = 0;

    int    iData[1];
    double dData[7];
    int numInt = 1;

    if (OPS_GetIntInput(&numInt, iData) != 0)  {
        opserr << "WARNING invalid uniaxialMaterial TakedaSlip tag" << endln;
        return 0;
    }

    int numDouble = 7;

    if (OPS_GetDoubleInput(&numDouble, dData) != 0)  {
        opserr << "Invalid Args want: uniaxialMaterial TakedaSlip tag?";
        opserr << "Crack_Disp? Yield_Disp?";
        opserr << "Crack_Stiffness? Yield_Stiffness? Plastic_Stiffness?";
        opserr << "unload_from_global_factor? unload_from_local_factor?";
        return 0;
    }


 // Parsing was successful, allocate the material
    theMaterial = new TakedaSlip(iData[0],
        dData[0], dData[1],
        dData[2], dData[3], dData[4],
        dData[5], dData[6]);

    if (theMaterial == 0)  {
        opserr << "WARNING could not create uniaxialMaterial of type TakedaSlip Material\n";
        return 0;
    }

    return theMaterial;
}

TakedaSlip::TakedaSlip(int tag,
    double p_Uc, double p_Uy,
    double p_Kc, double p_Ky, double p_Kp,
    double p_b0, double p_b1)
    :UniaxialMaterial(tag, 0),
    d_crack(p_Uc), d_yield(p_Uy),
    k_crack(p_Kc), k_yield(p_Ky), k_plastic(p_Kp),
    unload_from_global_factor(p_b0), unload_from_local_factor(p_b1)
{
    this->revertToStart();
}

TakedaSlip::TakedaSlip()
    :UniaxialMaterial(0, 0),
    d_crack(0), d_yield(0),
    k_crack(0), k_yield(0), k_plastic(0),
    unload_from_global_factor(0), unload_from_local_factor(0)
{
    this->revertToStart();
}

TakedaSlip::~TakedaSlip()
{
 // does nothing
}

int TakedaSlip::setTrialStrain(double strain, double strainRate)
{
 // all variables to the last commit
    this->revertToLastCommit();

 // state determination algorithm: defines the current force and tangent stiffness
    const double k_pinch_factor = 3;
    const double k_global_factor = 1;
    const double d_old = d_new;
    const double f_old = f_new;
    d_new = strain;
    const int ex_branch = branch;

// Unloading
// Positive
    if ((branch == 3 || branch == 4 || branch == 5 || branch == 6) && d_new < d_old) {
        const bool on_backbone = (branch == 5 || branch == 6);
        branch = 1;
        if (on_backbone) {
            pos_f_global = f_old;
            pos_d_global = d_old;
        }
        d_local = d_old;
        f_local = f_old;
        const double pos_f_unload = pos_d_global < pos_d_yield ? pos_f_global : pos_f_yield;
        const double pos_d_unload = pos_d_global < pos_d_yield ? pos_d_global : pos_d_yield;
        k_unload_global = (neg_f_crack - pos_f_unload) / (neg_d_crack - pos_d_unload) * pow(pos_d_unload / pos_d_global, unload_from_global_factor);
        k_tangent = k_unload_global * (on_backbone ? 1 : unload_from_local_factor);
        d_zero = d_local - f_local / k_tangent;
    }
// Negative
    if ((branch == 13 || branch == 14 || branch == 15 || branch == 16) && d_old < d_new) {
        const bool on_backbone = (branch == 15 || branch == 16);
        branch = 11;
        if (on_backbone) {
            neg_f_global = f_old;
            neg_d_global = d_old;
        }
        d_local = d_old;
        f_local = f_old;
        const double neg_f_unload = neg_d_yield < neg_d_global ? neg_f_global : neg_f_yield;
        const double neg_d_unload = neg_d_yield < neg_d_global ? neg_d_global : neg_d_yield;
        k_unload_global = (pos_f_crack - neg_f_unload) / (pos_d_crack - neg_d_unload) * pow(neg_d_unload / neg_d_global, unload_from_global_factor);
        k_tangent = k_unload_global * (on_backbone ? 1 : unload_from_local_factor);
        d_zero = d_local - f_local / k_tangent;
    }

// Forward to Reloading from Unloading
// Towards Positive
    if (branch == 11 && d_zero < d_new) {
        if (pos_d_global <= pos_d_crack) {
            branch = 2;
        } else if (pos_d_global < pos_d_yield) {
            branch = 4;
            const double k_to_global = pos_f_global / (pos_d_global - d_zero);
            const double k_to_yield = pos_f_yield / (pos_d_yield - d_zero);
            if (k_to_global <= k_to_yield) {
                pos_d_global = pos_d_yield;
                pos_f_global = pos_f_yield;
            }
            d_pinch = d_zero;
            f_pinch = 0;
            k_tangent = pos_f_global / (pos_d_global - d_zero);
        } else {
            const double k_from_global = pos_f_global / pos_d_global * k_global_factor;
            const double k_to_global = pos_f_global / (pos_d_global - d_zero);
            if (k_to_global > k_from_global) {
                branch = 4;
                d_pinch = d_zero;
                f_pinch = 0;
                k_tangent = k_to_global;
            } else {
                branch = 3;
                const double global_d_zero = neg_d_global - neg_f_global / k_unload_global;
                const double k_pinch = pos_f_global / (pos_d_global - global_d_zero) * pow(pos_d_global / (pos_d_global - global_d_zero), k_pinch_factor);
                d_pinch = (k_from_global * pos_d_global - k_pinch * d_zero - pos_f_global) / (k_from_global - k_pinch);
                f_pinch = (d_pinch - d_zero) * k_pinch;
                k_tangent = k_pinch;
            }
        }
        // d_zero remained
    }
    if (branch == 12 && d_zero < d_new) {
        branch = 1;
        // k_tangent remained
        // d_zero remained
    }
// Towards Negative
    if (branch == 1 && d_new < d_zero) {
        if (neg_d_crack <= neg_d_global) {
            branch = 12;
        } else if (neg_d_yield < neg_d_global) {
            branch = 14;
            const double k_to_global = neg_f_global / (neg_d_global - d_zero);
            const double k_to_yield = neg_f_yield / (neg_d_yield - d_zero);
            if (k_to_global <= k_to_yield) {
                neg_d_global = neg_d_yield;
                neg_f_global = neg_f_yield;
            }
            d_pinch = d_zero;
            f_pinch = 0;
            k_tangent = neg_f_global / (neg_d_global - d_zero);
        } else {
            const double k_from_global = neg_f_global / neg_d_global * k_global_factor;
            const double k_to_global = neg_f_global / (neg_d_global - d_zero);
            if (k_to_global > k_from_global) {
                branch = 14;
                d_pinch = d_zero;
                f_pinch = 0;
                k_tangent = k_to_global;
            } else {
                branch = 13;
                const double global_d_zero = pos_d_global - pos_f_global / k_unload_global;
                const double k_pinch = neg_f_global / (neg_d_global - global_d_zero) * pow(neg_d_global / (neg_d_global - global_d_zero), k_pinch_factor);
                d_pinch = (k_from_global * neg_d_global - k_pinch * d_zero - neg_f_global) / (k_from_global - k_pinch);
                f_pinch = (d_pinch - d_zero) * k_pinch;
                k_tangent = k_pinch;
            }
        }
        // d_zero remained
    }
    if (branch == 2 && d_new < d_zero) {
        branch = 11;
    }

// Reloading
// Positive
    if (branch == 1 && d_old < d_new && d_local < d_new) {
        if (pos_d_yield < pos_d_global) {
            branch = 3;
            k_tangent = (f_pinch - f_local) / (d_pinch - d_local);
        } else {
            branch = 4;
            k_tangent = (pos_f_global - f_local) / (pos_d_global - d_local);
        }
        d_zero = d_local - f_local / k_tangent;
    }
    if (branch == 2 && d_zero + pos_f_crack / k_tangent < d_new) {
        if (neg_d_yield < neg_d_global && pos_d_global < pos_d_crack) {
            branch = 5;
            k_tangent = k_yield;
            d_zero = pos_d_yield - pos_f_yield / k_tangent;
        } else {
            branch = 4;
            pos_f_global = pos_f_crack;
            pos_d_global = pos_d_crack;
            k_tangent = pos_f_global / (pos_d_global - d_zero);
            // d_zero remained
        }
    }
    if (branch == 3 && d_pinch < d_new) {
        branch = 4;
        k_tangent = (pos_f_global - f_pinch) / (pos_d_global - d_pinch);
        d_zero = pos_d_global - pos_f_global / k_tangent;
    }
    if ((branch == 4 && pos_d_global < d_new) || (branch == 0 && pos_d_crack < d_new)) {
        branch = 5;
        k_tangent = k_yield;
        d_zero = pos_d_yield - pos_f_yield / k_tangent;
    }
    if (branch == 5 && pos_d_yield < d_new) {
        branch = 6;
        k_tangent = k_plastic;
        d_zero = pos_d_yield - pos_f_yield / k_tangent;
    }
// Negative
    if (branch == 11 && d_new < d_old && d_new < d_local) {
        if (neg_d_global < neg_d_yield) {
            branch = 13;
            k_tangent = (f_pinch - f_local) / (d_pinch - d_local);
        } else {
            branch = 14;
            k_tangent = (neg_f_global - f_local) / (neg_d_global - d_local);
        }
        std::cout << d_zero << "\n";
        d_zero = d_local - f_local / k_tangent;
        std::cout << d_zero << "\n";
    }
    if (branch == 12 && d_new < d_zero + neg_f_crack / k_tangent) {
        if (neg_d_crack < neg_d_global && pos_d_global < pos_d_yield) {
            branch = 15;
            k_tangent = k_yield;
            d_zero = neg_d_yield - neg_f_yield / k_tangent;
        } else {
            branch = 14;
            neg_f_global = neg_f_crack;
            neg_d_global = neg_d_crack;
            k_tangent = neg_f_global / (neg_d_global - d_zero);
            // d_zero remained
        }
    }
    if (branch == 13 && d_new < d_pinch) {
        branch = 14;
        k_tangent = (neg_f_global - f_pinch) / (neg_d_global - d_pinch);
        d_zero = neg_d_global - neg_f_global / k_tangent;
    }
    if ((branch == 14 && d_new < neg_d_global) || (branch == 0 && d_new < neg_d_crack)) {
        branch = 15;
        k_tangent = k_yield;
        d_zero = neg_d_yield - neg_f_yield / k_tangent;
    }
    if (branch == 15 && d_new < neg_d_yield) {
        branch = 16;
        k_tangent = k_plastic;
        d_zero = neg_d_yield - neg_f_yield / k_tangent;
    }

// Calculate Force
    f_new = (d_new - d_zero) * k_tangent;

    // if (branch != ex_branch) {
    //     std::cout << ex_branch << " -> " << branch << "\n";
    // }

    return 0;
}

double TakedaSlip::getStress(void)
{
 // cout << " getStress" << endln;
    return (f_new);
}

double TakedaSlip::getTangent(void)
{
 // cout << " getTangent" << endln;
    return (k_tangent);
}

double TakedaSlip::getInitialTangent(void)
{
 // cout << " getInitialTangent" << endln;
    return (k_crack);
}

double TakedaSlip::getStrain(void)
{
 // cout << " getStrain" << endln;
    return (d_new);
}

int TakedaSlip::commitState(void)
{
    cbranch = branch;
    ck_tangent = k_tangent;
    k_unload_global = k_unload_global;
    cf_new = f_new;
    cd_new = d_new;
    cf_local = f_local;
    cd_local = d_local;
    cf_pinch = f_pinch;
    cd_pinch = d_pinch;
    cd_zero = d_zero;
    cpos_f_global = pos_f_global;
    cpos_d_global = pos_d_global;
    cneg_f_global = neg_f_global;
    cneg_d_global = neg_d_global;
    return 0;
}

int TakedaSlip::revertToLastCommit(void)
{
    branch = cbranch;
    k_tangent = ck_tangent;
    k_unload_global = k_unload_global;
    f_new = cf_new;
    d_new = cd_new;
    f_local = cf_local;
    d_local = cd_local;
    f_pinch = cf_pinch;
    d_pinch = cd_pinch;
    d_zero = cd_zero;
    pos_f_global = cpos_f_global;
    pos_d_global = cpos_d_global;
    neg_f_global = cneg_f_global;
    neg_d_global = cneg_d_global;
    return 0;
}

int TakedaSlip::revertToStart(void)
{
    branch = cbranch = 0;
    k_tangent = ck_tangent = k_crack;
    const double f_crack = d_crack * k_crack;
    const double f_yield = f_crack + k_yield * (d_yield - d_crack);

    pos_f_global = cpos_f_global = pos_f_crack = f_crack;
    pos_d_global = cpos_d_global = pos_d_crack = d_crack;
    neg_f_global = cneg_f_global = neg_f_crack = - f_crack;
    neg_d_global = cneg_d_global = neg_d_crack = - d_crack;
    pos_f_yield = f_yield;
    pos_d_yield = d_yield;
    neg_f_yield = - f_yield;
    neg_d_yield = - d_yield;

    d_zero = cd_zero = 0;
    f_pinch = cf_pinch = 0;

    return 0;
}

UniaxialMaterial *
TakedaSlip::getCopy(void)
{
    TakedaSlip * theCopy = new TakedaSlip(
        this->getTag(),
        d_crack, d_yield,
        k_crack, k_yield, k_plastic,
        unload_from_global_factor, unload_from_local_factor);

    theCopy->branch = branch;
    theCopy->k_tangent = k_tangent;
    theCopy->k_unload_global = k_unload_global;
    theCopy->d_new = d_new;
    theCopy->f_new = f_new;
    theCopy->f_local = f_local;
    theCopy->d_local = d_local;
    theCopy->f_pinch = f_pinch;
    theCopy->d_pinch = d_pinch;
    theCopy->d_zero = d_zero;
    theCopy->pos_f_global = pos_f_global;
    theCopy->pos_d_global = pos_d_global;
    theCopy->neg_f_global = neg_f_global;
    theCopy->neg_d_global = neg_d_global;
    theCopy->pos_f_crack = pos_f_crack;
    theCopy->pos_d_crack = pos_d_crack;
    theCopy->pos_f_yield = pos_f_yield;
    theCopy->pos_d_yield = pos_d_yield;
    theCopy->neg_f_crack = neg_f_crack;
    theCopy->neg_d_crack = neg_d_crack;
    theCopy->neg_f_yield = neg_f_yield;
    theCopy->neg_d_yield = neg_d_yield;

    theCopy->cbranch = cbranch;
    theCopy->ck_tangent = ck_tangent;
    theCopy->k_unload_global = k_unload_global;
    theCopy->cd_new = cd_new;
    theCopy->cf_new = cf_new;
    theCopy->cf_local = cf_local;
    theCopy->cd_local = cd_local;
    theCopy->cf_pinch = cf_pinch;
    theCopy->cd_pinch = cd_pinch;
    theCopy->cd_zero = cd_zero;
    theCopy->cpos_f_global = cpos_f_global;
    theCopy->cpos_d_global = cpos_d_global;
    theCopy->cneg_f_global = cneg_f_global;
    theCopy->cneg_d_global = cneg_d_global;
    return theCopy;
}

int TakedaSlip::sendSelf(int cTag, Channel &theChannel)
{
    int res = 0;
    cout << " sendSelf" << endln;

    static Vector data(137);
    data(0) = this->getTag();
    data(11) = branch;
    data(13) = d_new;
    data(15) = f_new;
    data(16) = k_tangent;
    data(18) = f_local;
    data(20) = d_local;
    data(25) = d_pinch;


    data(28) = k_unload_global;

    data(111) = cbranch;
    data(113) = cd_new;
    data(115) = cf_new;
    data(116) = ck_tangent;
    data(118) = cf_local;
    data(120) = cd_local;
    data(125) = cd_pinch;


    data(128) = k_unload_global;
    res = theChannel.sendVector(this->getDbTag(), cTag, data);
    if (res < 0)
        opserr << "TakedaSlip::sendSelf() - failed to send data\n";
    return res;
}

int TakedaSlip::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    int res = 0;
    static Vector data(137);
    res = theChannel.recvVector(this->getDbTag(), cTag, data);

    if (res < 0) {
        opserr << "TakedaSlip::recvSelf() - failed to receive data\n";
        this->setTag(0);
    }
    else {
        cout << " recvSelf" << endln;
        this->setTag((int)data(0));
        branch = data(11);
        d_new = data(13);
        f_new = data(15);
        k_tangent = data(16);
        f_local = data(18);
        d_local = data(20);
        d_pinch = data(25);


        k_unload_global = data(28);

        cbranch = data(111);
        cd_new = data(113);
        cf_new = data(115);
        ck_tangent = data(116);
        cf_local = data(118);
        cd_local = data(120);
        cd_pinch = data(125);


        k_unload_global = data(128);
    }

    return res;
}

void TakedaSlip::Print(OPS_Stream &s, int flag)
{
    cout << "TakedaSlip tag: " << this->getTag() << endln;
}