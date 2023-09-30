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

#include < math.h >
#include < TakedaSlip.h >
#include < elementAPI.h >
#include < Vector.h >
#include < Channel.h >
#include < OPS_Globals.h >
#include < algorithm >

#include < cstdlib >
#include < iomanip >
#include < iostream >
#include < stdio.h >
#include < stdlib.h >
#include < string >
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
    const double error = 0.0001;
    const double k_pinch_factor = 3;
    const double k_global_factor = 1;
    const double d_old = d_new;
    const double f_old = f_new;
    d_new = strain;
    int is = f_old > 0 ? 1 : 2;
    int sign = f_old > 0 ? 1 : -1;

// Reloading to Unloading
    if ((branch == 2 || branch == 3 || branch == 6 || branch == 7 || branch == 9) && (d_new - d_old) * sign <= 0) {
        if (branch == 2 || branch == 3) {
            f_global[is] = f_old;
            d_global[is] = d_old;
        }
        d_local = d_old;
        f_local = f_old;
        if (branch == 2) {
            k_local = (abs(f_global[is]) + f_crack) / (abs(d_global[is]) + d_crack);
        } else {
            k_local = (f_crack + f_yield) / (d_crack + d_yield) * pow((d_yield / abs(d_global[is])), unload_from_global_factor);
        }
        if (branch == 2 || branch == 3) {
            branch = 4;
        } else {
            branch = 8;
            k_local *= unload_from_local_factor;
        }
        k_tangent = k_local;
        d_zero = d_local - f_local / k_local;
    }

// 4 and 5 are same branch
    if (branch == 5 && (d_new - d_zero) * sign <= 0) {
        branch = 4;
        is = d_new > d_old ? 1 : 2;
        sign = d_new > d_old ? 1 : -1;
    }

// Reloading From Unloading
    if ((branch == 4 || branch == 8) && (d_new - d_zero) * sign <= 0) {
        is = d_new > d_old ? 1 : 2;
        sign = d_new > d_old ? 1 : -1;
        // -Crack: 5
        // -Yield: 9
        // -: 6
        if (abs(d_global[is]) <= d_crack) {
            branch = 5;
            d_unload = d_zero + sign * f_crack / k_local;
            f_unload = f_crack * sign;
        } else if (abs(d_global[is]) <= d_yield) {
            branch = 7;
            const double k_to_global = f_global[is] / (d_global[is] - d_zero);
            const double k_to_yield = f_yield * sign / (d_yield * sign - d_zero);
            if (k_to_global <= k_to_yield) {
                d_global[is] = d_yield * sign;
                f_global[is] = f_yield * sign;
            }
            k_tangent = f_global[is] / (d_global[is] - d_zero);
        } else {
            const double k_global = f_global[is] / d_global[is] * k_global_factor;
            const double d_zero_from_global = d_global[is] - f_global[is] / k_global;
            d_pinch = d_zero;
            if ((d_zero_from_global - d_zero) * sign <= 0) {
                branch = 7;
                k_tangent = f_global[is] / (d_global[is] - d_zero);
            } else {
                branch = 6;
                const double k_global = f_global[is] / d_global[is] * k_global_factor;
                const double k_temp = (f_crack + f_yield) / (d_crack + d_yield) * pow((d_yield / abs(d_global[3 - is])), unload_from_global_factor);
                const double d_zero_from_global = d_global[3 - is] - f_global[3 - is] / k_temp;
                const double k_pinch = abs(f_global[is] / (d_zero_from_global - d_global[is])) * pow( abs(d_global[is] / (d_zero_from_global - d_global[is])), k_pinch_factor);
                if (abs(k_pinch - k_global) <= error) {
                    d_pinch = d_zero;
                } else {
                    d_pinch = (k_global * d_global[is] - k_pinch * d_zero - f_global[is]) / (k_global - k_pinch);
                }
                k_tangent = k_pinch;
            }
        }
    }

// Unloading Branch to Reloading
// 4 -> 2
    if (branch == 4 && (d_global[is] - d_new) * sign <= 0) {
        branch = 2;
        k_tangent = k_yield;
    }
// 5 -> 2, 9
    if (branch == 5 && (d_unload - d_new) * sign <= 0) {
        if (abs(d_global[is]) <= d_crack && abs(d_global[3 - is]) <= d_yield) {
            branch = 2;
            k_tangent = k_yield;
        } else {
            branch = 7;
            f_global[is] = f_crack * sign;
            d_global[is] = d_crack * sign;
            k_tangent = f_global[is] / (d_global[is] - d_zero);
            // %d_zero = d_yield * sign - f_yield * sign * (d_yield * sign - d_unload) / (f_yield * sign - f_unload);%%%%�o������ψʂƙ��f�͂𐳕��̂ǂ�����Ђъ���_�𒴂���܂łЂъ���_���L������悤�ɂ���
        }
    }
// 8 -> 6, 9
    if (branch == 8 && (d_local - d_new) * sign <= 0) {
        if (d_yield <= abs(d_global[is])) {
            branch = 6;
            const double k_temp = (f_crack + f_yield) / (d_crack + d_yield) * pow((d_yield / abs(d_global[3 - is])), unload_from_global_factor);
            const double d_zero_from_global = d_global[3 - is] - f_global[3 - is] / k_temp;
            const double k_pinch = abs(f_global[is] / (d_zero_from_global - d_global[is])) * pow( abs(d_global[is] / (d_zero_from_global - d_global[is])), k_pinch_factor);
            k_tangent = k_pinch;
        } else {
            branch = 7;
            k_tangent = f_global[is] / (d_global[is] - d_zero);
        }
    }
// Pinching and Reloading
// 6 -> 7
    if (branch == 6 && (d_pinch - d_new) * sign <= 0) {
        branch = 7;
        const double k_global = f_global[is] / d_global[is] * k_global_factor;
        const double d_zero_from_global = d_global[is] - f_global[is] / k_global;
        k_tangent = k_global;
        d_zero = d_zero_from_global;
    }
// 7 -> 2
    if (branch == 7 && (d_global[is] - d_new) * sign <= 0) {
        branch = 2;
        k_tangent = k_yield;
    }
// Backbone
// 1 -> 2
    if (branch == 1 && d_crack <= abs(d_new))  {
        branch = 2;
        k_tangent = k_yield;
    }
// 2 -> 3
    if (branch == 2 && d_yield <= abs(d_new)) {
        branch = 3;
        k_tangent = k_plastic;
    }

// Calculate Force
    if (branch == 1) {
        f_new = k_crack * d_new;
    } else if (branch == 2) {
        f_new = sign * f_crack + (d_new - sign * d_crack) * k_yield;
    } else if (branch == 3) {
        f_new = f_yield * sign + (d_new - d_yield * sign) * k_plastic;
    } else if (branch == 4) {
        f_new = f_global[is] + (d_new - d_global[is]) * k_local;
    } else if (branch == 5 || branch == 6 || branch == 7) {
        f_new = (d_new - d_zero) * k_tangent;
    } else if (branch == 8) {
        f_new = f_local + (d_new - d_local) * k_local;
    }
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
    cd_new = d_new;
    cf_new = f_new;
    ck_tangent = k_tangent;
    cf_unload = f_unload;
    cf_local = f_local;
    cd_unload = d_unload;
    cd_local = d_local;
    cf_global = f_global;
    cd_global = d_global;
    cd_pinch = d_pinch;
    ck_local = k_local;
    return 0;
}

int TakedaSlip::revertToLastCommit(void)
{
    branch = cbranch;
    d_new = cd_new;
    f_new = cf_new;
    k_tangent = ck_tangent;
    f_unload = cf_unload;
    f_local = cf_local;
    d_unload = cd_unload;
    d_local = cd_local;
    f_global = cf_global;
    d_global = cd_global;
    d_pinch = cd_pinch;
    k_local = ck_local;
    return 0;
}

int TakedaSlip::revertToStart(void)
{
    branch = cbranch = 1;
    k_tangent = ck_tangent = k_crack;
    f_crack = d_crack * k_crack;
    f_yield = f_crack + k_yield * (d_yield - d_crack);

    f_global[1] = cf_global[1] = f_crack;
    f_global[2] = cf_global[2] = - f_crack;
    d_global[1] = cd_global[1] = d_crack;
    d_global[2] = cd_global[2] = - d_crack;


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
    theCopy->d_new = d_new;
    theCopy->f_new = f_new;
    theCopy->k_tangent = k_tangent;
    theCopy->f_unload = f_unload;
    theCopy->f_local = f_local;
    theCopy->d_unload = d_unload;
    theCopy->d_local = d_local;
    theCopy->f_global = f_global;
    theCopy->d_global = d_global;
    theCopy->d_pinch = d_pinch;

    theCopy->k_local = k_local;

    theCopy->cbranch = cbranch;
    theCopy->cd_new = cd_new;
    theCopy->cf_new = cf_new;
    theCopy->ck_tangent = ck_tangent;
    theCopy->cf_unload = cf_unload;
    theCopy->cf_local = cf_local;
    theCopy->cd_unload = cd_unload;
    theCopy->cd_local = cd_local;
    theCopy->cf_global = cf_global;
    theCopy->cd_global = cd_global;
    theCopy->cd_pinch = cd_pinch;

    theCopy->ck_local = ck_local;
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
    data(17) = f_unload;
    data(18) = f_local;
    data(19) = d_unload;
    data(20) = d_local;
    data(21) = f_global[1];
    data(22) = f_global[2];
    data(23) = d_global[1];
    data(24) = d_global[2];
    data(25) = d_pinch;


    data(28) = k_local;

    data(111) = cbranch;
    data(113) = cd_new;
    data(115) = cf_new;
    data(116) = ck_tangent;
    data(117) = cf_unload;
    data(118) = cf_local;
    data(119) = cd_unload;
    data(120) = cd_local;
    data(121) = cf_global[1];
    data(122) = cf_global[2];
    data(123) = cd_global[1];
    data(124) = cd_global[2];
    data(125) = cd_pinch;


    data(128) = ck_local;
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
        f_unload = data(17);
        f_local = data(18);
        d_unload = data(19);
        d_local = data(20);
        f_global[1] = data(21);
        f_global[2] = data(22);
        d_global[1] = data(23);
        d_global[2] = data(24);
        d_pinch = data(25);


        k_local = data(28);

        cbranch = data(111);
        cd_new = data(113);
        cf_new = data(115);
        ck_tangent = data(116);
        cf_unload = data(117);
        cf_local = data(118);
        cd_unload = data(119);
        cd_local = data(120);
        cf_global[1] = data(121);
        cf_global[2] = data(122);
        cd_global[1] = data(123);
        cd_global[2] = data(124);
        cd_pinch = data(125);


        ck_local = data(128);
    }

    return res;
}

void TakedaSlip::Print(OPS_Stream &s, int flag)
{
    cout << "TakedaSlip tag: " << this->getTag() << endln;
}