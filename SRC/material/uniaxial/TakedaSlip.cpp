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
    bool on_backbone  = (branch == 4 || branch == 5);
    if ((branch == 2 || branch == 3 || on_backbone) && (d_new - d_old) * sign <= 0) {
        branch = 1;
        if (on_backbone) {
            f_global[is] = f_old;
            d_global[is] = d_old;
        }
        d_local = d_old;
        f_local = f_old;
        const double f_reload = d_yield > abs(d_global[is]) ? abs(f_global[is]) : f_yield;
        const double d_reload = d_yield > abs(d_global[is]) ? abs(d_global[is]) : d_yield;
        k_unload = (f_reload + f_crack) / (d_reload + d_crack) * pow(d_reload / abs(d_global[is]), unload_from_global_factor);
        if (!on_backbone) {
            k_unload *= unload_from_local_factor;
        }
        k_tangent = k_unload;
        d_zero = d_local - f_local / k_unload;
    }

// Forward to Reloading from Unloading
    if (branch == 1 && (d_new - d_zero) * sign <= 0) {
        is = d_new > d_old ? 1 : 2;
        sign = d_new > d_old ? 1 : -1;
        if (abs(d_global[is]) <= d_crack) {
            // 降伏している場合、p_crackに向かわない
            branch = 15;
            d_unload = d_zero + sign * f_crack / k_unload;
            f_unload = f_crack * sign;
        } else if (abs(d_global[is]) <= d_yield) {
            branch = 3;
            const double k_to_global = f_global[is] / (d_global[is] - d_zero);
            const double k_to_yield = f_yield * sign / (d_yield * sign - d_zero);
            // 現在降伏ブランチに負勾配を設定しているので問題はないが、一般的には不等号が逆
            if (k_to_global <= k_to_yield) {
                d_global[is] = d_yield * sign;
                f_global[is] = f_yield * sign;
            }
            k_tangent = f_global[is] / (d_global[is] - d_zero);
        } else {
            // ここでは正側を見てる
            const double k_from_global = f_global[is] / d_global[is] * k_global_factor;
            const double k_to_global = f_global[is] / (d_global[is] - d_zero);
            if (k_to_global > k_from_global) {
                branch = 3;
                k_tangent = k_to_global;
            // }
            // const double k_global = f_global[is] / d_global[is] * k_global_factor;
            // const double d_zero_from_global = d_global[is] - f_global[is] / k_global;
            // d_pinch = d_zero;
            // if ((d_zero_from_global - d_zero) * sign <= 0) {
            //     branch = 3;
            //     k_tangent = f_global[is] / (d_global[is] - d_zero);
            } else {
                branch = 2;
                const double k_global = f_global[is] / d_global[is] * k_global_factor;
                // なんで反対側の値が必要なのか一向に謎
                const double k_temp = (f_crack + f_yield) / (d_crack + d_yield) * pow((d_yield / abs(d_global[3 - is])), unload_from_global_factor);
                const double d_zero_from_global = d_global[3 - is] - f_global[3 - is] / k_temp;
                const double k_pinch = abs(f_global[is] / (d_zero_from_global - d_global[is])) * pow( abs(d_global[is] / (d_zero_from_global - d_global[is])), k_pinch_factor);
                if (abs(k_pinch - k_global) <= error) {
                    d_pinch = d_zero;
                } else {
                    d_pinch = (k_global * d_global[is] - k_pinch * d_zero - f_global[is]) / (k_global - k_pinch);
                }
                f_pinch = (d_pinch - d_zero) * k_pinch;
                k_tangent = k_pinch;
            }
        }
    }

    if (branch == 15 && (d_unload - d_new) * sign <= 0) {
        if (abs(d_global[is]) <= d_crack && abs(d_global[3 - is]) <= d_yield) {
            branch = 4;
            k_tangent = k_yield;
        } else {
            // p_unloadからp_crackに向かうわけではないので、段差が生じる
            branch = 3;
            f_global[is] = f_crack * sign;
            d_global[is] = d_crack * sign;
            k_tangent = f_global[is] / (d_global[is] - d_zero);
            // %d_zero = d_yield * sign - f_yield * sign * (d_yield * sign - d_unload) / (f_yield * sign - f_unload);%%%%�o������ψʂƙ��f�͂𐳕��̂ǂ�����Ђъ���_�𒴂���܂łЂъ���_���L������悤�ɂ���
        }
    }

// Back to the Reloding from Unloading
    if (branch == 15 && (d_new - d_zero) * sign <= 0) {
        branch = 1;
        is = d_new > d_old ? 1 : 2;
        sign = d_new > d_old ? 1 : -1;
    }
    if (branch == 1 && (d_local - d_new) * sign <= 0) {
        if (d_yield <= abs(d_global[is])) {
            branch = 2;
            k_tangent = f_pinch / (d_pinch - d_zero);
        } else {
            branch = 3;
            k_tangent = f_global[is] / (d_global[is] - d_zero);
        }
    }
// Pinching and Reloading
// 2 -> 3
    if (branch == 2 && (d_pinch - d_new) * sign <= 0) {
        branch = 3;
        const double k_global = (f_global[is] - f_pinch) / (d_global[is] - d_pinch);
        k_tangent = k_global;
        d_zero = d_global[is] - f_global[is] / k_global;
    }
// 3 -> 4
    if (branch == 3 && (d_global[is] - d_new) * sign <= 0) {
        branch = 4;
        k_tangent = k_yield;
    }
// Backbone
// 0 -> 4
    if (branch == 0 && d_crack <= abs(d_new))  {
        branch = 4;
        k_tangent = k_yield;
    }
// 4 -> 5
    if (branch == 4 && d_yield <= abs(d_new)) {
        branch = 5;
        k_tangent = k_plastic;
    }

// Calculate Force
    if (branch == 0 || branch == 1 || branch == 2 || branch == 3 || branch == 15) {
        f_new = (d_new - d_zero) * k_tangent;
    } else if (branch == 4 || branch == 5) {
        f_new = sign * f_yield + (d_new - sign * d_yield) * k_tangent;
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
    ck_tangent = k_tangent;
    ck_unload = k_unload;
    cf_new = f_new;
    cd_new = d_new;
    cf_global = f_global;
    cd_global = d_global;
    cf_local = f_local;
    cd_local = d_local;
    cf_unload = f_unload;
    cd_unload = d_unload;
    cf_pinch = f_pinch;
    cd_pinch = d_pinch;
    cd_zero = d_zero;
    return 0;
}

int TakedaSlip::revertToLastCommit(void)
{
    branch = cbranch;
    k_tangent = ck_tangent;
    k_unload = ck_unload;
    f_new = cf_new;
    d_new = cd_new;
    f_global = cf_global;
    d_global = cd_global;
    f_local = cf_local;
    d_local = cd_local;
    f_unload = cf_unload;
    d_unload = cd_unload;
    f_pinch = cf_pinch;
    d_pinch = cd_pinch;
    d_zero = cd_zero;
    return 0;
}

int TakedaSlip::revertToStart(void)
{
    branch = cbranch = 0;
    k_tangent = ck_tangent = k_crack;
    f_crack = d_crack * k_crack;
    f_yield = f_crack + k_yield * (d_yield - d_crack);

    f_global[1] = cf_global[1] = f_crack;
    f_global[2] = cf_global[2] = - f_crack;
    d_global[1] = cd_global[1] = d_crack;
    d_global[2] = cd_global[2] = - d_crack;

    d_zero = cd_zero = 0;

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
    theCopy->k_unload = k_unload;
    theCopy->d_new = d_new;
    theCopy->f_new = f_new;
    theCopy->f_global = f_global;
    theCopy->d_global = d_global;
    theCopy->f_local = f_local;
    theCopy->d_local = d_local;
    theCopy->f_unload = f_unload;
    theCopy->d_unload = d_unload;
    theCopy->f_pinch = f_pinch;
    theCopy->d_pinch = d_pinch;
    theCopy->d_zero = d_zero;


    theCopy->cbranch = cbranch;
    theCopy->ck_tangent = ck_tangent;
    theCopy->ck_unload = ck_unload;
    theCopy->cd_new = cd_new;
    theCopy->cf_new = cf_new;
    theCopy->cf_global = cf_global;
    theCopy->cd_global = cd_global;
    theCopy->cf_local = cf_local;
    theCopy->cd_local = cd_local;
    theCopy->cf_unload = cf_unload;
    theCopy->cd_unload = cd_unload;
    theCopy->cf_pinch = cf_pinch;
    theCopy->cd_pinch = cd_pinch;
    theCopy->cd_zero = cd_zero;
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


    data(28) = k_unload;

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


    data(128) = ck_unload;
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


        k_unload = data(28);

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


        ck_unload = data(128);
    }

    return res;
}

void TakedaSlip::Print(OPS_Stream &s, int flag)
{
    cout << "TakedaSlip tag: " << this->getTag() << endln;
}