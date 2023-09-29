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

std::tuple<int, float, float>TakedaSlip::Tslip_120(double d_yield, double d_new, int sign, double f_crack, double d_crack, double k_yield, double k_plastic, double f_yield)
{
    int branch;
    double k_tangent;
    double f_new;
    if (d_yield - abs(d_new) > 0)  {
        branch = 2;
        k_tangent = k_yield;
        f_new = sign * f_crack + (d_new - sign * d_crack) * k_yield;
        return {branch, k_tangent, f_new};
    } else {
        branch = 3;
        k_tangent = k_plastic;
        f_new = f_yield * sign + (d_new - d_yield * sign) * k_plastic;
        return {branch, k_tangent, f_new};
    }
}

int TakedaSlip::setTrialStrain(double strain, double strainRate)
{
 // all variables to the last commit
    this->revertToLastCommit();

 // state determination algorithm: defines the current force and tangent stiffness
    double k_global_from_pinch, k_global, k_toYield, k_pinch, d_zero_from_pinch;
    const double error = 0.0001;
    const double k_pinch_factor = 3;
    const double k_global_factor = 1;
    const double d_old = d_new;
    const double f_old = f_new;
    d_new = strain;
    int is = f_old > 0 ? 1 : 2;
    int sign = f_old > 0 ? 1 : -1;
    if (branch == 1)  {
        if (d_crack - abs(d_new) > 0)  {
            f_new = k_crack * d_new;

        } else {
            f_global[1] = f_crack;
            f_global[2] = - f_crack;
            d_global[1] = d_crack;
            d_global[2] = - d_crack;
            k_unload[1] = k_crack;
            k_unload[2] = k_crack;
            is = d_new > 0 ? 1 : 2;
            sign = d_new > 0 ? 1 : -1;
            std::tuple<int, float, float>llssff = this->Tslip_120(d_yield, d_new, sign, f_crack, d_crack, k_yield, k_plastic, f_yield);
            branch = std::get<0>(llssff);
            k_tangent = std::get<1>(llssff);
            f_new = std::get<2>(llssff);
        }
    }

    if (branch == 2) {
        if ((d_new - d_old) * sign > 0) {
            std::tuple<int, float, float>llssff = this->Tslip_120(d_yield, d_new, sign, f_crack, d_crack, k_yield, k_plastic, f_yield);
            branch = std::get<0>(llssff);
            k_tangent = std::get<1>(llssff);
            f_new = std::get<2>(llssff);

        } else {
            f_global[is] = f_old;
            d_global[is] = d_old;
            k_unload[is] = (abs(f_global[is]) + f_crack) / (abs(d_global[is]) + d_crack);
            k_global_from_pinch = k_global_factor * f_global[is] / d_global[is];
            d_zero = d_old - f_old / k_unload[is];
            if ((d_new - d_zero) * sign >= 0) {
                branch = 4;
                k_tangent = k_unload[is];
                f_new = f_global[is] + (d_new - d_global[is]) * k_unload[is];

            } else {
 // % loop of 430 %%
                is = is == 1 ? 2 : 1;
                sign = sign == 1 ? -1 : 1;
                if (abs(d_global[is]) <= d_crack) {
                    d_reload = d_zero + sign * f_crack / k_unload[3 - is];
                    f_reload = f_crack * sign;
                    if ((d_reload - d_new) * sign > 0) {
                        branch = 5;
                        f_new = (d_new - d_zero) * k_tangent;

                    } else {
 // % loop of 530 %%
                        if ((d_global[is] * sign <= d_crack) && (abs(d_global[3 - is]) <= d_yield)) {
                            std::tuple<int, float, float>llssff = this->Tslip_120(d_yield, d_new, sign, f_crack, d_crack, k_yield, k_plastic, f_yield);
                            branch = std::get<0>(llssff);
                            k_tangent = std::get<1>(llssff);
                            f_new = std::get<2>(llssff);

                        } else {
                            f_global[is] = f_yield * sign;
                            d_global[is] = d_yield * sign;
                            d_zero = d_yield * sign - f_yield * sign * (d_yield * sign - d_reload) / (f_yield * sign - f_reload);
                            if ((d_new - d_global[is]) * sign >= 0) {
                                std::tuple<int, float, float>llssff = this->Tslip_120(d_yield, d_new, sign, f_crack, d_crack, k_yield, k_plastic, f_yield);
                                branch = std::get<0>(llssff);
                                k_tangent = std::get<1>(llssff);
                                f_new = std::get<2>(llssff);

                            } else {
                                branch = 9;
                                k_tangent = f_global[is] / (d_global[is] - d_zero);
                                f_new = k_tangent * (d_new - d_zero);

                            }
                        }
                    }
                } else if (abs(d_global[is]) <= d_yield) {
                    k_global = f_global[is] / (d_global[is] - d_zero);
                    k_toYield = f_yield * sign / (d_yield * sign - d_zero);
                    if (k_global < k_toYield) {
                        d_global[is] = d_yield * sign;
                        f_global[is] = f_yield * sign;
                    }
                    if ((d_new - d_global[is]) * sign >= 0) {
                        std::tuple<int, float, float>llssff = this->Tslip_120(d_yield, d_new, sign, f_crack, d_crack, k_yield, k_plastic, f_yield);
                    branch = std::get<0>(llssff);
                    k_tangent = std::get<1>(llssff);
                    f_new = std::get<2>(llssff);

                    } else {
                        branch = 9;
                        k_tangent = f_global[is] / (d_global[is] - d_zero);
                        f_new = k_tangent * (d_new - d_zero);

                    }
                } else {
 // % loop of 440 %%
                    k_global = d_global[3 - is] - f_global[3 - is] / k_unload[3 - is];
                    k_pinch = abs(f_global[is] / (k_global - d_global[is])) * pow( abs(d_global[is] / (k_global - d_global[is])), k_pinch_factor);
                    k_global_from_pinch = f_global[is] / d_global[is] * k_global_factor;
                    if (abs(k_pinch - k_global_from_pinch) < error) {
                        d_pinch = d_zero;
                    } else {
                        d_pinch = (k_global_from_pinch * d_global[is] - k_pinch * d_zero - f_global[is]) / (k_global_from_pinch - k_pinch);
                    }

                    if ((d_new - d_pinch) * sign > 0) {
                        if ((d_new - d_global[is]) * sign >= 0) {
                            std::tuple<int, float, float>llssff = this->Tslip_120(d_yield, d_new, sign, f_crack, d_crack, k_yield, k_plastic, f_yield);
                            branch = std::get<0>(llssff);
                            k_tangent = std::get<1>(llssff);
                            f_new = std::get<2>(llssff);

                        } else {
                            branch = 7;
                            k_global_from_pinch = f_global[is] / d_global[is] * k_global_factor;
                            k_tangent = k_global_from_pinch;
                            d_zero_from_pinch = d_global[is] - f_global[is] / k_global_from_pinch;
                            d_zero = d_zero_from_pinch;
                            f_new = k_global_from_pinch * (d_new - d_zero_from_pinch);

                        }
                    } else {
                        branch = 6;
                        k_global = d_global[3 - is] - f_global[3 - is] / k_unload[3 - is];
                        k_pinch = abs(f_global[is] / (k_global - d_global[is])) * pow( abs(d_global[is] / (k_global - d_global[is])), k_pinch_factor);
                        k_tangent = k_pinch;
                        f_new = k_pinch * (d_new - d_zero);

                    }
                }
            }
        }
    }

 // % rule 3 loading on post yielding primary curve %%
    if (branch == 3) {
        if ((d_new - d_old) * sign > 0) {
            std::tuple<int, float, float>llssff = this->Tslip_120(d_yield, d_new, sign, f_crack, d_crack, k_yield, k_plastic, f_yield);
            branch = std::get<0>(llssff);
            k_tangent = std::get<1>(llssff);
            f_new = std::get<2>(llssff);

        } else {
            k_unload[is] = pow((d_yield / abs(d_old)), unload_from_global_factor) * (f_crack + f_yield) / (d_crack + d_yield);
            f_global[is] = f_old;
            d_global[is] = d_old;
            k_global_from_pinch = k_global_factor * f_global[is] / d_global[is];
            d_zero = d_old - f_old / k_unload[is];
            if ((d_new - d_zero) * sign >= 0) {
                branch = 4;
                k_tangent = k_unload[is];
                f_new = f_global[is] + (d_new - d_global[is]) * k_unload[is];

            } else {
 // % loop of 430 %%
                is = is == 1 ? 2 : 1;
                sign = sign == 1 ? -1 : 1;
                if (abs(d_global[is]) <= d_crack) {
                    d_reload = d_zero + sign * f_crack / k_unload[3 - is];
                    f_reload = f_crack * sign;
                    if ((d_reload - d_new) * sign > 0) {
                        branch = 5;
                        f_new = (d_new - d_zero) * k_tangent;

                    } else {
 // % loop of 530 %%
                        if ((d_global[is] * sign <= d_crack) && (abs(d_global[3 - is]) <= d_yield)) {
                            std::tuple<int, float, float>llssff = this->Tslip_120(d_yield, d_new, sign, f_crack, d_crack, k_yield, k_plastic, f_yield);
                            branch = std::get<0>(llssff);
                            k_tangent = std::get<1>(llssff);
                            f_new = std::get<2>(llssff);

                        } else {
                            f_global[is] = f_yield * sign;
                            d_global[is] = d_yield * sign;
                            d_zero = d_yield * sign - f_yield * sign * (d_yield * sign - d_reload) / (f_yield * sign - f_reload);
                            if ((d_new - d_global[is]) * sign >= 0) {
                                std::tuple<int, float, float>llssff = this->Tslip_120(d_yield, d_new, sign, f_crack, d_crack, k_yield, k_plastic, f_yield);
                                branch = std::get<0>(llssff);
                                k_tangent = std::get<1>(llssff);
                                f_new = std::get<2>(llssff);

                            } else {
                                branch = 9;
                                k_tangent = f_global[is] / (d_global[is] - d_zero);
                                f_new = k_tangent * (d_new - d_zero);

                            }
                        }
                    }
                } else if (abs(d_global[is]) <= d_yield) {
                    k_global = f_global[is] / (d_global[is] - d_zero);
                    k_toYield = f_yield * sign / (d_yield * sign - d_zero);
                    if (k_global < k_toYield) {
                        d_global[is] = d_yield * sign;
                        f_global[is] = f_yield * sign;
                    }
                    if ((d_new - d_global[is]) * sign >= 0) {
                        std::tuple<int, float, float>llssff = this->Tslip_120(d_yield, d_new, sign, f_crack, d_crack, k_yield, k_plastic, f_yield);
                        branch = std::get<0>(llssff);
                        k_tangent = std::get<1>(llssff);
                        f_new = std::get<2>(llssff);

                    } else {
                        branch = 9;
                        k_tangent = f_global[is] / (d_global[is] - d_zero);
                        f_new = k_tangent * (d_new - d_zero);

                    }
                } else {
 // % loop of 440 %%
                    k_global = d_global[3 - is] - f_global[3 - is] / k_unload[3 - is];
                    k_pinch = abs(f_global[is] / (k_global - d_global[is])) * pow( abs(d_global[is] / (k_global - d_global[is])), k_pinch_factor);
                    k_global_from_pinch = f_global[is] / d_global[is] * k_global_factor;
                    if (abs(k_pinch - k_global_from_pinch) < error) {
                        d_pinch = d_zero;
                    } else {
                        d_pinch = (k_global_from_pinch * d_global[is] - k_pinch * d_zero - f_global[is]) / (k_global_from_pinch - k_pinch);
                    }

                    if ((d_new - d_pinch) * sign > 0) {
                        if ((d_new - d_global[is]) * sign >= 0) {
                            std::tuple<int, float, float>llssff = this->Tslip_120(d_yield, d_new, sign, f_crack, d_crack, k_yield, k_plastic, f_yield);
                            branch = std::get<0>(llssff);
                            k_tangent = std::get<1>(llssff);
                            f_new = std::get<2>(llssff);

                        } else {
                            branch = 7;
                            k_global_from_pinch = f_global[is] / d_global[is] * k_global_factor;
                            k_tangent = k_global_from_pinch;
                            d_zero_from_pinch = d_global[is] - f_global[is] / k_global_from_pinch;
                            d_zero = d_zero_from_pinch;
                            f_new = k_global_from_pinch * (d_new - d_zero_from_pinch);

                        }
                    } else {
                        branch = 6;
                        k_global = d_global[3 - is] - f_global[3 - is] / k_unload[3 - is];
                        k_pinch = abs(f_global[is] / (k_global - d_global[is])) * pow( abs(d_global[is] / (k_global - d_global[is])), k_pinch_factor);
                        k_tangent = k_pinch;
                        f_new = k_pinch * (d_new - d_zero);

                    }
                }
            }
        }
    }

 // % rule 4 unloading from peak (d_global,f_global) on the primary curve %%
    if (branch == 4) {
        if ((d_global[is] - d_new) * sign <= 0) {
            std::tuple<int, float, float>llssff = this->Tslip_120(d_yield, d_new, sign, f_crack, d_crack, k_yield, k_plastic, f_yield);
            branch = std::get<0>(llssff);
            k_tangent = std::get<1>(llssff);
            f_new = std::get<2>(llssff);

        } else {
            if ((d_new - d_zero) * sign >= 0) {
                branch = 4;
                k_tangent = k_unload[is];
                f_new = f_global[is] + (d_new - d_global[is]) * k_unload[is];

            } else {
 // % loop of 430 %%
                is = is == 1 ? 2 : 1;
                sign = sign == 1 ? -1 : 1;
                if (abs(d_global[is]) <= d_crack) {
                    d_reload = d_zero + sign * f_crack / k_unload[3 - is];
                    f_reload = f_crack * sign;
                    if ((d_reload - d_new) * sign > 0) {
                        branch = 5;
                        f_new = (d_new - d_zero) * k_tangent;

                    } else {
 // % loop of 530 %%
                        if ((d_global[is] * sign <= d_crack) && (abs(d_global[3 - is]) <= d_yield)) {
                            std::tuple<int, float, float>llssff = this->Tslip_120(d_yield, d_new, sign, f_crack, d_crack, k_yield, k_plastic, f_yield);
                            branch = std::get<0>(llssff);
                            k_tangent = std::get<1>(llssff);
                            f_new = std::get<2>(llssff);

                        } else {
                            f_global[is] = f_yield * sign;
                            d_global[is] = d_yield * sign;
                            d_zero = d_yield * sign - f_yield * sign * (d_yield * sign - d_reload) / (f_yield * sign - f_reload);
                            if ((d_new - d_global[is]) * sign >= 0) {
                                std::tuple<int, float, float>llssff = this->Tslip_120(d_yield, d_new, sign, f_crack, d_crack, k_yield, k_plastic, f_yield);
                                branch = std::get<0>(llssff);
                                k_tangent = std::get<1>(llssff);
                                f_new = std::get<2>(llssff);

                            } else {
                                branch = 9;
                                k_tangent = f_global[is] / (d_global[is] - d_zero);
                                f_new = k_tangent * (d_new - d_zero);

                            }
                        }
                    }
                } else if (abs(d_global[is]) <= d_yield) {
                    k_global = f_global[is] / (d_global[is] - d_zero);
                    k_toYield = f_yield * sign / (d_yield * sign - d_zero);
                    if (k_global < k_toYield) {
                        d_global[is] = d_yield * sign;
                        f_global[is] = f_yield * sign;
                    }
                    if ((d_new - d_global[is]) * sign >= 0) {
                        std::tuple<int, float, float>llssff = this->Tslip_120(d_yield, d_new, sign, f_crack, d_crack, k_yield, k_plastic, f_yield);
                        branch = std::get<0>(llssff);
                        k_tangent = std::get<1>(llssff);
                        f_new = std::get<2>(llssff);

                    } else {
                        branch = 9;
                        k_tangent = f_global[is] / (d_global[is] - d_zero);
                        f_new = k_tangent * (d_new - d_zero);

                    }
                } else {
 // % loop of 440 %%
                    k_global = d_global[3 - is] - f_global[3 - is] / k_unload[3 - is];
                    k_pinch = abs(f_global[is] / (k_global - d_global[is])) * pow( abs(d_global[is] / (k_global - d_global[is])), k_pinch_factor);
                    k_global_from_pinch = f_global[is] / d_global[is] * k_global_factor;
                    if (abs(k_pinch - k_global_from_pinch) < error) {
                        d_pinch = d_zero;
                    } else {
                        d_pinch = (k_global_from_pinch * d_global[is] - k_pinch * d_zero - f_global[is]) / (k_global_from_pinch - k_pinch);
                    }

                    if ((d_new - d_pinch) * sign > 0) {
                        if ((d_new - d_global[is]) * sign >= 0) {
                            std::tuple<int, float, float>llssff = this->Tslip_120(d_yield, d_new, sign, f_crack, d_crack, k_yield, k_plastic, f_yield);
                            branch = std::get<0>(llssff);
                            k_tangent = std::get<1>(llssff);
                            f_new = std::get<2>(llssff);

                        } else {
                            branch = 7;
                            k_global_from_pinch = f_global[is] / d_global[is] * k_global_factor;
                            k_tangent = k_global_from_pinch;
                            d_zero_from_pinch = d_global[is] - f_global[is] / k_global_from_pinch;
                            d_zero = d_zero_from_pinch;
                            f_new = k_global_from_pinch * (d_new - d_zero_from_pinch);

                        }
                    } else {
                        branch = 6;
                        k_global = d_global[3 - is] - f_global[3 - is] / k_unload[3 - is];
                        k_pinch = abs(f_global[is] / (k_global - d_global[is])) * pow( abs(d_global[is] / (k_global - d_global[is])), k_pinch_factor);
                        k_tangent = k_pinch;
                        f_new = k_pinch * (d_new - d_zero);

                    }
                }
            }
        }
    }

 // % rule 5 load reversed at zeor - crossing paint (d_zero,0) withyout %%
 // % previous cracking in a new direction                        %%
    if (branch == 5) {
        if ((d_new - d_zero) * sign > 0) {
            if ((d_reload - d_new) * sign > 0) {
                        branch = 5;
                f_new = (d_new - d_zero) * k_tangent;

            } else {
 // % loop of 530 %%
                if ((d_global[is] * sign <= d_crack) && (abs(d_global[3 - is]) <= d_yield)) {
                    std::tuple<int, float, float>llssff = this->Tslip_120(d_yield, d_new, sign, f_crack, d_crack, k_yield, k_plastic, f_yield);
                    branch = std::get<0>(llssff);
                    k_tangent = std::get<1>(llssff);
                    f_new = std::get<2>(llssff);

                } else {
                    f_global[is] = f_crack * sign;
                    d_global[is] = d_crack * sign;
 // %d_zero = d_yield * sign - f_yield * sign * (d_yield * sign - d_reload) / (f_yield * sign - f_reload);%%%%�o������ψʂƙ��f�͂𐳕��̂ǂ�����Ђъ���_�𒴂���܂łЂъ���_���L������悤�ɂ���
                    if ((d_new - d_global[is]) * sign >= 0) {
                        std::tuple<int, float, float>llssff = this->Tslip_120(d_yield, d_new, sign, f_crack, d_crack, k_yield, k_plastic, f_yield);
                        branch = std::get<0>(llssff);
                        k_tangent = std::get<1>(llssff);
                        f_new = std::get<2>(llssff);

                    } else {
                        branch = 9;
                        k_tangent = f_global[is] / (d_global[is] - d_zero);
                        f_new = k_tangent * (d_new - d_zero);

                    }
                }
            }
        } else {
            is = is == 1 ? 2 : 1;
            sign = sign == 1 ? -1 : 1;
            if ((d_global[is] - d_new) * sign > 0) {
                branch = 4;
                k_tangent = k_unload[is];
                f_new = f_global[is] + (d_new - d_global[is]) * k_unload[is];

            } else {
                std::tuple<int, float, float>llssff = this->Tslip_120(d_yield, d_new, sign, f_crack, d_crack, k_yield, k_plastic, f_yield);
                branch = std::get<0>(llssff);
                k_tangent = std::get<1>(llssff);
                f_new = std::get<2>(llssff);

            }
        }
    }

 // % rule 6 reloading with soft spring toward stiffness %%
 // % changing point d_pinch                                  %%
    if (branch == 6) {
        if ((d_new - d_old) * sign > 0) {
            if ((d_pinch - d_new) * sign > 0) {
                k_global = d_global[3 - is] - f_global[3 - is] / k_unload[3 - is];
                k_pinch = abs(f_global[is] / (k_global - d_global[is])) * pow( abs(d_global[is] / (k_global - d_global[is])), k_pinch_factor);
                f_new = k_pinch * (d_new - d_zero);

            } else {
                if ((d_new - d_global[is]) * sign > 0) {
                    branch = 3;
                    k_tangent = k_plastic;
                    f_new = f_yield * sign + (d_new - d_yield * sign) * k_plastic;

                } else {
                    branch = 7;
                    k_global_from_pinch = f_global[is] / d_global[is] * k_global_factor;
                    k_tangent = k_global_from_pinch;
                    d_zero_from_pinch = d_global[is] - f_global[is] / k_global_from_pinch;
                    d_zero = d_zero_from_pinch;
                    f_new = k_global_from_pinch * (d_new - d_zero_from_pinch);

                }
            }
        } else {
 // % loop of 630 %%
            d_local = d_old;
            f_local = f_old;
            k_local = k_unload[is] * unload_from_local_factor;
            d_zero_from_local = d_local - f_local / k_local;
            if ((d_zero_from_local - d_new) * sign < 0) {
                branch = 8;
                k_tangent = k_local;
                f_new = f_local + (d_new - d_local) * k_local;

            } else {
 // % loop of 840 %%
                d_zero = d_zero_from_local;
                is = is == 1 ? 2 : 1;
                sign = sign == 1 ? -1 : 1;
                if (abs(d_global[is]) >= d_yield) {
                    k_global_from_pinch = f_global[is] / d_global[is] * k_global_factor;
                    d_zero_from_pinch = d_global[is] - f_global[is] / k_global_from_pinch;
                    d_pinch = d_zero;
                    if ((d_zero_from_pinch - d_zero) * sign <= 0) {
                        if ((d_new - d_global[is]) * sign >= 0) {
                            std::tuple<int, float, float>llssff = this->Tslip_120(d_yield, d_new, sign, f_crack, d_crack, k_yield, k_plastic, f_yield);
                            branch = std::get<0>(llssff);
                            k_tangent = std::get<1>(llssff);
                            f_new = std::get<2>(llssff);

                        } else {
                            branch = 9;
                            k_tangent = f_global[is] / (d_global[is] - d_zero);
                            f_new = k_tangent * (d_new - d_zero);

                        }
                    } else {
 // % loop of 440 %%
                        k_global = d_global[3 - is] - f_global[3 - is] / k_unload[3 - is];
                        k_pinch = abs(f_global[is] / (k_global - d_global[is])) * pow( abs(d_global[is] / (k_global - d_global[is])), k_pinch_factor);
                        k_global_from_pinch = f_global[is] / d_global[is] * k_global_factor;
                        if (abs(k_pinch - k_global_from_pinch) < error) {
                            d_pinch = d_zero;
                        } else {
                            d_pinch = (k_global_from_pinch * d_global[is] - k_pinch * d_zero - f_global[is]) / (k_global_from_pinch - k_pinch);
                        }

                        if ((d_new - d_pinch) * sign > 0) {
                            if ((d_new - d_global[is]) * sign >= 0) {
                                std::tuple<int, float, float>llssff = this->Tslip_120(d_yield, d_new, sign, f_crack, d_crack, k_yield, k_plastic, f_yield);
                                branch = std::get<0>(llssff);
                                k_tangent = std::get<1>(llssff);
                                f_new = std::get<2>(llssff);

                            } else {
                                branch = 7;
                                k_global_from_pinch = f_global[is] / d_global[is] * k_global_factor;
                                k_tangent = k_global_from_pinch;
                                d_zero_from_pinch = d_global[is] - f_global[is] / k_global_from_pinch;
                                d_zero = d_zero_from_pinch;
                                f_new = k_global_from_pinch * (d_new - d_zero_from_pinch);

                            }
                        } else {
                            branch = 6;
                            k_global = d_global[3 - is] - f_global[3 - is] / k_unload[3 - is];
                            k_pinch = abs(f_global[is] / (k_global - d_global[is])) * pow( abs(d_global[is] / (k_global - d_global[is])), k_pinch_factor);
                            k_tangent = k_pinch;
                            f_new = k_pinch * (d_new - d_zero);

                        }
                    }
                } else {
                    if ((d_new - d_global[is]) * sign >= 0) {
                        std::tuple<int, float, float>llssff = this->Tslip_120(d_yield, d_new, sign, f_crack, d_crack, k_yield, k_plastic, f_yield);
                        branch = std::get<0>(llssff);
                        k_tangent = std::get<1>(llssff);
                        f_new = std::get<2>(llssff);

                    } else {
                        branch = 9;
                        k_tangent = f_global[is] / (d_global[is] - d_zero);
                        f_new = k_tangent * (d_new - d_zero);

                    }
                }
            }
        }
    }

 // % rule 7 loading with hard spring toward previous yield point %%
    if (branch == 7) {
        if ((d_new - d_old) * sign > 0) {
            if ((d_new - d_global[is]) * sign > 0) {
                branch = 3;
                k_tangent = k_plastic;
                f_new = f_yield * sign + (d_new - d_yield * sign) * k_plastic;

            } else {
                branch = 7;
                k_global_from_pinch = f_global[is] / d_global[is] * k_global_factor;
                d_zero_from_pinch = d_global[is] - f_global[is] / k_global_from_pinch;
                f_new = k_global_from_pinch * (d_new - d_zero_from_pinch);
// Important Addition
                k_tangent = k_global_from_pinch;
                d_zero = d_zero_from_pinch;

            }
        } else {
 // % loop of 630 %%
            d_local = d_old;
            f_local = f_old;
            k_local = k_unload[is] * unload_from_local_factor;
            d_zero_from_local = d_local - f_local / k_local;
            if ((d_zero_from_local - d_new) * sign < 0) {
                branch = 8;
                k_tangent = k_local;
                f_new = f_local + (d_new - d_local) * k_local;

            } else {
 // % loop of 840 %%
                d_zero = d_zero_from_local;
                is = is == 1 ? 2 : 1;
                sign = sign == 1 ? -1 : 1;
                if (abs(d_global[is]) >= d_yield) {
                    k_global_from_pinch = f_global[is] / d_global[is] * k_global_factor;
                    d_zero_from_pinch = d_global[is] - f_global[is] / k_global_from_pinch;
                    d_pinch = d_zero;
                    if ((d_zero_from_pinch - d_zero) * sign <= 0) {
                        if ((d_new - d_global[is]) * sign >= 0) {
                            std::tuple<int, float, float>llssff = this->Tslip_120(d_yield, d_new, sign, f_crack, d_crack, k_yield, k_plastic, f_yield);
                            branch = std::get<0>(llssff);
                            k_tangent = std::get<1>(llssff);
                            f_new = std::get<2>(llssff);

                        } else {
                            branch = 9;
                            k_tangent = f_global[is] / (d_global[is] - d_zero);
                            f_new = k_tangent * (d_new - d_zero);

                        }
                    } else {
 // % loop of 440 %%
                        k_global = d_global[3 - is] - f_global[3 - is] / k_unload[3 - is];
                        k_pinch = abs(f_global[is] / (k_global - d_global[is])) * pow( abs(d_global[is] / (k_global - d_global[is])), k_pinch_factor);
                        k_global_from_pinch = f_global[is] / d_global[is] * k_global_factor;
                        if (abs(k_pinch - k_global_from_pinch) < error) {
                            d_pinch = d_zero;
                        } else {
                            d_pinch = (k_global_from_pinch * d_global[is] - k_pinch * d_zero - f_global[is]) / (k_global_from_pinch - k_pinch);
                        }

                        if ((d_new - d_pinch) * sign > 0) {
                            if ((d_new - d_global[is]) * sign >= 0) {
                                std::tuple<int, float, float>llssff = this->Tslip_120(d_yield, d_new, sign, f_crack, d_crack, k_yield, k_plastic, f_yield);
                                branch = std::get<0>(llssff);
                                k_tangent = std::get<1>(llssff);
                                f_new = std::get<2>(llssff);

                            } else {
                                branch = 7;
                                k_global_from_pinch = f_global[is] / d_global[is] * k_global_factor;
                                k_tangent = k_global_from_pinch;
                                d_zero_from_pinch = d_global[is] - f_global[is] / k_global_from_pinch;
                                d_zero = d_zero_from_pinch;
                                f_new = k_global_from_pinch * (d_new - d_zero_from_pinch);

                            }
                        } else {
                            branch = 6;
                            k_global = d_global[3 - is] - f_global[3 - is] / k_unload[3 - is];
                            k_pinch = abs(f_global[is] / (k_global - d_global[is])) * pow( abs(d_global[is] / (k_global - d_global[is])), k_pinch_factor);
                            k_tangent = k_pinch;
                            f_new = k_pinch * (d_new - d_zero);

                        }
                    }
                } else {
                    if ((d_new - d_global[is]) * sign >= 0) {
                        std::tuple<int, float, float>llssff = this->Tslip_120(d_yield, d_new, sign, f_crack, d_crack, k_yield, k_plastic, f_yield);
                        branch = std::get<0>(llssff);
                        k_tangent = std::get<1>(llssff);
                        f_new = std::get<2>(llssff);

                    } else {
                        branch = 9;
                        k_tangent = f_global[is] / (d_global[is] - d_zero);
                        f_new = k_tangent * (d_new - d_zero);

                    }
                }
            }
        }
    }

 // % rule 8 unloading from inner peak (d_local,f_local) %%
    if (branch == 8) {
        if ((d_zero_from_local - d_new) * sign < 0) {
            if ((d_local - d_new) * sign > 0) {
                f_new = f_local + (d_new - d_local) * k_local;

            } else {
// % loop of 830 %%
                if (abs(d_global[is]) > d_yield) {
                    if ((d_new - d_pinch) * sign <= 0) {
                        branch = 6;
                        k_global = d_global[3 - is] - f_global[3 - is] / k_unload[3 - is];
                        k_pinch = abs(f_global[is] / (k_global - d_global[is])) * pow( abs(d_global[is] / (k_global - d_global[is])), k_pinch_factor);
                        k_tangent = k_pinch;
                        f_new = k_pinch * (d_new - d_zero);

                    } else {
                        if ((d_new - d_global[is]) * sign >= 0) {
                            std::tuple<int, float, float>llssff = this->Tslip_120(d_yield, d_new, sign, f_crack, d_crack, k_yield, k_plastic, f_yield);
                            branch = std::get<0>(llssff);
                            k_tangent = std::get<1>(llssff);
                            f_new = std::get<2>(llssff);

                        } else {
                            branch = 9;
                            k_tangent = f_global[is] / (d_global[is] - d_zero);
                            f_new = k_tangent * (d_new - d_zero);

                        }
                    }
                } else {
                    if ((d_new - d_global[is]) * sign >= 0) {
                            std::tuple<int, float, float>llssff = this->Tslip_120(d_yield, d_new, sign, f_crack, d_crack, k_yield, k_plastic, f_yield);
                            branch = std::get<0>(llssff);
                            k_tangent = std::get<1>(llssff);
                            f_new = std::get<2>(llssff);

                    } else {
                            branch = 9;
                            k_tangent = f_global[is] / (d_global[is] - d_zero);
                            f_new = k_tangent * (d_new - d_zero);

                        }
                }
            }
        } else {
 // % loop of 840 %%
            d_zero = d_zero_from_local;
            is = is == 1 ? 2 : 1;
            sign = sign == 1 ? -1 : 1;
            if (abs(d_global[is]) >= d_yield) {
                k_global_from_pinch = f_global[is] / d_global[is] * k_global_factor;
                d_zero_from_pinch = d_global[is] - f_global[is] / k_global_from_pinch;
                d_pinch = d_zero;
                if ((d_zero_from_pinch - d_zero) * sign <= 0) {
                    if ((d_new - d_global[is]) * sign >= 0) {
                        std::tuple<int, float, float>llssff = this->Tslip_120(d_yield, d_new, sign, f_crack, d_crack, k_yield, k_plastic, f_yield);
                        branch = std::get<0>(llssff);
                        k_tangent = std::get<1>(llssff);
                        f_new = std::get<2>(llssff);

                    } else {
                        branch = 9;
                        k_tangent = f_global[is] / (d_global[is] - d_zero);
                        f_new = k_tangent * (d_new - d_zero);

                    }
                } else {
 // % loop of 440 %%
                    k_global = d_global[3 - is] - f_global[3 - is] / k_unload[3 - is];
                    k_pinch = abs(f_global[is] / (k_global - d_global[is])) * pow( abs(d_global[is] / (k_global - d_global[is])), k_pinch_factor);
                    k_global_from_pinch = f_global[is] / d_global[is] * k_global_factor;
                    if (abs(k_pinch - k_global_from_pinch) < error) {
                        d_pinch = d_zero;
                    } else {
                        d_pinch = (k_global_from_pinch * d_global[is] - k_pinch * d_zero - f_global[is]) / (k_global_from_pinch - k_pinch);
                    }

                    if ((d_new - d_pinch) * sign > 0) {
                        if ((d_new - d_global[is]) * sign >= 0) {
                            std::tuple<int, float, float>llssff = this->Tslip_120(d_yield, d_new, sign, f_crack, d_crack, k_yield, k_plastic, f_yield);
                            branch = std::get<0>(llssff);
                            k_tangent = std::get<1>(llssff);
                            f_new = std::get<2>(llssff);

                        } else {
                            branch = 7;
                            k_global_from_pinch = f_global[is] / d_global[is] * k_global_factor;
                            k_tangent = k_global_from_pinch;
                            d_zero_from_pinch = d_global[is] - f_global[is] / k_global_from_pinch;
                            d_zero = d_zero_from_pinch;
                            f_new = k_global_from_pinch * (d_new - d_zero_from_pinch);

                        }
                    } else {
                        branch = 6;
                        k_global = d_global[3 - is] - f_global[3 - is] / k_unload[3 - is];
                        k_pinch = abs(f_global[is] / (k_global - d_global[is])) * pow( abs(d_global[is] / (k_global - d_global[is])), k_pinch_factor);
                        k_tangent = k_pinch;
                        f_new = k_pinch * (d_new - d_zero);

                    }
                }
            } else {
                if ((d_new - d_global[is]) * sign >= 0) {
                    std::tuple<int, float, float>llssff = this->Tslip_120(d_yield, d_new, sign, f_crack, d_crack, k_yield, k_plastic, f_yield);
                    branch = std::get<0>(llssff);
                    k_tangent = std::get<1>(llssff);
                    f_new = std::get<2>(llssff);

                } else {
                    branch = 9;
                    k_tangent = f_global[is] / (d_global[is] - d_zero);
                    f_new = k_tangent * (d_new - d_zero);

                }
            }
        }
    }

 // % rule 9 reloading toward peak (d2,f2) directly %%
    if (branch == 9) {
    if ((d_new - d_old) * sign > 0) {
            if ((d_new - d_global[is]) * sign > 0) {
                std::tuple<int, float, float>llssff = this->Tslip_120(d_yield, d_new, sign, f_crack, d_crack, k_yield, k_plastic, f_yield);
                branch = std::get<0>(llssff);
                k_tangent = std::get<1>(llssff);
                f_new = std::get<2>(llssff);

            } else {
                f_new = k_tangent * (d_new - d_zero);

            }
    } else {
 // % loop of 630 %%
            d_local = d_old;
            f_local = f_old;
            k_local = k_unload[is] * unload_from_local_factor;
            d_zero_from_local = d_local - f_local / k_local;
            if ((d_zero_from_local - d_new) * sign < 0) {
                branch = 8;
                k_tangent = k_local;
                f_new = f_local + (d_new - d_local) * k_local;

            } else {
 // % loop of 840 %%
                d_zero = d_zero_from_local;
                is = is == 1 ? 2 : 1;
                sign = sign == 1 ? -1 : 1;
                if (abs(d_global[is]) >= d_yield) {
                    k_global_from_pinch = f_global[is] / d_global[is] * k_global_factor;
                    d_zero_from_pinch = d_global[is] - f_global[is] / k_global_from_pinch;
                    d_pinch = d_zero;
                    if ((d_zero_from_pinch - d_zero) * sign <= 0) {
                        if ((d_new - d_global[is]) * sign >= 0) {
                            std::tuple<int, float, float>llssff = this->Tslip_120(d_yield, d_new, sign, f_crack, d_crack, k_yield, k_plastic, f_yield);
                            branch = std::get<0>(llssff);
                            k_tangent = std::get<1>(llssff);
                            f_new = std::get<2>(llssff);

                        } else {
                            branch = 9;
                            k_tangent = f_global[is] / (d_global[is] - d_zero);
                            f_new = k_tangent * (d_new - d_zero);

                        }
                    } else {
 // % loop of 440 %%
                        k_global = d_global[3 - is] - f_global[3 - is] / k_unload[3 - is];
                        k_pinch = abs(f_global[is] / (k_global - d_global[is])) * pow( abs(d_global[is] / (k_global - d_global[is])), k_pinch_factor);
                        k_global_from_pinch = f_global[is] / d_global[is] * k_global_factor;
                        if (abs(k_pinch - k_global_from_pinch) < error) {
                            d_pinch = d_zero;
                        } else {
                            d_pinch = (k_global_from_pinch * d_global[is] - k_pinch * d_zero - f_global[is]) / (k_global_from_pinch - k_pinch);
                        }

                        if ((d_new - d_pinch) * sign > 0) {
                            if ((d_new - d_global[is]) * sign >= 0) {
                                std::tuple<int, float, float>llssff = this->Tslip_120(d_yield, d_new, sign, f_crack, d_crack, k_yield, k_plastic, f_yield);
                                branch = std::get<0>(llssff);
                                k_tangent = std::get<1>(llssff);
                                f_new = std::get<2>(llssff);

                            } else {
                                branch = 7;
                                k_global_from_pinch = f_global[is] / d_global[is] * k_global_factor;
                                k_tangent = k_global_from_pinch;
                                d_zero_from_pinch = d_global[is] - f_global[is] / k_global_from_pinch;
                                d_zero = d_zero_from_pinch;
                                f_new = k_global_from_pinch * (d_new - d_zero_from_pinch);

                            }
                        } else {
                            branch = 6;
                            k_global = d_global[3 - is] - f_global[3 - is] / k_unload[3 - is];
                            k_pinch = abs(f_global[is] / (k_global - d_global[is])) * pow( abs(d_global[is] / (k_global - d_global[is])), k_pinch_factor);
                            k_tangent = k_pinch;
                            f_new = k_pinch * (d_new - d_zero);

                        }
                    }
                } else {
                    if ((d_new - d_global[is]) * sign >= 0) {
                        std::tuple<int, float, float>llssff = this->Tslip_120(d_yield, d_new, sign, f_crack, d_crack, k_yield, k_plastic, f_yield);
                        branch = std::get<0>(llssff);
                        k_tangent = std::get<1>(llssff);
                        f_new = std::get<2>(llssff);

                    } else {
                        branch = 9;
                        k_tangent = f_global[is] / (d_global[is] - d_zero);
                        f_new = k_tangent * (d_new - d_zero);

                    }
                }
            }
        }
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
    cf_reload = f_reload;
    cf_local = f_local;
    cd_reload = d_reload;
    cd_local = d_local;
    cf_global = f_global;
    cd_global = d_global;
    cd_pinch = d_pinch;
    ck_unload = k_unload;
    ck_local = k_local;
    return 0;
}

int TakedaSlip::revertToLastCommit(void)
{
    branch = cbranch;
    d_new = cd_new;
    f_new = cf_new;
    k_tangent = ck_tangent;
    f_reload = cf_reload;
    f_local = cf_local;
    d_reload = cd_reload;
    d_local = cd_local;
    f_global = cf_global;
    d_global = cd_global;
    d_pinch = cd_pinch;
    k_unload = ck_unload;
    k_local = ck_local;
    return 0;
}

int TakedaSlip::revertToStart(void)
{
    branch = cbranch = 1;
    k_tangent = ck_tangent = k_crack;
    f_crack = d_crack * k_crack;
    f_yield = f_crack + k_yield * (d_yield - d_crack);
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
    theCopy->f_reload = f_reload;
    theCopy->f_local = f_local;
    theCopy->d_reload = d_reload;
    theCopy->d_local = d_local;
    theCopy->f_global = f_global;
    theCopy->d_global = d_global;
    theCopy->d_pinch = d_pinch;
    theCopy->k_unload = k_unload;
    theCopy->k_local = k_local;

    theCopy->cbranch = cbranch;
    theCopy->cd_new = cd_new;
    theCopy->cf_new = cf_new;
    theCopy->ck_tangent = ck_tangent;
    theCopy->cf_reload = cf_reload;
    theCopy->cf_local = cf_local;
    theCopy->cd_reload = cd_reload;
    theCopy->cd_local = cd_local;
    theCopy->cf_global = cf_global;
    theCopy->cd_global = cd_global;
    theCopy->cd_pinch = cd_pinch;
    theCopy->ck_unload = ck_unload;
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
    data(17) = f_reload;
    data(18) = f_local;
    data(19) = d_reload;
    data(20) = d_local;
    data(21) = f_global[1];
    data(22) = f_global[2];
    data(23) = d_global[1];
    data(24) = d_global[2];
    data(25) = d_pinch;
    data(26) = k_unload[1];
    data(27) = k_unload[2];
    data(28) = k_local;

    data(111) = cbranch;
    data(113) = cd_new;
    data(115) = cf_new;
    data(116) = ck_tangent;
    data(117) = cf_reload;
    data(118) = cf_local;
    data(119) = cd_reload;
    data(120) = cd_local;
    data(121) = cf_global[1];
    data(122) = cf_global[2];
    data(123) = cd_global[1];
    data(124) = cd_global[2];
    data(125) = cd_pinch;
    data(126) = ck_unload[1];
    data(127) = ck_unload[2];
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
        f_reload = data(17);
        f_local = data(18);
        d_reload = data(19);
        d_local = data(20);
        f_global[1] = data(21);
        f_global[2] = data(22);
        d_global[1] = data(23);
        d_global[2] = data(24);
        d_pinch = data(25);
        k_unload[1] = data(26);
        k_unload[2] = data(27);
        k_local = data(28);

        cbranch = data(111);
        cd_new = data(113);
        cf_new = data(115);
        ck_tangent = data(116);
        cf_reload = data(117);
        cf_local = data(118);
        cd_reload = data(119);
        cd_local = data(120);
        cf_global[1] = data(121);
        cf_global[2] = data(122);
        cd_global[1] = data(123);
        cd_global[2] = data(124);
        cd_pinch = data(125);
        ck_unload[1] = data(126);
        ck_unload[2] = data(127);
        ck_local = data(128);
    }

    return res;
}

void TakedaSlip::Print(OPS_Stream &s, int flag)
{
    cout << "TakedaSlip tag: " << this->getTag() << endln;
}