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

#include <math.h>
#include <IMKBilin.h>
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

static int numIMKBilinMaterials	= 0;

void *
OPS_IMKBilin(void)
{
    if (numIMKBilinMaterials == 0) {
        numIMKBilinMaterials++;
        OPS_Error("Mod. IMK Bilinear Model - AE-Aug22\n", 1);
    }

    // Pointer to a uniaxial material that will be returned
    UniaxialMaterial *theMaterial	= 0;

    int    iData[1];
    double	dData[21];
    int numData	= 1;

    if (OPS_GetIntInput(&numData, iData) != 0) {
        opserr << "WARNING invalid uniaxialMaterial IMKBilin tag" << endln;
        return 0;
    }

    numData	= 21;

    if (OPS_GetDoubleInput(&numData, dData) != 0) {
        opserr << "Invalid Args want: uniaxialMaterial IMKBilin tag? Ke? ";
        opserr << "Up_pos? Upc_pos? Uu_pos? Mpe_pos? MmaxMpe_pos? ResM_pos? ";
        opserr << "Up_neg? Upc_neg? Uu_neg? Mpe_neg? MmaxMpe_neg? ResM_neg? ";
        opserr << "LamdaS?  LamdaC? LamdaK? Cs? Cc? Ck? D_pos? D_neg? ";
        return 0;
    }


    // Parsing was successful, allocate the material
    theMaterial	= new IMKBilin(iData[0],
        dData[0],
        dData[1], dData[2], dData[3], dData[4], dData[5], dData[6],
        dData[7], dData[8], dData[9], dData[10], dData[11], dData[12],
        dData[13], dData[14], dData[15], dData[16], dData[17], dData[18], dData[19], dData[20]);

    if (theMaterial == 0) {
        opserr << "WARNING could not create uniaxialMaterial of type IMKBilin Material\n";
        return 0;
    }

    return theMaterial;
}

IMKBilin::IMKBilin(int tag, double	p_Ke,
    double	p_posUp_0, double	p_posUpc_0, double	p_posUu_0, double	p_posMpe_0, double	p_posMmaxMpe_0, double	p_posResM_0,
    double	p_negUp_0, double	p_negUpc_0, double	p_negUu_0, double	p_negMpe_0, double	p_negMmaxMpe_0, double	p_negResM_0,
    double	p_LAMBDA_S, double	p_LAMBDA_C, double	p_LAMBDA_K, double	p_c_S, double	p_c_C, double	p_c_K, double	p_D_pos, double	p_D_neg)
    :UniaxialMaterial(tag, 0), Ke(p_Ke),
    posUp_0(p_posUp_0), posUpc_0(p_posUpc_0), posUu_0(p_posUu_0), posMpe_0(p_posMpe_0), posMmaxMpe_0(p_posMmaxMpe_0), posResM_0(p_posResM_0),
    negUp_0(p_negUp_0), negUpc_0(p_negUpc_0), negUu_0(p_negUu_0), negMpe_0(p_negMpe_0), negMmaxMpe_0(p_negMmaxMpe_0), negResM_0(p_negResM_0),
    LAMBDA_S(p_LAMBDA_S), LAMBDA_C(p_LAMBDA_C), LAMBDA_K(p_LAMBDA_K), c_S(p_c_S), c_C(p_c_C), c_K(p_c_K), D_pos(p_D_pos), D_neg(p_D_neg)
{
    this->revertToStart();
}

IMKBilin::IMKBilin()
    :UniaxialMaterial(0, 0), Ke(0),
    posUp_0(0), posUpc_0(0), posUu_0(0), posMpe_0(0), posMmaxMpe_0(0), posResM_0(0),
    negUp_0(0), negUpc_0(0), negUu_0(0), negMpe_0(0), negMmaxMpe_0(0), negResM_0(0),
    LAMBDA_S(0), LAMBDA_C(0), LAMBDA_K(0), c_S(0), c_C(0), c_K(0), D_pos(0), D_neg(0)
{
    this->revertToStart();
}

IMKBilin::~IMKBilin()
{
    // does nothing
}

int IMKBilin::setTrialStrain(double	strain, double	strainRate)
{
    //all variables to the last commit
    this->revertToLastCommit();

    //state determination algorithm: defines the current force and tangent stiffness
    U	    = strain; //set trial displacement
    Ui_1	= Ui;
    Fi_1	= Fi;
    Di_1	= Di;
    Ui	    = U;

    double dU   = Ui - Ui_1;
    double dEi;
    if (dU==0) {
        Fi  = Fi_1;
        dEi = 0;
    }
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////  MAIN CODE //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // %%%%%%%% INITIALIZE CURRENT BACKBONE VALUES AS PREVIOUS %%%%%%%%%%%%
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    double	posMmax	    = posMmax_1;
    double	posMpe	    = posMpe_1;
    double	posMpeProj	= posMpeProj_1;
    double	posMmaxProj	= posMmaxProj_1;
    double	posKp	    = posKp_1;
    double	posKpc	    = posKpc_1;
    double	posUy	    = posUy_1;
    double	posUmax	    = posUmax_1;

    double	negMmax	    = negMmax_1;
    double	negMpe	    = negMpe_1;
    double	negMpeProj	= negMpeProj_1;
    double	negMmaxProj	= negMmaxProj_1;
    double	negKp	    = negKp_1;
    double	negKpc	    = negKpc_1;
    double	negUy	    = negUy_1;
    double	negUmax	    = negUmax_1;

    double	Mi_boundary	= 0.0;

    double	QuarterFlag, Rintrsct_K, DISP_Rev;
    double	betaS, betaC, betaK;
    double	K_j;
    double	Ki, Kpi, Kpci, Umaxi, MpeProji, MmaxProji;

    double	Mi_temp;

    QuarterFlag	    = 0;
    Reversal_Flag	= 0;

        ///////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////
        ////////////////// BRANCH DETERMINATION AND FLAG RAISE ////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////
    //  Branch
    //      0:  Elastic
    //      1:  Unloading Branch

    //      5:  Towards Capping Point   +
    //      6:  Towards Residual Point  +
    //      7:  Residual Branch         +

    //      15: Towards Capping Point   -
    //      16: Towards Residual Point  -
    //      17: Residual Branch         -

    //  Flag
    //      Yield_Flag:     Preserved.      When the deformation exceeds yield capacity for the first time.
    //      Excursion_Flag: Not preserved.  When crossing X-axis. Evokes re-considering of the deteriorations and which peak to go for.
    //      Reversal_Flag:  Not preserved.  When unloading starts. Evokes re-condiersing of the stiffness deterioration and peak point registration.

    // Find the direction of the current increment "Di": 1:if Ui is moving right    and -1: if Ui is moving left
    if (Ui >= Ui_1) {
        Di	= 1;
    }
    else {
        Di	= -1;
    }

    //  Simple Notation for current step parameters
    Fi	= Fi_1 + K_j_1 * dU;

    //  Check for Fail Flag
    if ((Ui >= posUu_0)) {
        Fail_FlagPos	= 1;
    }
    if ((Ui <= -negUu_0)) {
        Fail_FlagNeg	= 1;
    }

    //  Get Information before first Yield
    if (((Fi >= posMpe_0) || (Fi <= -negMpe_0)) && (Yield_Flag == 0)) {
        Yield_Flag	= 1;
    }

    //  Check if previous point was a reversal point
    if (Di_1 / Di < 0) {
        Reversal_Flag	= 1;
        Ulocal	= Ui_1;
        Flocal	= Fi_1;
    }

    // Update loading / unloading stiffness at load reversals
    if (Reversal_Flag == 1) {
        Rintrsct_K  = Ulocal - Flocal / K_j_1;
        DISP_Rev    = engTotl - engExcr_1 - 0.5*Flocal *(Rintrsct_K - Ulocal);
        betaK       = pow((DISP_Rev / (2 * refEnergyK - engTotl + 0.5*Flocal * (Rintrsct_K - Ulocal))), c_K);

        K_j	        = K_j_1 * (1 - betaK);
        if ((Mrpos_Flag == 1) || (Mrneg_Flag == 1)) {
            K_j	    = 0.5*Ke;
        }
    }
    else {
        betaK	= betaK_1;
        K_j	        = K_j_1;
    }

        ///////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////
        /////////////////// UPDATE BACKBONE CURVE /////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////
        // Update Positive Backbone and Target Peak Point

    //  Calculate Backbone parameters at current excursion based on Energy Dissipated in the previous Excursion
    if (Excursion_Flag == 1.0) {
        betaS	= pow((engExcr / (refEnergyS - engTotl)), c_S);
        betaC	= pow((engExcr / (refEnergyC - engTotl)), c_C);

        if (Ui > Ui_1) {
            //  Update Mpe, Mmax Projion, Kp, and Kpc for Current Step
            posMpe	    = posMpe_1	    * (1.0 - betaS * D_pos);
            posKp	    = posKp_1	    * (1.0 - betaS * D_pos);
            posMmaxProj	= posMmaxProj_1	* (1.0 - betaC * D_pos);
            if (posMr_0 == 0.0) {
                posKpc	= posKpc_0 * (posMmaxProj - posMr_0) / posMmaxProj;
            }
            else {
                posKpc	= posKpc_0 * (posMpe - posMr_0) / (posMpe_0 - posMr_0);
            }

            //  Calculate Rotation at Capping Point (Intersection of the two Slopes)
            posUy	        = posMpe / K_j;
            posMpeProj	    = posMpe - posKp * posUy;
            posUmax	        = fabs((posMmaxProj - posMpeProj) / (posKpc + posKp));
            posMmax	        = posMpeProj + posUmax * posKp;

            if ((posMmax - posMr_0) / (posUmax + fabs(Ui)- posMr_0/K_j) < posKp) {
                posKp	    = (posMmax - posMr_0) / (posUmax + fabs(Ui) - posMr_0 / K_j);
                posMpeProj	= posMpe - posKp * posUy;
                posUmax	    = fabs((posMmaxProj - posMpeProj) / (posKpc + posKp));
                posMmax	    = posMpeProj + posUmax * posKp;
            }
        }
        else {
            //  Update Mpe, Mmax Projion, Kp, and Kpc for Current Step
            negMpe	    = negMpe_1 * (1.0 - betaS * D_neg);
            negMmaxProj	= negMmaxProj_1 * (1.0 - betaC * D_neg);
            negKp	= negKp_1 * (1.0 - betaS * D_neg);
            if (negMr_0 == 0.0) {
                negKpc	= negKpc_0 * (negMmaxProj - negMr_0) / negMmaxProj;
            }
            else {
                negKpc	= negKpc_0 * (negMpe - negMr_0) / (negMpe_0 - negMr_0);
            }

            //  Calculate Rotation at Capping Point (Intersection of the two Slopes)
            negUy       = negMpe / K_j;
            negMpeProj  = negMpe - negKp * negUy;
            negUmax     = fabs((negMmaxProj - negMpeProj) / (negKpc + negKp));
            negMmax     = negMpeProj + negUmax * negKp;

            if ((negMmax - negMr_0) / (negUmax + fabs(Ui) - negMr_0 / K_j) < negKp) {
                negKp       = (negMmax - negMr_0) / (negUmax + fabs(Ui) - negMr_0 / K_j);
                negMpeProj  = negMpe - negKp * negUy;
                negUmax     = fabs((negMmaxProj - negMpeProj) / (negKpc + negKp));
                negMmax     = negMpeProj + negUmax * negKp;
            }
        }

    }
    else {

        betaS   = betaS_1;
        betaC   = betaC_1;

        if (Di >= 0.0) {
            posMpe      = posMpe_1;
            posMpeProj  = posMpeProj_1;
            posMmaxProj = posMmaxProj_1;
            posMmax     = posMmax_1;
            posKp       = posKp_1;
            posKpc      = posKpc_1;
            posUy       = posUy_1;
            posUmax     = posUmax_1;
        }
        else {
            negMpe      = negMpe_1;
            negMpeProj  = negMpeProj_1;
            negMmaxProj = negMmaxProj_1;
            negMmax     = negMmax_1;
            negKp       = negKp_1;
            negKpc      = negKpc_1;
            negUy       = negUy_1;
            negUmax     = negUmax_1;
        }
    }

        ///////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////// COMPUTE FORCE INCREMENT /////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////

    // If the residual moment is reached in a given direction, Override the values of Mmax, Umax, Kp and Kpc
    if (Di > 0.0) {
        if (posMmax < posMr_0) {
            posMmax = posMr_0;
            posKpc  = pow(10., -6);
            posKp   = pow(10., -6);
            posUmax = pow(10., -6);
        }
    }
    else {
        if (negMmax < negMr_0) {
            negMmax = negMr_0;
            negKpc  = pow(10., -6);
            negKp   = pow(10., -6);
            negUmax = pow(10., -6);
        }
    }

    //  Simple and unified notation for current bacbone parameters
    if (Fi_1 + K_j * dU > 0.0) {
        Ki          = K_j;
        Kpi         = posKp;
        Kpci        = posKpc;
        Umaxi       = posUmax;
        MpeProji    = posMpeProj;
        MmaxProji   = posMmaxProj;
    }
    else {
        Ki          = K_j;
        Kpi         = negKp;
        Kpci        = negKpc;
        Umaxi       = negUmax;
        MpeProji    = negMpeProj;
        MmaxProji   = negMmaxProj;
    }


    //  Moment Calculation Based on unloading/reloading stiffeness
    Fi              = Fi_1 + Ki * dU;

    //  Location Flags
    if ((Ui >= 0.0) && (Fi >= 0.0)) {
        QuarterFlag = 1;
    }
    else if ((Ui >= 0.0) && (Fi < 0.0)) {
        QuarterFlag = 2;
    }
    else if ((Ui <= 0.0) && (Fi < 0.0)) {
        QuarterFlag = 3;
    }
    else if ((Ui <= 0.0) && (Fi > 0.0)) {
        QuarterFlag = 4;
    }

    // Get Boundary Moment at Current Step Based on Current BackBone Curve
    if (QuarterFlag == 1) {
        if (fabs(Ui) <= Umaxi) {
            Mi_boundary = MpeProji + Kpi * Ui;
        }
        else if (fabs(Ui) > Umaxi) {
            Mi_boundary = max(posMr_0, MmaxProji - Kpci * Ui);
        }
        if (Mi_boundary <= posMr_0) {
            Mrpos_Flag  = 1;
        }
    }
    else if (QuarterFlag == 3) {
        if (fabs(Ui) <= Umaxi) {
            Mi_boundary = -MpeProji + Kpi * Ui;
            //cout << "        MpeProji=" << MpeProji << " Kpi=" << Kpi << " Mbound=" << Mi_boundary << endln;

        }
        else if (fabs(Ui) > Umaxi) {
            Mi_boundary = min(-negMr_0, -MmaxProji - Kpci * Ui);
        }
        if (Mi_boundary >= -negMr_0) {
            Mrneg_Flag  = 1;
        }
    }
    else if (QuarterFlag == 2) {
        Mi_boundary     = min(-negMr_0, -MpeProji + Kpi * fabs(Ui));
        if (Mi_boundary == -negMr_0 && TangentK==1.e-6) {
            Mrneg_Flag  = 1;
        }
    }
    else if (QuarterFlag == 4) {
        Mi_boundary     = max(posMr_0, MpeProji - Kpi * fabs(Ui));
        if (Mi_boundary == posMr_0 && TangentK == 1.e-6) {
            Mrneg_Flag  = 1;
        }
    }

    //cout << "                Fi_1=" << Fi_1 << " Fi=" << Fi << " TangentK=" << TangentK << " Mbound=" << Mi_boundary << " Q=" << QuarterFlag << endln;

    // If Failure took place in a given direction (Fail_Flag_dir=1), Set the Boundary Moment in the opposite direction to Mr
    if ((Ui <= 0.0) && (Di > 0.0) && (Fail_FlagNeg == 1)) {
        Mi_boundary     = posMr_0;
    }
    else if ((Ui >= 0.0) && (Di < 0.0) && (Fail_FlagPos == 1)) {
        Mi_boundary     = -negMr_0;
    }


    // %%%%%%% Current Step Moment Calculation %%%%%%%
    // If current moment based on unloading/reloading Ki is larger than the boundary moment, set it equal to the boundary moment
    if (QuarterFlag == 1 && Di >= 0.0 && Fi >= Mi_boundary) {
        Fi  = Mi_boundary;
    }
    else if (QuarterFlag == 3 && Di <= 0.0 && Fi <= Mi_boundary) {
        Fi  = Mi_boundary;
    }
    else if (QuarterFlag == 2 && Fi <= Mi_boundary) {
        Fi  = Mi_boundary;
    }
    else if (QuarterFlag == 4 && Fi >= Mi_boundary) {
        Fi  = Mi_boundary;
    }

    if ((Mrneg_Flag == 1) || (Mrpos_Flag == 1)) {
        if (QuarterFlag == 1 && Di > 0 && Fi_1 == posMr_0)  {
            Fi  = posMr_0;
        }
        if  (QuarterFlag == 3 && Di < 0 && Fi_1 == -negMr_0) {
            Fi  = -negMr_0;
        }
    }

    // if fail flag is reached in any loading direction, set current moment equal to zero
    if ((Fail_FlagPos == 1.0) || (Fail_FlagNeg == 1.0) || (Energy_Flag == 1)) {
        Fi  = 0.0;
    }

    if (Yield_Flag != 1) {
        if (Ui >= posUy_0) {
            Fi  = posMpe_0 + posKp_0 * (Ui - posUy_0);
        }
        else {
            Fi  = Ke * (Ui);
        }

        if (Ui <= -negUy_0) {
            Fi  = -negMpe_0 - negKp_0 * fabs(Ui - negUy_0);
        }
        else {
            Fi  = Ke * (Ui);
        }
    }

    //cout << "                Fi_1=" << Fi_1 << " Fi=" << Fi << " TangentK=" << TangentK << " Mbound=" << Mi_boundary << " Q=" << QuarterFlag << endln;

    // %%%%%%%%%%%%% Energy Calculation %%%%%%%%%%%%%
    dEi         = (Fi + Fi_1) * 0.5 * dU;
    engTotl   += dEi; //  total energy dissipated till current incremental step

    // Energy calculation at each new excursion
    if (Fi / Fi_1 <= 0.0) {
        engExcr       = engTotl - engExcr_1;  // total energy dissipated in current excursion
        engExcr_1     = engTotl;                // total energy dissipated in previous excursion
        Excursion_Flag  = 1;
    }
    else {
        Excursion_Flag  = 0.0;
    }

    // Check if the Component inherit Reference Energy is Consumed
    if (Excursion_Flag == 1) {
        if ((engTotl >= refEnergyS) || (engTotl >= refEnergyC)) {
            Energy_Flag = 1;
        }
        if ((betaS > 1) || (betaC > 1)) {
            Energy_Flag = 1;
        }
    }
    else if (Reversal_Flag == 1) {
        if (engTotl >= refEnergyK) {
            Energy_Flag = 1;
        }
        if (betaK > 1) {
            Energy_Flag = 1;
        }
    }

    // %%%%%%%%%% PREPARE RETURN VALUES %%%%%%%%%%%%%

    if (Ui >= Ui_1) {
        posUy_1         = posUy;
        posUmax_1       = posUmax;
        posKp_1         = posKp;
        posKpc_1        = posKpc;
        posMpe_1        = posMpe;
        posMpeProj_1    = posMpeProj;
        posMmax_1       = posMmax;
        posMmaxProj_1   = posMmaxProj;
    }
    else {
        negUy_1         = negUy;
        negUmax_1       = negUmax;
        negKp_1         = negKp;
        negKpc_1        = negKpc;
        negMpe_1        = negMpe;
        negMpeProj_1    = negMpeProj;
        negMmax_1       = negMmax;
        negMmaxProj_1   = negMmaxProj;
    }

    K_j_1	    = K_j;
    betaS_1	= betaS;
    betaC_1	= betaC;
    betaK_1	= betaK;

    // Tangent Stiffeness Calculation
    if (Fi == posMr_0 || Fi == -negMr_0) {
        TangentK	= pow(10., -6);
    }

    if (Ui == Ui_1) {
        TangentK	= Ke;
        Fi	= Fi_1;
    }
    else {
        TangentK	= (Fi - Fi_1) / dU;
        if (TangentK == 0) {
            TangentK	= pow(10., -6);
        }
    }

    //cout << "                Fi_1=" << Fi_1 << " Fi=" << Fi << " Ke=" << Ke << " TangentK=" << TangentK << " Mbound=" << Mi_boundary << endln;

    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // %%%%%%%%%%%%%%%%%%%%%% END OF MAIN CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    return 0;
}

double	IMKBilin::getStress(void)
{
    //cout << " getStress" << endln;
    return (Fi);
}

double	IMKBilin::getTangent(void)
{
    //cout << " getTangent" << endln;
    return (TangentK);
}

double	IMKBilin::getInitialTangent(void)
{
    //cout << " getInitialTangent" << endln;
    return (Ke);
}

double	IMKBilin::getStrain(void)
{
    //cout << " getStrain" << endln;
    return (U);
}

int IMKBilin::commitState(void)
{
    //cout << " commitState" << endln;

    //commit trial  variables

    cU	            = U;

    cUi	            = Ui;
    cFi	            = Fi;
    cDi	            = Di;
    cUi_1	        = Ui_1;
    cFi_1	        = Fi_1;
    cDi_1	        = Di_1;

    cUlocal	        = Ulocal;
    cFlocal	        = Flocal;
    cTangentK	    = TangentK;

    cK_j_1	        = K_j_1;
    cPosUy_1	    = posUy_1;
    cPosUmax_1	    = posUmax_1;
    cPosKp_1	    = posKp_1;
    cPosKpc_1	    = posKpc_1;
    cPosMpe_1	    = posMpe_1;
    cPosMpeProj_1	= posMpeProj_1;
    cPosMmax_1	    = posMmax_1;
    cPosMmaxProj_1	= posMmaxProj_1;

    cNegUy_1	    = negUy_1;
    cNegUmax_1	    = negUmax_1;
    cNegKp_1	= negKp_1;
    cNegKpc_1	= negKpc_1;
    cNegMpe_1	    = negMpe_1;
    cNegMpeProj_1	= negMpeProj_1;
    cNegMmax_1	    = negMmax_1;
    cNegMmaxProj_1	= negMmaxProj_1;

    cBetaS_1	    = betaS_1;
    cBetaC_1	    = betaC_1;
    cBetaK_1	    = betaK_1;

    cExcursion_Flag	= Excursion_Flag;
    cReversal_Flag	= Reversal_Flag;
    cYield_Flag	    = Yield_Flag;
    cFail_FlagPos	= Fail_FlagPos;
    cFail_FlagNeg	= Fail_FlagNeg;
    cMrpos_Flag	    = Mrpos_Flag;
    cMrneg_Flag	    = Mrneg_Flag;
    cEnergy_Flag	= Energy_Flag;

    cEngExcr_1	= engExcr_1;
    cEngExcr	    = engExcr;
    cEngRvrs	    = engRvrs;
    cEngTotl	    = engTotl;

    return 0;
}

int IMKBilin::revertToLastCommit(void)
{
    //cout << " revertToLastCommit" << endln;

    //the opposite of commit trial history variables
    U	            = cU;

    Ui	            = cUi;
    Fi	            = cFi;
    Di	            = cDi;
    Ui_1	        = cUi_1;
    Fi_1	        = cFi_1;
    Di_1	        = cDi_1;

    Ulocal	    = cUlocal;
    Flocal	    = cFlocal;
    TangentK	    = cTangentK;

    K_j_1	        = cK_j_1;
    posUy_1	        = cPosUy_1;
    posUmax_1	    = cPosUmax_1;
    posKp_1	= cPosKp_1;
    posKpc_1	= cPosKpc_1;
    posMpe_1	    = cPosMpe_1;
    posMpeProj_1	= cPosMpeProj_1;
    posMmax_1	    = cPosMmax_1;
    posMmaxProj_1	= cPosMmaxProj_1;

    negUy_1	        = cNegUy_1;
    negUmax_1	    = cNegUmax_1;
    negKp_1	= cNegKp_1;
    negKpc_1	= cNegKpc_1;
    negMpe_1	    = cNegMpe_1;
    negMpeProj_1	= cNegMpeProj_1;
    negMmax_1	    = cNegMmax_1;
    negMmaxProj_1	= cNegMmaxProj_1;

    betaS_1	    = cBetaS_1;
    betaC_1	    = cBetaC_1;
    betaK_1	    = cBetaK_1;

    Excursion_Flag	= cExcursion_Flag;
    Reversal_Flag	= cReversal_Flag;
    Yield_Flag	    = cYield_Flag;
    Fail_FlagPos	= cFail_FlagPos;
    Fail_FlagNeg	= cFail_FlagNeg;
    Mrpos_Flag	    = cMrpos_Flag;
    Mrneg_Flag	    = cMrneg_Flag;
    Energy_Flag	    = cEnergy_Flag;

    engExcr_1	    = cEngExcr_1;
    engExcr	    = cEngExcr;
    engRvrs	    = cEngRvrs;
    engTotl	    = cEngTotl;

    return 0;
}

int IMKBilin::revertToStart(void)
{
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\\
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ONE TIME CALCULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\\
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

    if (posResM_0 == 0.0) {
        posResM_0	= 0.01;
    }
    if (negResM_0 == 0.0) {
        negResM_0	= 0.01;
    }

    posUy_0	                            = posMpe_0 / Ke;
    posUmax_0	                        = posUy_0 + posUp_0;
    posKp_0	                    = (posMmax_0 - posMpe_0) / (posUp_0);
    posKpc_0	                    = posMmax_0 / (posUpc_0);
    // posMpe_0	                        = posMpe_0;
    posMmax_0	                        = posMmaxMpe_0 * posMpe_0;
    posMpeProj_0	                    = posMmax_0 - posKp_0  * posUmax_0;
    posMmaxProj_0	                    = posMmax_0 + posKpc_0 * posUmax_0;

    negUy_0	                            = negMpe_0 / Ke;
    negUmax_0	                        = negUy_0 + negUp_0;
    negKp_0	                    = (negMmax_0 - negMpe_0) / (negUp_0);
    negKpc_0	                    = negMmax_0 / (negUpc_0);
    // negMpe_0	                        = negMpe_0;
    negMmax_0	                        = negMmaxMpe_0 * negMpe_0;
    negMpeProj_0	                    = negMmax_0 - negKp_0  * negUmax_0;
    negMmaxProj_0	                    = negMmax_0 + negKpc_0 * negUmax_0;

    posMr_0	                            = posResM_0*posMpe_0;
    negMr_0	                            = negResM_0*negMpe_0;

    refEnergyS	                    = LAMBDA_S*posMpe_0;
    refEnergyC	                    = LAMBDA_C*posMpe_0;
    refEnergyK	                    = LAMBDA_K*posMpe_0;

    K_j_1	        = cK_j_1	        = Ke;
    posUy_1	        = cPosUy_1	        = posMpe_0 / Ke;
    posUmax_1	    = cPosUmax_1	    = posUy_0 + posUp_0;
    posKp_1	= cPosKp_1	    = (posMmax_0 - posMpe_0) / (posUp_0);
    posKpc_1	= cPosKpc_1	= posMmax_0 / (posUpc_0);
    posMpe_1	    = cPosMpe_1	        = posMpe_0;
    posMmax_1	    = cPosMmax_1	    = posMmaxMpe_0*posMpe_0;
    posMpeProj_1	= cPosMpeProj_1	    = posMmax_1 - posKp_1 *posUmax_1;
    posMmaxProj_1	= cPosMmaxProj_1	= posMmax_1 + posKpc_1*posUmax_1;

    negUy_1	        = cNegUy_1	        = negMpe_0 / Ke;
    negUmax_1	    = cNegUmax_1	    = negUy_0 + negUp_0;
    negKp_1	= cNegKp_1	    = (negMmax_0 - negMpe_0) / (negUp_0);
    negKpc_1	= cNegKpc_1	= negMmax_0 / (negUpc_0);
    negMpe_1	    = cNegMpe_1	        = negMpe_0;
    negMmax_1	    = cNegMmax_1	    = negMmaxMpe_0*negMpe_0;
    negMpeProj_1	= cNegMpeProj_1	    = negMmax_1 - negKp_1 *negUmax_1;
    negMmaxProj_1	= cNegMmaxProj_1	= negMmax_1 + negKpc_1*negUmax_1;

    //initially I zero everything   
    U	            = cU	            = 0;

    Ui	            = cUi	            = 0;
    Fi	            = cFi	            = 0;
    Di	            = cDi	            = 0;
    Ui_1	        = cUi_1	            = 0;
    Fi_1	        = cFi_1	            = 0;
    Di_1	        = cDi_1	            = 0;
    TangentK	    = cTangentK	        = Ke;

    Ulocal	    = cUlocal	    = 0;
    Flocal	    = cFlocal	    = 0;

    betaS_1	        = cBetaS_1	        = 0;
    betaC_1	        = cBetaC_1	        = 0;
    betaK_1	        = cBetaK_1	        = 0;

    Excursion_Flag	= cExcursion_Flag	= 0;
    Reversal_Flag	= cReversal_Flag	= 0;
    Yield_Flag	    = cYield_Flag	    = 0;
    Fail_FlagPos	= cFail_FlagPos	    = 0;
    Fail_FlagNeg	= cFail_FlagNeg	    = 0;
    Mrpos_Flag	    = cMrpos_Flag	    = 0;
    Mrneg_Flag	    = cMrneg_Flag	    = 0;
    Energy_Flag	    = cEnergy_Flag	    = 0;

    engExcr_1	    = cEngExcr_1	    = 0;
    engExcr	    = cEngExcr	    = 0;
    engRvrs	    = cEngRvrs	    = 0;
    engTotl	    = cEngTotl	    = 0;

    return 0;
}

UniaxialMaterial *
IMKBilin::getCopy(void)
{
    IMKBilin *theCopy	= new IMKBilin(this->getTag(), Ke,
        posUp_0, posUpc_0, posUu_0, posMpe_0, posMmaxMpe_0, posResM_0,
        negUp_0, negUpc_0, negUu_0, negMpe_0, negMmaxMpe_0, negResM_0,
        LAMBDA_S, LAMBDA_C, LAMBDA_K, c_S, c_C, c_K, D_pos, D_neg);

    //cout << " getCopy" << endln;

    theCopy->posMr_0	        = posMr_0;
    theCopy->negMr_0	        = negMr_0;

    theCopy->U	                = U;
    theCopy->cU	                = cU;

    theCopy->Ui	                = Ui;
    theCopy->Fi	                = Fi;
    theCopy->Di	                = Di;
    theCopy->Ui_1	            = Ui_1;
    theCopy->Fi_1	            = Fi_1;
    theCopy->Di_1	            = Di_1;

    theCopy->Ulocal	        = Ulocal;
    theCopy->Flocal	        = Flocal;
    theCopy->TangentK	        = TangentK;

    theCopy->K_j_1	            = K_j_1;
    theCopy->posUy_1	        = posUy_1;
    theCopy->posUmax_1	        = posUmax_1;
    theCopy->posKp_1	    = posKp_1;
    theCopy->posKpc_1	    = posKpc_1;
    theCopy->posMpe_1	        = posMpe_1;
    theCopy->posMpeProj_1	    = posMpeProj_1;
    theCopy->posMmax_1	        = posMmax_1;
    theCopy->posMmaxProj_1	    = posMmaxProj_1;

    theCopy->negUy_1	        = negUy_1;
    theCopy->negUmax_1	        = negUmax_1;
    theCopy->negKp_1	    = negKp_1;
    theCopy->negKpc_1	    = negKpc_1;
    theCopy->negMpe_1	        = negMpe_1;
    theCopy->negMpeProj_1	    = negMpeProj_1;
    theCopy->negMmax_1	        = negMmax_1;
    theCopy->negMmaxProj_1	    = negMmaxProj_1;

    theCopy->betaS_1	        = betaS_1;
    theCopy->betaC_1	        = betaC_1;
    theCopy->betaK_1	        = betaK_1;

    theCopy->Excursion_Flag	    = Excursion_Flag;
    theCopy->Reversal_Flag	    = Reversal_Flag;
    theCopy->Yield_Flag	        = Yield_Flag;
    theCopy->Fail_FlagPos	    = Fail_FlagPos;
    theCopy->Fail_FlagNeg	    = Fail_FlagNeg;
    theCopy->Mrpos_Flag	        = Mrpos_Flag;
    theCopy->Mrneg_Flag	        = Mrneg_Flag;
    theCopy->Energy_Flag	    = Energy_Flag;

    theCopy->engExcr_1	    = engExcr_1;
    theCopy->engExcr	        = engExcr;
    theCopy->engRvrs	        = engRvrs;
    theCopy->engTotl	        = engTotl;


    theCopy->cPosMr_0	        = cPosMr_0;
    theCopy->cNegMr_0	        = cNegMr_0;

    theCopy->cUi	            = cUi;
    theCopy->cFi	            = cFi;
    theCopy->cDi	            = cDi;
    theCopy->cUi_1	            = cUi_1;
    theCopy->cFi_1	            = cFi_1;
    theCopy->cDi_1	            = cDi_1;

    theCopy->cUlocal	        = cUlocal;
    theCopy->cFlocal	        = cFlocal;
    theCopy->cTangentK	        = cTangentK;

    theCopy->cK_j_1	            = cK_j_1;
    theCopy->cPosUy_1	        = cPosUy_1;
    theCopy->cPosUmax_1	        = cPosUmax_1;
    theCopy->cPosKp_1	    = cPosKp_1;
    theCopy->cPosKpc_1	    = cPosKpc_1;
    theCopy->cPosMpe_1	        = cPosMpe_1;
    theCopy->cPosMpeProj_1	    = cPosMpeProj_1;
    theCopy->cPosMmax_1	        = cPosMmax_1;
    theCopy->cPosMmaxProj_1	    = cPosMmaxProj_1;

    theCopy->cNegUy_1	        = cNegUy_1;
    theCopy->cNegUmax_1	        = cNegUmax_1;
    theCopy->cNegKp_1	    = cNegKp_1;
    theCopy->cNegKpc_1	    = cNegKpc_1;
    theCopy->cNegMpe_1	        = cNegMpe_1;
    theCopy->cNegMpeProj_1	    = cNegMpeProj_1;
    theCopy->cNegMmax_1	        = cNegMmax_1;
    theCopy->cNegMmaxProj_1	    = cNegMmaxProj_1;

    theCopy->cBetaS_1	        = cBetaS_1;
    theCopy->cBetaC_1	        = cBetaC_1;
    theCopy->cBetaK_1	        = cBetaK_1;

    theCopy->cExcursion_Flag	= cExcursion_Flag;
    theCopy->cReversal_Flag	    = cReversal_Flag;
    theCopy->cYield_Flag	    = cYield_Flag;
    theCopy->cFail_FlagPos	    = cFail_FlagPos;
    theCopy->cFail_FlagNeg	    = cFail_FlagNeg;
    theCopy->cMrpos_Flag	    = cMrpos_Flag;
    theCopy->cMrneg_Flag	    = cMrneg_Flag;
    theCopy->cEnergy_Flag	    = cEnergy_Flag;

    theCopy->cEngExcr_1	    = cEngExcr_1;
    theCopy->cEngExcr	        = cEngExcr;
    theCopy->cEngRvrs	        = cEngRvrs;
    theCopy->cEngTotl	        = cEngTotl;

    return theCopy;
}

int IMKBilin::sendSelf(int cTag, Channel &theChannel)
{
    int res	= 0;
    cout << " sendSelf" << endln;

    static Vector data(113);
    data(0)	    = this->getTag();
    data(1)	    = Ke;
    data(2)	    = posUp_0;
    data(3)	    = posUpc_0;
    data(4)	    = posUu_0;
    data(5)	    = posMpe_0;
    data(6)	    = posMmaxMpe_0;
    data(7)	    = posResM_0;
    data(8)	    = negUp_0;
    data(9)	    = negUpc_0;
    data(10)	= negUu_0;
    data(11)	= negMpe_0;
    data(12)	= negMmaxMpe_0;
    data(13)	= negResM_0;
    data(14)	= LAMBDA_S;
    data(15)	= LAMBDA_C;
    data(16)	= LAMBDA_K;
    data(17)	= c_S;
    data(18)	= c_C;
    data(19)	= c_K;
    data(20)	= D_pos;
    data(21)	= D_neg;

    data(22)	= Ui;
    data(23)	= Fi;
    data(24)	= Di;
    data(25)	= Ui_1;
    data(26)	= Fi_1;
    data(27)	= Di_1;
    data(28)	= Ulocal;
    data(29)	= Flocal;
    data(30)	= TangentK;

    data(31)	= K_j_1;
    data(32)	= posUy_1;
    data(33)	= posUmax_1;
    data(34)	= posKp_1;
    data(35)	= posKpc_1;
    data(36)	= posMpe_1;
    data(37)	= posMpeProj_1;
    data(38)	= posMmax_1;
    data(39)	= posMmaxProj_1;

    data(40)	= negUy_1;
    data(41)	= negUmax_1;
    data(42)	= negKp_1;
    data(43)	= negKpc_1;
    data(44)	= negMpe_1;
    data(45)	= negMpeProj_1;
    data(46)	= negMmax_1;
    data(47)	= negMmaxProj_1;

    data(48)	= betaS_1;
    data(49)	= betaC_1;
    data(50)	= betaK_1;

    data(51)	= refEnergyS;
    data(52)	= refEnergyC;
    data(53)	= refEnergyK;

    data(54)	= Excursion_Flag;
    data(55)	= Reversal_Flag;
    data(56)	= Yield_Flag;
    data(57)	= Fail_FlagPos;
    data(58)	= Fail_FlagNeg;
    data(59)	= Mrpos_Flag;
    data(60)	= Mrneg_Flag;
    data(61)	= Energy_Flag;

    data(62)	= engExcr_1;
    data(63)	= engExcr;
    data(64)	= engRvrs;
    data(65)	= engTotl;

    data(66)	= cUi;
    data(67)	= cFi;
    data(68)	= cDi;
    data(69)	= cUi_1;
    data(70)	= cFi_1;
    data(71)	= cDi_1;
    data(72)	= cUlocal;
    data(73)	= cFlocal;
    data(74)	= cTangentK;

    data(75)	= cK_j_1;
    data(76)	= cPosUy_1;
    data(77)	= cPosUmax_1;
    data(78)	= cPosKp_1;
    data(79)	= cPosKpc_1;
    data(80)	= cPosMpe_1;
    data(81)	= cPosMpeProj_1;
    data(82)	= cPosMmax_1;
    data(83)	= cPosMmaxProj_1;

    data(84)	= cNegUy_1;
    data(85)	= cNegUmax_1;
    data(86)	= cNegKp_1;
    data(87)	= cNegKpc_1;
    data(88)	= cNegMpe_1;
    data(89)	= cNegMpeProj_1;
    data(90)	= cNegMmax_1;
    data(91)	= cNegMmaxProj_1;

    data(92)	= cBetaS_1;
    data(93)	= cBetaC_1;
    data(94)	= cBetaK_1;

    data(95)	= cExcursion_Flag;
    data(96)	= cReversal_Flag;
    data(97)	= cYield_Flag;
    data(98)	= cFail_FlagPos;
    data(99)	= cFail_FlagNeg;
    data(100)	= cMrpos_Flag;
    data(101)	= cMrneg_Flag;
    data(102)	= cEnergy_Flag;

    data(103)	= cEngExcr_1;
    data(104)	= cEngExcr;
    data(105)	= cEngRvrs;
    data(106)	= cEngTotl;

    data(107)	= posMr_0;
    data(108)	= negMr_0;
    data(109)	= cPosMr_0;
    data(110)	= cNegMr_0;

    data(111)	= U;
    data(112)	= cU;

    res	= theChannel.sendVector(this->getDbTag(), cTag, data);
    if (res < 0)
        opserr << "IMKBilin::sendSelf() - failed to send data\n";

    return res;
}

int IMKBilin::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    int res	= 0;
    static Vector data(113);
    res	= theChannel.recvVector(this->getDbTag(), cTag, data);

    if (res < 0) {
        opserr << "IMKBilin::recvSelf() - failed to receive data\n";
        this->setTag(0);
    }
    else {
        cout << " recvSelf" << endln;
        this->setTag((int)data(0));
        Ke	            = data(1);
        posUp_0	        = data(2);
        posUpc_0	    = data(3);
        posUu_0	        = data(4);
        posMpe_0	    = data(5);
        posMmaxMpe_0	= data(6);
        posResM_0	    = data(7);
        negUp_0	        = data(8);
        negUpc_0	    = data(9);
        negUu_0	        = data(10);
        negMpe_0	    = data(11);
        negMmaxMpe_0	= data(12);
        negResM_0	    = data(13);
        LAMBDA_S	    = data(14);
        LAMBDA_C	    = data(15);
        LAMBDA_K	    = data(16);
        c_S	            = data(17);
        c_C	            = data(18);
        c_K	            = data(19);
        D_pos	        = data(20);
        D_neg	        = data(21);

        Ui	            = data(22);
        Fi	            = data(23);
        Di	            = data(24);
        Ui_1	        = data(25);
        Fi_1	        = data(26);
        Di_1	        = data(27);
        Ulocal	    = data(28);
        Flocal	    = data(29);
        TangentK	    = data(30);

        K_j_1	        = data(31);
        posUy_1	        = data(32);
        posUmax_1	    = data(33);
        posKp_1	= data(34);
        posKpc_1	= data(35);
        posMpe_1	    = data(36);
        posMpeProj_1	= data(37);
        posMmax_1	    = data(38);
        posMmaxProj_1	= data(39);

        negUy_1	        = data(40);
        negUmax_1	    = data(41);
        negKp_1	= data(42);
        negKpc_1	= data(43);
        negMpe_1	    = data(44);
        negMpeProj_1	= data(45);
        negMmax_1	    = data(46);
        negMmaxProj_1	= data(47);

        betaS_1	        = data(48);
        betaC_1	        = data(49);
        betaK_1	        = data(50);

        refEnergyS	= data(51);
        refEnergyC	= data(52);
        refEnergyK	= data(53);

        Excursion_Flag	= data(54);
        Reversal_Flag	= data(55);
        Yield_Flag	    = data(56);
        Fail_FlagPos	= data(57);
        Fail_FlagNeg	= data(58);
        Mrpos_Flag	    = data(59);
        Mrneg_Flag	    = data(60);
        Energy_Flag	    = data(61);

        engExcr_1	    = data(62);
        engExcr	    = data(63);
        engRvrs	    = data(64);
        engTotl	    = data(65);

        cUi	            = data(66);
        cFi	            = data(67);
        cDi	            = data(68);
        cUi_1	        = data(69);
        cFi_1	        = data(70);
        cDi_1	        = data(71);

        cUlocal	    = data(72);
        cFlocal	    = data(73);
        cTangentK	    = data(74);

        cK_j_1	        = data(75);
        cPosUy_1	    = data(76);
        cPosUmax_1	    = data(77);
        cPosKp_1	= data(78);
        cPosKpc_1	= data(79);
        cPosMpe_1	    = data(80);
        cPosMpeProj_1	= data(81);
        cPosMmax_1	    = data(82);
        cPosMmaxProj_1	= data(83);

        cNegUy_1	    = data(84);
        cNegUmax_1	    = data(85);
        cNegKp_1	= data(86);
        cNegKpc_1	= data(87);
        cNegMpe_1	    = data(88);
        cNegMpeProj_1	= data(89);
        cNegMmax_1	    = data(90);
        cNegMmaxProj_1	= data(91);

        cBetaS_1	    = data(92);
        cBetaC_1	    = data(93);
        cBetaK_1	    = data(94);

        cExcursion_Flag	= data(95);
        cReversal_Flag	= data(96);
        cYield_Flag	    = data(97);
        cFail_FlagPos	= data(98);
        cFail_FlagNeg	= data(99);
        cMrpos_Flag	    = data(100);
        cMrneg_Flag	    = data(101);
        cEnergy_Flag	= data(102);

        cEngExcr_1	= data(103);
        cEngExcr	    = data(104);
        cEngRvrs	    = data(105);
        cEngTotl	    = data(106);

        posMr_0	        = data(107);
        negMr_0	        = data(108);
        cPosMr_0	    = data(109);
        cNegMr_0	    = data(110);

        U	            = data(111);
        cU	            = data(112);

    }

    return res;
}

void IMKBilin::Print(OPS_Stream &s, int flag)
{
    cout << "IMKBilin tag: " << this->getTag() << endln;
}
