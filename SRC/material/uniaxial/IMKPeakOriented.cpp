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
#include <IMKPeakOriented.h>
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

static int numIMKPeakOrientedMaterials	= 0;

void *
OPS_IMKPeakOriented()
{
    if (numIMKPeakOrientedMaterials == 0) {
        numIMKPeakOrientedMaterials++;
        OPS_Error("IMK with Peak-Oriented Response - Code by Elkady & Eljisr (July22)\n", 1);
    }

    // Pointer to a uniaxial material that will be returned
    UniaxialMaterial *theMaterial	= 0;

    int    iData[1];
    double dData[23];
    int numData	= 1;

    if (OPS_GetIntInput(&numData, iData) != 0) {
        opserr << "WARNING invalid uniaxialMaterial IMKPeakOriented tag" << endln;
        return 0;
    }

    numData	= 23;


    if (OPS_GetDoubleInput(&numData, dData) != 0) {
        opserr << "Invalid Args want: uniaxialMaterial IMKPeakOriented tag? Ke? ";
        opserr << "Up_pos? Upc_pos? Uu_pos? Fy_pos? FcapFy_pos? ResF_pos? ";
        opserr << "Up_neg? Upc_neg? Uu_neg? Fy_neg? FcapFy_neg? ResF_neg? ";
        opserr << "LamdaS? LamdaC? LamdaA? LamdaK? Cs? Cc? Ca? Ck? D_pos? D_neg? ";
        return 0;
    }



    // Parsing was successful, allocate the material
    theMaterial	= new IMKPeakOriented(iData[0],
        dData[0],
        dData[1], dData[2], dData[3], dData[4], dData[5], dData[6],
        dData[7], dData[8], dData[9], dData[10], dData[11], dData[12],
        dData[13], dData[14], dData[15], dData[16], dData[17], dData[18], dData[19], dData[20],
        dData[21], dData[22]);

    if (theMaterial == 0) {
        opserr << "WARNING could not create uniaxialMaterial of type IMKPeakOriented Material\n";
        return 0;
    }

    return theMaterial;
}

IMKPeakOriented::IMKPeakOriented(int tag, double p_Ke,
    double p_Up_pos, double p_Upc_pos, double p_Uu_pos, double p_Fy_pos, double p_FcapFy_pos, double p_ResF_pos,
    double p_Up_neg, double p_Upc_neg, double p_Uu_neg, double p_Fy_neg, double p_FcapFy_neg, double p_ResF_neg,
    double p_LAMBDA_S, double p_LAMBDA_C, double p_LAMBDA_A, double p_LAMBDA_K, double p_c_S, double p_c_C, double p_c_A, double p_c_K, double p_D_pos, double p_D_neg)
    : UniaxialMaterial(tag, 0), Ke(p_Ke),
    Up_pos(p_Up_pos), Upc_pos(p_Upc_pos), Uu_pos(p_Uu_pos), Fy_pos(p_Fy_pos), FcapFy_pos(p_FcapFy_pos), ResF_pos(p_ResF_pos),
    Up_neg(p_Up_neg), Upc_neg(p_Upc_neg), Uu_neg(p_Uu_neg), Fy_neg(p_Fy_neg), FcapFy_neg(p_FcapFy_neg), ResF_neg(p_ResF_neg),
    LAMBDA_S(p_LAMBDA_S), LAMBDA_C(p_LAMBDA_C), LAMBDA_A(p_LAMBDA_A), LAMBDA_K(p_LAMBDA_K), c_S(p_c_S), c_C(p_c_C), c_A(p_c_A), c_K(p_c_K), D_pos(p_D_pos), D_neg(p_D_neg)
{
    this->revertToStart();
}

IMKPeakOriented::IMKPeakOriented()
    :UniaxialMaterial(0, 0), Ke(0),
    Up_pos(0), Upc_pos(0), Uu_pos(0), Fy_pos(0), FcapFy_pos(0), ResF_pos(0),
    Up_neg(0), Upc_neg(0), Uu_neg(0), Fy_neg(0), FcapFy_neg(0), ResF_neg(0),
    LAMBDA_S(0), LAMBDA_C(0), LAMBDA_A(0), LAMBDA_K(0), c_S(0), c_C(0), c_A(0), c_K(0), D_pos(0), D_neg(0)
{
    this->revertToStart();
}

IMKPeakOriented::~IMKPeakOriented()
{
    // does nothing
}

int IMKPeakOriented::setTrialStrain(double strain, double strainRate)
{
    //all variables to the last commit
    this->revertToLastCommit();

    //state determination algorithm: defines the current force and tangent stiffness
    U		= strain; //set trial displacement
    Ui_1	= Ui;
    Fi_1	= Fi;
    Ui		= U;

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////  MAIN CODE //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // Incremental deformation at current step
    dU	= Ui - Ui_1;


    if (Failure_Flag) {     // When a failure has already occured
        Fi 	= 0;
        dEi	= 0;
    } else if (dU == 0) {   // When deformation doesn't change from the last
        Fi 	= Fi_1;
        dEi	= 0;
    } else {
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
    //      3:  Towards local Peak      +
    //      4:  Towards global Peak     +
    //      5:  Towards Capping Point   +
    //      6:  Towards Residual Point  +
    //      7:  Residual Branch         +
    //      13: Towards local Peak      -
    //      14: Towards global Peak     -
    //      15: Towards Capping Point   -
    //      16: Towards Residual Point  -
    //      17: Residual Branch         -
    //  Flag
    //      Yield_Flag:     Preserved.      When the deformation exceeds yield capacity for the first time.
    //      Excursion_Flag: Not preserved.  When crossing X-axis. Evokes re-considering of the deteriorations and which peak to go for.
    //      Reversal_Flag:  Not preserved.  When unloading starts. Evokes re-condiersing of the stiffness deterioration and peak point registration.
        exBranch       	= Branch;
        Excursion_Flag 	= false;
        Reversal_Flag  	= false;
        if (Branch == 0) {
            // CHECK FOR YIELDING
            if (Ui > posUy) {
                Yield_Flag 	= true;
                Branch 	= 5;
            } else if (Ui < negUy) {
                Yield_Flag 	= true;
                Branch 	= 15;
            }
        } else if (Branch == 1) {
            if (Fi_1*(Fi_1+dU*K_unload) <= 0) {
            // CHECK FOR NEW EXCURSION
                Excursion_Flag 	= true;
            } else if (Ui > posUlocal) {
                Branch 	= 4;
            } else if (Ui < negUlocal) {
                Branch 	= 14;
            }
        } else if (Fi_1*dU < 0) {
            Reversal_Flag  	= true;
            Branch 	= 1;
        }
    // Branch shifting from 3 -> 4 -> 5 -> 6 -> 7 can be considered.
        if (Branch == 3 && Ui > posUlocal) {
            Branch 	= 4;
        }
        if (Branch == 4 && Ui > posUglobal) {
            Branch 	= 5;
        }
        if (Branch == 5 && Ui > posUcap) {
            Branch 	= 6;
        }
        if (Branch == 6 && Ui > posUres) {
            Branch 	= 7;
        }

        if (Branch == 13 && Ui < negUlocal) {
            Branch 	= 14;
        }
        if (Branch == 14 && Ui < negUglobal) {
            Branch 	= 15;
        }
        if (Branch == 15 && Ui < negUcap) {
            Branch 	= 16;
        }
        if (Branch == 16 && Ui < negUres) {
            Branch 	= 17;
        }
    // UPDATE PEAK POINTS
        if (Reversal_Flag) {
            if ( Fi_1 > 0 ){
                posUlocal	= Ui_1;             // UPDATE local
                posFlocal	= Fi_1;
                if ( Ui_1 > posUglobal ) {    // UPDATE GLOBAL
                    posUglobal   	= Ui_1;
                    posFglobal   	= Fi_1;
                }
            } else {
                negUlocal	= Ui_1;             // UPDATE local
                negFlocal	= Fi_1;
                if ( Ui_1 < negUglobal ) {    // UPDATE GLOBAL
                    negUglobal   	= Ui_1;
                    negFglobal   	= Fi_1;
                }
            }
        }

        ///////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////
        /////////////////// UPDATE DETERIORATION PARAMETERS ///////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////

        // UPDATE DETERIORATION PARAMETERS AT EACH NEW EXCURSION

        if (Reversal_Flag) {
            EpjK        = engAcml            - 0.5*(Fi_1 / K_unload)*Fi_1;
            EiK         = engAcml - engDspt  - 0.5*(Fi_1 / K_unload)*Fi_1;
            betaK       = pow( (EiK / (EtK - EpjK)), c_K );
            K_unload    = K_unload * (1 - betaK);
            TangentK    = K_unload;
        // Detect unloading completed in a step.
            if (Fi_1*(Fi_1+dU*K_unload) <= 0) {
                Excursion_Flag  = true;
                Reversal_Flag   = false;
            }
        }
        else {
            betaK   = 0;
        }
        if (Excursion_Flag) {
            //Epj   = engAcml + dEi;
            Ei      = fmax(0, engAcml - engDspt);
            betaS   = pow((Ei / (EtS - engAcml)), c_S);
            betaC   = pow((Ei / (EtC - engAcml)), c_C);
            betaA   = pow((Ei / (EtA - engAcml)), c_A);
            engDspt = engAcml;
        }
        else {
            //Epj   = engDspt;
            betaS   = 0;
            betaC   = 0;
            betaA   = 0;
        }
        ///////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////
        /////////////////// UPDATE BACKBONE CURVE /////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////
        // Update Positive Backbone and Target Peak Point
        if ( Excursion_Flag && Yield_Flag ) {
            // Positive loading backbone
            if (Fi_1 < 0) {
                // Basic strength deterioration: Yield point
                // Basic strength deterioration: Post-yield Stiffness
                posFy   *= (1 - betaS * D_pos);
                posKp   *= (1 - betaS * D_pos);
                if (posFy < posFres) {
                    posFy   = posFres;
                    posKp   = 0;
                }
                posUy       = posFy / Ke;
                // Basic strength deterioration: Capping Point
                sPCsp   = (posFy - posUy * posKp - posFcap + posKpc * posUcap) / (posKpc - posKp);
                posFcap = posFcap + (sPCsp - posUcap)*posKpc;
                posUcap = sPCsp;
                // Post-capping strength deterioration: Capping point
                sPCpcp  = max(posUcap + betaC * D_pos*(posFcap - posKpc * posUcap) / (posKpc - posKp), posUy);
                posFcap = posFcap + (sPCpcp - posUcap)*posKp;
                posUcap = sPCpcp;
                // Accelerated reloading stiffness deterioration: Target peak deformation point
                posUglobal  = (1 + betaA * D_pos)*posUglobal;
                if (posUglobal < posUy) {
                    posFglobal  = Ke * posUglobal;
                    // Target peak deformation in post-yield branch of the updated backbone
                }
                else if (posUglobal < posUcap) {
                    posFglobal  = posKp * (posUglobal - posUy) + posFy;
                    // Target peak deformation in post-capping branch of the updated backbone
                }
                else {
                    posFglobal  = max(posKpc*(posUglobal - posUcap) + posFcap, posFres);
                }
                posUres = (posFres - posFcap + posKpc * posUcap) / posKpc;
            }
            else {
                // Update Negative Backbone and Target Peak Point
                // Basic strength deterioration: Yield point
                // Basic strength deterioration: Post-yield stiffness
                negFy	*= (1 - betaS * D_neg);
                negKp	*= (1 - betaS * D_neg);
                if (negFy > negFres) {
                    negFy	= negFres;
                    negKp	= 0;
                }
                negUy		= negFy / Ke;
                // Basic strength deterioration: Capping point
                sPCsn		= (negFy - negUy * negKp - negFcap + negKpc * negUcap) / (negKpc - negKp);
                negFcap	= negFcap + (sPCsn - negUcap)*negKpc;
                negUcap	= sPCsn;
                // Post-capping strength deterioration: Capping point
                sPCpcn		= min(negUcap + betaC * D_neg*(negFcap - negKpc * negUcap) / (negKpc - negKp), negUy);
                negFcap	= negFcap + (sPCpcn - negUcap)*negKp;
                negUcap	= sPCpcn;
                // Accelerated reloading stiffness deterioration: Target peak deformation point
                negUglobal	= (1 + betaA * D_neg)*negUglobal;
                // Target peak deformation in reloading branch of the updated backbone
                if (negUglobal > negUy) {
                    negFglobal	= Ke * negUglobal;
                    // Target peak deformation in post-yield branch of the updated backbone
                }
                else if (negUglobal > negUcap) {
                    negFglobal	= negKp * (negUglobal - negUy) + negFy;
                    // Target peak deformation in post-capping branch of the updated backbone
                }
                else {
                    negFglobal	= min(negKpc*(negUglobal - negUcap) + negFcap, negFres);
                }
                negUres  	= (negFres - negFcap + negKpc * negUcap) / negKpc;
            }
        }

        ///////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////// COMPUTE FORCE INCREMENT /////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////

        if (Excursion_Flag) {
            // Detection of reloading completed in a step might be needed, while it's not as severe as a one step unloading.
            u0 	= Ui_1 - (Fi_1 / K_unload);
            if (dU > 0) {
                K_local    	= posFlocal   / (posUlocal  - u0);
                K_global   	= posFglobal  / (posUglobal - u0);
                if ( (posFlocal < posFglobal) && (K_local > K_global)) {
                    Branch     	= 3;
                    K_reload   	= K_local;
                }
                else {
                    Branch     	= 4;
                    K_reload   	= K_global;
                }
            }
            else {
                K_local    	= negFlocal   / (negUlocal  - u0);
                K_global   	= negFglobal  / (negUglobal - u0);
                if ( (negFlocal > negFglobal) && (K_local > K_global)) {
                    Branch     	= 13;
                    K_reload   	= K_local;
                }
                else {
                    Branch     	= 14;
                    K_reload   	= K_global;
                }
            }
            dF 	= 0             - Fi_1 + K_reload*  (Ui - u0);
            TangentK			= K_reload;
// With Branch Change
    // Positive Force
        }
        // CASE 4: WHEN RELOADING BUT BETWEEN LAST CYCLE PEAK POINT AND GLOBAL PEAK POINT
        // CASE 5: WHEN LOADING IN GENERAL TOWARDS THE TARGET PEAK
        // CASE 6: WHEN LOADING IN GENERAL TOWARDS THE LAST CYCLE PEAK POINT BUT BEYOND IT
        else if (Branch == 4 && exBranch != 4) {
            K_reload   	= (posFglobal - posFlocal) / (posUglobal - posUlocal);
            TangentK	= K_reload;
            dF 	= posFlocal   - Fi_1 + K_reload*  (Ui - posUlocal);
        }
        // CASE 7: WHEN LOADING BEYOND THE TARGET PEAK BUT BEFORE THE CAPPING POINT
        else if (Branch == 5 && exBranch == 0) {
            TangentK	= posKp;
            dF 	= posFy       - Fi_1 + posKp*   (Ui - posUy);
        }
        else if (Branch == 5 && exBranch != 5) {
            TangentK	= posKp;
            dF 	= posFglobal  - Fi_1 + posKp*   (Ui - posUglobal);
        }
        // CASE 8: WHEN LOADING AND BETWEEN THE CAPPING POINT AND THE RESIDUAL POINT
        else if (Branch == 6 && exBranch == 5) {
            TangentK	= posKpc;
            dF 	= posFcap     - Fi_1 + posKpc*  (Ui - posUcap);
        }
        else if (Branch == 6 && exBranch != 6) {
            TangentK	= posKpc;
            dF 	= posFglobal  - Fi_1 + posKpc*  (Ui - posUglobal);
        }
        // CASE 9: WHEN LOADING AND BEYOND THE RESIDUAL POINT
        else if (Branch == 7 && exBranch != 7) {
            TangentK	= 0;
            dF 	= posFres     - Fi_1;
    // Negative Force
        }
        else if (Branch == 14 && exBranch != 14) {
        // CASE 4: WHEN RELOADING BUT BETWEEN LAST CYCLE PEAK POINT AND GLOBAL PEAK POINT
        // CASE 5: WHEN LOADING IN GENERAL TOWARDS THE TARGET PEAK
        // CASE 6: WHEN LOADING IN GENERAL TOWARDS THE LAST CYCLE PEAK POINT BUT BEYOND IT
            K_reload   	= (negFglobal - negFlocal) / (negUglobal - negUlocal);
            TangentK	= K_reload;
            dF         	= negFlocal - Fi_1 + K_reload*(Ui - negUlocal);
        }
        // CASE 7: WHEN LOADING BEYOND THE TARGET PEAK BUT BEFORE THE CAPPING POINT
        else if (Branch == 15 && exBranch == 0) {
            TangentK	= negKp;
            dF 	= negFy - Fi_1 + negKp*(Ui - negUy);
        }
        else if (Branch == 15 && exBranch != 15) {
            TangentK	= negKp;
            dF 	= negFglobal - Fi_1 + negKp*(Ui - negUglobal);
        }
        // CASE 8: WHEN LOADING AND BETWEEN THE CAPPING POINT AND THE RESIDUAL POINT
        else if (Branch == 16 && exBranch == 15) {
            TangentK	= negKpc;
            dF 	= negFcap - Fi_1 + negKpc*(Ui - negUcap);
        }
        else if (Branch == 16 && exBranch != 16) {
            TangentK	= negKpc;
            dF 	= negFglobal - Fi_1 + negKpc*(Ui - negUglobal);
        }
        // CASE 9: WHEN LOADING AND BEYOND THE RESIDUAL POINT
        else if (Branch == 17 && exBranch != 17) {
            TangentK	= 0;
            dF 	= negFres - Fi_1;
// Without Branch Change
        } else {
            dF	= dU*TangentK;
        }
    // Branch Change check
        // if (Branch!=exBranch) {
        //  std::cout << exBranch << " -> " << Branch << "\n";
        // }

        // Force
        Fi	= Fi_1 + dF;

        ///////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////
        // CHECK FOR FAILURE
        ///////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////

        // Failure criteria (Tolerance	= 1//)
    // I have no idea about why it can' t be 0 nor 1.
        FailS	= ( betaS < -0.01 || betaS > 1.01 );
        FailC	= ( betaC < -0.01 || betaC > 1.01 );
        FailA	= ( betaA < -0.01 || betaA > 1.01 );
        FailK	= ( betaK < -0.01 || betaK > 1.01 );
        FailPp 	= ( posFglobal == 0               );
        FailPn 	= ( negFglobal == 0               );
        FailDp 	= ( Ui >  Uu_pos                  );
        FailDn 	= ( Ui < -Uu_neg                  );
        FailRp 	= ( Branch ==  7 && posFres == 0  );
        FailRn 	= ( Branch == 17 && negFres == 0  );
        Failure_Flag    = (FailS||FailC||FailA||FailK||FailPp||FailPn||FailRp||FailRn||FailDp||FailDn);
        if (Failure_Flag) {
            Fi 	= 0;
        }
        dEi	= 0.5*(Fi + Fi_1)*dU; // Internal energy increment
    }
    //// Energy
    engAcml	+= dEi; 	

    //// Update Variables
    du_i_1	= dU;

    // Tangent Stiffeness Calculation
    ki	= TangentK;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////// END OF MAIN CODE ///////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    return 0;
}

double IMKPeakOriented::getStress(void)
{
    //cout << " getStress" << endln;
    return (Fi);
}

double IMKPeakOriented::getTangent(void)
{
    //cout << " getTangent" << endln;
    return (TangentK);
}

double IMKPeakOriented::getInitialTangent(void)
{
    //cout << " getInitialTangent" << endln;
    return (Ke);
}

double IMKPeakOriented::getStrain(void)
{
    //cout << " getStrain" << endln;
    return (U);
}

int IMKPeakOriented::commitState(void)
{
    //cout << " commitState" << endln;

    //commit trial  variables

    cU		= U;

    cUi		= Ui;
    cFi		= Fi;
    cUi_1	= Ui_1;
    cFi_1	= Fi_1;

    cTangentK	= TangentK;

    cdu_i_1	= du_i_1;

    cPosUy	= posUy;
    cPosUcap	= posUcap;
    cPosFy	= posFy;
    cPosFcap	= posFcap;
    cPosUglobal	= posUglobal;
    cPosFglobal	= posFglobal;

    cPosUres	= posUres;
    cPosFres	= posFres;
    cPosKp	= posKp;
    cPosKpc	= posKpc;

    cNegUy	= negUy;
    cNegUcap	= negUcap;
    cNegFy	= negFy;
    cNegFcap	= negFcap;
    cNegUglobal	= negUglobal;
    cNegFglobal	= negFglobal;

    cNegUres	= negUres;
    cNegFres	= negFres;
    cNegKp	= negKp;
    cNegKpc	= negKpc;

    cK_unload	= K_unload;

    cEngAcml	= engAcml;
    cEngDspt	= engDspt;

    // cu0	= u0;

    cPosUlocal	= posUlocal;
    cPosFlocal	= posFlocal;
    cNegUlocal	= negUlocal;
    cNegFlocal	= negFlocal;

    cFailure_Flag		= Failure_Flag;
    // cExcursion_Flag	= Excursion_Flag;
    cExBranch			= exBranch;
    cBranch				= Branch;
    // cTargetPeak_Flag= TargetPeak_Flag;
    cYield_Flag			= Yield_Flag;
    // cReversal_Flag	= Reversal_Flag;

    cK_reload	= K_reload;

    return 0;
}
int IMKPeakOriented::revertToLastCommit(void)
{
    //cout << " revertToLastCommit" << endln;
    //the opposite of commit trial history variables
    U	= cU;
    Ui	= cUi;
    Fi	= cFi;
    Ui_1	= cUi_1;
    Fi_1	= cFi_1;
    TangentK	= cTangentK;
    du_i_1	= cdu_i_1;

    posUy	= cPosUy;
    posUcap	= cPosUcap;
    posFy	= cPosFy;
    posFcap	= cPosFcap;
    posUglobal	= cPosUglobal;
    posFglobal	= cPosFglobal;

    posUres	= cPosUres;
    posFres	= cPosFres;
    posKp	= cPosKp;
    posKpc	= cPosKpc;


    negUy	= cNegUy;
    negUcap	= cNegUcap;
    negFy	= cNegFy;
    negFcap	= cNegFcap;
    negUglobal	= cNegUglobal;
    negFglobal	= cNegFglobal;

    negUres	= cNegUres;
    negFres	= cNegFres;
    negKp	= cNegKp;
    negKpc	= cNegKpc;


    K_unload	= cK_unload;



    engAcml	= cEngAcml;
    engDspt	= cEngDspt;

    posUlocal	= cPosUlocal;
    posFlocal	= cPosFlocal;
    negUlocal	= cNegUlocal;
    negFlocal	= cNegFlocal;

    Failure_Flag	= cFailure_Flag;
    // Excursion_Flag	= cExcursion_Flag;
    exBranch	= cExBranch;
    Branch		= cBranch;
    // TargetPeak_Flag	= cTargetPeak_Flag;
    Yield_Flag	= cYield_Flag;
    // Reversal_Flag	= cReversal_Flag;

    // u0	= cu0;

    K_reload	= cK_reload;

    return 0;
}

int IMKPeakOriented::revertToStart(void)
{
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\\
    //////////////////////////////////////////////////////////////////// ONE TIME CALCULATIONS ////////////////////////////////////////////////////////////////////\\
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/



    betaS	= 0;
    betaC	= 0;
    betaK	= 0;
    betaA	= 0;

    Uy_pos  	= Fy_pos / Ke;
    Ucap_pos	= Uy_pos + Up_pos;
    Fcap_pos	= FcapFy_pos*Fy_pos;
    Kp_pos 		= (Fcap_pos - Fy_pos) / Up_pos;
    Kpc_pos 	= Fcap_pos / Upc_pos;

    Uy_neg 		= Fy_neg / Ke;
    Ucap_neg	= Uy_neg + Up_neg;
    Fcap_neg	= FcapFy_neg*Fy_neg;
    Kp_neg 		= (Fcap_neg - Fy_neg) / Up_neg;
    Kpc_neg 	= Fcap_neg / Upc_neg;

    posUglobal	= cPosUglobal	= Uy_pos;
    posFglobal	= cPosFglobal	= Fy_pos;
    negUglobal	= cNegUglobal	= -Uy_neg;
    negFglobal	= cNegFglobal	= -Fy_neg;

    posUy		= cPosUy	=  Uy_pos;
    posFy    	= cPosFy 	=  Fy_pos;
    posKp    	= cPosKp 	=  Kp_pos;
    posKpc   	= cPosKpc	= -Kpc_pos;
    negUy    	= cNegUy	= -Uy_neg;
    negFy    	= cNegFy	= -Fy_neg;
    negKp    	= cNegKp	=  Kp_neg;
    negKpc   	= cNegKpc	= -Kpc_neg;

    posUcap  	= cPosUcap	= Ucap_pos;
    posFcap  	= cPosFcap	= Fcap_pos;
    posFres  	= cPosFres	= Fy_pos*ResF_pos;
    negUcap  	= cNegUcap	= -Ucap_neg;
    negFcap  	= cNegFcap	= -Fcap_neg;
    negFres  	= cNegFres	= -Fy_neg*ResF_neg;

    K_unload	= cK_unload	= Ke;

    posUres	= cPosUres	= (posFres - posFcap) / posKpc + posUcap;
    negUres	= cNegUres	= (negFres - negFcap) / negKpc + negUcap;

    engAcml 	= cEngAcml	= 0.0;
    engDspt	= cEngDspt	= 0.0;

    u0	= 0.0;

    EtS	= LAMBDA_S *Fy_pos;
    EtC	= LAMBDA_C *Fy_pos;
    EtA	= LAMBDA_A *Fy_pos;
    EtK	= LAMBDA_K *Fy_pos;

    Failure_Flag 	= cFailure_Flag	  	= false;
    Excursion_Flag 	= cExcursion_Flag 	= false;
    Branch      	= cBranch         	= 0;
    exBranch    	= cExBranch   	  	= false;
    // TargetPeak_Flag 	= cTargetPeak_Flag 	= false;
    Yield_Flag		= cYield_Flag	  	= false;
    Reversal_Flag	= cReversal_Flag  	= false;

    posUlocal	= cPosUlocal	=  Uy_pos;
    posFlocal	= cPosFlocal	=  Fy_pos;
    negUlocal	= cNegUlocal	= -Uy_neg;
    negFlocal	= cNegFlocal	= -Fy_neg;

    K_reload		= cK_reload 	= Ke;

    cdu_i_1	= 0;

    //initially I zero everything   
    U	= cU	= 0;
    Ui 	= cUi 	= 0;
    Fi 	= cFi 	= 0;
    Ui_1  	= cUi_1 	= 0;
    Fi_1	= cFi_1	= 0;

    TangentK	= cTangentK	= Ke;
    //cout << " revertToStart:" << endln; //<< " U=" << U << " Ui=" << Ui << " TanK=" << TangentK << endln;

    return 0;
}

UniaxialMaterial *
IMKPeakOriented::getCopy(void)
{
    IMKPeakOriented *theCopy	= new IMKPeakOriented(this->getTag(), Ke,
        Uy_pos, Ucap_pos, Uu_pos, Fy_pos, FcapFy_pos, ResF_pos,
        Uy_neg, Ucap_neg, Uu_neg, Fy_neg, FcapFy_neg, ResF_neg,
        LAMBDA_S, LAMBDA_C, LAMBDA_A, LAMBDA_K, c_S, c_C, c_A, c_K, D_pos, D_neg);

    //cout << " getCopy" << endln;

    theCopy->U	= U;
    theCopy->cU	= cU;

    theCopy->TangentK	= TangentK;

    theCopy->Ui		= Ui;
    theCopy->Fi		= Fi;
    theCopy->Ui_1	= Ui_1;
    theCopy->Fi_1	= Fi_1;
    theCopy->du_i_1	= du_i_1;

    theCopy->posUy		= posUy;
    theCopy->posUcap	= posUcap;
    theCopy->posFy		= posFy;
    theCopy->posFcap	= posFcap;
    theCopy->posUglobal	= posUglobal;
    theCopy->posFglobal	= posFglobal;
    theCopy->posUres	= posUres;
    theCopy->posFres	= posFres;
    theCopy->posKp		= posKp;
    theCopy->posKpc	= posKpc;

    theCopy->negUy		= negUy;
    theCopy->negUcap	= negUcap;
    theCopy->negFy		= negFy;
    theCopy->negFcap	= negFcap;
    theCopy->negUglobal	= negUglobal;
    theCopy->negFglobal	= negFglobal;
    theCopy->negUres	= negUres;
    theCopy->negFres	= negFres;
    theCopy->negKp		= negKp;
    theCopy->negKpc	= negKpc;

    theCopy->K_unload	= K_unload;

    theCopy->engAcml		= engAcml;
    theCopy->engDspt	= engDspt;

    theCopy->u0	= u0;

    theCopy->posUlocal	= posUlocal;
    theCopy->posFlocal	= posFlocal;
    theCopy->negUlocal	= negUlocal;
    theCopy->negFlocal	= negFlocal;

    theCopy->Failure_Flag 		= Failure_Flag;
    theCopy->Excursion_Flag 	= Excursion_Flag;
    theCopy->exBranch 	= exBranch;
    // theCopy->TargetPeak_Flag= TargetPeak_Flag;
    theCopy->Yield_Flag 		= Yield_Flag;
    theCopy->Reversal_Flag		= Reversal_Flag;

    theCopy->K_reload	= K_reload;

    theCopy->cTangentK	= cTangentK;

    theCopy->cUi	= cUi;
    theCopy->cFi	= cFi;
    theCopy->cUi_1	= cUi_1;
    theCopy->cFi_1	= cFi_1;
    theCopy->cdu_i_1	= cdu_i_1;

    theCopy->cPosUy	= cPosUy;
    theCopy->cPosUcap	= cPosUcap;
    theCopy->cPosFy	= cPosFy;
    theCopy->cPosFcap	= cPosFcap;
    theCopy->cPosUglobal	= cPosUglobal;
    theCopy->cPosFglobal	= cPosFglobal;
    theCopy->cPosUres	= cPosUres;
    theCopy->cPosFres	= cPosFres;
    theCopy->cPosKp	= cPosKp;
    theCopy->cPosKpc	= cPosKpc;

    theCopy->cNegUy	= cNegUy;
    theCopy->cNegUcap	= cNegUcap;
    theCopy->cNegFy	= cNegFy;
    theCopy->cNegFcap	= cNegFcap;
    theCopy->cNegUglobal	= cNegUglobal;
    theCopy->cNegFglobal	= cNegFglobal;
    theCopy->cNegUres	= cNegUres;
    theCopy->cNegFres	= cNegFres;
    theCopy->cNegKp	= cNegKp;
    theCopy->cNegKpc	= cNegKpc;

    theCopy->cK_unload	= cK_unload;

    theCopy->cEngAcml 	= cEngAcml;
    theCopy->cEngDspt	= cEngDspt;

    theCopy->cu0	= cu0;

    theCopy->cPosUlocal	= cPosUlocal;
    theCopy->cPosFlocal	= cPosFlocal;
    theCopy->cNegUlocal	= cNegUlocal;
    theCopy->cNegFlocal	= cNegFlocal;

    theCopy->cFailure_Flag		= cFailure_Flag;
    theCopy->cExcursion_Flag	= cExcursion_Flag;
    theCopy->cExBranch			= cExBranch;
    // theCopy->cTargetPeak_Flag	= cTargetPeak_Flag;
    theCopy->cYield_Flag 		= cYield_Flag;
    theCopy->cReversal_Flag		= cReversal_Flag;

    theCopy->cK_reload	= cK_reload;

    return theCopy;
}

int IMKPeakOriented::sendSelf(int cTag, Channel &theChannel)
{
    int res	= 0;
    cout << " sendSelf" << endln;

    static Vector data(137);
    data(0)	= this->getTag();
    data(1)  	= Ke;
    data(2)  	= Uy_pos;
    data(3)  	= Ucap_pos;
    data(4)  	= Uu_pos;
    data(5)  	= Fy_pos;
    data(6)  	= FcapFy_pos;
    data(7)  	= ResF_pos;
    data(8)  	= Uy_neg;
    data(9)  	= Ucap_neg;
    data(10) 	= Uu_neg;
    data(11) 	= Fy_neg;
    data(12) 	= FcapFy_neg;
    data(13) 	= ResF_neg;
    data(14) 	= LAMBDA_S;
    data(15) 	= LAMBDA_C;
    data(16) 	= LAMBDA_A;
    data(17) 	= LAMBDA_K;
    data(18) 	= c_S;
    data(19) 	= c_C;
    data(20) 	= c_A;
    data(21) 	= c_K;
    data(22) 	= D_pos;
    data(23) 	= D_neg;
    data(24) 	= Ui;
    data(25) 	= Fi;
    data(26) 	= Ui_1;
    data(27) 	= Fi_1;
    data(28) 	= du_i_1;

    data(29) 	= posUy;
    data(30) 	= posUcap;
    data(31) 	= posFy;
    data(32) 	= posFcap;
    data(33) 	= posUglobal;
    data(34) 	= posFglobal;
    data(35) 	= posUres;
    data(36) 	= posFres;
    data(37) 	= posKp;
    data(38) 	= posKpc;

    data(39) 	= negUy;
    data(40) 	= negUcap;
    data(41) 	= negFy;
    data(42) 	= negFcap;
    data(43) 	= negUglobal;
    data(44) 	= negFglobal;
    data(45) 	= negUres;
    data(46) 	= negFres;
    data(47) 	= negKp;
    data(48) 	= negKpc;

    data(49) 	= K_unload;

    data(50) 	= Failure_Flag;
    data(51) 	= Excursion_Flag;
    data(52) 	= Branch;
    data(53) 	= exBranch;
    // data(54) 	= TargetPeak_Flag;
    data(55) 	= Yield_Flag;

    data(56) 	= engAcml;
    data(57) 	= engDspt;

    data(58) 	= u0;
    data(59) 	= dU;
    data(60) 	= dF;

    data(61) 	= FailS;
    data(62) 	= FailC;
    data(63) 	= FailA;
    data(64) 	= FailK;

    data(65)	= Ei;
    data(66)	= dEi;
    data(67)	= Epj;
    data(68)	= EpjK;
    data(69)	= EiK;
    data(70)	= c_S;
    data(71)	= c_C;
    data(72)	= c_A;
    data(73)	= c_K;
    data(74)	= EtS;
    data(75)	= EtC;
    data(76)	= EtA;
    data(77)	= EtK;
    data(78)	= betaS;
    data(79)	= betaC;
    data(80)	= betaA;
    data(81)	= betaK;
    data(82)	= sPCsp;
    data(83)	= sPCpcp;

    data(84)	= TangentK;

    data(85)	= Uy_pos;
    data(86)	= Ucap_pos;
    data(87)	= Fcap_pos;
    data(88)	= Kp_pos;
    data(89)	= Kpc_pos;

    data(90)	= Uy_neg;
    data(91)	= Ucap_neg;
    data(92)	= Fcap_neg;
    data(93)	= Kp_neg;
    data(94)	= Kpc_neg;

    data(95)	= cUi;
    data(96)	= cFi;
    data(97)	= cUi_1;
    data(98)	= cFi_1;
    data(99)	= cdu_i_1;

    data(100)	= cPosUy;
    data(101)	= cPosUcap;
    data(102)	= cPosFy;
    data(103)	= cPosFcap;
    data(104)	= cPosUglobal;
    data(105)	= cPosFglobal;
    data(106)	= cPosUres;
    data(107)	= cPosFres;
    data(108)	= cPosKp;
    data(109)	= cPosKpc;

    data(110)	= cNegUy;
    data(111)	= cNegUcap;
    data(112)	= cNegFy;
    data(113)	= cNegFcap;
    data(114)	= cNegUglobal;
    data(115)	= cNegFglobal;
    data(116)	= cNegUres;
    data(117)	= cNegFres;
    data(118)	= cNegKp;
    data(119)	= cNegKpc;

    data(120)	= cK_unload;

    data(121)	= cPosUlocal;
    data(122)	= cPosFlocal;
    data(123)	= cNegUlocal;
    data(124)	= cNegFlocal;

    data(125)	= cFailure_Flag;
    data(126)	= cExcursion_Flag;
    data(127)	= cExBranch;
    data(128)	= cBranch;
    // data(129)	= cTargetPeak_Flag;
    data(130)	= cYield_Flag;

    data(131)	= cK_reload;

    data(132)	= K_local;
    data(133)	= K_global;
    // data(134)	= K_check;

    data(135)	= cReversal_Flag;
    data(136)	= Reversal_Flag;

    res	= theChannel.sendVector(this->getDbTag(), cTag, data);
    if (res < 0)
        opserr << "IMKPeakOriented::sendSelf() - failed to send data\n";

    return res;
}

int IMKPeakOriented::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    int res	= 0;
    static Vector data(137);
    res	= theChannel.recvVector(this->getDbTag(), cTag, data);

    if (res < 0) {
        opserr << "IMKPeakOriented::recvSelf() - failed to receive data\n";
        this->setTag(0);
    }
    else {
        cout << " recvSelf" << endln;
        this->setTag((int)data(0));
        Ke					= data(1);
        Up_pos				= data(2);
        Upc_pos				= data(3);
        Uu_pos				= data(4);
        Fy_pos				= data(5);
        FcapFy_pos			= data(6);
        ResF_pos			= data(7);
        Up_neg				= data(8);
        Upc_neg				= data(9);
        Uu_neg				= data(10);
        Fy_neg				= data(11);
        FcapFy_neg			= data(12);
        ResF_neg			= data(13);
        LAMBDA_S			= data(14);
        LAMBDA_C			= data(15);
        LAMBDA_A			= data(16);
        LAMBDA_K			= data(17);
        c_S					= data(18);
        c_C					= data(19);
        c_A					= data(20);
        c_K					= data(21);
        D_pos				= data(22);
        D_neg				= data(23);
        Ui					= data(24);
        Fi					= data(25);
        Ui_1				= data(26);
        Fi_1				= data(27);
        du_i_1				= data(28);
        posUy				= data(29);
        posUcap			= data(30);
        posFy				= data(31);
        posFcap			= data(32);
        posUglobal		= data(33);
        posFglobal		= data(34);
        posUres			= data(35);
        posFres			= data(36);
        posKp				= data(37);
        posKpc			= data(38);
        negUy				= data(39);
        negUcap			= data(40);
        negFy				= data(41);
        negFcap			= data(42);
        negUglobal		= data(43);
        negFglobal		= data(44);
        negUres			= data(45);
        negFres			= data(46);
        negKp				= data(47);
        negKpc			= data(48);
        Failure_Flag		= data(49);
        Excursion_Flag		= data(50);
        exBranch			= data(51);
        Branch				= data(52);
        // TargetPeak_Flag		= data(53);
        Yield_Flag			= data(54);
        K_unload			= data(55);
        engAcml			= data(56);
        engDspt			= data(57);
        u0					= data(58);
        dU					= data(59);
        dF					= data(60);
        FailS				= data(61);
        FailC				= data(62);
        FailA				= data(63);
        FailK				= data(64);
        Ei					= data(65);
        dEi					= data(66);
        Epj					= data(67);
        EpjK				= data(68);
        EiK					= data(79);
        c_S					= data(70);
        c_C					= data(71);
        c_A					= data(72);
        c_K					= data(73);
        EtS					= data(74);
        EtC					= data(75);
        EtA					= data(76);
        EtK					= data(77);
        betaS				= data(78);
        betaC				= data(79);
        betaA				= data(80);
        betaK				= data(81);
        sPCsp				= data(82);
        sPCpcp				= data(83);
        TangentK			= data(84);
        Uy_pos				= data(85);
        Ucap_pos			= data(86);
        Fcap_pos			= data(87);
        Kp_pos				= data(88);
        Kpc_pos				= data(89);
        Uy_neg				= data(90);
        Ucap_neg			= data(91);
        Fcap_neg			= data(92);
        Kp_neg				= data(93);
        Kpc_neg				= data(94);
        cUi					= data(95);
        cFi					= data(96);
        cUi_1				= data(97);
        cFi_1				= data(98);
        cdu_i_1				= data(99);
        cPosUy			= data(100);
        cPosUcap			= data(101);
        cPosFy			= data(102);
        cPosFcap			= data(103);
        cPosUglobal		= data(104);
        cPosFglobal		= data(105);
        cPosUres			= data(106);
        cPosFres			= data(107);
        cPosKp			= data(108);
        cPosKpc			= data(109);
        cNegUy			= data(110);
        cNegUcap			= data(111);
        cNegFy			= data(112);
        cNegFcap			= data(113);
        cNegUglobal		= data(114);
        cNegFglobal		= data(115);
        cNegUres			= data(116);
        cNegFres			= data(117);
        cNegKp			= data(118);
        cNegKpc			= data(119);
        cK_unload			= data(120);
        cPosUlocal		= data(121);
        cPosFlocal		= data(122);
        cNegUlocal		= data(123);
        cNegFlocal		= data(124);
        cFailure_Flag		= data(125);
        cExcursion_Flag		= data(126);
        cExBranch			= data(127);
        cBranch				= data(128);
        // cTargetPeak_Flag	= data(129);
        cYield_Flag   		= data(130);
        cK_reload   		= data(131);
        K_local				= data(132);
        K_global			= data(133);
        // K_check			= data(134);
        cReversal_Flag	= data(135);
        Reversal_Flag	= data(136);
    }

    return res;
}

void IMKPeakOriented::Print(OPS_Stream &s, int flag)
{
    cout << "IMKPeakOriented tag: " << this->getTag() << endln;
}
