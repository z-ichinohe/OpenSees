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
    UniaxialMaterial *theMaterial = 0;

    int    iData[1];
    double dData[7];
    int numInt = 1;

    if (OPS_GetIntInput(&numInt, iData) != 0) {
        opserr << "WARNING invalid uniaxialMaterial TakedaSlip tag" << endln;
        return 0;
    }

    int numDouble = 7;

    if (OPS_GetDoubleInput(&numDouble, dData) != 0) {
        opserr << "Invalid Args want: uniaxialMaterial TakedaSlip tag?";
        opserr << "Crack_Disp? Yield_Disp?";
        opserr << "Crack_Stiffness? Yield_Stiffness? Plastic_Stiffness?";
        opserr << "b0? b1?";
        return 0;
    }


    // Parsing was successful, allocate the material
    theMaterial = new TakedaSlip(iData[0],
        dData[0], dData[1],
        dData[2], dData[3], dData[4],
        dData[5], dData[6]);

    if (theMaterial == 0) {
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
    Uc(p_Uc), Uy(p_Uy),
    Kc(p_Kc), Ky(p_Ky), Kp(p_Kp),
    b0(p_b0), b1(p_b1)
{
    this->revertToStart();
}

TakedaSlip::TakedaSlip()
    :UniaxialMaterial(0, 0),
    Uc(0), Uy(0),
    Kc(0), Ky(0), Kp(0),
    b0(0), b1(0)
{
    this->revertToStart();
}

TakedaSlip::~TakedaSlip()
{
    // does nothing
}

int TakedaSlip::setTrialStrain(double strain, double strainRate)
{
    //all variables to the last commit
    this->revertToLastCommit();

    //state determination algorithm: defines the current force and tangent stiffness
    const double b2 = 3;
    const double b3 = 1;
    const double kappaD = 0.5;
    const double kappaF = 0.5;
    const double Ui_1 = Ui;
    const double Fi_1 = Fi;
    U = strain; //set trial displacement
    Ui = U;
    const double dU = Ui - Ui_1;    // Incremental deformation at current step
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////  MAIN CODE //////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if (dU == 0) {   // When deformation doesn't change from the last
        Fi  = Fi_1;
    } else {
        const bool onBackbone = (Branch > 1);
    ///////////////////////////////////////////////////////////////////////////////////////////
    /////////////////// WHEN REVERSAL /////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
        if ( (onBackbone && Fi_1*dU < 0) || (onBackbone && Fi_1==0 && Ui_1*dU<=0) ) {
            Branch          = 1;
    /////////////////////////// UPDATE PEAK POINTS ////////////////////////////////////////////
            if ( Fi_1 > 0 ){
                posUlocal	= Ui_1;           // UPDATE LOCAL
                posFlocal	= Fi_1;
                if ( Ui_1 > posUglobal ) {    // UPDATE GLOBAL
                    posUglobal   	= Ui_1;
                    posFglobal   	= Fi_1;
                }
            } else {
                negUlocal	= Ui_1;           // UPDATE LOCAL
                negFlocal	= Fi_1;
                if ( Ui_1 < negUglobal ) {    // UPDATE GLOBAL
                    negUglobal   	= Ui_1;
                    negFglobal   	= Fi_1;
                }
            }
    /////////////////// UPDATE UNLOADING STIFFNESS ////////////////////////////////////////////
        }
        Fi	    = Fi_1 + Kunload * dU;
    ///////////////////////////////////////////////////////////////////////////////////////////
    /////////////////// WHEN NEW EXCURSION /////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
        if (Branch == 1 && Fi_1*Fi <= 0.0) {
    ////////////////////////// RELOADING TARGET DETERMINATION /////////////////////////////////
            if (dU > 0) {
                double u0 	    = Ui_1 - (Fi_1 / Kunload);
                double Uplstc   = posUglobal    - (posFglobal / Kunload);
                Upinch          = (1-kappaD)*Uplstc;
                Fpinch          = kappaF*posFglobal*(Upinch-u0)/(posUglobal-u0);
                double Kpinch   = Fpinch        / (Upinch       - u0);
                double Kglobal  = posFglobal    / (posUglobal - u0);
                double Klocal   = posFlocal     / (posUlocal  - u0);
                if (u0 < Upinch) {
                    Branch  = 2;
                    Kreload = Kpinch;
                } else  if ( u0 < posUlocal && posFlocal < posFglobal && Klocal > Kglobal) {
                    Branch  = 3;
                    Kreload = Klocal;
                } else {
                    Branch  = 4;
                    Kreload = Kglobal;
                }
            }
            else {
                double u0 	    = Ui_1 - (Fi_1 / Kunload);
                double Uplstc   = negUglobal    - (negFglobal / Kunload);
                Upinch          = (1-kappaD)*Uplstc;
                Fpinch          = kappaF*negFglobal*(Upinch-u0)/(negUglobal-u0);
                double Kpinch   = Fpinch        / (Upinch       - u0);
                double Kglobal  = negFglobal    / (negUglobal - u0);
                double Klocal   = negFlocal     / (negUlocal  - u0);
                if (u0 > Upinch) {
                    Branch 	= 12;
                    Kreload = Kpinch;
                } else  if ( u0 > negUlocal && negFlocal > negFglobal && Klocal > Kglobal) {
                    Branch  = 13;
                    Kreload = Klocal;
                } else {
                    Branch  = 14;
                    Kreload = Kglobal;
                }
            }
         }
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    ////////////////// BRANCH SHIFT CHECK AND TANGENT STIFFNESS UPDATE ////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    //  Branch
    //      0:  Elastic
    //      1:  Unloading Branch
    //      2:  Towards Pinching Point  +
    //      3:  Towards Local Peak      +
    //      4:  Towards Global Peak     +
    //      5:  Towards Yield Point  +
    //      6:  Plastic Branch         +
    //      12: Towards Pinching Point  -
    //      13: Towards Local Peak      -
    //      14: Towards Global Peak     -
    //      15: Towards Yield Point  -
    //      16: Plastic Branch         -
    // Branch shifting from 2 -> 3 -> 4 -> 5 -> 6 -> 7
        // int exBranch = Branch;
        if (Branch == 0 && Ui > posUc) {            // Yield in Positive
            Branch  = 5;
        } else if (Branch == 0 && Ui < negUc) {     // Yield in Negative
            Branch  = 15;
        } else if (Branch == 1 && Fi_1 > 0 && Ui > posUlocal) {
            const double Kpinch  = (Fpinch       - posFlocal) / (Upinch      - posUlocal);
            const double Kglobal = (posFglobal   - posFlocal) / (posUglobal  - posUlocal);
            if (posUlocal < Upinch && posFlocal < Fpinch && Upinch < posUglobal && Fpinch < posFglobal && Kpinch < Kglobal) {
                Kreload = Kpinch;
                Branch  = 2;                        // Pinching Branch
            } else {
                Kreload = Kglobal;
                Branch  = 4;                        // Towards Global Peak
            }
        } else if (Branch == 1 && Fi_1 < 0 && Ui < negUlocal) {    // Back to Reloading (Negative)
            const double Kpinch  = (Fpinch       - negFlocal) / (Upinch      - negUlocal);
            const double Kglobal = (negFglobal   - negFlocal) / (negUglobal  - negUlocal);
            if (negUglobal < Upinch && negFglobal < Fpinch && Upinch < negUlocal && Fpinch < negFlocal && Kpinch < Kglobal) {
                Kreload = Kpinch;
                Branch  = 12;                   // Pinching Branch
            } else {
                Kreload = Kglobal;
                Branch  = 14;                   // Towards Global Peak
            }
        }
    // Positive
        if (Branch == 2 && Ui > Upinch) {
            const double Kglobal  = (posFglobal   - Fpinch) / (posUglobal - Upinch);
            const double Klocal   = (posFlocal    - Fpinch) / (posUlocal  - Upinch);
            if (Upinch < posUlocal && Fpinch < posFlocal && posFlocal < posFglobal && Klocal > Kglobal) {
                Kreload = Klocal;
                Branch  = 3;
            } else {
                Kreload = Kglobal;
                Branch 	= 4;
            }
        }
        if (Branch == 3 && Ui > posUlocal) {
            Kreload     = (posFglobal - posFlocal) / (posUglobal - posUlocal);
            Branch      = 4;
        }
        if (Branch == 4 && Ui > posUglobal) {
            Branch 	    = 5;
        }
        if (Branch == 5 && Ui > posUy) {
            Branch 	    = 6;
        }
    // Negative
        if (Branch == 12 && Ui < Upinch) {
            double Kglobal  = (negFglobal   - Fpinch) / (negUglobal - Upinch);
            double Klocal   = (negFlocal    - Fpinch) / (negUlocal  - Upinch);
            if (negFglobal < negFlocal && negUlocal < Upinch && negFlocal < Fpinch && Klocal > Kglobal) {
                Kreload = Klocal;
                Branch  = 13;
            } else {
                Kreload = Kglobal;
                Branch 	= 14;
            }
        }
        if (Branch == 13 && Ui < negUlocal) {
            Kreload     = (negFglobal - negFlocal) / (negUglobal - negUlocal);
            Branch      = 14;
        }
        if (Branch == 14 && Ui < negUglobal) {
            Branch 	    = 15;
        }
        if (Branch == 15 && Ui < negUy) {
            Branch 	    = 16;
        }
    // Branch Change check
        // if (Branch!=exBranch) {
        //     std::cout << exBranch << " -> " << Branch << "\n";
        // }
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////// COMPUTE FORCE BASED ON BRANCH /////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
        if (Branch == 0) {
            Fi  = Kc*Ui;
        } else if (Branch == 1) {
            Fi  = Fi_1 + Kunload*dU;
    // Positive
        } else if (Branch == 2) {
            Fi  = Fpinch    + Kreload*(Ui - Upinch);
        } else if (Branch == 3) {
            Fi  = posFlocal + Kreload*(Ui - posUlocal);
        } else if (Branch == 4) {
            Fi  = posFglobal + Kreload*(Ui - posUglobal);
        } else if (Branch == 5) {
            Fi  = posFy   + Ky*(Ui - posUy);
        } else if (Branch == 6) {
            Fi  = posFy   + Kp*(Ui - posUy);
    // Negative
        } else if (Branch == 12) {
            Fi  = Fpinch    + Kreload*(Ui - Upinch);
        } else if (Branch == 13) {
            Fi  = negFlocal + Kreload*(Ui - negUlocal);
        } else if (Branch == 14) {
            Fi  = negFglobal + Kreload*(Ui - negUglobal);
        } else if (Branch == 15) {
            Fi  = negFy   + Ky*(Ui - negUy);
        } else if (Branch == 16) {
            Fi  = negFy   + Kp*(Ui - negUy);
        }
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    // CHECK FOR FAILURE
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
        KgetTangent  = (Fi - Fi_1) / dU;
    }
    if (KgetTangent==0) {
        KgetTangent  = 1e-6;
    }
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////// END OF MAIN CODE ///////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    return 0;
}

double TakedaSlip::getStress(void)
{
    //cout << " getStress" << endln;
    return (Fi);
}

double TakedaSlip::getTangent(void)
{
    //cout << " getTangent" << endln;
    return (KgetTangent);
}

double TakedaSlip::getInitialTangent(void)
{
    //cout << " getInitialTangent" << endln;
    return (Kc);
}

double TakedaSlip::getStrain(void)
{
    //cout << " getStrain" << endln;
    return (U);
}

int TakedaSlip::commitState(void)
{
    //cout << " commitState" << endln;

    //commit trial  variables
// 10 Pos U and F
    cPosUlocal	= posUlocal;
    cPosFlocal	= posFlocal;
    cPosUglobal	= posUglobal;
    cPosFglobal	= posFglobal;
    cNegUlocal	= negUlocal;
    cNegFlocal	= negFlocal;
    cNegUglobal	= negUglobal;
    cNegFglobal	= negFglobal;
    cFpinch     = Fpinch;
    cUpinch     = Upinch;
// 3 State
    cU		    = U;
    cUi	    	= Ui;
    cFi	        = Fi;
// 2 Stiffness
    cKreload	= Kreload;
    cKunload	= Kunload;
// 1 Flag
    cBranch         = Branch;
    return 0;
}

int TakedaSlip::revertToLastCommit(void)
{
    //cout << " revertToLastCommit" << endln;
    //the opposite of commit trial history variables
// 10 Positive U and F
    posUlocal	    = cPosUlocal;
    posFlocal	    = cPosFlocal;
    posUglobal	    = cPosUglobal;
    posFglobal	    = cPosFglobal;
    negUlocal	    = cNegUlocal;
    negFlocal	    = cNegFlocal;
    negUglobal	    = cNegUglobal;
    negFglobal	    = cNegFglobal;
    Fpinch          = cFpinch;
    Upinch          = cUpinch;
// 3 State Variables
    U	            = cU;
    Ui	            = cUi;
    Fi	            = cFi;
// 2 Stiffness
    Kreload	    = cKreload;
    Kunload	        = cKunload;
// 1 Flag
    Branch          = cBranch;
    return 0;
}

int TakedaSlip::revertToStart(void)
{
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\\
    //////////////////////////////////////////////////////////////////// ONE TIME CALCULATIONS ////////////////////////////////////////////////////////////////////\\
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
// 14 Initial Values
    posUc = Uc;
    posUy = Uy;
    posFc = Uc * Kc;
    posFy = Uc * Kc + (Uy - Uc) * Ky;
    negUc = -Uc;
    negUy = -Uy;
    negFc = - Uc * Kc;
    negFy = - Uc * Kc - (Uy - Uc) * Ky;
// 10 U and F
    posUlocal	= cPosUlocal	= posUc;
    posFlocal	= cPosFlocal	= posFc;
    posUglobal	= cPosUglobal	= posUc;
    posFglobal	= cPosFglobal	= posFc;
    negUlocal	= cNegUlocal	= -negUc;
    negFlocal	= cNegFlocal	= -negFc;
    negUglobal	= cNegUglobal	= -negUc;
    negFglobal	= cNegFglobal	= -negFc;
    Fpinch      = cFpinch       = 0.0;
    Upinch      = cUpinch       = 0.0;
// 3 State Values
    U	        = cU	        = 0;
    Ui      	= cUi 	        = 0;
    Fi 	        = cFi 	        = 0;
// 2 Stiffness
    Kreload	    = cKreload	    = Kc;
    Kunload	    = cKunload	    = Kc;
    KgetTangent = Kc;
// 1 Flag
    Branch      	= cBranch       = 0;
    return 0;
}

UniaxialMaterial *
TakedaSlip::getCopy(void)
{
    TakedaSlip *theCopy = new TakedaSlip(
        this->getTag(),
        Uc, Uy,
        Kc, Ky, Kp,
        b0, b1);
// 10 U and F
    theCopy->posUlocal = posUlocal;
    theCopy->posFlocal = posFlocal;
    theCopy->posUglobal = posUglobal;
    theCopy->posFglobal = posFglobal;
    theCopy->negUlocal = negUlocal;
    theCopy->negFlocal = negFlocal;
    theCopy->negUglobal = negUglobal;
    theCopy->negFglobal = negFglobal;
// 2 Pinching
    theCopy->Fpinch = Fpinch;
    theCopy->Upinch = Upinch;
// 3 State Values
    theCopy->U = U;
    theCopy->Ui = Ui;
    theCopy->Fi = Fi;
// 2 Stiffness
    theCopy->Kreload = Kreload;
    theCopy->Kunload = Kunload;
// 2 Flag
    theCopy->Branch = Branch;
// 10 U and F
    theCopy->cPosUlocal = cPosUlocal;
    theCopy->cPosFlocal = cPosFlocal;
    theCopy->cPosUglobal = cPosUglobal;
    theCopy->cPosFglobal = cPosFglobal;
    theCopy->cNegUglobal = cNegUglobal;
    theCopy->cNegFglobal = cNegFglobal;
    theCopy->cNegUlocal = cNegUlocal;
    theCopy->cNegFlocal = cNegFlocal;
// 3 State
    theCopy->cU = cU;
    theCopy->cUi = cUi;
    theCopy->cFi = cFi;
// 2 Stiffness
    theCopy->cKreload = cKreload;
    theCopy->cKunload = cKunload;
// 2 Pinching
    theCopy->cFpinch = cFpinch;
    theCopy->cUpinch = cUpinch;
// 2 Flag
    theCopy->cBranch = cBranch;
    return theCopy;
}

int TakedaSlip::sendSelf(int cTag, Channel &theChannel)
{
    int res = 0;
    cout << " sendSelf" << endln;

    static Vector data(137);
    data(0) = this->getTag();
// 25 Fixed Input Material Parameters 1-25
    data(1)  	= Uc;
    data(2)  	= Uy;
    data(3)  	= Kc;
    data(4)  	= Ky;
    data(5)  	= Kp;
    data(6)  	= b0;
    data(7)  	= b1;
// 14 Initial Values 31-44
    data(31)	= posUc;
    data(32)	= posUy;
    data(33)	= posFc;
    data(34)	= posFy;
    data(36)	= negUc;
    data(37)	= negUy;
    data(38)	= negFc;
    data(39)	= negFy;
// 12 Positive U and F 51-62
    data(55) 	= posUlocal;
    data(56) 	= posFlocal;
    data(57) 	= posUglobal;
    data(58) 	= posFglobal;
// 3 State Variables 63-65
    data(63)    = U;
    data(64) 	= Ui;
    data(65) 	= Fi;
// 2 Stiffness 66 67
    data(66)	= Kreload;
    data(67) 	= Kunload;
// 12 Negative U and F 71-82
    data(75) 	= negUlocal;
    data(76) 	= negFlocal;
    data(77) 	= negUglobal;
    data(78) 	= negFglobal;
// 2 Pinching 83 84
    data(83)    = Fpinch;
    data(84)    = Upinch;
// 2 Flag 85 86
    data(86) 	= Branch;
// 12 Positive U and F 101-112
    data(105)	= cPosUlocal;
    data(106)	= cPosFlocal;
    data(107)	= cPosUglobal;
    data(108)	= cPosFglobal;
// 3 State Variables 113-115
    data(113)   = cU;
    data(114)	= cUi;
    data(115)	= cFi;
// 2 Stiffness 116 117
    data(116)   = cKreload;
    data(117)	= cKunload;
// 12 Negative U and F 121-132
    data(125)	= cNegUlocal;
    data(126)	= cNegFlocal;
    data(127)	= cNegUglobal;
    data(128)	= cNegFglobal;
// 2 Pinching 133 134
    data(133)   = cFpinch;
    data(134)   = cUpinch;
// 2 Flag 135 136
    data(136)	= cBranch;
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
    // 25 Fixed Input Material Parameters
        Uc				= data(1);
        Uy			= data(2);
        Kc		= data(3);
        Ky			= data(4);
        Kp			= data(5);
        b0		= data(6);
        b1		= data(7);
    // 14 Initial Values
        posUc			= data(31);
        posUy		= data(32);
        posFc		= data(33);
        posFy			= data(34);
        negUy			= data(36);
        negUc		= data(37);
        negFy		= data(38);
        negFc			= data(39);
    // 12 Positive U and F
        posUlocal	    = data(55);
        posFlocal	    = data(56);
        posUglobal	    = data(57);
        posFglobal	    = data(58);
    // 3 State Variables
        U               = data(63);
        Ui				= data(64);
        Fi				= data(65);
    // 2 Stiffness
        Kreload		= data(66);
        Kunload	    	= data(67);
    // 12 Negative U and F
        negUlocal	    = data(75);
        negFlocal  	    = data(76);
        negUglobal	    = data(77);
        negFglobal  	= data(78);
    // 2 Pinching
        Fpinch          = data(83);
        Upinch          = data(84);
    // 2 Flag
        Branch			= data(86);
    // 12 Positive U and F
        cPosUlocal		= data(105);
        cPosFlocal		= data(106);
        cPosUglobal		= data(107);
        cPosFglobal		= data(108);
    // 3 State Variables
        cU              = data(113);
        cUi				= data(114);
        cFi				= data(115);
    // 2 Stiffness
        cKreload       = data(116);
        cKunload		= data(117);
    // 12 Negative U and F
        cNegUlocal		= data(125);
        cNegFlocal		= data(126);
        cNegUglobal		= data(127);
        cNegFglobal		= data(128);
    // 2 Pinching
        cFpinch         = data(133);
        cUpinch         = data(134);
    // 2 Flag
        cBranch			= data(136);
    }

    return res;
}

void TakedaSlip::Print(OPS_Stream &s, int flag)
{
    cout << "TakedaSlip tag: " << this->getTag() << endln;
}