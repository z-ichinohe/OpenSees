/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999,	The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California,	Berkeley,	is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,	 and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

//Modified Ibarra-Medina-Krawinkler with Bilin Hysteretic Response

//**********************************************************************                                                                     
// Code Developed by: Ahmed Elkady
// Lecturer,	University of Southampton
// Last Updated: October 25th 2021
//**********************************************************************

#ifndef IMKBilin_h
#define IMKBilin_h

#include <UniaxialMaterial.h>

class IMKBilin : public UniaxialMaterial
{
public:
	IMKBilin(int tag,	double	Ke,
		double	posUp_0,	double	posUpc_0,	double	posUu_0,	double	posMpe_0,	double	posMmaxMpe_0,	double	posResM_0,
		double	negUp_0,	double	negUpc_0,	double	negUu_0,	double	negMpe_0,	double	negMmaxMpe_0,	double	negResM_0,
		double	LAMBDA_S,	double	LAMBDA_C,	double	LAMBDA_K,	double	c_S,	double	c_C,	double	c_K,	double	D_pos,	double	D_neg);
	IMKBilin();
	~IMKBilin();
	const char *getClassType(void) const { return "IMKBilin"; };
	int setTrialStrain(double	strain,	double	strainRate = 0.0);
	double	getStrain(void);
	double	getStress(void);
	double	getTangent(void);
	double	getInitialTangent(void);
	int commitState(void);
	int revertToLastCommit(void);
	int revertToStart(void);
	UniaxialMaterial *getCopy(void);
	int sendSelf(int commitTag,	Channel &theChannel);
	int recvSelf(int commitTag,	Channel &theChannel,	FEM_ObjectBroker &theBroker);
	void Print(OPS_Stream &s,	int flag = 0);


protected:

private:
	//my functions


	//Fixed input material parameters 
	double	Ke;
	double	posUp_0;
	double	posUpc_0;
	double	posUu_0;
	double	posMpe_0;
	double	posMmaxMpe_0;
	double	posResM_0;
	double	negUp_0;
	double	negUpc_0;
	double	negUu_0;
	double	negMpe_0;
	double	negMmaxMpe_0;
	double	negResM_0;
	double	LAMBDA_S;
	double	LAMBDA_C;
	double	LAMBDA_K;
	double	c_S;
	double	c_C;
	double	c_K;
	double	D_pos;
	double	D_neg;

	//State variables 
	double	U,	cU;

	//History variables 
	double	posMr_0,		cPosMr_0;
	double	negMr_0,		cNegMr_0;
	double	posUy_0;
	double	posUmax_0;
	double	posKp_0;
	double	posKpc_0;
	double	posMmax_0;
	double	posMpeProj_0;
	double	posMmaxProj_0;

	double	negUy_0;
	double	negUmax_0;
	double	negKp_0;
	double	negKpc_0;
	double	negMmax_0;
	double	negMpeProj_0;
	double	negMmaxProj_0;

	double	refEnergyS;
	double	refEnergyC;
	double	refEnergyK;

	double	K_j_1,			cK_j_1;
	double	posUy_1,		cPosUy_1;
	double	posUmax_1,		cPosUmax_1;
	double	posKp_1,		cPosKp_1;
	double	posKpc_1,		cPosKpc_1;
	double	posMpe_1,		cPosMpe_1;
	double	posMpeProj_1,	cPosMpeProj_1;
	double	posMmax_1,		cPosMmax_1;
	double	posMmaxProj_1,	cPosMmaxProj_1;

	double	negUy_1,		cNegUy_1;
	double	negUmax_1,		cNegUmax_1;
	double	negKp_1,		cNegKp_1;
	double	negKpc_1,		cNegKpc_1;
	double	negMpe_1,		cNegMpe_1;
	double	negMpeProj_1,	cNegMpeProj_1;
	double	negMmax_1,		cNegMmax_1;
	double	negMmaxProj_1,	cNegMmaxProj_1;

	double	Ui,				cUi;
	double	Fi,				cFi;
	double	Di,				cDi;
	double	Ui_1,			cUi_1;
	double	Fi_1,			cFi_1;
	double	Di_1,			cDi_1;

	double	betaS_1,		cBetaS_1;
	double	betaC_1,		cBetaC_1;
	double	betaK_1,		cBetaK_1;

	double	Excursion_Flag,	cExcursion_Flag;
	double	Reversal_Flag,	cReversal_Flag;
	double	Yield_Flag,		cYield_Flag;
	double	Fail_FlagPos,	cFail_FlagPos;
	double	Fail_FlagNeg,	cFail_FlagNeg;
	double	Mrpos_Flag,		cMrpos_Flag;
	double	Mrneg_Flag,		cMrneg_Flag;
	double	Energy_Flag,	cEnergy_Flag;

	double	engExcr_1,		cEngExcr_1;
	double	engExcr,		cEngExcr;
	double	engRvrs,		cEngRvrs;
	double	engTotl,		cEngTotl;

	double	Ulocal,		cUlocal;
	double	Flocal,		cFlocal;
	double	TangentK,		cTangentK;


};


#endif

//////////////////////////// Variables Definitions /////////////////////////////
/*
Ke 					Initial elastic stiffness
posUp_0 		Initial pre-capping plastic rotation in the +ve loading direction
posUpc_0   	Initial post-capping plastic rotation in the +ve loading direction
posUu_0 		Ultimate rotation in the +ve loading direction
posMpe_0 			Initial effective plastic moment in the +ve loading direction
posMmaxMpe_0 		Initial maximum-to-effective plastic moment ratio in the +ve loading direction
posResM_0 			Residual moment in the +ve loading direction
negUp_0 		Initial pre-capping plastic rotation in the -ve loading direction
negUpc_0   	Initial post-capping plastic rotation in the -ve loading direction
negUu_0    	Ultimate rotation in the -ve loading direction
negMpe_0        	Initial effective plastic moment in the -ve loading direction
negMmaxMpe_0    	Initial maximum-to-effective plastic moment ratio in the -ve loading direction
negResM_0       	Residual moment in the -ve loading direction
LAMBDA_S 			Cyclic deterioration parameter for strength deterioration 
LAMBDA_C 			Cyclic deterioration parameter for post-capping strength deterioration
LAMBDA_K 			Cyclic deterioration parameter for unloading stiffness deterioration 
c_S 				Rate of strength deterioration.
c_C 				Rate of post-capping strength deterioration.
c_K 				Rate of unloading stiffness deterioration.
D_pos 				Rate of cyclic deterioration in the +ve loading direction
D_neg 				Rate of cyclic deterioration in the -ve loading direction
n					Parameter identifying the offset rotation on the unloading side
Roffset				Offset rotation identifying the Smooth Transition region
LAMBDA_F			Cyclic deterioration parameter for Smooth Transition deterioration 
c_F					Cyclic deterioration parameter for Smooth Transition deterioration
Ui 					Rotation at current step
Fi 					Moment at current step 
Di 					Rotation Direction at current step 
Ui_1  				Rotation at previous step
Fi_1            	Moment at previous step 
Di_1            	Rotation Direction at previous step
Ulocal 			Rotation at direction reversal points
Flocal 			Moment at direction reversal points
TangentK 			Tangent stiffness
K_j_1 				Unloading stiffness in the previous excursion
posUy_1 	Yielding rotation in the previous +ve excursion
posUmax_1 	Capping point rotation in the previous +ve excursion
posKp_1 	Pre-capping slope in the previous +ve excursion
posKpc_1 	Post-capping slope in the previous +ve excursion
posMpe_1 		Effective plastic moment in the previous +ve excursion
posMpeProj_1  Projed effective plastic moment in the previous +ve excursion
posMmax_1        Maximum moment in the previous +ve excursion
posMmaxProj_1 Projed maximum  moment in the previous +ve excursion
negUy_1     Yielding rotation in the previous -ve excursion
negUmax_1   Capping point rotation in the previous -ve excursion
negKp_1     Pre-capping slope in the previous -ve excursion
negKpc_1    Post-capping slope in the previous -ve excursion
negMpe_1         Effective plastic moment in the previous -ve excursion
negMpeProj_1  Projed effective plastic moment in the previous -ve 
negMmax_1        Maximum moment in the previous -ve excursion
negMmaxProj_1 Projed maximum  moment in the previous -ve excursion
betaS_1    
betaC_1    
betaK_1    
betaF_1 	  
refEnergyS    	Refernence energy for strength deterioration
refEnergyC    	Refernence energy for post-capping strength deterioration
refEnergyK    	Refernence energy for unloading stiffness deterioration
refEnergyF    	Refernence energy for Smooth Transition deterioration
Excursion_Flag 		Flag for Excursion occurrence (i.e.,	crossing the x-axis)
Reversal_Flag 		Flag for Loading direction reversal occurrence
Yield_Flag 			Flag for Yielding occurrence
Fail_FlagPos 		Flag for reaching the ultimate rotation in the +ve loading direction
Fail_FlagNeg 		Flag for reaching the ultimate rotation in the -ve loading direction
Mrpos_Flag 			Flag for reaching the residual moment in the +ve loading direction
Mrneg_Flag 			Flag for reaching the residual moment in the -ve loading direction
Energy_Flag 		Flag for reaching the reference energy
engExcr_1	Dissipated energy in previous excursion
engExcr 		Dissipated energy in current excursion
engRvrs 			Total dissipated energy till previous load reversal point
engTotl 		Total dissipated energy till current step
*/
