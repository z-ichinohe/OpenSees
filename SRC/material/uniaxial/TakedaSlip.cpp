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

    if (OPS_GetIntInput(&numInt, iData) != 0)  {
        opserr << "WARNING invalid uniaxialMaterial TakedaSlip tag" << endln;
        return 0;
    }

    int numDouble = 7;

    if (OPS_GetDoubleInput(&numDouble, dData) != 0)  {
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
    dc(p_Uc), dy(p_Uy),
    sc(p_Kc), sy(p_Ky), su(p_Kp),
    b0(p_b0), b1(p_b1)
{
    this->revertToStart();
}

TakedaSlip::TakedaSlip()
    :UniaxialMaterial(0, 0),
    dc(0), dy(0),
    sc(0), sy(0), su(0),
    b0(0), b1(0)
{
    this->revertToStart();
}

TakedaSlip::~TakedaSlip()
{
    // does nothing
}

std::tuple<int, float, float> TakedaSlip::Tslip_120(double dy, double dd, int sn, double fc, double dc, double sy, double su, double fy)
{
    int ll;
    double ss;
    double ff;
    if (dy - abs(dd) > 0)  {
        ll = 2;
        ss = sy;
        ff = sn * fc + (dd - sn * dc) * sy;
        return {ll, ss, ff};
    } else {
        ll = 3;
        ss = su;
        ff = fy * sn + (dd - dy * sn) * su;
        return {ll, ss, ff};
    }
}

int TakedaSlip::setTrialStrain(double strain, double strainRate)
{
    //all variables to the last commit
    this->revertToLastCommit();

    //state determination algorithm: defines the current force and tangent stiffness
    const double err = 0.0001;
    double eu, x, y, et, xf;
    b2 = 3;
    b3 = 1;
    ds = dd;
    fs = ff;
    dd = strain;
    fc = dc * sc;
    fy = fc + sy * (dy - dc);
    int is = 1;
    if (fs < 0)  {
        is = 2;
    }
    int sn = 3 - 2 * is;
    if (ll == 1)  {
        if (dc - abs(dd) > 0)  {
            ff = sc * dd;
            return 0;
        } else {
            fm[1] = fc;
            fm[2] = -fc;
            dm[1] = dc;
            dm[2] = -dc;
            s1[1] = sc;
            s1[2] = sc;
            is = 1;
            if (dd < 0)  {
                is = 2;
            }
            sn = 3 - 2 * is;
            std::tuple<int, float, float> llssff = this->Tslip_120(dy, dd, sn, fc, dc, sy, su, fy);
            ll = std::get<0>(llssff);
            ss = std::get<1>(llssff);
            ff = std::get<2>(llssff);
        }
    }

    if (ll==2) {
        if ((dd-ds)*sn >0) {
            if ( dy-abs(dd) >0) {
                ff=sn*fc+(dd-sn*dc)*sy;
                return 0;
            } else {
                ll=3;
                ss=su;
                ff=fy*sn+(dd-dy*sn)*su;
                return 0;
            }
        } else {
            fm[is]=fs;
            dm[is]=ds;
            s1[is]=(abs(fm[is])+fc)/(abs(dm[is])+dc);
            eu=b3*fm[is]/dm[is];
            x0=ds-fs/s1[is];
            if ( (dd-x0)*sn >=0) {
                ll=4;
                ss=s1[is];
                ff=fm[is]+(dd-dm[is])*s1[is];
                return 0;
            } else {
                // % loop of 430 %%
                is=3-is;
                sn=3-2*is;
                if ( abs(dm[is]) <= dc) {
                    d0=x0+sn*fc/s1[3-is];
                    f0=fc*sn;
                    if ( (d0-dd)*sn >0) {
                        ll=5;
                        ff=(dd-x0)*ss;
                        return 0;
                    } else {
                        // % loop of 530 %%
                        if ((dm[is]*sn <=dc) && (abs(dm[3-is]) <=dy)) {
                            std::tuple<int, float, float> llssff = this->Tslip_120(dy, dd, sn, fc, dc, sy, su, fy);
                            ll = std::get<0>(llssff);
                            ss = std::get<1>(llssff);
                            ff = std::get<2>(llssff);
                            return 0;
                        } else {
                            fm[is]=fy*sn;
                            dm[is]=dy*sn;
                            x0=dy*sn-fy*sn*(dy*sn-d0)/(fy*sn-f0);
                            if ( (dd-dm[is])*sn >=0) {
                                std::tuple<int, float, float> llssff = this->Tslip_120(dy, dd, sn, fc, dc, sy, su, fy);
                                ll = std::get<0>(llssff);
                                ss = std::get<1>(llssff);
                                ff = std::get<2>(llssff);
                                return 0;
                            } else {
                                ll=9;
                                ss=fm[is]/(dm[is]-x0);
                                ff=ss*(dd-x0);
                                return 0;
                            }
                        }
                    }
                } else if ( abs(dm[is]) <= dy) {
                    x=fm[is]/(dm[is]-x0);
                    y=fy*sn/(dy*sn-x0);
                    if (x < y) {
                        dm[is]=dy*sn;
                        fm[is]=fy*sn;
                    }
                    if ((dd-dm[is])*sn >=0) {
                        std::tuple<int, float, float> llssff = this->Tslip_120(dy, dd, sn, fc, dc, sy, su, fy);
                    ll = std::get<0>(llssff);
                    ss = std::get<1>(llssff);
                    ff = std::get<2>(llssff);
                        return 0;
                    } else {
                        ll=9;
                        ss=fm[is]/(dm[is]-x0);
                        ff=ss*(dd-x0);
                        return 0;
                    }
                } else {
                    // % loop of 440 %%
                    x=dm[3-is]-fm[3-is]/s1[3-is];
                    et=abs(fm[is]/(x-dm[is]))*pow( abs(dm[is]/(x-dm[is])), b2);
                    eu=fm[is]/dm[is]*b3;
                    if ( abs(et-eu) < err) {
                        xm=x0;
                    } else {
                        xm=(eu*dm[is]-et*x0-fm[is])/(eu-et);
                    }

                    if ( (dd-xm)*sn >0) {
                        if ((dd-dm[is])*sn >=0) {
                            std::tuple<int, float, float> llssff = this->Tslip_120(dy, dd, sn, fc, dc, sy, su, fy);
                            ll = std::get<0>(llssff);
                            ss = std::get<1>(llssff);
                            ff = std::get<2>(llssff);
                            return 0;
                        } else {
                            ll=7;
                            eu=fm[is]/dm[is]*b3;
                            ss=eu;
                            xf=dm[is]-fm[is]/eu;
                            x0=xf;
                            ff=eu*(dd-xf);
                            return 0;
                        }
                    } else {
                        ll=6;
                        x=dm[3-is]-fm[3-is]/s1[3-is];
                        et=abs(fm[is]/(x-dm[is]))*pow( abs(dm[is]/(x-dm[is])), b2);
                        ss=et;
                        ff=et*(dd-x0);
                        return 0;
                    }
                }
            }
        }
    }

    // % rule 3 loading on post yielding primary curve %%
    if ( ll==3) {
        if ( (dd-ds)*sn >0) {
            ff=fy*sn+(dd-dy*sn)*su;
            return 0;
        } else {
            s1[is]=pow( (dy/abs(ds)), b0)*(fc+fy)/(dc+dy);
            fm[is]=fs;
            dm[is]=ds;
            eu=b3*fm[is]/dm[is];
            x0=ds-fs/s1[is];
            if ( (dd-x0)*sn >=0) {
                ll=4;
                ss=s1[is];
                ff=fm[is]+(dd-dm[is])*s1[is];
                return 0;
            } else {
                // % loop of 430 %%
                is=3-is;
                sn=3-2*is;
                if ( abs(dm[is]) <= dc) {
                    d0=x0+sn*fc/s1[3-is];
                    f0=fc*sn;
                    if ( (d0-dd)*sn >0) {
                        ll=5;
                        ff=(dd-x0)*ss;
                        return 0;
                    } else {
                        // % loop of 530 %%
                        if ( ((dm[is]*sn <=dc) && (abs(dm[3-is]) <=dy))) {
                            std::tuple<int, float, float> llssff = this->Tslip_120(dy, dd, sn, fc, dc, sy, su, fy);
                            ll = std::get<0>(llssff);
                            ss = std::get<1>(llssff);
                            ff = std::get<2>(llssff);
                            return 0;
                        } else {
                            fm[is]=fy*sn;
                            dm[is]=dy*sn;
                            x0=dy*sn-fy*sn*(dy*sn-d0)/(fy*sn-f0);
                            if ( (dd-dm[is])*sn >=0) {
                                std::tuple<int, float, float> llssff = this->Tslip_120(dy, dd, sn, fc, dc, sy, su, fy);
                                ll = std::get<0>(llssff);
                                ss = std::get<1>(llssff);
                                ff = std::get<2>(llssff);
                                return 0;
                            } else {
                                ll=9;
                                ss=fm[is]/(dm[is]-x0);
                                ff=ss*(dd-x0);
                                return 0;
                            }
                        }
                    }
                } else if ( abs(dm[is]) <= dy) {
                    x=fm[is]/(dm[is]-x0);
                    y=fy*sn/(dy*sn-x0);
                    if (x < y) {
                        dm[is]=dy*sn;
                        fm[is]=fy*sn;
                    }
                    if ((dd-dm[is])*sn >=0) {
                        std::tuple<int, float, float> llssff = this->Tslip_120(dy, dd, sn, fc, dc, sy, su, fy);
                        ll = std::get<0>(llssff);
                        ss = std::get<1>(llssff);
                        ff = std::get<2>(llssff);
                        return 0;
                    } else {
                        ll=9;
                        ss=fm[is]/(dm[is]-x0);
                        ff=ss*(dd-x0);
                        return 0;
                    }
                } else {
                    // % loop of 440 %%
                    x=dm[3-is]-fm[3-is]/s1[3-is];
                    et=abs(fm[is]/(x-dm[is]))*pow( abs(dm[is]/(x-dm[is])), b2);
                    eu=fm[is]/dm[is]*b3;
                    if ( abs(et-eu) < err) {
                        xm=x0;
                    } else {
                        xm=(eu*dm[is]-et*x0-fm[is])/(eu-et);
                    }

                    if ( (dd-xm)*sn >0) {
                        if ((dd-dm[is])*sn >=0) {
                            std::tuple<int, float, float> llssff = this->Tslip_120(dy, dd, sn, fc, dc, sy, su, fy);
                            ll = std::get<0>(llssff);
                            ss = std::get<1>(llssff);
                            ff = std::get<2>(llssff);
                            return 0;
                        } else {
                            ll=7;
                            eu=fm[is]/dm[is]*b3;
                            ss=eu;
                            xf=dm[is]-fm[is]/eu;
                            x0=xf;
                            ff=eu*(dd-xf);
                            return 0;
                        }
                    } else {
                        ll=6;
                        x=dm[3-is]-fm[3-is]/s1[3-is];
                        et=abs(fm[is]/(x-dm[is]))*pow( abs(dm[is]/(x-dm[is])), b2);
                        ss=et;
                        ff=et*(dd-x0);
                        return 0;
                    }
                }
            }
        }
    }

    // % rule 4 unloading from peak (dm,fm) on the primary curve %%
    if ( ll==4) {
        if ( (dm[is]-dd)*sn <=0) {
            std::tuple<int, float, float> llssff = this->Tslip_120(dy, dd, sn, fc, dc, sy, su, fy);
            ll = std::get<0>(llssff);
            ss = std::get<1>(llssff);
            ff = std::get<2>(llssff);
            return 0;
        } else {
            if ( (dd-x0)*sn >0) {
                ff=fm[is]+(dd-dm[is])*s1[is];
                return 0;
            } else {
                // % loop of 430 %%
                is=3-is;
                sn=3-2*is;
                if ( abs(dm[is]) <= dc) {
                    d0=x0+sn*fc/s1[3-is];
                    f0=fc*sn;
                    if ( (d0-dd)*sn >0) {
                        ll=5;
                        ff=(dd-x0)*ss;
                        return 0;
                    } else {
                        // % loop of 530 %%
                        if ( ((dm[is]*sn <=dc) && (abs(dm[3-is]) <=dy))) {
                            std::tuple<int, float, float> llssff = this->Tslip_120(dy, dd, sn, fc, dc, sy, su, fy);
                            ll = std::get<0>(llssff);
                            ss = std::get<1>(llssff);
                            ff = std::get<2>(llssff);
                            return 0;
                        } else {
                            fm[is]=fy*sn;
                            dm[is]=dy*sn;
                            x0=dy*sn-fy*sn*(dy*sn-d0)/(fy*sn-f0);
                            if ( (dd-dm[is])*sn >=0) {
                                std::tuple<int, float, float> llssff = this->Tslip_120(dy, dd, sn, fc, dc, sy, su, fy);
                                ll = std::get<0>(llssff);
                                ss = std::get<1>(llssff);
                                ff = std::get<2>(llssff);
                                return 0;
                            } else {
                                ll=9;
                                ss=fm[is]/(dm[is]-x0);
                                ff=ss*(dd-x0);
                                return 0;
                            }
                        }
                    }
                } else if ( abs(dm[is]) <= dy) {
                    x=fm[is]/(dm[is]-x0);
                    y=fy*sn/(dy*sn-x0);
                    if (x < y) {
                        dm[is]=dy*sn;
                        fm[is]=fy*sn;
                    }
                    if ((dd-dm[is])*sn >=0) {
                        std::tuple<int, float, float> llssff = this->Tslip_120(dy, dd, sn, fc, dc, sy, su, fy);
                        ll = std::get<0>(llssff);
                        ss = std::get<1>(llssff);
                        ff = std::get<2>(llssff);
                        return 0;
                    } else {
                        ll=9;
                        ss=fm[is]/(dm[is]-x0);
                        ff=ss*(dd-x0);
                        return 0;
                    }
                } else {
                    // % loop of 440 %%
                    x=dm[3-is]-fm[3-is]/s1[3-is];
                    et=abs(fm[is]/(x-dm[is]))*pow( abs(dm[is]/(x-dm[is])), b2);
                    eu=fm[is]/dm[is]*b3;
                    if ( abs(et-eu) < err) {
                        xm=x0;
                    } else {
                        xm=(eu*dm[is]-et*x0-fm[is])/(eu-et);
                    }

                    if ( (dd-xm)*sn >0) {
                        if ((dd-dm[is])*sn >=0) {
                            std::tuple<int, float, float> llssff = this->Tslip_120(dy, dd, sn, fc, dc, sy, su, fy);
                            ll = std::get<0>(llssff);
                            ss = std::get<1>(llssff);
                            ff = std::get<2>(llssff);
                            return 0;
                        } else {
                            ll=7;
                            eu=fm[is]/dm[is]*b3;
                            ss=eu;
                            xf=dm[is]-fm[is]/eu;
                            x0=xf;
                            ff=eu*(dd-xf);
                            return 0;
                        }
                    } else {
                        ll=6;
                        x=dm[3-is]-fm[3-is]/s1[3-is];
                        et=abs(fm[is]/(x-dm[is]))*pow( abs(dm[is]/(x-dm[is])), b2);
                        ss=et;
                        ff=et*(dd-x0);
                        return 0;
                    }
                }
            }
        }
    }

    // % rule 5 load reversed at zeor-crossing paint (x0,0) withyout %%
    // % previous cracking in a new direction                        %%
    if (ll==5) {
        if ((dd-x0)*sn >0) {
            if ((d0-dd)*sn >0) {
                ff=(dd-x0)*ss;
                return 0;
            } else {
                // % loop of 530 %%
                if ( ((dm[is]*sn <=dc) && (abs(dm[3-is]) <=dy))) {
                    std::tuple<int, float, float> llssff = this->Tslip_120(dy, dd, sn, fc, dc, sy, su, fy);
                    ll = std::get<0>(llssff);
                    ss = std::get<1>(llssff);
                    ff = std::get<2>(llssff);
                    return 0;
                } else {
                    fm[is]=fc*sn;
                    dm[is]=dc*sn;
                    // %x0=dy*sn-fy*sn*(dy*sn-d0)/(fy*sn-f0);%%%%�o������ψʂƙ��f�͂𐳕��̂ǂ�����Ђъ���_�𒴂���܂łЂъ���_���L������悤�ɂ���
                    if ( (dd-dm[is])*sn >=0) {
                        std::tuple<int, float, float> llssff = this->Tslip_120(dy, dd, sn, fc, dc, sy, su, fy);
                        ll = std::get<0>(llssff);
                        ss = std::get<1>(llssff);
                        ff = std::get<2>(llssff);
                        return 0;
                    } else {
                        ll=9;
                        ss=fm[is]/(dm[is]-x0);
                        ff=ss*(dd-x0);
                        return 0;
                    }
                }
            }
        } else {
            is=3-is;
            sn=3-2*is;
            if ( (dm[is]-dd)*sn >0) {
                ll=4;
                ss=s1[is];
                ff=fm[is]+(dd-dm[is])*s1[is];
                return 0;
            } else {
                std::tuple<int, float, float> llssff = this->Tslip_120(dy, dd, sn, fc, dc, sy, su, fy);
                ll = std::get<0>(llssff);
                ss = std::get<1>(llssff);
                ff = std::get<2>(llssff);
                return 0;
            }
        }
    }

    // % rule 6 reloading with soft spring toward stiffness %%
    // % changing point xm                                  %%
    if ( ll==6) {
        if ( (dd-ds)*sn >0) {
            if ( (xm-dd)*sn >0) {
                x=dm[3-is]-fm[3-is]/s1[3-is];
                et=abs(fm[is]/(x-dm[is]))*pow( abs(dm[is]/(x-dm[is])), b2);
                ff=et*(dd-x0);
                return 0;
            } else {
                if ( (dd-dm[is])*sn >0) {
                    ll=3;
                    ss=su;
                    ff=fy*sn+(dd-dy*sn)*su;
                    return 0;
                } else {
                    ll=7;
                    eu=fm[is]/dm[is]*b3;
                    ss=eu;
                    xf=dm[is]-fm[is]/eu;
                    x0=xf;
                    ff=eu*(dd-xf);
                    return 0;
                }
            }
        } else {
            // % loop of 630 %%
            d1=ds;
            f1=fs;
            s2=s1[is]*b1;
            x1=d1-f1/s2;
            if ((x1-dd)*sn <0) {
                ll=8;
                ss=s2;
                ff=f1+(dd-d1)*s2;
                return 0;
            } else {
                // % loop of 840 %%
                x0=x1;
                is=3-is;
                sn=3-2*is;
                if ( abs(dm[is]) >=dy) {
                    eu=fm[is]/dm[is]*b3;
                    xf=dm[is]-fm[is]/eu;
                    xm=x0;
                    if ((xf-x0)*sn <=0) {
                        if ( (dd-dm[is])*sn >=0) {
                            std::tuple<int, float, float> llssff = this->Tslip_120(dy, dd, sn, fc, dc, sy, su, fy);
                            ll = std::get<0>(llssff);
                            ss = std::get<1>(llssff);
                            ff = std::get<2>(llssff);
                            return 0;
                        } else {
                            ll=9;
                            ss=fm[is]/(dm[is]-x0);
                            ff=ss*(dd-x0);
                            return 0;
                        }
                    } else {
                        // % loop of 440 %%
                        x=dm[3-is]-fm[3-is]/s1[3-is];
                        et=abs(fm[is]/(x-dm[is]))*pow( abs(dm[is]/(x-dm[is])), b2);
                        eu=fm[is]/dm[is]*b3;
                        if ( abs(et-eu) < err) {
                            xm=x0;
                        } else {
                            xm=(eu*dm[is]-et*x0-fm[is])/(eu-et);
                        }

                        if ( (dd-xm)*sn >0) {
                            if ((dd-dm[is])*sn >=0) {
                                std::tuple<int, float, float> llssff = this->Tslip_120(dy, dd, sn, fc, dc, sy, su, fy);
                                ll = std::get<0>(llssff);
                                ss = std::get<1>(llssff);
                                ff = std::get<2>(llssff);
                                return 0;
                            } else {
                                ll=7;
                                eu=fm[is]/dm[is]*b3;
                                ss=eu;
                                xf=dm[is]-fm[is]/eu;
                                x0=xf;
                                ff=eu*(dd-xf);
                                return 0;
                            }
                        } else {
                            ll=6;
                            x=dm[3-is]-fm[3-is]/s1[3-is];
                            et=abs(fm[is]/(x-dm[is]))*pow( abs(dm[is]/(x-dm[is])), b2);
                            ss=et;
                            ff=et*(dd-x0);
                            return 0;
                        }
                    }
                } else {
                    if ( (dd-dm[is])*sn >=0) {
                        std::tuple<int, float, float> llssff = this->Tslip_120(dy, dd, sn, fc, dc, sy, su, fy);
                        ll = std::get<0>(llssff);
                        ss = std::get<1>(llssff);
                        ff = std::get<2>(llssff);
                        return 0;
                    } else {
                        ll=9;
                        ss=fm[is]/(dm[is]-x0);
                        ff=ss*(dd-x0);
                        return 0;
                    }
                }
            }
        }
    }

    // % rule 7 loading with hard spring toward previous yield point %%
    if (ll==7) {
        if ((dd-ds)*sn >0) {
            if ( (dd-dm[is])*sn >=0) {
                ll=3;
                ss=su;
                ff=fy*sn+(dd-dy*sn)*su;
                return 0;
            } else {
                eu=fm[is]/dm[is]*b3;
                xf=dm[is]-fm[is]/eu;
                ff=eu*(dd-xf);
                return 0;
            }
        } else {
            // % loop of 630 %%
            d1=ds;
            f1=fs;
            s2=s1[is]*b1;
            x1=d1-f1/s2;
            if ((x1-dd)*sn <0) {
                ll=8;
                ss=s2;
                ff=f1+(dd-d1)*s2;
                return 0;
            } else {
                // % loop of 840 %%
                x0=x1;
                is=3-is;
                sn=3-2*is;
                if ( abs(dm[is]) >dy) {
                    eu=fm[is]/dm[is]*b3;
                    xf=dm[is]-fm[is]/eu;
                    xm=x0;
                    if ((xf-x0)*sn <=0) {
                        if ( (dd-dm[is])*sn >=0) {
                            std::tuple<int, float, float> llssff = this->Tslip_120(dy, dd, sn, fc, dc, sy, su, fy);
                            ll = std::get<0>(llssff);
                            ss = std::get<1>(llssff);
                            ff = std::get<2>(llssff);
                            return 0;
                        } else {
                            ll=9;
                            ss=fm[is]/(dm[is]-x0);
                            ff=ss*(dd-x0);
                            return 0;
                        }
                    } else {
                        // % loop of 440 %%
                        x=dm[3-is]-fm[3-is]/s1[3-is];
                        et=abs(fm[is]/(x-dm[is]))*pow( abs(dm[is]/(x-dm[is])), b2);
                        eu=fm[is]/dm[is]*b3;
                        if ( abs(et-eu) < err) {
                            xm=x0;
                        } else {
                            xm=(eu*dm[is]-et*x0-fm[is])/(eu-et);
                        }

                        if ( (dd-xm)*sn >0) {
                            if ((dd-dm[is])*sn >=0) {
                                std::tuple<int, float, float> llssff = this->Tslip_120(dy, dd, sn, fc, dc, sy, su, fy);
                                ll = std::get<0>(llssff);
                                ss = std::get<1>(llssff);
                                ff = std::get<2>(llssff);
                                return 0;
                            } else {
                                ll=7;
                                eu=fm[is]/dm[is]*b3;
                                ss=eu;
                                xf=dm[is]-fm[is]/eu;
                                x0=xf;
                                ff=eu*(dd-xf);
                                return 0;
                            }
                        } else {
                            ll=6;
                            x=dm[3-is]-fm[3-is]/s1[3-is];
                            et=abs(fm[is]/(x-dm[is]))*pow( abs(dm[is]/(x-dm[is])), b2);
                            ss=et;
                            ff=et*(dd-x0);
                            return 0;
                        }
                    }
                } else {
                    if ( (dd-dm[is])*sn >=0) {
                        std::tuple<int, float, float> llssff = this->Tslip_120(dy, dd, sn, fc, dc, sy, su, fy);
                        ll = std::get<0>(llssff);
                        ss = std::get<1>(llssff);
                        ff = std::get<2>(llssff);
                        return 0;
                    } else {
                        ll=9;
                        ss=fm[is]/(dm[is]-x0);
                        ff=ss*(dd-x0);
                        return 0;
                    }
                }
            }
        }
    }

    // % rule 8 unloading from inner peak (d1,f1) %%
    if ( ll==8) {
        if ((x1-dd)*sn >=0) {
            // % loop of 840 %%
            x0=x1;
            is=3-is;
            sn=3-2*is;
            if ( abs(dm[is]) >dy) {
                eu=fm[is]/dm[is]*b3;
                xf=dm[is]-fm[is]/eu;
                xm=x0;
                if ((xf-x0)*sn <=0) {
                    if ( (dd-dm[is])*sn >=0) {
                        std::tuple<int, float, float> llssff = this->Tslip_120(dy, dd, sn, fc, dc, sy, su, fy);
                ll = std::get<0>(llssff);
                ss = std::get<1>(llssff);
                ff = std::get<2>(llssff);
                        return 0;
                    } else {
                        ll=9;
                        ss=fm[is]/(dm[is]-x0);
                        ff=ss*(dd-x0);
                        return 0;
                    }
                } else {
                    // % loop of 440 %%
                    x=dm[3-is]-fm[3-is]/s1[3-is];
                    et=abs(fm[is]/(x-dm[is]))*pow( abs(dm[is]/(x-dm[is])), b2);
                    eu=fm[is]/dm[is]*b3;
                    if ( abs(et-eu) < err) {
                        xm=x0;
                    } else {
                        xm=(eu*dm[is]-et*x0-fm[is])/(eu-et);
                    }

                    if ( (dd-xm)*sn >0) {
                        if ((dd-dm[is])*sn >=0) {
                            std::tuple<int, float, float> llssff = this->Tslip_120(dy, dd, sn, fc, dc, sy, su, fy);
                            ll = std::get<0>(llssff);
                            ss = std::get<1>(llssff);
                            ff = std::get<2>(llssff);
                            return 0;
                        } else {
                            ll=7;
                            eu=fm[is]/dm[is]*b3;
                            ss=eu;
                            xf=dm[is]-fm[is]/eu;
                            x0=xf;
                            ff=eu*(dd-xf);
                            return 0;
                        }
                    } else {
                            ll=6;
                            x=dm[3-is]-fm[3-is]/s1[3-is];
                            et=abs(fm[is]/(x-dm[is]))*pow( abs(dm[is]/(x-dm[is])), b2);
                            ss=et;
                            ff=et*(dd-x0);
                            return 0;
                        }
                    }
                } else {
                    if ( (dd-dm[is])*sn >=0) {
                        std::tuple<int, float, float> llssff = this->Tslip_120(dy, dd, sn, fc, dc, sy, su, fy);
                        ll = std::get<0>(llssff);
                        ss = std::get<1>(llssff);
                        ff = std::get<2>(llssff);
                        return 0;
                    } else {
                        ll=9;
                        ss=fm[is]/(dm[is]-x0);
                        ff=ss*(dd-x0);
                        return 0;
                    }
                }
            } else {
                if ( (d1-dd)*sn >0) {
                    ff=f1+(dd-d1)*s2;
                    return 0;
                } else {
                    // % loop of 830 %%
                    if (abs(dm[is]) >dy) {
                        if ( (dd-xm)*sn <=0) {
                            ll=6;
                            x=dm[3-is]-fm[3-is]/s1[3-is];
                            et=abs(fm[is]/(x-dm[is]))*pow( abs(dm[is]/(x-dm[is])), b2);
                            ss=et;
                            ff=et*(dd-x0);
                            return 0;
                        } else {
                            if ( (dd-dm[is])*sn >=0) {
                                std::tuple<int, float, float> llssff = this->Tslip_120(dy, dd, sn, fc, dc, sy, su, fy);
                                ll = std::get<0>(llssff);
                                ss = std::get<1>(llssff);
                                ff = std::get<2>(llssff);
                                return 0;
                            } else {
                                ll=9;
                                ss=fm[is]/(dm[is]-x0);
                                ff=ss*(dd-x0);
                                return 0;
                            }
                        }
                    } else {
                        if ( (dd-dm[is])*sn >=0) {
                                std::tuple<int, float, float> llssff = this->Tslip_120(dy, dd, sn, fc, dc, sy, su, fy);
                                ll = std::get<0>(llssff);
                                ss = std::get<1>(llssff);
                                ff = std::get<2>(llssff);
                                return 0;
                        } else {
                                ll=9;
                                ss=fm[is]/(dm[is]-x0);
                                ff=ss*(dd-x0);
                                return 0;
                            }
                    }
                }
        }
    }

    // % rule 9 reloading toward peak (d2,f2) directly %%
    if (ll==9) {
    if ( (dd-ds)*sn >0) {
            if ( (dd-dm[is])*sn >0) {
                std::tuple<int, float, float> llssff = this->Tslip_120(dy, dd, sn, fc, dc, sy, su, fy);
                ll = std::get<0>(llssff);
                ss = std::get<1>(llssff);
                ff = std::get<2>(llssff);
                return 0;
            } else {
                ff=ss*(dd-x0);
                return 0;
            }
    } else {
            d1=ds;
            f1=fs;
            s2=s1[is]*b1;
            x1=d1-f1/s2;
            if ( (x1-dd)*sn >=0) {
                // % loop of 840 %%
                x0=x1;
                is=3-is;
                sn=3-2*is;
                if ( abs(dm[is]) >dy) {
                    eu=fm[is]/dm[is]*b3;
                    xf=dm[is]-fm[is]/eu;
                    xm=x0;
                    if ((xf-x0)*sn <=0) {
                        if ( (dd-dm[is])*sn >=0) {
                            std::tuple<int, float, float> llssff = this->Tslip_120(dy, dd, sn, fc, dc, sy, su, fy);
                            ll = std::get<0>(llssff);
                            ss = std::get<1>(llssff);
                            ff = std::get<2>(llssff);
                            return 0;
                        } else {
                            ll=9;
                            ss=fm[is]/(dm[is]-x0);
                            ff=ss*(dd-x0);
                            return 0;
                        }
                    } else {
                        // % loop of 440 %%
                        x=dm[3-is]-fm[3-is]/s1[3-is];
                        et=abs(fm[is]/(x-dm[is]))*pow( abs(dm[is]/(x-dm[is])), b2);
                        eu=fm[is]/dm[is]*b3;
                        if ( abs(et-eu) < err) {
                            xm=x0;
                        } else {
                            xm=(eu*dm[is]-et*x0-fm[is])/(eu-et);
                        }

                        if ( (dd-xm)*sn >0) {
                            if ((dd-dm[is])*sn >=0) {
                                std::tuple<int, float, float> llssff = this->Tslip_120(dy, dd, sn, fc, dc, sy, su, fy);
                                ll = std::get<0>(llssff);
                                ss = std::get<1>(llssff);
                                ff = std::get<2>(llssff);
                                return 0;
                            } else {
                                ll=7;
                                eu=fm[is]/dm[is]*b3;
                                ss=eu;
                                xf=dm[is]-fm[is]/eu;
                                x0=xf;
                                ff=eu*(dd-xf);
                                return 0;
                            }
                        } else {
                            ll=6;
                            x=dm[3-is]-fm[3-is]/s1[3-is];
                            et=abs(fm[is]/(x-dm[is]))*pow( abs(dm[is]/(x-dm[is])), b2);
                            ss=et;
                            ff=et*(dd-x0);
                            return 0;
                        }
                    }
                } else {
                    if ( (dd-dm[is])*sn >=0) {
                        std::tuple<int, float, float> llssff = this->Tslip_120(dy, dd, sn, fc, dc, sy, su, fy);
                        ll = std::get<0>(llssff);
                        ss = std::get<1>(llssff);
                        ff = std::get<2>(llssff);
                        return 0;
                    } else {
                        ll=9;
                        ss=fm[is]/(dm[is]-x0);
                        ff=ss*(dd-x0);
                        return 0;
                    }
                }
            } else {
                ll=8;
                ss=s2;
                ff=f1+(dd-d1)*s2;
                return 0;
            }
        }
    }
}

double TakedaSlip::getStress(void)
{
    //cout << " getStress" << endln;
    return (ff);
}

double TakedaSlip::getTangent(void)
{
    //cout << " getTangent" << endln;
    return (ss);
}

double TakedaSlip::getInitialTangent(void)
{
    //cout << " getInitialTangent" << endln;
    return (sc);
}

double TakedaSlip::getStrain(void)
{
    //cout << " getStrain" << endln;
    return (dd);
}

int TakedaSlip::commitState(void)
{
    cll = ll;
    cds = ds;
    cdd = dd;
    cfs = fs;
    cff = ff;
    css = ss;
    cf0 = f0;
    cf1 = f1;
    cd0 = d0;
    cd1 = d1;
    cfm = fm;
    cdm = dm;
    cxm = xm;
    cs1 = s1;
    cs2 = s2;
    return 0;
}

int TakedaSlip::revertToLastCommit(void)
{
    ll = cll;
    ds = cds;
    dd = cdd;
    fs = cfs;
    ff = cff;
    ss = css;
    f0 = cf0;
    f1 = cf1;
    d0 = cd0;
    d1 = cd1;
    fm = cfm;
    dm = cdm;
    xm = cxm;
    s1 = cs1;
    s2 = cs2;
    return 0;
}

int TakedaSlip::revertToStart(void)
{
    ll = cll = 1;
    ss = css = sc;
    return 0;
}

UniaxialMaterial *
TakedaSlip::getCopy(void)
{
    TakedaSlip *theCopy = new TakedaSlip(
        this->getTag(),
        dc, dy,
        sc, sy, su,
        b0, b1);
    theCopy->ll = ll;
    theCopy->ds = ds;
    theCopy->dd = dd;
    theCopy->fs = fs;
    theCopy->ff = ff;
    theCopy->ss = ss;
    theCopy->f0 = f0;
    theCopy->f1 = f1;
    theCopy->d0 = d0;
    theCopy->d1 = d1;
    theCopy->fm = fm;
    theCopy->dm = dm;
    theCopy->xm = xm;
    theCopy->s1 = s1;
    theCopy->s2 = s2;

    theCopy->cll = cll;
    theCopy->cds = cds;
    theCopy->cdd = cdd;
    theCopy->cfs = cfs;
    theCopy->cff = cff;
    theCopy->css = css;
    theCopy->cf0 = cf0;
    theCopy->cf1 = cf1;
    theCopy->cd0 = cd0;
    theCopy->cd1 = cd1;
    theCopy->cfm = cfm;
    theCopy->cdm = cdm;
    theCopy->cxm = cxm;
    theCopy->cs1 = cs1;
    theCopy->cs2 = cs2;
    return theCopy;
}

int TakedaSlip::sendSelf(int cTag, Channel &theChannel)
{
    int res = 0;
    cout << " sendSelf" << endln;

    static Vector data(137);
    data(0) = this->getTag();
    data(11) = ll;
    data(12) = ds;
    data(13) = dd;
    data(14) = fs;
    data(15) = ff;
    data(16) = ss;
    data(17) = f0;
    data(18) = f1;
    data(19) = d0;
    data(20) = d1;
    data(21) = fm[1];
    data(22) = fm[2];
    data(23) = dm[1];
    data(24) = dm[2];
    data(25) = xm;
    data(26) = s1[1];
    data(27) = s1[2];
    data(28) = s2;

    data(111) = cll;
    data(112) = cds;
    data(113) = cdd;
    data(114) = cfs;
    data(115) = cff;
    data(116) = css;
    data(117) = cf0;
    data(118) = cf1;
    data(119) = cd0;
    data(120) = cd1;
    data(121) = cfm[1];
    data(122) = cfm[2];
    data(123) = cdm[1];
    data(124) = cdm[2];
    data(125) = cxm;
    data(126) = cs1[1];
    data(127) = cs1[2];
    data(128) = cs2;
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
        ll = data(11);
        ds = data(12);
        dd = data(13);
        fs = data(14);
        ff = data(15);
        ss = data(16);
        f0 = data(17);
        f1 = data(18);
        d0 = data(19);
        d1 = data(20);
        fm[1] = data(21);
        fm[2] = data(22);
        dm[1] = data(23);
        dm[2] = data(24);
        xm = data(25);
        s1[1] = data(26);
        s1[2] = data(27);
        s2 = data(28);

        cll = data(111);
        cds = data(112);
        cdd = data(113);
        cfs = data(114);
        cff = data(115);
        css = data(116);
        cf0 = data(117);
        cf1 = data(118);
        cd0 = data(119);
        cd1 = data(120);
        cfm[1] = data(121);
        cfm[2] = data(122);
        cdm[1] = data(123);
        cdm[2] = data(124);
        cxm = data(125);
        cs1[1] = data(126);
        cs1[2] = data(127);
        cs2 = data(128);
    }

    return res;
}

void TakedaSlip::Print(OPS_Stream &s, int flag)
{
    cout << "TakedaSlip tag: " << this->getTag() << endln;
}