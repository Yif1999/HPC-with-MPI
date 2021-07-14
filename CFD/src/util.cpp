#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "util.h"

/*输入进程数返回分块维度结果数组*/
void Block_Divide(int n,int *a){
    float ref=sqrt(n);
    int refINT=int(ceil(ref));
    for (int i=refINT;i>0;i--){
        if (1.0*n/i==int(1.0*n/i)){
            if (i<n/i){
                a[0]=i;
                a[1]=n/i;
            }else{
                a[1]=i;
                a[0]=n/i;
            }
            return;
        }
    }
}

/*x轴向的SW分裂并返回更新数据的结构u*/
unit Steger_Warming_X(unit u){  
    float lmd[4],lmdp[4],lmdn[4]; 
    float ratio=u.param.rho/(2.0*gama) ;
    lmd[0]=u.param.vel.u;
    lmd[1]=lmd[0];
    lmd[2]=lmd[0]-u.param.c;
    lmd[3]=lmd[0]+u.param.c;

    for (int k=0;k<4;k++){
        lmdp[k]=0.5*(lmd[k]+sqrt(lmd[k]*lmd[k]+eps*eps));
        lmdn[k]=0.5*(lmd[k]-sqrt(lmd[k]*lmd[k]+eps*eps));
    }
    u.param.fp[0] = ratio* (2.0 * (gama- 1.0) * lmdp[0] + lmdp[2] + lmdp[3]);
    u.param.fp[1] = ratio* (2.0 * (gama- 1.0) *u.param.vel.u* lmdp[0] + (u.param.vel.u - u.param.c) * lmdp[2] + (u.param.vel.u + u.param.c) * lmdp[3]);
    u.param.fp[2] = ratio* (2.0 * (gama- 1.0) * u.param.vel.v * lmdp[0] +u.param.vel.v * lmdp[2] + u.param.vel.v * lmdp[3]);
    u.param.fp[3] = ratio* ((gama- 1.0) * (u.param.vel.u * u.param.vel.u + u.param.vel.v * u.param.vel.v) * lmdp[0] + 0.5 * ((u.param.vel.u - u.param.c) * (u.param.vel.u - u.param.c) +u.param.vel.v * u.param.vel.v) * lmdp[2] + 0.5 * ((u.param.vel.u + u.param.c) * (u.param.vel.u + u.param.c) + u.param.vel.v * u.param.vel.v) * lmdp[3] + ((3.0 - gama) / (2.0 * (gama- 1)) * u.param.c * u.param.c * (lmdp[2] + lmdp[3])));
    u.param.fn[0] = ratio* (2.0 * (gama- 1.0) * lmdn[0] + lmdn[2] + lmdn[3]);
    u.param.fn[1] = ratio* (2.0 * (gama- 1.0) *u.param.vel.u* lmdn[0] + (u.param.vel.u - u.param.c) * lmdn[2] + (u.param.vel.u + u.param.c) * lmdn[3]);
    u.param.fn[2] = ratio* (2.0 * (gama- 1.0) * u.param.vel.v * lmdn[0] +u.param.vel.v * lmdn[2] + u.param.vel.v * lmdn[3]);
    u.param.fn[3] = ratio* ((gama- 1.0) * (u.param.vel.u * u.param.vel.u + u.param.vel.v * u.param.vel.v) * lmdn[0] + 0.5 * ((u.param.vel.u - u.param.c) * (u.param.vel.u - u.param.c) +u.param.vel.v * u.param.vel.v) * lmdn[2] + 0.5 * ((u.param.vel.u + u.param.c) * (u.param.vel.u + u.param.c) + u.param.vel.v * u.param.vel.v) * lmdn[3] + ((3.0 - gama) / (2.0 * (gama- 1)) * u.param.c * u.param.c * (lmdn[2] + lmdn[3])));

    return u;
}

/*y轴向的SW分裂并返回更新数据的结构u*/
unit Steger_Warming_Y(unit u){  
    float miu[4],miup[4],miun[4]; 
    float ratio=u.param.rho/(2.0*gama) ;
    miu[0]=u.param.vel.v;
    miu[1]=miu[0];
    miu[2]=miu[0]-u.param.c;
    miu[3]=miu[0]+u.param.c;
    for (int k=0;k<4;k++){
        miup[k]=0.5*(miu[k]+sqrt(miu[k]*miu[k]+eps*eps));
        miun[k]=0.5*(miu[k]-sqrt(miu[k]*miu[k]+eps*eps));
    }
    u.param.gp[0] = ratio * (2.0 * (gama- 1.0) * miup[0] + miup[2] + miup[3]);
    u.param.gp[1] = ratio * (2.0 * (gama- 1.0) *u.param.vel.u* miup[0] + u.param.vel.u * miup[2] + u.param.vel.u * miup[3]);
    u.param.gp[2] = ratio * (2.0 * (gama- 1.0) * u.param.vel.v * miup[0] +(u.param.vel.v-u.param.c)* miup[2] + (u.param.vel.v+u.param.c)* miup[3]);
    u.param.gp[3] = ratio * ((gama- 1.0) * (u.param.vel.u * u.param.vel.u + u.param.vel.v * u.param.vel.v) * miup[0] + 0.5 * ((u.param.vel.v - u.param.c) * (u.param.vel.v - u.param.c) +u.param.vel.u * u.param.vel.u) * miup[2] + 0.5 * ((u.param.vel.v + u.param.c) * (u.param.vel.v + u.param.c) + u.param.vel.u * u.param.vel.u) * miup[3] + ((3.0 - gama) / (2.0 * (gama- 1)) * u.param.c * u.param.c * (miup[2] + miup[3])));
    u.param.gn[0] = ratio * (2.0 * (gama- 1.0) * miun[0] + miun[2] + miun[3]);
    u.param.gn[1] = ratio * (2.0 * (gama- 1.0) *u.param.vel.u* miun[0] + u.param.vel.u* miun[2] + u.param.vel.u * miun[3]);
    u.param.gn[2] = ratio * (2.0 * (gama- 1.0) * u.param.vel.v * miun[0] +(u.param.vel.v-u.param.c) * miun[2] + (u.param.vel.v+u.param.c) * miun[3]);
    u.param.gn[3] = ratio * ((gama- 1.0) * (u.param.vel.u * u.param.vel.u + u.param.vel.v * u.param.vel.v) * miun[0] + 0.5 * ((u.param.vel.v - u.param.c) * (u.param.vel.v- u.param.c) +u.param.vel.u * u.param.vel.u) * miun[2] + 0.5 * ((u.param.vel.v + u.param.c) * (u.param.vel.v + u.param.c) + u.param.vel.u * u.param.vel.u) * miun[3] + ((3.0 - gama) / (2.0 * (gama- 1)) * u.param.c * u.param.c * (miun[2] + miun[3])));

    return u;
}

/*x轴向的WENO差分计算*/
void WENO_X(float (*fp)[4],float (*fn)[4],float *fx){
    int k;
    float fp_WENO[4],fn_WENO[4];
    float f1p[4],f2p[4],f3p[4],f1n[4],f2n[4],f3n[4];
    float C1=0.1,C2=0.6,C3=0.3;
    float omega1[4],omega2[4],omega3[4];
    float alpha1[4],alpha2[4],alpha3[4];
    float IS1[4],IS2[4],IS3[4];

    //正通量计算
    for (k=0;k<4;k++){
        IS1[k]=0.25*pow((fp[1][k]-4.0*fp[2][k]+3.0*fp[3][k]),2.0)+13.0/12.0*pow((fp[1][k]-2.0*fp[2][k]+fp[3][k]),2.0);
        IS2[k]=0.25*pow((fp[2][k]-fp[4][k]),2.0)+13.0/12.0*pow((fp[2][k]-2.0*fp[3][k]+fp[4][k]),2.0);
        IS3[k]=0.25*pow((3.0*fp[3][k]-4.0*fp[4][k]+fp[5][k]),2.0)+13.0/12.0*pow((fp[3][k]-2.0*fp[4][k]+fp[5][k]),2.0);
        alpha1[k]=C1*1.0/pow((eps+IS1[k]),2.0);
        alpha2[k]=C2*1.0/pow((eps+IS2[k]),2.0);
        alpha3[k]=C3*1.0/pow((eps+IS3[k]),2.0);
        omega1[k]=alpha1[k]/(alpha1[k]+alpha2[k]+alpha3[k]);
        omega2[k]=alpha2[k]/(alpha1[k]+alpha2[k]+alpha3[k]);
        omega3[k]=alpha3[k]/(alpha1[k]+alpha2[k]+alpha3[k]);
        f1p[k]=1.0/3*fp[1][k]-7.0/6*fp[2][k]+11.0/6*fp[3][k];
        f2p[k]=-1.0/6*fp[2][k]+5.0/6*fp[3][k]+1.0/3*fp[4][k];
        f3p[k]=1.0/3*fp[3][k]+5.0/6*fp[4][k]-1.0/6*fp[5][k];
        fp_WENO[k]=omega1[k]*f1p[k]+omega2[k]*f2p[k]+omega3[k]*f3p[k];

        IS1[k]=0.25*pow((fp[0][k]-4.0*fp[1][k]+3.0*fp[2][k]),2.0)+13.0/12.0*pow((fp[0][k]-2.0*fp[1][k]+fp[2][k]),2.0);
        IS2[k]=0.25*pow((fp[1][k]-fp[3][k]),2.0)+13.0/12.0*pow((fp[1][k]-2.0*fp[2][k]+fp[3][k]),2.0);
        IS3[k]=0.25*pow((3.0*fp[2][k]-4.0*fp[3][k]+fp[4][k]),2.0)+13.0/12.0*pow((fp[2][k]-2.0*fp[3][k]+fp[4][k]),2.0);
        alpha1[k]=C1*1.0/pow((eps+IS1[k]),2.0);
        alpha2[k]=C2*1.0/pow((eps+IS2[k]),2.0);
        alpha3[k]=C3*1.0/pow((eps+IS3[k]),2.0);
        omega1[k]=alpha1[k]/(alpha1[k]+alpha2[k]+alpha3[k]);
        omega2[k]=alpha2[k]/(alpha1[k]+alpha2[k]+alpha3[k]);
        omega3[k]=alpha3[k]/(alpha1[k]+alpha2[k]+alpha3[k]);
        f1n[k]=1.0/3*fp[0][k]-7.0/6*fp[1][k]+11.0/6*fp[2][k];
        f2n[k]=-1.0/6*fp[1][k]+5.0/6*fp[2][k]+1.0/3*fp[3][k];
        f3n[k]=1.0/3*fp[2][k]+5.0/6*fp[3][k]-1.0/6*fp[4][k];     
        fn_WENO[k]=omega1[k]*f1n[k]+omega2[k]*f2n[k]+omega3[k]*f3n[k];
        
        fx[k]+=(fp_WENO[k]-fn_WENO[k])/dx;
    }
    //负通量计算
    for (k=0;k<4;k++){
        IS1[k]=0.25*pow((fn[5][k]-4.0*fn[4][k]+3.0*fn[3][k]),2.0)+13.0/12.0*pow((fn[5][k]-2.0*fn[4][k]+fn[3][k]),2.0);
        IS2[k]=0.25*pow((fn[4][k]-fn[2][k]),2.0)+13.0/12.0*pow((fn[4][k]-2.0*fn[3][k]+fn[2][k]),2.0);
        IS3[k]=0.25*pow((3.0*fn[3][k]-4.0*fn[2][k]+fn[1][k]),2.0)+13.0/12.0*pow((fn[3][k]-2.0*fn[2][k]+fn[1][k]),2.0);
        alpha1[k]=C1*1.0/pow((eps+IS1[k]),2.0);
        alpha2[k]=C2*1.0/pow((eps+IS2[k]),2.0);
        alpha3[k]=C3*1.0/pow((eps+IS3[k]),2.0);
        omega1[k]=alpha1[k]/(alpha1[k]+alpha2[k]+alpha3[k]);
        omega2[k]=alpha2[k]/(alpha1[k]+alpha2[k]+alpha3[k]);
        omega3[k]=alpha3[k]/(alpha1[k]+alpha2[k]+alpha3[k]);
        f1n[k]=1.0/3*fn[5][k]-7.0/6*fn[4][k]+11.0/6*fn[3][k];
        f2n[k]=-1.0/6*fn[4][k]+5.0/6*fn[3][k]+1.0/3*fn[2][k];
        f3n[k]=1.0/3*fn[3][k]+5.0/6*fn[2][k]-1.0/6*fn[1][k];
        fn_WENO[k]=omega1[k]*f1n[k]+omega2[k]*f2n[k]+omega3[k]*f3n[k];

        IS1[k]=0.25*pow((fn[6][k]-4.0*fn[5][k]+3.0*fn[4][k]),2.0)+13.0/12.0*pow((fn[6][k]-2.0*fn[5][k]+fn[4][k]),2.0);
        IS2[k]=0.25*pow((fn[5][k]-fn[3][k]),2.0)+13.0/12.0*pow((fn[5][k]-2.0*fn[4][k]+fn[3][k]),2.0);
        IS3[k]=0.25*pow((3.0*fn[4][k]-4.0*fn[3][k]+fn[2][k]),2.0)+13.0/12.0*pow((fn[4][k]-2.0*fn[3][k]+fn[2][k]),2.0);
        alpha1[k]=C1*1.0/pow((eps+IS1[k]),2.0);
        alpha2[k]=C2*1.0/pow((eps+IS2[k]),2.0);
        alpha3[k]=C3*1.0/pow((eps+IS3[k]),2.0);
        omega1[k]=alpha1[k]/(alpha1[k]+alpha2[k]+alpha3[k]);
        omega2[k]=alpha2[k]/(alpha1[k]+alpha2[k]+alpha3[k]);
        omega3[k]=alpha3[k]/(alpha1[k]+alpha2[k]+alpha3[k]);
        f1p[k]=1.0/3*fn[6][k]-7.0/6*fn[5][k]+11.0/6*fn[4][k];
        f2p[k]=-1.0/6*fn[5][k]+5.0/6*fn[4][k]+1.0/3*fn[3][k];
        f3p[k]=1.0/3*fn[4][k]+5.0/6*fn[3][k]-1.0/6*fn[2][k];  
        fp_WENO[k]=omega1[k]*f1p[k]+omega2[k]*f2p[k]+omega3[k]*f3p[k];

        fx[k]+=(fp_WENO[k]-fn_WENO[k])/dx;
    }

    return ;
}

/*y轴向的WENO差分计算*/
void WENO_Y(float (*gp)[4],float (*gn)[4],float *gy){
    int k;
    float gp_WENO[4],gn_WENO[4];
    float g1p[4],g2p[4],g3p[4],g1n[4],g2n[4],g3n[4];
    float C1=0.1,C2=0.6,C3=0.3;
    float omega1[4],omega2[4],omega3[4];
    float alpha1[4],alpha2[4],alpha3[4];
    float IS1[4],IS2[4],IS3[4];

    //正通量计算
    for (k=0;k<4;k++){
        IS1[k]=0.25*pow((gp[1][k]-4.0*gp[2][k]+3.0*gp[3][k]),2.0)+13.0/12.0*pow((gp[1][k]-2.0*gp[2][k]+gp[3][k]),2.0);
        IS2[k]=0.25*pow((gp[2][k]-gp[4][k]),2.0)+13.0/12.0*pow((gp[2][k]-2.0*gp[3][k]+gp[4][k]),2.0);
        IS3[k]=0.25*pow((3.0*gp[3][k]-4.0*gp[4][k]+gp[5][k]),2.0)+13.0/12.0*pow((gp[3][k]-2.0*gp[4][k]+gp[5][k]),2.0);
        alpha1[k]=C1*1.0/pow((eps+IS1[k]),2.0);
        alpha2[k]=C2*1.0/pow((eps+IS2[k]),2.0);
        alpha3[k]=C3*1.0/pow((eps+IS3[k]),2.0);
        omega1[k]=alpha1[k]/(alpha1[k]+alpha2[k]+alpha3[k]);
        omega2[k]=alpha2[k]/(alpha1[k]+alpha2[k]+alpha3[k]);
        omega3[k]=alpha3[k]/(alpha1[k]+alpha2[k]+alpha3[k]);
        g1p[k]=1.0/3*gp[1][k]-7.0/6*gp[2][k]+11.0/6*gp[3][k];
        g2p[k]=-1.0/6*gp[2][k]+5.0/6*gp[3][k]+1.0/3*gp[4][k];
        g3p[k]=1.0/3*gp[3][k]+5.0/6*gp[4][k]-1.0/6*gp[5][k];
        gp_WENO[k]=omega1[k]*g1p[k]+omega2[k]*g2p[k]+omega3[k]*g3p[k];

        IS1[k]=0.25*pow((gp[0][k]-4.0*gp[1][k]+3.0*gp[2][k]),2.0)+13.0/12.0*pow((gp[0][k]-2.0*gp[1][k]+gp[2][k]),2.0);
        IS2[k]=0.25*pow((gp[1][k]-gp[3][k]),2.0)+13.0/12.0*pow((gp[1][k]-2.0*gp[2][k]+gp[3][k]),2.0);
        IS3[k]=0.25*pow((3.0*gp[2][k]-4.0*gp[3][k]+gp[4][k]),2.0)+13.0/12.0*pow((gp[2][k]-2.0*gp[3][k]+gp[4][k]),2.0);
        alpha1[k]=C1*1.0/pow((eps+IS1[k]),2.0);
        alpha2[k]=C2*1.0/pow((eps+IS2[k]),2.0);
        alpha3[k]=C3*1.0/pow((eps+IS3[k]),2.0);
        omega1[k]=alpha1[k]/(alpha1[k]+alpha2[k]+alpha3[k]);
        omega2[k]=alpha2[k]/(alpha1[k]+alpha2[k]+alpha3[k]);
        omega3[k]=alpha3[k]/(alpha1[k]+alpha2[k]+alpha3[k]);        
        g1n[k]=1.0/3*gp[0][k]-7.0/6*gp[1][k]+11.0/6*gp[2][k];
        g2n[k]=-1.0/6*gp[1][k]+5.0/6*gp[2][k]+1.0/3*gp[3][k];
        g3n[k]=1.0/3*gp[2][k]+5.0/6*gp[3][k]-1.0/6*gp[4][k];     
        gn_WENO[k]=omega1[k]*g1n[k]+omega2[k]*g2n[k]+omega3[k]*g3n[k];

        gy[k]+=(gp_WENO[k]-gn_WENO[k])/dy;
    }
    //负通量计算
    for (k=0;k<4;k++){
        IS1[k]=0.25*pow((gn[5][k]-4.0*gn[4][k]+3.0*gn[3][k]),2.0)+13.0/12.0*pow((gn[5][k]-2.0*gn[4][k]+gn[3][k]),2.0);
        IS2[k]=0.25*pow((gn[4][k]-gn[2][k]),2.0)+13.0/12.0*pow((gn[4][k]-2.0*gn[3][k]+gn[2][k]),2.0);
        IS3[k]=0.25*pow((3.0*gn[3][k]-4.0*gn[2][k]+gn[1][k]),2.0)+13.0/12.0*pow((gn[3][k]-2.0*gn[2][k]+gn[1][k]),2.0);
        alpha1[k]=C1*1.0/pow((eps+IS1[k]),2.0);
        alpha2[k]=C2*1.0/pow((eps+IS2[k]),2.0);
        alpha3[k]=C3*1.0/pow((eps+IS3[k]),2.0);
        omega1[k]=alpha1[k]/(alpha1[k]+alpha2[k]+alpha3[k]);
        omega2[k]=alpha2[k]/(alpha1[k]+alpha2[k]+alpha3[k]);
        omega3[k]=alpha3[k]/(alpha1[k]+alpha2[k]+alpha3[k]);
        g1n[k]=1.0/3*gn[5][k]-7.0/6*gn[4][k]+11.0/6*gn[3][k];
        g2n[k]=-1.0/6*gn[4][k]+5.0/6*gn[3][k]+1.0/3*gn[2][k];
        g3n[k]=1.0/3*gn[3][k]+5.0/6*gn[2][k]-1.0/6*gn[1][k];
        gn_WENO[k]=omega1[k]*g1n[k]+omega2[k]*g2n[k]+omega3[k]*g3n[k];

        IS1[k]=0.25*pow((gn[6][k]-4.0*gn[5][k]+3.0*gn[4][k]),2.0)+13.0/12.0*pow((gn[6][k]-2.0*gn[5][k]+gn[4][k]),2.0);
        IS2[k]=0.25*pow((gn[5][k]-gn[3][k]),2.0)+13.0/12.0*pow((gn[5][k]-2.0*gn[4][k]+gn[3][k]),2.0);
        IS3[k]=0.25*pow((3.0*gn[4][k]-4.0*gn[3][k]+gn[2][k]),2.0)+13.0/12.0*pow((gn[4][k]-2.0*gn[3][k]+gn[2][k]),2.0);
        alpha1[k]=C1*1.0/pow((eps+IS1[k]),2.0);
        alpha2[k]=C2*1.0/pow((eps+IS2[k]),2.0);
        alpha3[k]=C3*1.0/pow((eps+IS3[k]),2.0);
        omega1[k]=alpha1[k]/(alpha1[k]+alpha2[k]+alpha3[k]);
        omega2[k]=alpha2[k]/(alpha1[k]+alpha2[k]+alpha3[k]);
        omega3[k]=alpha3[k]/(alpha1[k]+alpha2[k]+alpha3[k]);
        g1p[k]=1.0/3*gn[6][k]-7.0/6*gn[5][k]+11.0/6*gn[4][k];
        g2p[k]=-1.0/6*gn[5][k]+5.0/6*gn[4][k]+1.0/3*gn[3][k];
        g3p[k]=1.0/3*gn[4][k]+5.0/6*gn[3][k]-1.0/6*gn[2][k];  
        gp_WENO[k]=omega1[k]*g1p[k]+omega2[k]*g2p[k]+omega3[k]*g3p[k];

        gy[k]+=(gp_WENO[k]-gn_WENO[k])/dy;
    }

    return ;
}

/*输入寻址i、j并返回对应一维结果数组索引地址*/
int Find_Adrs(int i,int j,int blockL, int blockH,int *dims){
    int index,overhead;
    overhead=(i/blockH*dims[1]+j/blockL)*blockH*blockL;
    index=(i%blockH)*blockL+j%blockL+overhead;
    return index;
}