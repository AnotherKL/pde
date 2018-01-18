/*
*-d/dx(p*du/dx)+ru=f
*u(0)=0;u(1)=0;0<=x<=1
*u(x)=sin(pi*x)
*p(x)=e^x
*r(x)=1
*-d/dx(p*du/dx)+ru=PI*e^x*(PI*sin(PI*x)-cos(PI*x))+r*sin(PI*x)
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "nsde.h"

//f(x)的值
double fi(double x);
//u(xi)精确解
double ui(double xi);
//p(x)
double px(double xi);
//u(x)
double rx(double xi);
//打印矩阵
//void prtmatrix();

int main(void) {

    //初始化
    double a[N];    //系数矩阵左下非零元集合，从左到右，从上到下
    double b[N];    //系数矩阵对角线元素集合,从左到右，从上到下
    double c[N];    //系数矩阵右上非零元集合，从左到右，从上到下
    double F[N];    //方程组常数项集合
    double x[N];    //网格节点
    double xx[N];   //网格节点(对偶点) x[1]为x[1+1/2],x[2-1/2]为x[1+1/2]即x[1]
    double e[N];    //误差
    double maxe;    //最大误差
    double nom;     //2范数 norm of matrice
    double ux[N];   //U(xi)近似解
    double h;       //网格步长
    int n;          //n等分
    int r;          //系数矩阵阶数
    double *p_maxe; //指向 maxe
    double *p_nom;  //指向 nom

    maxe = 0;
    nom = 0;
    p_maxe = &maxe;
    p_nom = &nom;

    //输入：n等分
    printf("Please enter an integer for n (n >= 3): \n");
    scanf("%d", &n);
    if (n >= 3) {

    //网格剖分
    h = 1.0 / n;
    for (int ix = 0; ix <= n ; ix++) {
        x[ix] = ix * h;
        xx[ix] = (ix + 1.0 / 2) * h;
    }

    //计算矩阵的秩
    r = n - 1;

    //Ax=F
    //构造系数矩阵A以及F
    for (int ib = 0; ib < r; ib++) {
        b[ib] = px(xx[ib])*(1/(h*h)) + px(xx[ib+1])*(1/(h*h)) + rx(x[ib+1]);
    }
    for (int ia = 0; ia < (r-1); ia++) {
        a[ia] = - px(xx[ia+1]) * (1 / (h * h));
    }
    for (int ic = 0; ic < (r-1); ic++) {
        c[ic] = - px(xx[ic+1]) * (1 / (h * h));
    }
    F[0] = fi(x[1]) - (- px(xx[0]) * (1 / (h * h)));
    F[r-1] = fi(x[r]) - (- px(xx[r]) * (1 / (h * h)));
    for (int ir = 1; ir < (r-1); ir++) {
        F[ir] = fi(x[ir+1]);
    }

    //求解三对角矩阵
    solvematrix(a, b, c, F, ux, r);

    //计算误差
    err(x, ux, e, p_maxe, p_nom, r);

    //打印结果
    res(ux, e, maxe, nom, r);
    }
    else {
        printf("Sorry, n is an illegal value.\n");
        exit(1);
    }

    return 0;
}

double ui(double x) {

    double uxi;

    uxi = sin(PI * x);

    return uxi;
}

double px(double xi) {

    double pxi;

    pxi = exp(xi);

    return pxi;
}

double rx(double xi) {

    double rxi;

    rxi = 1.0;

    return rxi;
}

double fi(double x) {

    double fx;

    fx = PI*exp(x)*(PI*sin(PI*x)-cos(PI*x))+rx(x)*sin(PI*x);

    return fx;
}
