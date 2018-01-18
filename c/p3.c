/*
*-d/dx(p*du/dx)+ru=f
*u(0)=0;u(1)=0;0<=x<=1
*u(x)=sin(pi*x)/p(x)
*p(x)=1     x<1/2
*p(x)=1000  x>1/2
*r(x)=1
*-d/dx(p*du/dx)+ru=PI^2*sin(PI*x)+rx(x)*sin(PI*x)/px(x)
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "nsde.h"

//f(x)的值
double fi(double x);
//u(xi)精确解
double ui(double x);
//p(x)
double px(double xi);
//r(x)
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
    double dx;      //积分函数横坐标间隔数dx：越小越精确
    double maxe;    //最大误差
    double nom;     //2范数 norm of matrice
    double ux[N];   //U(xi)近似解
    double h;       //网格步长
    int n;          //n等分
    int r;          //系数矩阵阶数
    double *p_maxe; //指向 maxe
    double *p_nom;  //指向 nom

    dx = 0.000001;
    maxe = 0;
    nom = 0;
    p_maxe = &maxe;
    p_nom = &nom;

    //输入：n等分
    printf("Please enter an integer for n (n >= 3): \n");
    scanf("%d", &n);
    if (n >= 3) {

    //录入dx
    //printf("dx = \n");
    //scanf("%d", &dx);

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
        b[ib] = (1/((1/h)*integrate(px,x[ib],x[ib+1],dx)))*(1/(h*h));
        b[ib] += (1/((1/h)*integrate(px,x[ib+1],x[ib+2],dx)))*(1/(h*h));
        b[ib] += ((1/h)*integrate(rx,xx[ib],xx[ib+1],dx));
    }
    for (int ia = 0; ia < (r-1); ia++) {
        a[ia] = - (1/((1/h) * integrate(px,x[ia+1],x[ia+2],dx))) * (1/(h*h));
    }
    for (int ic = 0; ic < (r-1); ic++) {
        c[ic] = - (1/((1/h) * integrate(px,x[ic+1],x[ic+2],dx))) * (1/(h*h));
    }
    F[0] = (1/h)*integrate(fi,xx[0],xx[1],dx);
    F[0] -= -(1/((1/h)*integrate(px,x[0],x[1],dx)))*(1/(h*h));
    F[r-1] = (1/h)*integrate(fi,xx[r-1],xx[r],dx);
    F[r-1] -= -(1/((1/h)*integrate(px,x[r],x[r+1],dx))) * (1/(h*h));
    for (int ir = 1; ir < (r-1); ir++) {
        F[ir] = (1/h)*integrate(fi,xx[ir],xx[ir+1],dx);
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

    uxi = sin(PI*x) / px(x);

    return uxi;
}

double px(double xi) {

    if (xi < (1.0/2)) {
        return 1.0;
    }
    else {
        return 1000.0;
    }

}

double rx(double xi) {

    double rxi;

    rxi = 1.0;

    return rxi;
}

double fi(double x) {

    double fx;

    fx = PI*PI*sin(PI*x)+rx(x)*sin(PI*x)/px(x);

    return fx;
}
