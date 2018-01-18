/*
-u"+pu'+qu=f
u(0)=0;u(1)=0;0<=x<=1
u(x)=sin(pi*x)
*/
/*-u"(x)+u'(x)+u(x)=PI^2*sin(PI*x)+PI*cos(PI*x)+sin(PI*x)*/
/*-u"(x)+pu'(x)+qu(x)=PI^2*sin(PI*x)+p*PI*cos(PI*x)+q*sin(PI*x)*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "nsde.h"

//f(x)的值
double fi(double x, double p, double q);
//u(xi)精确解
double ui(double xi);
//打印矩阵
//void prtmatrix();

int main(void) {
    //初始化
    double a[N];    //系数矩阵左下非零元集合，从左到右，从上到下
    double b[N];    //系数矩阵对角线元素集合,从左到右，从上到下
    double c[N];    //系数矩阵右上非零元集合，从左到右，从上到下
    double F[N];    //方程组常数项集合
    double x[N];    //网格节点
    double e[N];    //误差
    double maxe;    //最大误差
	double nom;     //2范数 norm of matrice
    double p, q;    //-u"+pu'+qu=f
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

    //录入系数o,p,q
    printf("-u\"+pu\'+qu=f\np, q (eg: p q)= \n");
    scanf("%d %d", &p, &q);

    //网格剖分
    h = 1.0 / n;
    for (int ix = 0; ix <= n ; ix++)
        x[ix] = ix * h;

    //计算矩阵的秩
    r = n - 1;

    //Ax=F
    //构造系数矩阵A以及F
    for (int ia = 0; ia < (r-1); ia++) {
        a[ia] = - p / (2 * h) - 1 / (h * h);
        c[ia] = p / (2 * h) - 1 / (h * h);
    }
    for (int ib = 0; ib < r; ib++) {
        b[ib] = 2 / (h * h) + q;
        F[ib] = fi(x[ib+1], p, q);
    }

    //求解三对角矩阵
    solvematrix(a, b, c, F, ux, r);

    //计算误差
    err(x, ux, e, p_maxe, p_nom, r)

    //打印结果
    res(ux, e, maxe, nom, r);
    }
    else {
        printf("Sorry, n is an illegal value.\n");
        exit(1);
    }

    return 0;
}

double fi(double x, double p, double q) {

    double fx;

    fx = PI*PI*sin(PI*x)+p*PI*cos(PI*x)+q*sin(PI*x);

    return fx;
}

double ui(double x) {

    double uxi;

    uxi = sin(PI * x);

    return uxi;
}
