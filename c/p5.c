/*******************************************************************************
-d/dx(du/dx)-d/dy(du/dy)+a(du/dx)+b(du/dy)+cu=f
u|D = g(x,y); D = (0,1)x(0,1)
u(x)=x*(1-x)*y*(1-y)
f=-(-2y+2y^2)-(-2x+2x^2)+a(y-y^2-2xy+2xy^2)+b(x-2xy-x^2+2yx^2)+c(x(1-x)y(1-y))
*******************************************************************************/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "nsde.h"

//f(x)的值
double fi(double x, double y, double xsa, double xsb, double xsc);
//u(xi)精确解
double uii(double x, double y);
//误差分析
void err2(double x[],double y[],double ux[],double e[],double *maxe,double *nom,int r,int ot);
//ui(不使用) 防止编译nsde.c时出错
double ui(double x);
void ress(double ux[], int r);

int main(void) {

    //初始化
    double x[N];    //网格节点
    double y[N];    //网格节点：y轴方向
    double e[N];    //误差
    double maxe;    //最大误差
    double nom;     //2范数 norm of matrice
    double h;       //网格步长
    int n;          //n等分
    int r;          //系数矩阵阶数
    double *p_maxe; //指向 maxe
    double *p_nom;  //指向 nom
    int ot;         //n-1

    int ta1;
    int tf;

    double xsa,xsb,xsc; //a,b,c
    double mk1,mk2,mk3,mk4,mk5;//5点系数

    ta1 = 0;
    tf = 0;

    maxe = 0;
    nom = 0;
    p_maxe = &maxe;
    p_nom = &nom;

    //输入：n等分
    printf("Please enter an integer for n (n in [4,31]): \n");
    scanf("%d", &n);
    if (n >= 4 && n <= 31) {

    //录入a,b,c
    printf("a,b,c: \n");
    scanf("%lf,%lf,%lf\n", &xsa, &xsb, &xsc);

    printf("000\n");

    //网格剖分
    h = 1.0 / n;
    for (int ix = 0; ix <= n ; ix++) {
        x[ix] = ix * h;
        y[ix] = ix * h;
    }

    //计算矩阵的秩
    r = (n-1)*(n-1);
    ot = n-1;

    printf("001\n");
    printf("001\n");
    

    //方程组常数项集合
    double F[r];
    printf("0011\n");

    //U(xi)近似解
    double ux[r];
    for (int tq1 = 0; tq1 < r; tq1++) {
        ux[tq1] = 0.0;
    }

    printf("0012\n");

    //系数矩阵A
    double ax[r][1000];
    for (int iaxi = 0; iaxi < r; iaxi++) {
        for (int ixaj = 0; ixaj < r; ixaj++) {
            ax[iaxi][ixaj] = 0.0;
        }
    }

    printf("002\n");

    //计算5点系数
    mk1 = -(1/(h*h))-(xsb/(2*h));
    mk2 = -(1/(h*h))-(xsa/(2*h));
    mk3 =  (4/(h*h))+xsc;
    mk4 = -(1/(h*h))+(xsa/(2*h));
    mk5 = -(1/(h*h))+(xsb/(2*h));

    printf("003\n");

    //赋值系数矩阵A(Ax=b)
    for (int ii3 = 0; ii3 < r; ii3++) {
        ax[ii3][ii3] = mk3;
    }
    for (int ii1 = ot; ii1 < r; ii1++) {
        ax[ii1][ta1] = mk1;
        ta1++;
    }
    for (int ii5 = 0; ii5 < (r-ot); ii5++) {
        ax[ii5][ot+ii5] = mk5;
    }
    for (int ii2 = 1; ii2 < r; ii2++) {
        ax[ii2][ii2-1] = mk2;
    }
    for (int ij2 = ot; ij2 < r; ij2 += ot) {
        ax[ij2][ij2-1] = 0.0;
    }
    for (int ii4 = 0; ii4 < (r-1); ii4++) {
        ax[ii4][ii4+1] = mk4;
    }
    for (int ij4 = ot; ij4 < r; ij4 += ot) {
        ax[ij4-1][ij4] = 0.0;
    }

    printf("004\n");

    //赋值系数矩阵A(Ax=b)
    for (int ibfj = 0; ibfj < ot; ibfj++) {
        for (int ibfi = 0; ibfi < ot; ibfi++) {
            F[ibfi+tf] = fi(x[ibfi+1],y[ibfj+1],xsa,xsb,xsc);
        }
        tf += ot;
    }

    printf("005\n");

    //高斯-赛德尔迭代法求解线性方程组
    gs(ax, ux, F, r);

    printf("006\n");

    //计算误差
    err2(x, y, ux, e, p_maxe, p_nom, r, ot);

    printf("007\n");

    //打印结果
    res(ux, e, maxe, nom, r);
    ress(ux, r);


    printf("008\n");

    }
    else {
        printf("Sorry, n is an illegal value.\n");
        exit(1);
    }
    return 0;
}

double uii(double x, double y) {

    double uxi;

    uxi = x*(1-x)*y*(1-y);

    return uxi;
}

double fi(double x, double y, double xsa, double xsb, double xsc) {

    double fx;
    double i1, i2, i3, i4, i5;

    i1 = -(-2*y+2*y*y);
    i2 = -(-2*x+2*x*x);
    i3 = xsa*(y-y*y-2*x*y+2*x*y*y);
    i4 = xsb*(x-2*x*y-x*x+2*y*x*x);
    i5 = xsc*(x*(1-x)*y*(1-y));

    fx = i1 + i2 + i3 + i4 + i5;

    return fx;
}

void err2(double x[],double y[],double ux[],double e[],double *maxe,double *nom,int r,int ot){

    double sum = 0;
    int tp = 0;

    for (int i = 0; i < ot; i++) {
        for (int j = 0; j < ot; j++) {
            e[j+ot] = fabs(uii(x[j+1],y[i+1]) - ux[i]);
            if (e[i] > *maxe)
                *maxe = e[i];
        }
        tp += ot;
    }

    //2-范数
    //||u-un||=o(h^2)
    for (int k = 0; k < r; k++) {
        sum += e[k] * e[k];
    }
    *nom = sqrt(sum);
}

double ui(double x) {
    double uix;
    uix = x;
    return uix;
}

void ress(double ux[], int r) {

    FILE *fp;
    int i;
    if ((fp = fopen("matres.txt", "a")) == NULL)
    {
        printf("Connot open this file.");
        exit(3);
    }
    fprintf(fp, "n=%d时的数值解：\n", sqrt(r)+1);
    for (i = 0; i < r; i++)
    {
        fprintf(fp, "%f\n", ux[i]);
    }
    fclose(fp);
}
