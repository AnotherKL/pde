/*******************************************************************************
nsde.c
偏微分方程有关的函数
*******************************************************************************/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "nsde.h"

//局部函数原型

//函数接口
//误差分析
void err(double x[], double ux[], double e[], double *maxe, double *nom, int r){

    double sum = 0;

    for (int i = 0; i < r; i++) {
        e[i] = fabs(ui(x[i+1]) - ux[i]);
        if (e[i] > *maxe)
            *maxe = e[i];
        sum += e[i] * e[i];
    }

    //2-范数
    //||u-un||=o(h^2)
    *nom = sqrt(sum);
}

//追赶法求解三对角矩阵
void solvematrix(double a[],double b[],double c[],double F[],double ux[],int r){

        double beta[N], gamma[N];

        //beta递推
        beta[0] = c[0] / b[0];
        for (int i = 1; i < (r - 1); i++)
            beta[i] = c[i] / (b[i] - a[i] * beta[i-1]);

        //求Ly=f
        gamma[0] = F[0] / b[0];
        for (int j = 1; j < r; j++) {
            gamma[j] = (F[j] - a[j] * gamma[j-1]) / (b[j] - a[j] * beta[j-1]);
        }

        //求Ux=y
        ux[r-1] = gamma[r-1];
        for (int k = (r - 2); k >= 0; k--)
            ux[k] = gamma[k] - beta[k] * ux[k+1];
}

//记录结果
void res(double ux[], double e[], double maxe, double nom, int r) {

    FILE *fp;
    int i;
    if ((fp = fopen("Result.txt", "a")) == NULL)
    {
        printf("Connot open this file.");
        exit(3);
    }
    fprintf(fp, "n=%d时的数值解：\n", r + 1);
    for (i = 0; i < r; i++)
    {
        fprintf(fp, "u%d=%f\n", i + 1, ux[i]);
    }
    fprintf(fp, "\n对应数值解的误差的绝对值：\n");
    for (i = 0; i < r; i++)
    {
        fprintf(fp, "u%d:%10.3e\n", i + 1, e[i]);
    }
    fprintf(fp, "\n数值解的最大误差：%10.3e\n", maxe);
    fprintf(fp, "2范数：%f\n\n\n", nom);
    fclose(fp);
}

//求一重定积分 积分区间[a,b]，横坐标间隔数dx：越小越精确
double integrate(double (*funcallback)(double x),double a,double b,double dx) {

    double result = 0;  //保存最后结果

    for (double i = a; i <= b; i += dx)
        result += funcallback(i) * dx;

    return result;
}

//高斯-赛德尔迭代法求线性方程组
void gs(double ax[][1000], double ux[], double F[], int r) {

	double tmp1 = 0.0;
	double tmp2 = 0.0;
	double tmp3 = 0.0;
    double eps = 1.0e-5;    //迭代精度
	double dif = 1.0;		//相邻迭代的结果差
	double dis = 0.0;
	double uxc[r-1];		//上一次的迭代结果
	int k = 0;				//迭代次数
	int max = 1000;				//最大迭代次数 undo：求收敛速度

	while ((k < max) && (dif > eps)) {
		for (int i2 = 0; i2 < r; i2++) {
			uxc[i2] = ux[i2];
		}
		for (int i = 0; i < r; i++) {
            tmp1 = 0.0;
        	tmp2 = 0.0;
			for (int j = 0; j <= (i-1); j++) {
				tmp1 += ax[i][j] * ux[j];
			}
			for (int t = (i+1); t < r; t++) {
				tmp2 += ax[i][t] * ux[t];
			}
            ux[i] = (F[i] - tmp1 - tmp2)/ax[i][i];
		}
		//判断是否继续迭代
		tmp3 = 0.0;
		for (int i3 = 0; i3 < r; i3++) {
			dis = fabs(ux[i3] - uxc[i3]);
			if (dis > tmp3) {
				tmp3 = dis;
			}
		}
		dif = tmp3;
        k++;
	}
}
