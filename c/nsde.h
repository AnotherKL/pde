/*******************************************************************************
nsde.h
偏微分方程数值解头文件
*******************************************************************************/
#ifndef NSDE_H_
#define NSDE_H_
//#include <stdbool.h>   //C99特性
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define N 10000
#define PI 3.1415926

/*一些变量符号
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
*/

/*误差分析
 *需要函数: ui()
 *eg: err(x, ux, e, p_maxe, p_nom, r);
 */
void err(double x[], double ux[], double e[], double *maxe, double *nom, int r);

/*追赶法求解三对角矩阵
 *eg: solvematrix(a, b, c, F, ux, r);
 */
void solvematrix(double a[],double b[],double c[],double F[],double ux[],int r);

/*记录结果
 *eg: res(ux, e, maxe, nom, r);
 */
void res(double ux[], double e[], double maxe, double nom, int r);

/*矩形公式求一重定积分
 *积分区间[a,b]，横坐标间隔数dx：越小越精确
 *eg: integrate(function, a, b, dx);
 */
double integrate(double (*funcallback)(double x),double a,double b,double dx);

/*e高斯-赛德尔迭代法求解矩阵（最多31阶）
 *eg：gs(ax, ux, F, r);
 */
void gs(double ax[][1000], double ux[], double F[], int r);

//引用
extern double ui(double x);

#endif
