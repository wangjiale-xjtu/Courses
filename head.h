#include <iostream>
#include <fstream>
#include <math.h>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
// #include "windows.h"
#include <cstring>

#define pi 3.1415926

using namespace std;
struct Type_bran                       //支路(修改)
{
	int i,j; 
	double R,X,YK;
};   

struct Type_gene                   //发电机节点
{
	int i,j;
	double P,Q,V;
};

struct Type_load                        //负荷节点
{      
	int i,j;
	double P,Q,V;
};

struct Type_PV                     //PV节点数组
{  
	double V;
	int i;
};


struct Type_Yii           //自导纳
{
	   double G,B;
};
struct Type_Yij          //互导纳
{
	   double G,B;
	   int j;
};
struct Typet_Yij      //临时互导纳结构体：解决排序和并联问题  
{
	   double G,B,B1;
	   int i,j;
};

struct Type_U      //因子表上三角
{
	double value;
	int j;
}; 

struct Type_volt             //节点电压
{
	double V,theta;
	int i;
};






Type_Yii   *Y_ii,*Y_Bii;  //加入Y_Bii在形成B'时应用，和YBij=Y_ij不用另外占用空间
Type_Yii   *Y_ii_tem;
Type_Yij   *Y_ij;
Type_Yij	*Y_ij_tem;
Typet_Yij  *Y_ijt;
Type_U *U,*U1;

Type_bran     *bran;
Type_gene  *gene;
Type_load       *Load;
Type_PV     *PVNode;
Type_volt *NodalVotage;

int *NYseq,*NYsum;
int Maxtm,dtime;                 //最大迭代次数Maxtm
int	N,Nb,Nt,Ng,Nl,Npv,KP,KQ,ErrorNode;   //Nt变压器支路数
double epsilon,V0,MaxError;          
double **Power_node,**Power_gene;
double *DI;                          //计算各节点的功率误差
double *D,*B,*D1,*B1;                //因子表参数
int *NUsum,*sum_nu_t;
int *SOA,*SOB,*SOD,**SOX;              //节点优化 
                
