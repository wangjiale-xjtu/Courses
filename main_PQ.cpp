#include "head.h"
#include<stdlib.h>
/*========================���ݵĶ���===========================*/	
void Input_Data(char *filename)
{
	
	int i; 
	
	ifstream in(filename);
	if(!in) cout<<"the input-file is not found!"<<endl;
	in>>Maxtm>>epsilon>>V0;
	in>>N>>Nb>>Nt>>Ng>>Nl;
	
	Power_node=new double *[N+1];
    for(int k=1;k<=N;k++)
		Power_node[k]=new double[3]; 
	Power_gene=new double *[N+1];
	for(int k=1;k<=N;k++)
		Power_gene[k]=new double[3]; 
	
    bran=new Type_bran[Nb+Nt+1];
    gene=new Type_gene[Ng+1];
    Load=new Type_load[Nl+1];
    PVNode=new Type_PV[Ng];
	
	for(i=1;i<=Nb+Nt;i++)
		in>>bran[i].i>>bran[i].j>>bran[i].R>>bran[i].X>>bran[i].YK;
	for(i=1;i<=Ng;i++)
		in>>gene[i].i>>gene[i].j>>gene[i].P>>gene[i].Q>>gene[i].V;
	for(i=1;i<=Nl;i++)
		in>>Load[i].i>>Load[i].j>>Load[i].P>>Load[i].Q>>Load[i].V;    
    in.close();
  
}
/*=======================�ڵ��Ż����=============================*/


void Sort_OptNode()
{
SOA=new int[N+1];
SOB=new int[N+1];
SOD=new int[N+1];
SOX=new int *[N+1];
for(int k=1;k<=N;k++)
SOX[k]=new int[10];                 //10�ڵ����֧·��(��̬����)
int I,I1,I2,I3,I4;
int J,J1,J2,J3,J4;
int K,K1,K11,K2,K3,K4,Num_bala;
/*=================================================================*/
for(K=1;K<=N;K++)                  //�Զ˽ڵ�ų�ֵΪ0
   for(K2=1;K2<=10;K2++)
    SOX[K][K2]=0;
for(I=1;I<=N;I=I+1)                //�ڵ������ֵ��0
    SOD[I]=0;
for(K=1;K<=Ng;K++)                 //Ѱ��ƽ��ڵ��
{
	if(gene[K].j==0) 
	{
		Num_bala=gene[K].i;
		break;
	}
}

for(K=1;K<=Nb+Nt;K++)
{   
	I=bran[K].i;J=bran[K].j;          //ȡ��֧·�Ľڵ��
	if(I==J) continue;                  //�Ƿ�Ϊ�ӵ�֧·
	if(I==Num_bala||J==Num_bala)
	{
		if(I==Num_bala) SOD[I]=100;          //ƽ��ڵ㶨��Ϊ���ȡ����
		if(J==Num_bala) SOD[J]=100;
		continue;
	}
    I1=SOD[I];J1=SOD[J];                   //���˵�����ӵ�֧·��
	for(K1=1;K1<=I1;K1++)
		if (SOX[I][K1]==J)    break;       //��I�����ӵ�֧·����J����ѭ��
	for(K11=1;K11<=J1;K11++)
		if (SOX[J][K11]==I)   break;       //��J�����ӵ�֧·����I����ѭ��
	if (K1>I1||K11>J1)
	{
		I1=I1+1;J1=J1+1; 
		SOD[I]=I1;SOD[J]=J1;
		SOX[I][I1]=J;
		SOX[J][J1]=I;
	}
}
/*===========================���Ը��ڵ��========================================*/

/*===============================================================================*/
    I4=1;
	do 
	{
    I1=101;                  //��ֵҪ����ƽ��ڵ�Ķ�
	for(I=1;I<=N;I++)
	{
		I2=SOD[I];
		if(I2<I1)
		{I1=I2;J=I;} 
	}
 //	  cout<<"ԭ�ڵ�: "<<J<<"    ";
	/*===============================================================================*/
	SOA[I4]=J;SOB[J]=I4;SOD[J]=101; //SOA��žɽڵ���SOB����½ڵ��ţ�ɾ��J
	for(J1=1;J1<=I1;J1++)           //ɾ����J�ڵ����ӵ�����֧·
	{   if(J==Num_bala) break;
		  I3=SOX[J][J1];J2=SOD[I3];
		  for(J3=1;J3<=J2;J3++)
		  {
			  J4=SOX[I3][J3];
			  if (J4==J) break;
		  }
		  SOX[I3][J3]=SOX[I3][J2];
		  SOD[I3]=J2-1;
	}
	/*===========================================================*/
	for(J1=1;J1<I1;J1++)              //Ѱ��ע��֧·
	{   if(J==Num_bala) break;
		  K1=SOX[J][J1];K3=SOD[K1];
		  for(J2=J1+1;J2<=I1;J2++)
		  {
			  K2=SOX[J][J2];
			  for(J3=1;J3<=K3;J3++)
				  if(SOX[K1][J3]==K2) break;
				  if(J3>K3)
				  {K3=K3+1;SOX[K1][K3]=K2;SOD[K1]=K3;
				  K4=SOD[K2]+1;SOX[K2][K4]=K1;SOD[K2]=K4;}
		  }
	}

	I4=I4+1;
	 
	} while (I4<=N);


	  	  for(K=1;K<=Nb+Nt;K++)
	  {  
		  I=bran[K].i;J=bran[K].j;
		  bran[K].i=SOB[I];
		  if(K>Nb)
		      bran[K].j=-SOB[J];       //��¼��ѹ��
		  else
			  bran[K].j=SOB[J];
	  }
	  for(K=1;K<=Nb+Nt;K++)                       //����С����ǰ
		 if(abs(bran[K].i)>abs(bran[K].j))
		  {
			  I=bran[K].i;
		      bran[K].i=bran[K].j;
			  bran[K].j=I;
		  }
      
	 for(K=1;K<=Ng;K++)                           //�ڵ���������
	 { I=gene[K].i;gene[K].i=SOB[I];}
	 for(K=1;K<=Nl;K++)                          
	 { I=Load[K].i;Load[K].i=SOB[I];}

}
/*==========================�γɽڵ㵼�ɾ���=============================*/
 

void Mat_Y()
{
	Y_ii=new Type_Yii[N+1];
	Y_Bii=new Type_Yii[N+1];
	Y_ii_tem=new Type_Yii[N+1];
	Y_ij=new Type_Yij[Nb+Nt+1];
	Y_ij_tem=new Type_Yij[Nb+Nt+1];
	Y_ijt=new Typet_Yij[Nb+Nt+1]; //add
	NYseq=new int[N+1];
	NYsum=new int[N+1];
    int     i,j,L,n;
	double	R,X,YK,Gij,Bij,Zmag2,b_ij;
	int     K,K0,K1,K2,I,J,J1,J3;

	L=0;                                 //���㻥���ɵĸ���
	for(i=1;i<=N;i++)                    //��ʼ��
	{ Y_ii[i].G=0.0  ;  Y_ii[i].B=0.0;
      Y_Bii[i].G=0.0  ;  Y_Bii[i].B=0.0;
	  Y_ii_tem[i].G=0.0 ;  Y_ii_tem[i].B=0.0;}
	for(i=1;i<=Nb+Nt;i++)
	{
		Y_ijt[i].B=0.0;Y_ijt[i].G=0.0;Y_ijt[i].B1=0.0;
	    Y_ij[i].B=0.0;Y_ij[i].G=0.0;
	    Y_ij_tem[i].B=0.0;Y_ij_tem[i].G=0.0;
	}

	for(n=1;n<=Nb+Nt;n++)               //֧·���ݴ���
	{
		i=bran[n].i;
		j=bran[n].j;
		R=bran[n].R;
		X=bran[n].X;
		YK=bran[n].YK;
		if(i==j)                       //�ӵ�֧·����XΪ�ӵص翹
		{
		    Y_ii[i].B=Y_ii[i].B-1.0/X;
            Y_ii_tem[i].B=Y_ii_tem[i].B-1.0/X;
			continue; 
		}                             
		Zmag2=R*R+X*X;
		Gij=R/Zmag2;
		Bij=-X/Zmag2;
		b_ij=-1.0/X;
        if(n<=Nb)
		{
			Y_ii[i].G=Y_ii[i].G+Gij;            //��·֧·�Ե���
			Y_ii[i].B=Y_ii[i].B+Bij+YK;
			Y_ii[j].G=Y_ii[j].G+Gij;
			Y_ii[j].B=Y_ii[j].B+Bij+YK;
			Y_ii_tem[i].B=Y_ii_tem[i].B+b_ij+YK;
			Y_ii_tem[j].B=Y_ii_tem[j].B+b_ij+YK;
            
			L=L+1;                            //��·֧·�Ļ�����
			Y_ijt[L].G=-Gij;                  
	    	Y_ijt[L].B=-Bij;
			Y_ijt[L].B1=-b_ij;
		}
		else
		{   
			if(j<0)                          //��ѹ��֧·�Ե���
			{
			Y_ii[i].G=Y_ii[i].G+Gij;         
			Y_ii[i].B=Y_ii[i].B+Bij;
			Y_ii_tem[i].B=Y_ii_tem[i].B+b_ij;
			Y_ii[abs(j)].G=Y_ii[abs(j)].G+Gij/YK/YK;
			Y_ii[abs(j)].B=Y_ii[abs(j)].B+Bij/YK/YK;
			Y_ii_tem[abs(j)].B=Y_ii_tem[abs(j)].B+b_ij/YK/YK;
			}
			else
			{
			Y_ii[j].G=Y_ii[j].G+Gij;            
			Y_ii[j].B=Y_ii[j].B+Bij;
			Y_ii_tem[j].B=Y_ii_tem[j].B+b_ij;
			Y_ii[abs(i)].G=Y_ii[abs(i)].G+Gij/YK/YK;
			Y_ii[abs(i)].B=Y_ii[abs(i)].B+Bij/YK/YK;
			Y_ii_tem[abs(i)].B=Y_ii_tem[abs(i)].B+b_ij/YK/YK;
			}

           	L=L+1;                          //��ѹ��֧·�Ļ�����
			Y_ijt[L].G=-Gij/YK;              
			Y_ijt[L].B=-Bij/YK;
			Y_ijt[L].B1=-b_ij/YK;
		}
                            
		Y_ijt[L].i=i;                       //��·�ͱ�ѹ���Ķ˵�Ŵ洢
		Y_ijt[L].j=j;		
	}

/*================��˳����д�����ɾ���===========================*/
	J=0;K0=0;
	for(I=1;I<=N;I++)
	{   J1=0;
		for(K=1;K<=L;K++)
		{
			if (abs(Y_ijt[K].i)!=I) continue;
			J3=abs(Y_ijt[K].j);
			for(K1=1;K1<=J1;K1++)
			{
		    	K2=K0+K1;
				if(abs(Y_ij[K2].j)==J3)
					break;
			}
			if(K1>J1)
			{
				J=J+1;
				Y_ij[J].G=Y_ijt[K].G;
				Y_ij[J].B=Y_ijt[K].B;
                Y_ij_tem[J].B=Y_ijt[K].B1;
				Y_ij[J].j=J3;
                Y_ij_tem[J].j=J3;
                J1=J1+1;
			}
			else
			{
				Y_ij[K2].G=Y_ij[K2].G+Y_ijt[K].G;
				Y_ij[K2].B=Y_ij[K2].B+Y_ijt[K].B;
				Y_ij_tem[K2].B=Y_ij_tem[K2].B+Y_ijt[K].B1;
			}
		}
        NYsum[I]=J1;                     //Y Y'��ͬ��һ��NYsum;
        K0=K0+J1;
	}
	NYseq[1]=1;
    for(I=1;I<=N-1;I++)
	NYseq[I+1]=NYseq[I]+NYsum[I]; 	

	for(n=1;n<=N;n++)                   //���ݴ����γ�B'����
	{
		Y_Bii[n].B=Y_ii[n].B;
		//cout<<Y_Bii[n].B<<endl;
	}
	for(n=1;n<=Nb+Nt;n++)               
	{
		i=bran[n].i;
		j=bran[n].j;
		R=bran[n].R;
		X=bran[n].X;
		YK=bran[n].YK;

		Zmag2=R*R+X*X;
		Gij=R/Zmag2;
		Bij=-X/Zmag2;

		if(i==j)                        //�Ե�֧·����
		{
			Y_Bii[i].B=Y_ii[i].B+1.0/X;   //���벻�ӵ�֧·��Ӱ���Ե��ɣ�XΪ�ӵص���
			continue; 
		}

        if(n<=Nb)
		{
			Y_Bii[i].B=Y_Bii[i].B-YK;
			Y_Bii[j].B=Y_Bii[j].B-YK;
		}
		else
		{
			if(j<0)
			{
               Y_Bii[i].B=Y_Bii[i].B-Bij+Bij/YK;
			   Y_Bii[abs(j)].B=Y_Bii[abs(j)].B-Bij/YK/YK+Bij/YK;

			}
			else
			{
               Y_Bii[j].B=Y_Bii[j].B-Bij+Bij/YK;
			   Y_Bii[abs(i)].B=Y_Bii[abs(i)].B-Bij/YK/YK+Bij/YK;

			}
		}

	}

}
/*==========================B' B''�Լ����ӱ�=================================*/
  

void FactorForm(int flag)
{  	
	int i,j,count,n_u;
	int i_pv=0,i_above,m;    
	double	Btemp,n,p;              //m,nΪ��ʱ����
    Npv=Ng-1;
    for(i=2;i<=Ng;i++)              //��¼PV�ڵ����
	{PVNode[i-1].V=gene[i].V;
	PVNode[i-1].i=gene[i].i;}                                     
	U=new Type_U[2*(Nb+Nt)+1];      //U�ĸ�����ȷ��(Nb�ı�������ȡ2)
	U1=new Type_U[2*(Nb+Nt)+1];
	B=new double[N+1];
	B1=new double[N+1];
	D=new double[N+1];
	D1=new double[N+1];
    NUsum=new int[N+1];
	sum_nu_t=new int[N+1];
	                    
	

	do 
	{
		for(i=1;i<N;i++)
	{   		
		if(flag==2)
		{
			for(count=1;count<=Npv;count++)
			{
				j=PVNode[count].i;
				if(i==j) {i_pv=i;break;}
			}
		}                                     //Ѱ��PV�ڵ��
        if(flag==2&&i==i_pv) 
		{
			sum_nu_t[i]=0;D1[i]=0.0;
			continue;
		}
		for(count=1;count<=N;count++)        //�幤������
		{B[count]=0.0;B1[count]=0.0;}

		if(flag==2)
			B1[i]=Y_ii_tem[i].B;
		else B[i]=Y_Bii[i].B;

		for(count=NYseq[i];count<NYseq[i+1];count++)
		{   
			if(flag==2)
			{j=Y_ij_tem[count].j;B1[j]=Y_ij_tem[count].B;}
			else
			{j=Y_ij[count].j;B[j]=Y_ij[count].B;}

		}                                    //����A��B��F����i�еĳ�ʼ��
		/*======================================================*/
		if(flag==1) B[N]=0.0;
		if(flag==2)
		{
			//for(count=1;count<=Npv;count++)
			//{j=PVNode[count].i;B1[j]=0.0;}
			for(count=1;count<=Ng;count++)
			{j=gene[count].i;B1[j]=0.0;}
		}                                   //�������C����
		/*======================================================*/
	


		n_u=1;
		for(i_above=1;i_above<=i-1;i_above++)
		{   count=1;
            if(flag==2)  m=sum_nu_t[i_above];
			else m=NUsum[i_above];
			while(count<=m)
			{   if(flag==2) n=U1[n_u].j;
			    else n=U[n_u].j;
				if(n==i) break;
				n_u=n_u+1;
                count=count+1;
			}
            if(count>m) continue;
            if(flag==2) Btemp=U1[n_u].value/D1[i_above];
			else Btemp=U[n_u].value/D[i_above];
			while(count<=m)
			{
				if(flag==2)
				{j=U1[n_u].j;B1[j]=B1[j]-Btemp*U1[n_u].value;}
				else {j=U[n_u].j;B[j]=B[j]-Btemp*U[n_u].value;}
				count=count+1;
				n_u=n_u+1;
			}

		}                                     //��С��D����ȥ����
    

        if(flag==2)
		{Btemp=1.0/B1[i];D1[i]=Btemp;}
		else 
		{Btemp=1.0/B[i];D[i]=Btemp;}
		count=0;
		for(j=i+1;j<N;j++)
		{  
			if(flag==2)  p=B1[j];
		    else p=B[j];
			if(p==0) continue;
			if(flag==2)
			{
				U1[n_u].value=B1[j]*Btemp;
				U1[n_u].j=j;}
			else 
			{
				U[n_u].value=B[j]*Btemp;
				U[n_u].j=j;
			}
		    count=count+1;
			n_u=n_u+1;
		}
		if(flag==2) sum_nu_t[i]=count;                   
		else NUsum[i]=count;                //��С��E���洢����
	  /*======================================================*/		
	}
	flag=flag+1;
	} while (flag<=2);
	
	
}

/*========================ϵͳ��ѹ��ֵ����=============================*/
void Initialization()
{
	int i,n;
    NodalVotage=new Type_volt[N+1];
	KP=1;KQ=1;
	NodalVotage[N].V=gene[1].V;
    NodalVotage[N].theta=0.0;
	for(i=1;i<N;i++)
	{
    	NodalVotage[i].V=V0;
    	NodalVotage[i].theta=0.0;
	}
    for(n=1;n<=Npv;n++)
	{
		i=PVNode[n].i;
        NodalVotage[i].V=PVNode[n].V;
	}
	
}
/*===================�ڵ㹦�ʵļ���======================*/
	
void PQ_Calcu(int flag)
{   
    int i,j,n;
	double A,B,VV,theta,Vi;
   
	for(i=1;i<=N;i++)
	{
		Power_node[i][flag]=0.0;
		Power_gene[i][flag]=0.0;
	}
	for(i=1;i<=N;i++)
	{
		Vi=NodalVotage[i].V;
		if(flag==1)
			 A=Y_ii[i].G;
		else A=-Y_ii[i].B;
        Power_node[i][flag]=Power_node[i][flag]+Vi*Vi*A;
		if(i==N) break;
		for(n=NYseq[i];n<=NYseq[i+1]-1;n++)
		{
		    if(flag==1)
			{A=Y_ij[n].G;B=Y_ij[n].B;}
			else
            {A=-Y_ij[n].B;B=Y_ij[n].G;}
			j=Y_ij[n].j;
			VV=Vi*NodalVotage[j].V;
			theta=NodalVotage[i].theta-NodalVotage[j].theta;
			A=A*VV*cos(theta);
			B=B*VV*sin(theta);
			Power_node[i][flag]=Power_node[i][flag]+A+B;
			Power_node[j][flag]=Power_node[j][flag]+A-B;

		}

	}


}
//====================������ڵ�Ĺ�������������=============================//

void Max_err(int flag)
{   DI=new double[N+1]; 
	int i,j,count,n_l,n_g;
	int i_l,i_g;
	int n_pv,i_pv;
    double Wi,Vi,Wtemp;

	MaxError=0.0;
	i=1;n_g=1;n_l=1;n_pv=1;
	i_g=gene[2].i;
	i_l=Load[1].i;
	i_pv=PVNode[1].i;
	for(count=1;count<=N;count++)
    DI[count]=0.0;
	while(1)
	{ 
		Vi=NodalVotage[i].V;
        for(count=1;count<=Nl;count++)
		{
			j=Load[count].i;
			if(i==j) {i_l=i;break;}
		}

		if(i==i_l)
		{
			if(flag==1) Wi=-Load[count].P;
			else        Wi=-Load[count].Q;
		}
        else Wi=0.0;
		Wtemp=Wi;
		Wi=Wi-Power_node[i][flag];

		for(count=1;count<=Ng;count++)
		{
			j=gene[count].i;
			if(i==j) {i_g=i;break;}
		}
		if(i==i_g)
		{
           Power_node[i][flag]=Wtemp;
           Power_gene[i][flag]=-Wi;
		   if(flag==1) Wi=Wi+gene[count].P;
		   else        Wi=Wi+gene[count].Q;
		}

		if(i==N) break;

		if(flag==2)
		{
			for(count=1;count<=Ng;count++)
			{
				j=gene[count].i;
				if(i==j) {i_pv=i;break;}
			}
		}                                     //Ѱ��PV�ڵ��
		if(flag==2&&i==i_pv)
		{
			DI[i]=0.0;
		}
        else
		{
			if(fabs(Wi)>MaxError)
			{
                MaxError=fabs(Wi);
			    ErrorNode=i;
			}
			DI[i]=Wi/Vi;
		}
		DI[i]=Wi/Vi;
		i=i+1;

	}	
}
//====================�������������===============================//
void Solve_Equ(int flag)
{ 
	int i,j,n_u,count,m,n;
    double DItemp;
    n_u=1;
    for(i=1;i<=N-1;i++)
	{
 	   DItemp=DI[i];
	   if(flag==1) m=NUsum[i];
	   else m=sum_nu_t[i];
	   for(count=1;count<=m;count++)
	   {
		   if(flag==1)
		   { 
			   j=U[n_u].j;
		       DI[j]=DI[j]-DItemp*U[n_u].value;
		   }
		   else
		   { 
			   j=U1[n_u].j;
		       DI[j]=DI[j]-DItemp*U1[n_u].value;
		   }
		    n_u=n_u+1;
	   }
	   if(flag==1)
		   DI[i]=DItemp*D[i];
	   else
		   DI[i]=DItemp*D1[i];

	}

    for(i=N-1;i>=1;i--)
	{
		DItemp=DI[i];
		if(flag==1) n=NUsum[i];
		else n=sum_nu_t[i];
		for(count=1;count<=n;count++)
		{
		 n_u=n_u-1;
		 if(flag==1)
		 {
			 j=U[n_u].j;
			 DItemp=DItemp-DI[j]*U[n_u].value;
		 }
		 else
		 {
			 j=U1[n_u].j;
			 DItemp=DItemp-DI[j]*U1[n_u].value;
		 }
		}
	    DI[i]=DItemp;
	}
}
/*======================== ������  =============================*/
void Output_Results()
{   
	double Vmin,V,theta,P,Q;
    int VminNode,i,j,n_g,i_g,count;
	
	double PLoss,QLoss,Pl,Pg,Ql,Qg,R,X,YK,Zmag2;  //����������֧·�����й�
	double Vi,Vj;
	double Ei,Fi,Ej,Fj,DE,DF;
	double Pij,Pji,Qij,Qji,Ir,Ii;
	int n;
	PLoss=0.0;QLoss=0.0;Pl=0;Pg=0;Ql=0;Qg=0;//add


	Vmin=NodalVotage[1].V;
	i_g=gene[1].i;
	VminNode=1;
	n_g=1;
	P=0.0;Q=0.0;



	dtime=0;
	KP=1;
	KQ=1;
	while(dtime<=Maxtm)
	{
		PQ_Calcu(1);      //�й�����
		Max_err(1);
		if(MaxError<epsilon) 
		{
			KP=0;
			if(KQ==0) break;
			else goto L1;
		}
		Solve_Equ(1);
		for(i=1;i<=N-1;i++)
			NodalVotage[i].theta=NodalVotage[i].theta-DI[i]/V0;
	
L1:	
		PQ_Calcu(2);      //�޹�����
		Max_err(2);
		if(MaxError<epsilon)
		{
			KQ=0;
			if(KP==0) break;
			else goto L2;
		}
		Solve_Equ(2);
		for(i=1;i<=N-1;i++)
			NodalVotage[i].V=NodalVotage[i].V-DI[i];
		
L2:
		dtime=dtime+1; 		
	}

	ofstream out("output.txt");                         
	ofstream outcsvNode;
	if(!out) cout<<"Can't open the input file"<<endl;
	out<<endl<<endl<<endl;
    out<<"***********PQ�ֽⷨ����������(�ڵ���֧·�������ѱ�����.csv�ļ�)**********"<<endl;
    out<<endl;
    out<<"�ڵ���:   "<<N<<endl;
	out<<"֧·��:   "<<Nb+Nt<<endl;
	out<<"�������: "<<Ng<<endl;
    out<<"������: "<<Nl<<endl;
    out<<"��������: "<<dtime<<endl;
	out<<"��������: "<<epsilon<<endl;

	out.setf(ios::left);
	
	outcsvNode.open("outputNode.csv", ios::out | ios::trunc); //���csv�ļ�
	outcsvNode << "�ڵ�" << "," << "��ѹ��ֵ" << "," << "���" << "," << "�����й�" << "," << "�����޹�"<<"," << "������й�"<<"," << "������޹�" << endl;

	for(i=1;i<=N;i++)
	{
       theta=NodalVotage[i].theta/pi*180;
	   V=NodalVotage[i].V;
	   if(V<Vmin)
	   {
		   Vmin=V;
		   VminNode=i;
	   }
	   for(count=1;count<=Ng;count++)
	   {
			j=gene[count].i;
			if(i==j) {i_g=i;break;}
	   }
	   if(i==i_g)
	   {
         P=Power_gene[i][1];
		 Q=Power_gene[i][2];
	   }
	   else
	   {
         P=0.0;
		 Q=0.0;
	   }
	   Pl+=Power_node[i][1];Ql+=Power_node[i][2];  //����ͳ��
	   Pg+=P;Qg+=Q;

	   outcsvNode << SOA[i]<<","<<V << ","<<theta << ","<< Power_node[i][1]<<","<< Power_node[i][2]<<","<<P<<","<<Q<< endl;  //���CSV�ڵ���Ϣ
	}
	outcsvNode.close();
	
	out<<setw(20)<<"ϵͳ��С��ѹ:     "<<setw(20)<<Vmin<<endl;
	out<<setw(20)<<"ϵͳ��С��ѹ�ڵ�: "<<setw(20)<<SOA[VminNode]<<endl;
	out<<setw(20)<<"�ܸ����ù�����:   "<<setw(20)<<Pl << endl;
	
	
	out<<setw(20)<<"�ܸ����޹�����:   "<<setw(20)<<Ql<<endl;
	out <<setw(20)<<"�ܷ������ù�����: "<<setw(20)<<Pg<< endl;
	out<<setw(20)<<"�ܷ������޹�����: "<<setw(20)<<Qg<<endl;
	
     
  

    out.setf(ios::left);
	ofstream outcsvLine;
	outcsvLine.open("outputLine.csv", ios::out | ios::trunc); //���csv�ļ�
	outcsvLine <<"���"<<  "," << "�ڵ�i" <<  "," << "�ڵ�j" <<  "," << "�й�Pij" <<  "," << "�޹�Qij" << "," << "�й�Pji" <<  "," << "�޹�Qji" << endl;
	for(n=1;n<=Nb+Nt;n++)
	{
		i=abs(bran[n].i);j=abs(bran[n].j);
		R=bran[n].R;X=bran[n].X;YK=bran[n].YK;
		
		Vi=NodalVotage[i].V;
		theta=NodalVotage[i].theta;
		Ei=Vi*cos(theta);
		Fi=Vi*sin(theta);

     	Vj=NodalVotage[j].V;
    	theta=NodalVotage[j].theta;
	    Ej=Vj*cos(theta);
	    Fj=Vj*sin(theta);

		if(bran[n].i<0||bran[n].j<0)
		{
			if(bran[n].i<0)
			{
				Ei=Ei/YK;Fi=Fi/YK;
			}
		    else
			{
				Ej=Ej/YK;Fj=Fj/YK;
			}
			YK=0.0;
		}

		DE=Ei-Ej;DF=Fi-Fj;
		Zmag2=R*R+X*X;
		Ir=(DE*R+DF*X)/Zmag2;
		Ii=(DF*R-DE*X)/Zmag2;

		Pij=Ir*Ei+Ii*Fi;
		Qij=Ir*Fi-Ii*Ei;
		Pji=-Ir*Ej-Ii*Fj;
		Qji=-Ir*Fj+Ii*Ej;

		Qij=Qij-Vi*Vi*YK;
		Qji=Qji-Vj*Vj*YK;

		if(i==j) Qij=Vi*Vi*(1.0/X);

		PLoss=PLoss+Pij+Pji;
		QLoss=QLoss+Qij+Qji;
      	

		outcsvLine << n<< "," << SOA[i] << "," << SOA[j] << "," << Pij << "," << Qij << "," << Pji << "," << Qji << endl;
	}                   
	outcsvLine.close();

	    out<<setw(20)<<"ϵͳ���й����PLoss: "<<setw(20)<<PLoss<<endl;
	    out<<setw(20)<<"ϵͳ���޹����QLoss: "<<setw(20)<<QLoss<<endl;
        out.close();
		
}
/*===================================================================*/
int main()
{                              
	char filename[10],ch;
	while (1)
	{
	
	system("cls");
	cout<<endl
		<< "                    *                                           *        " << endl
		<< "                  **                                          **  " << endl
		<< " *------------ *** --------------------------------------- *** -------------*" << endl
		<< " |           ****                                        ****               |" << endl
		<< " |         *********        PQ�ֽⷨ�����������       *********            |" << endl
		<< " |           *****                                       *****              |" << endl
		<< " *--------- **** -------------------------------------- **** ---------------*" << endl
		<< "            ***                                         ***  " << endl
		<< "           **                                          **  " << endl
		<< "          *                                           * " << endl
		
		
		<<endl
		<<"��ѡ����Ҫʹ�õĲ���ϵͳ"<<endl
		<<endl
		<<"                               1:  IEEE 5 �ڵ�ϵͳ"<<endl
		<<"                               2:  IEEE 14�ڵ�ϵͳ"<<endl
		<<"                               3:  IEEE 30�ڵ�ϵͳ"<<endl
		<<"                               4:  IEEE 57�ڵ�ϵͳ"<<endl
		<<"                               5:  IEEE118�ڵ�ϵͳ"<<endl
		<<"                               0:  �˳�       "<<endl
		   <<endl;	
	cin>>ch;
	if(ch=='0') break;
	switch (ch)
	{
		case '1' : strcpy_s(filename,"5.txt");
			break;
		case '2' : strcpy_s(filename,"14.txt");
			break;
		case '3' : strcpy_s(filename,"30.txt");
			break;
		case '4' : strcpy_s(filename,"57.txt");
			break;
		case '5' : strcpy_s(filename,"118.txt");
			break;
		default: cout<<"��������������ѡ��"<<endl;
			system("pause");
			continue;
		}
	
	Input_Data(filename);        //��������

 	Sort_OptNode();          //�ڵ��Ż�

	Mat_Y();	                 //�ڿ���ϡ�����洢������ǰ�����γɵ��ɾ���

	FactorForm(1);                  //�γ����ӱ�	

	Initialization();                //ϵͳ����ֵ 
	
    Output_Results();               //���������
	
	
	cout<<"����ѱ�����.txt�ļ���.csv�ļ���"<<endl;
	cout<<"���뷵�ز���ϵͳѡ��׶Σ��밴�»س���"<<endl;
	getchar();
	if(getchar()=='0') break;
	}
	 	
}
