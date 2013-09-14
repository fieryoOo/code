#include <iostream>

using namespace std;

int main(int na, char *arg[])
{ 
  if(na!=2)
    {
      cout<<"usage:lf_curve_fit input_curve"<<endl;
      return 0;
    }
  //input five cubic B-spline basis
  double z_max,z_min;
  double small_step=0.01;
  int i,j,nbasis;
  int npts=100;
  nbasis=5;
  char filename[300], buff[300];
  double basis[nbasis][npts+1];
  double tvector[npts+1],temp_t;
  FILE *fbasis;
  for(i=0;i<nbasis;i++)
    {
      sprintf(filename,"B_spline_%d.txt",i);
      if((fbasis=fopen(filename,"r"))==NULL)
	{
	  cout<<filename<<" not found!!"<<endl;
	  return 1;
	}
      for(j=0;j<=npts;j++)
	{
	  if(fscanf(fbasis,"%lf %lf",&temp_t,&basis[i][j])==EOF)
	    {
	      cout<<"erro in npts!!"<<endl;
	      return 1;
	    }
	  if(i==0)
	    tvector[j]=temp_t;
	  else
	    {
	      if(tvector[j]!=temp_t)
		{
		  cout<<"error in tvector!!"<<endl;
		  return 1;
		}
	    }
	}
      fclose(fbasis);
    }
  if(basis[0][0]!=1 || basis[nbasis-1][npts]!=1)
    {
      cout<<"error on boundary!!"<<endl;
      return 1;
    }
  z_max=tvector[npts];
  z_min=tvector[0];
  cout<<"done input basis z_max = "<<z_max<<" z_min = "<<z_min<<endl;
  FILE *f1,*f_fit;
  int n_input;
  double input_z[npts+1],input[npts+1];
  if((f1=fopen(arg[1],"r"))==NULL)
    {
      cout<<arg[1]<<" not found!!"<<endl;
      return 1;
    }
  j=0;
  while(fgets(buff, 300, f1)!=NULL) {
     sscanf(buff, "%lf %lf", &input_z[j], &input[j]);
     if(input_z[j]<z_max && input_z[j]>z_min) j++;
  }
  n_input = j;
  cout<<"n_input ="<<n_input<<endl;
    
  fclose(f1);
  cout<<"done input input n_input = "<<n_input<<endl;
  int input_index[n_input];
  for(j=0;j<npts;j++)
    {
      for(i=0;i<n_input;i++)
	{
	  if(tvector[j]<=input_z[i] && tvector[j+1]>input_z[i])
	    input_index[i]=j;
	}
    }
  //cout<<"pass1!!"<<endl;
  double fit_parameter[nbasis],old_parameter[nbasis];
  //make trial parameter
  for(i=0;i<nbasis;i++)
    {
      fit_parameter[i]=1.0/nbasis;
    }
  //calculate misfit  
  double mis_sum,old_mis_sum,temp1,temp2,temp,mis_cri,mis[nbasis],mis_square_sum,damp_fact;
  damp_fact=1.0;
  mis_cri=0.001;
  old_mis_sum=10000000;
  int k,ii,max_try;
  max_try=100;
  for(k=0;k<max_try;k++)
    {
      mis_sum=0;
      for(j=0;j<n_input;j++)
      	{
	  temp1=0;
	  temp2=0;
	  //cout<<"pass1.5!!"<<endl;
	  for(i=0;i<nbasis;i++)
	    {
	      temp1+=fit_parameter[i]*basis[i][input_index[j]];
	      temp2+=fit_parameter[i]*basis[i][input_index[j]+1];
	    }
	  //cout<<"pass2!!"<<endl;
	  temp=(temp2-temp1)/(tvector[input_index[j]+1]-tvector[input_index[j]])*(input_z[j]-tvector[input_index[j]])+temp1;
	  mis_sum+=(input[j]-temp)*(input[j]-temp);      
	}
      cout<<"parameter : ";
      for(i=0;i<nbasis;i++)
	cout<<fit_parameter[i]<<" ";
      cout<<endl<<"misfit : "<<mis_sum<<endl;
      if(mis_sum<mis_cri)
	{
	  f_fit=fopen("best_fit.txt","w");
	  for(j=0;j<=npts;j++)
	    {
	      temp=0;
	      for(i=0;i<nbasis;i++)
		temp+=fit_parameter[i]*basis[i][j];
	      fprintf(f_fit,"%lf %lf\n",tvector[j],temp);
	    }
	  fclose(f_fit);
	  cout<<"done fit"<<endl;
	  return 0;
	}
      else
	{
	  //patial derivative
	  if(mis_sum>old_mis_sum)
	    {
	      k-=1;
	      damp_fact=damp_fact/2.0;
	      mis_sum=old_mis_sum;
	      for(i=0;i<nbasis;i++)
		fit_parameter[i]=old_parameter[i];
	    }
	  mis_square_sum=0;
	  for(ii=0;ii<nbasis;ii++)
	    {
	      //mis_sum=0;
	      //for(j=0;j<n_input;j++)
	      //	{
	      //  temp1=0;
	      //  temp2=0;
	      //  for(i=0;i<nbasis;i++)
	      //    {
	      //      temp1+=fit_parameter[i]*basis[i][input_index[j]];
	      //      temp2+=fit_parameter[i]*basis[i][input_index[j]+1];
	      //    }
	      //  temp=(temp2-temp1)/(tvector[input_index[j]+1]-tvector[input_index[j]])*(input_z[j]-tvector[input_index[j]])+temp1;
	      //				      mis_sum+=(input[j]-temp)*(input[j]-temp);
	      //		      }
	      //cout<<"pass3!!"<<endl;
	      mis[ii]=0;
	      for(j=0;j<n_input;j++)
		{
		  temp1=0;
		  temp2=0;
		  for(i=0;i<nbasis;i++)
		    {
		      if(i==ii)
			{
			  temp1+=(fit_parameter[i]+small_step)*basis[i][input_index[j]];
			  temp2+=(fit_parameter[i]+small_step)*basis[i][input_index[j]+1];
			}
		      else
			{
			  temp1+=fit_parameter[i]*basis[i][input_index[j]];
			  temp2+=fit_parameter[i]*basis[i][input_index[j]+1];
			}
		    }
		  temp=(temp2-temp1)/(tvector[input_index[j]+1]-tvector[input_index[j]])*(input_z[j]-tvector[input_index[j]])+temp1;
		  mis[ii]+=(input[j]-temp)*(input[j]-temp);
		}
	      mis_square_sum+=(mis_sum-mis[ii])*(mis_sum-mis[ii]);
	    }
	  for(i=0;i<nbasis;i++)
	    {
	      old_parameter[i]=fit_parameter[i];
	      fit_parameter[i]+=mis_sum*(mis_sum-mis[i])/mis_square_sum*small_step*damp_fact;
	    }
	  old_mis_sum=mis_sum;
	}      
    }
  cout<<"max try("<<max_try<<") reach !!"<<endl;
  f_fit=fopen("best_fit.txt","w");
  for(j=0;j<=npts;j++)
    {
      temp=0;
      for(i=0;i<nbasis;i++)
	temp+=fit_parameter[i]*basis[i][j];
      fprintf(f_fit,"%lf %lf\n",tvector[j],temp);
    }
  fclose(f_fit);
  cout<<"done fit"<<endl;
  return 0;
}
