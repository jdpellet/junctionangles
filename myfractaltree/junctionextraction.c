#include<malloc.h>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define FREE_ARG char*
#define NR_END 1
#define stacksize 10000000000
#define PI 3.14159265359
#define sqrt2 1.414213562373
#define oneoversqrt2 0.707106781186
#define large 1.e^12
#define fillincrement 0.05

int **mask,**draindir,*lakeis,*lakejs,*iup,*idown,*jup,*jdown,lattice_size_x,lattice_size_y,ic,jc;
unsigned long count,*topovecind;
float drop,thresharea,threshhillslopecontrib,**topo,**topo2,**slope,**area,*topovec,flow1,flow2,flow3,flow4,flow5,flow6,flow7,flow8,deltax;

void free_lvector(unsigned long *v, long nl, long nh)
/* free an unsigned long vector allocated with lvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

unsigned long *lvector(long nl, long nh)
/* allocate an unsigned long vector with subscript range v[nl..nh] */
{
	unsigned long *v;

	v=(unsigned long *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(long)));
	return v-nl+NR_END;
}

int *ivector(nl,nh)
long nh,nl;
/* allocate an int vector with subscript range v[nl..nh] */
{
        int *v;

        v=(int *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(int)));
        return v-nl+NR_END;
}

float *vector(nl,nh)
long nh,nl;
/* allocate an int vector with subscript range v[nl..nh] */
{
        float *v;

        v=(float *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(float)));
        return v-nl+NR_END;
}

float **matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float **m;

	/* allocate pointers to rows */
	m=(float **) malloc((size_t)((nrow+NR_END)*sizeof(float*)));
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(float *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float)));
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

int **imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	int **m;

	/* allocate pointers to rows */
	m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;
#define M 7
#define NSTACK 100000

void indexx(unsigned long n, float arr[], unsigned long indx[])
{
	unsigned long i,indxt,ir=n,itemp,j,k,l=1;
	long jstack=0,*istack;
	float a;
	
        istack=lvector(1,NSTACK);
        for (j=1;j<=n;j++) indx[j]=j;
        for (;;) {
                if (ir-l < M) {
                        for (j=l+1;j<=ir;j++) {
                                indxt=indx[j];
                                a=arr[indxt];
                                for (i=j-1;i>=1;i--) {
                                        if (arr[indx[i]] <= a) break;
                                        indx[i+1]=indx[i];
                                }
                                indx[i+1]=indxt;
                        }
                        if (jstack == 0) break;
                        ir=istack[jstack--];
                        l=istack[jstack--];
                } else {
                        k=(l+ir) >> 1;
                        SWAP(indx[k],indx[l+1]);
                        if (arr[indx[l+1]] > arr[indx[ir]]) {
                                SWAP(indx[l+1],indx[ir])
                        }
                        if (arr[indx[l]] > arr[indx[ir]]) {
                                SWAP(indx[l],indx[ir])
                        }
                        if (arr[indx[l+1]] > arr[indx[l]]) {
                                SWAP(indx[l+1],indx[l])
                        }
                        i=l+1;
                        j=ir;
                        indxt=indx[l];
                        a=arr[indxt];
                        for (;;) {
                                do i++; while (arr[indx[i]] < a);
                                do j--; while (arr[indx[j]] > a);
                                if (j < i) break;
                                SWAP(indx[i],indx[j])
                        }
                        indx[l]=indx[j];
                        indx[j]=indxt;
                        jstack += 2;
                        if (ir-i+1 >= j-l) {
                                istack[jstack]=ir;
                                istack[jstack-1]=i;
                                ir=j-1;
                        } else {
                                istack[jstack]=j-1;
                                istack[jstack-1]=l;
                                l=i;
                        }
                }
        }
        free_lvector(istack,1,NSTACK);
}
#undef M
#undef NSTACK
#undef SWAP

void setupgridneighbors()
{    int i,j;

     idown=ivector(1,lattice_size_x);
     iup=ivector(1,lattice_size_x);
     jup=ivector(1,lattice_size_y);
     jdown=ivector(1,lattice_size_y);
     for (i=1;i<=lattice_size_x;i++)
      {idown[i]=i-1;
       iup[i]=i+1;}
     idown[1]=1;
     iup[lattice_size_x]=lattice_size_x;
     for (j=1;j<=lattice_size_y;j++)
      {jdown[j]=j-1;
       jup[j]=j+1;}
     jdown[1]=1;
     jup[lattice_size_y]=lattice_size_y;
}

void calculatedrainagedirections(int i,int j)
{
     float down;

     down=0;draindir[i][j]=0;
	 if (topo[i][j]-topo[iup[i]][j]>down)
	  {down=topo[i][j]-topo[iup[i]][j];draindir[i][j]=6;}
     if (topo[i][j]-topo[idown[i]][j]>down)
	  {down=topo[i][j]-topo[idown[i]][j];draindir[i][j]=2;} 
     if (topo[i][j]-topo[i][jup[j]]>down)
	  {down=topo[i][j]-topo[i][jup[j]];draindir[i][j]=4;}
     if (topo[i][j]-topo[i][jdown[j]]>down)
	  {down=topo[i][j]-topo[i][jdown[j]];draindir[i][j]=8;}
     if ((topo[i][j]-topo[iup[i]][jup[j]])*oneoversqrt2>down)
	  {down=(topo[i][j]-topo[iup[i]][jup[j]])*oneoversqrt2;draindir[i][j]=5;}
     if ((topo[i][j]-topo[idown[i]][jup[j]])*oneoversqrt2>down)
	  {down=(topo[i][j]-topo[idown[i]][jup[j]])*oneoversqrt2;draindir[i][j]=3;}
     if ((topo[i][j]-topo[iup[i]][jdown[j]])*oneoversqrt2>down)
	  {down=(topo[i][j]-topo[iup[i]][jdown[j]])*oneoversqrt2;draindir[i][j]=7;}
     if ((topo[i][j]-topo[idown[i]][jdown[j]])*oneoversqrt2>down)
	  {down=(topo[i][j]-topo[idown[i]][jdown[j]])*oneoversqrt2;draindir[i][j]=1;}
}

void d8route(int i,int j)
{
     switch (draindir[i][j])
	   {case 1 : area[idown[i]][jdown[j]]+=area[i][j];break;
		case 2 : area[idown[i]][j]+=area[i][j];break;
		case 3 : area[idown[i]][jup[j]]+=area[i][j];break;
		case 4 : area[i][jup[j]]+=area[i][j];break;
        case 5 : area[iup[i]][jup[j]]+=area[i][j];break;
		case 6 : area[iup[i]][j]+=area[i][j];break;
		case 7 : area[iup[i]][jdown[j]]+=area[i][j];break;
		case 8 : area[i][jdown[j]]+=area[i][j];break;}
}

void push(int i,int j)
{  
	 count++;
     lakeis[count]=i;
     lakejs[count]=j;
}

void pop()
{  
     ic=lakeis[count];
     jc=lakejs[count];
	 count--;
}

void hydrologiccorrection()
{    int i,j;
     float max;

     count=0;
     for (j=1;j<=lattice_size_y;j++)
	  for (i=1;i<=lattice_size_x;i++)
	   {push(i,j);
	    while (count>0) 
		 {pop();
		  max=topo[ic][jc];
          if (topo[iup[ic]][jc]<max) max=topo[iup[ic]][jc];
          if (topo[idown[ic]][jc]<max) max=topo[idown[ic]][jc];
          if (topo[ic][jup[jc]]<max) max=topo[ic][jup[jc]];
          if (topo[ic][jdown[jc]]<max) max=topo[ic][jdown[jc]];
	      if (topo[iup[ic]][jup[jc]]<max) max=topo[iup[ic]][jup[jc]];
          if (topo[idown[ic]][jdown[jc]]<max) max=topo[idown[ic]][jdown[jc]];
          if (topo[idown[ic]][jup[jc]]<max) max=topo[idown[ic]][jup[jc]];
          if (topo[iup[ic]][jdown[jc]]<max) max=topo[iup[ic]][jdown[jc]];
          if ((topo[ic][jc]>0)&&(topo[ic][jc]<=max)&&(ic>1)&&(jc>1)&&(ic<lattice_size_x)&&(jc<lattice_size_y)&&(count<stacksize))
		   {topo[ic][jc]=max+fillincrement;
			push(ic,jc);
			push(iup[ic],jc);
			push(idown[ic],jc);
			push(ic,jup[jc]);
			push(ic,jdown[jc]);
	        push(iup[ic],jup[jc]);
			push(idown[ic],jdown[jc]);
			push(idown[ic],jup[jc]);
	        push(iup[ic],jdown[jc]);}}}
}

void setupgrids()
{
	 topo=matrix(1,lattice_size_x,1,lattice_size_y);
	 topo2=matrix(1,lattice_size_x,1,lattice_size_y);
	 area=matrix(1,lattice_size_x,1,lattice_size_y);
	 mask=imatrix(1,lattice_size_x,1,lattice_size_y);
	 draindir=imatrix(1,lattice_size_x,1,lattice_size_y);
	 topovec=vector(1,lattice_size_x*lattice_size_y);
	 topovecind=lvector(1,lattice_size_x*lattice_size_y);
	 lakeis=ivector(1,stacksize);
	 lakejs=ivector(1,stacksize);
}

void computecontributingarea()
{   int i,j,m;

	for (j=1;j<=lattice_size_y;j++)
     for (i=1;i<=lattice_size_x;i++)
	  {calculatedrainagedirections(i,j);
       topovec[(j-1)*lattice_size_x+i]=topo[i][j];
	   area[i][j]=deltax*deltax;}
	indexx(lattice_size_x*lattice_size_y,topovec,topovecind);
	m=lattice_size_x*lattice_size_y+1;
	while (m>1)
	 {m--;
      i=(topovecind[m])%lattice_size_x;
      if (i==0) i=lattice_size_x;
      j=(topovecind[m])/lattice_size_x+1;
      if (i==lattice_size_x) j--;
	  if (topo[i][j]>0) d8route(i,j);}
	for (j=1;j<=lattice_size_y;j++)
     for (i=1;i<=lattice_size_x;i++)
	  if (area[i][j]<thresharea) area[i][j]=0;
		  
}

void identifyjunctions()
{   int i,j,n,ilast,jlast,il,jl,i1,j1,i2,j2,i3,j3,i2k,j2k,i1k,j1k;
    float max,max1,slope1,slope2,slope3,dropl,dist,diag,angle1,angle2,angle3,angle13,angle23;
    FILE *fp3;

	fp3=fopen("./myfractaltreejunctions.txt","w");
	for (j=3;j<=lattice_size_y-2;j++)
     for (i=3;i<=lattice_size_x-2;i++)
      if (area[i][j]>thresharea)
	   {//identify pixel in direction of upstream trib 1 (largest trib)
        max=0;
        if ((area[iup[i]][j]>max)&&(area[i][j]>area[iup[i]][j])) {max=area[iup[i]][j];i1=iup[i];j1=j;}
        if ((area[idown[i]][j]>max)&&(area[i][j]>area[idown[i]][j])) {max=area[idown[i]][j];i1=idown[i];j1=j;}
		if ((area[i][jup[j]]>max)&&(area[i][j]>area[i][jup[j]])) {max=area[i][jup[j]];i1=i;j1=jup[j];}
		if ((area[i][jdown[j]]>max)&&(area[i][j]>area[i][jdown[j]])) {max=area[i][jdown[j]];i1=i;j1=jdown[j];}
		if ((area[iup[i]][jup[j]]>max)&&(area[i][j]>area[iup[i]][jup[j]])) {max=area[iup[i]][jup[j]];i1=iup[i];j1=jup[j];}
		if ((area[iup[i]][jdown[j]]>max)&&(area[i][j]>area[iup[i]][jdown[j]])) {max=area[iup[i]][jdown[j]];i1=iup[i];j1=jdown[j];}
		if ((area[idown[i]][jup[j]]>max)&&(area[i][j]>area[idown[i]][jup[j]])) {max=area[idown[i]][jup[j]];i1=idown[i];j1=jup[j];}
        if ((area[idown[i]][jdown[j]]>max)&&(area[i][j]>area[idown[i]][jdown[j]])) {max=area[idown[i]][jdown[j]];i1=idown[i];j1=jdown[j];}
		i1k=i1;j1k=j1;
		//identify pixel in direction of upstream trib 2
		max1=max;max=0;
		if ((area[iup[i]][j]>max)&&(area[iup[i]][j]<max1)) {max=area[iup[i]][j];i2=iup[i];j2=j;diag=1;}
		if ((area[idown[i]][j]>max)&&(area[idown[i]][j]<max1)) {max=area[idown[i]][j];i2=idown[i];j2=j;diag=1;}
		if ((area[i][jup[j]]>max)&&(area[i][jup[j]]<max1)) {max=area[i][jup[j]];i2=i;j2=jup[j];diag=1;}
		if ((area[i][jdown[j]]>max)&&(area[i][jdown[j]]<max1)) {max=area[i][jdown[j]];i2=i;j2=jdown[j];diag=1;}
		if ((area[iup[i]][jup[j]]>max)&&(area[iup[i]][jup[j]]<max1)) {max=area[iup[i]][jup[j]];i2=iup[i];j2=jup[j];diag=sqrt2;}
		if ((area[iup[i]][jdown[j]]>max)&&(area[iup[i]][jdown[j]]<max1)) {max=area[iup[i]][jdown[j]];i2=iup[i];j2=jdown[j];diag=sqrt2;}
		if ((area[idown[i]][jup[j]]>max)&&(area[idown[i]][jup[j]]<max1)) {max=area[idown[i]][jup[j]];i2=idown[i];j2=jup[j];diag=sqrt2;}
        if ((area[idown[i]][jdown[j]]>max)&&(area[idown[i]][jdown[j]]<max1)) {max=area[idown[i]][jdown[j]];i2=idown[i];j2=jdown[j];diag=sqrt2;}
        i2k=i2;j2k=j2;
	    //if area of upstream trib 2 is above a threshold, try to extract junction
		if ((max>thresharea)&&(fabs((area[i1k][j1k]+area[i2k][j2k])/area[i][j])-1)<threshhillslopecontrib)
		 {//search downslope until elevation difference is greater than "jump" to compute slope and direction 
		  dropl=0;
		  i3=i;j3=j;
		  while ((dropl<drop)&&(i3>1)&&(i3<lattice_size_x)&&(j3>1)&&(j3<lattice_size_y)&&(area[i3][j3]>0))
		   {max=0;
	        if (area[iup[i3]][j3]>max) {max=area[iup[i3]][j3];il=iup[i3];jl=j3;diag=1;}
	        if (area[idown[i3]][j3]>max) {max=area[idown[i3]][j3];il=idown[i3];jl=j3;diag=1;}
		    if (area[i3][jup[j3]]>max) {max=area[i3][jup[j3]];il=i3;jl=jup[j3];diag=1;}
		    if (area[i3][jdown[j3]]>max) {max=area[i3][jdown[j3]];il=i3;jl=jdown[j3];diag=1;}
		    if (area[iup[i3]][jup[j3]]>max) {max=area[iup[i3]][jup[j3]];il=iup[i3];jl=jup[j3];diag=sqrt2;}
		    if (area[iup[i3]][jdown[j3]]>max) {max=area[iup[i3]][jdown[j3]];il=iup[i3];jl=jdown[j3];diag=sqrt2;}
		    if (area[idown[i3]][jup[j3]]>max) {max=area[idown[i3]][jup[j3]];il=idown[i3];jl=jup[j3];diag=sqrt2;}
		    if (area[idown[i3]][jdown[j3]]>max) {max=area[idown[i3]][jdown[j3]];il=idown[i3];jl=jdown[j3];diag=sqrt2;}
		    dist+=diag*deltax;
			i3=il;j3=jl;
			mask[i3][j3]=3;
			dropl=topo[i][j]-topo[i3][j3];}
		   slope3=dropl/dist; 
		   //angles are defined with positive values increasing counterclockwise from E, -180 to 180 degrees
	       if (i3!=i) angle3=180/PI*atan2(1.0*j3-j,1.0*i-i3)-180; else {if (j3>j) angle3=-90; else angle3=90;}
		   if (angle3>180) angle3=-(360-angle3);
		   //search upslope along tributary 2 to compute slope and direction
	       dist=0;
		   dropl=0;
		   i2=i2k;j2=j2k;
		   ilast=i2;jlast=j2;
		   while ((dropl<drop)&&(i2>1)&&(i2<lattice_size_x)&&(j2>1)&&(j2<lattice_size_y)&&(area[i2][j2]>0))
		    {max=0;
			 if ((area[iup[i2]][j2]>max)&&(area[iup[i2]][j2]<area[ilast][jlast])) {max=area[iup[i2]][j2];il=iup[i2];jl=j2;diag=1;}
			 if ((area[idown[i2]][j2]>max)&&(area[idown[i2]][j2]<area[ilast][jlast])) {max=area[idown[i2]][j2];il=idown[i2];jl=j2;diag=1;}
			 if ((area[i2][jup[j2]]>max)&&(area[i2][jup[j2]]<area[ilast][jlast])) {max=area[i2][jup[j2]];il=i2;jl=jup[j2];diag=1;}
		     if ((area[i2][jdown[j2]]>max)&&(area[i2][jdown[j2]]<area[ilast][jlast])) {max=area[i2][jdown[j2]];il=i2;jl=jdown[j2];diag=1;}
		     if ((area[iup[i2]][jup[j2]]>max)&&(area[iup[i2]][jup[j2]]<area[ilast][jlast])) {max=area[iup[i2]][jup[j2]];il=iup[i2];jl=jup[j2];diag=sqrt2;}
		     if ((area[iup[i2]][jdown[j2]]>max)&&(area[iup[i2]][jdown[j2]]<area[ilast][jlast])) {max=area[iup[i2]][jdown[j2]];il=iup[i2];jl=jdown[j2];diag=sqrt2;}
		     if ((area[idown[i2]][jup[j2]]>max)&&(area[idown[i2]][jup[j2]]<area[ilast][jlast])) {max=area[idown[i2]][jup[j2]];il=idown[i2];jl=jup[j2];diag=sqrt2;}
             if ((area[idown[i2]][jdown[j2]]>max)&&(area[idown[i2]][jdown[j2]]<area[ilast][jlast])) {max=area[idown[i2]][jdown[j2]];il=idown[i2];jl=jdown[j2];diag=sqrt2;}
		  	 dist+=diag*deltax;
			 ilast=i2;jlast=j2;
			 i2=il;j2=jl;
			 mask[i2][j2]=2;
			 dropl=topo[i2][j2]-topo[i][j];}
			slope2=dropl/dist; 
			if (i2!=i) angle2=-180+180/PI*atan2(1.0*j-j2,1.0*i2-i); else {if (j2>j) angle2=-90; else angle2=90;}
			if (angle2>180) angle2=-(360-angle2);
			dist=0;
			dropl=0;
			i1k=i1;j1k=j1;
			ilast=i1;jlast=j1;
			while ((dropl<drop)&&(i1>1)&&(i1<lattice_size_x)&&(j1>1)&&(j1<lattice_size_y)&&(area[i1][j1]>0))
		     {max=0;
			  if ((area[iup[i1]][j1]>max)&&(area[iup[i1]][j1]<area[ilast][jlast])) {max=area[iup[i1]][j1];il=iup[i1];jl=j1;diag=1;}
			  if ((area[idown[i1]][j1]>max)&&(area[idown[i1]][j1]<area[ilast][jlast])) {max=area[idown[i1]][j1];il=idown[i1];jl=j1;diag=1;}
		      if ((area[i1][jup[j1]]>max)&&(area[i1][jup[j1]]<area[ilast][jlast])) {max=area[i1][jup[j1]];il=i1;jl=jup[j1];diag=1;}
		      if ((area[i1][jdown[j1]]>max)&&(area[i1][jdown[j1]]<area[ilast][jlast])) {max=area[i1][jdown[j1]];il=i1;jl=jdown[j1];diag=1;}
		      if ((area[iup[i1]][jup[j1]]>max)&&(area[iup[i1]][jup[j1]]<area[ilast][jlast])) {max=area[iup[i1]][jup[j1]];il=iup[i1];jl=jup[j1];diag=sqrt2;}
		      if ((area[iup[i1]][jdown[j1]]>max)&&(area[iup[i1]][jdown[j1]]<area[ilast][jlast])) {max=area[iup[i1]][jdown[j1]];il=iup[i1];jl=jdown[j1];diag=sqrt2;}
		      if ((area[idown[i1]][jup[j1]]>max)&&(area[idown[i1]][jup[j1]]<area[ilast][jlast])) {max=area[idown[i1]][jup[j1]];il=idown[i1];jl=jup[j1];diag=sqrt2;}
              if ((area[idown[i1]][jdown[j1]]>max)&&(area[idown[i1]][jdown[j1]]<area[ilast][jlast])) {max=area[idown[i1]][jdown[j1]];il=idown[i1];jl=jdown[j1];diag=sqrt2;}
		  	  dist+=diag*deltax;
			  ilast=i1;jlast=j1;
			  i1=il;j1=jl;
			  mask[i1][j1]=1;
			  dropl=topo[i1][j1]-topo[i][j];}
		    slope1=dropl/dist; 
			if (i1!=i) angle1=180+180/PI*atan2(1.0*j-j1,1.0*i1-i); else {if (j1>j) angle1=-90; else angle1=90;}
			if (angle1>180) angle1=-(360-angle1);
			angle13=fabs(angle1-angle3);
			angle23=fabs(angle2-angle3);
			if (angle13>180) angle13=fabs(360-angle13);
			if (angle23>180) angle23=fabs(360-angle23);
			fprintf(fp3,"%d %d %f %f %f %f %f\n",i,j,angle1,angle2,angle3,angle13,angle23);}}
	fclose(fp3);
}

int main()
{
     FILE *fp0,*fp1,*fp2;
	 int i,j;
     
	 fp0=fopen("./myfractaltreelandscapetopo.txt","r");
	 fp1=fopen("./myfractaltreelandscapearea.txt","w");
	 fp2=fopen("./myfractaltreelandscapemask.txt","w");
	 lattice_size_x=1201;
	 lattice_size_y=800;
	 deltax=1;
	 thresharea=10000;              // m^2
	 threshhillslopecontrib=0.01;   // fraction
	 drop=0.05;                     // m
	 setupgridneighbors();
	 setupgrids();
     for (j=1;j<=lattice_size_y;j++)
	  for (i=1;i<=lattice_size_x;i++)
	   {fscanf(fp0,"%f",&topo[i][j]);
		mask[i][j]=0;} 
     fclose(fp0);	
	 computecontributingarea();
	 identifyjunctions();
	 for (j=1;j<=lattice_size_y;j++)
	  for (i=1;i<=lattice_size_x;i++)
	   {fprintf(fp2,"%d\n",mask[i][j]);
	    fprintf(fp1,"%f\n",area[i][j]);}
     fclose(fp1);fclose(fp2);		
}		 