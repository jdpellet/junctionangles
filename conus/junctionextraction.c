#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define FREE_ARG char*
#define NR_END 1
#define stacksize 1000000000
#define PI 3.14159265359
#define sqrt2 1.414213562373
#define oneoversqrt2 0.707106781186
#define large 10000000
#define small 0.000001
#define fillincrement 0.05
#define mfdweight 1.1

int **mask,**draindir,**intern,*avcount,*flatis,*flatjs,*avcount2,*avcount3,*avcount4,ind,numbins,threshmaskarea,threshfilled;
unsigned long count,*topovecind,*iup,*idown,*jup,*jdown,lattice_size_x,lattice_size_y,ic,jc,countermask;
float filled,max,threshdrop,threshslope,threshmfdtod8,maxarea,thresharea,threshhillslopecontrib,oneoverdeltax,elevc,*avangle,*avangle2,*avangle3,*avslope,**angle,**slopex,**slopey,**topoold,**topo,**slope,**area,*topovec,flow1,flow2,flow3,flow4,flow5,flow6,flow7,flow8,deltaxinkm,deltax;

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
{    unsigned long i,j;

     idown=lvector(1,lattice_size_x);
     iup=lvector(1,lattice_size_x);
     jup=lvector(1,lattice_size_y);
     jdown=lvector(1,lattice_size_y);
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

void setupgrids()
{
	 topo=matrix(1,lattice_size_x,1,lattice_size_y);
	 topoold=matrix(1,lattice_size_x,1,lattice_size_y);
	 area=matrix(1,lattice_size_x,1,lattice_size_y);
	 mask=imatrix(1,lattice_size_x,1,lattice_size_y);
	 intern=imatrix(1,lattice_size_x,1,lattice_size_y);
	 draindir=imatrix(1,lattice_size_x,1,lattice_size_y);
	 slope=matrix(1,lattice_size_x,1,lattice_size_y);
	 slopex=matrix(1,lattice_size_x,1,lattice_size_y);
	 slopey=matrix(1,lattice_size_x,1,lattice_size_y);
	 angle=matrix(1,lattice_size_x,1,lattice_size_y);
	 topovec=vector(1,lattice_size_x*lattice_size_y);
	 topovecind=lvector(1,lattice_size_x*lattice_size_y);
	 avcount=ivector(1,numbins);
	 avcount2=ivector(1,numbins);
	 avangle=vector(1,numbins);
	 avangle2=vector(1,numbins);
	 avslope=vector(1,numbins);
	 flatis=ivector(1,stacksize);
	 flatjs=ivector(1,stacksize);
}

void push(i,j)
unsigned long i,j;
{
	 count++;
	 flatis[count]=i;
     flatjs[count]=j;
}

void pop()
{
     ic=flatis[count];
     jc=flatjs[count];
	 count--;
}

void hydrologiccorrection()
{    unsigned long i,j;

     for (j=1;j<=lattice_size_y;j++)
	  for (i=1;i<=lattice_size_x;i++)
	   {filled=0;
		count=0;
	    if (intern[i][j]==0) push(i,j);
	    while ((count>0)&&(filled<threshfilled))
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
          if ((topo[ic][jc]<=max)&&(count<stacksize)&&(topo[ic][jc]>-99)&&(intern[ic][jc]==0))
		   {topo[ic][jc]=max+fillincrement;
		    if (topo[ic][jc]-topoold[ic][jc]>filled) filled=topo[ic][jc]-topoold[ic][jc];
		    push(ic,jc);
			push(iup[ic],jc);
			push(idown[ic],jc);
			push(ic,jup[jc]);
			push(ic,jdown[jc]);
	        push(iup[ic],jup[jc]);
			push(idown[ic],jdown[jc]);
			push(idown[ic],jup[jc]);
			push(iup[ic],jdown[jc]);}}
		count=0;
		if (filled>=threshfilled)
		 {push(i,j);
		  while (count>0)
		   {pop();
			if ((intern[ic][jc]==0)&&(topo[ic][jc]>(topoold[ic][jc]+0.001)))
			 {topo[ic][jc]=topoold[ic][jc];
			  intern[ic][jc]=1;
			  push(iup[ic],jc);
			  push(idown[ic],jc);
			  push(ic,jup[jc]);
			  push(ic,jdown[jc]);
	          push(iup[ic],jup[jc]);
			  push(idown[ic],jdown[jc]);
			  push(idown[ic],jup[jc]);
			  push(iup[ic],jdown[jc]);}}}}
	 count=0;
	 for (j=1;j<=lattice_size_y;j++)
	 {for (i=1;i<=lattice_size_x;i++)
	   {filled=0;
	    if (intern[i][j]==1) push(i,j);
	    while ((count>0)&&(filled<threshfilled))
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
          if ((intern[ic][jc]==1)&&(topo[ic][jc]<=max)&&(count<stacksize)&&(topo[ic][jc]>-99))
		   {topo[ic][jc]=max+fillincrement;
		    if (topo[ic][jc]-topoold[ic][jc]>filled) filled=topo[ic][jc]-topoold[ic][jc];
		    push(ic,jc);
			push(iup[ic],jc);
			push(idown[ic],jc);
			push(ic,jup[jc]);
			push(ic,jdown[jc]);
	        push(iup[ic],jup[jc]);
			push(idown[ic],jdown[jc]);
			push(idown[ic],jup[jc]);
			push(iup[ic],jdown[jc]);}}}}
}

void calculatedrainagedirections(i,j)
int i,j;
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
     slopex[i][j]=0.5*oneoverdeltax*(topo[idown[i]][j]-topo[iup[i]][j]);
	 slopey[i][j]=0.5*oneoverdeltax*(topo[i][jdown[j]]-topo[i][jup[j]]);
}

void d8route(i,j)
int i,j;
{
     switch (draindir[i][j])
	   {case 1 : area[idown[i]][jdown[j]]+=area[i][j];slopex[idown[i]][jdown[j]]+=slopex[i][j];slopey[idown[i]][jdown[j]]+=slopey[i][j];break;
		case 2 : area[idown[i]][j]+=area[i][j];slopex[idown[i]][j]+=slopex[i][j];slopey[idown[i]][j]+=slopey[i][j];break;
		case 3 : area[idown[i]][jup[j]]+=area[i][j];slopex[idown[i]][jup[j]]+=slopex[i][j];slopey[idown[i]][jup[j]]+=slopey[i][j];break;
		case 4 : area[i][jup[j]]+=area[i][j];slopex[i][jup[j]]+=slopex[i][j];slopey[i][jup[j]]+=slopey[i][j];break;
        case 5 : area[iup[i]][jup[j]]+=area[i][j];slopex[iup[i]][jup[j]]+=slopex[i][j];slopey[iup[i]][jup[j]]+=slopey[i][j];break;
		case 6 : area[iup[i]][j]+=area[i][j];slopex[iup[i]][j]+=slopex[i][j];slopey[iup[i]][j]+=slopey[i][j];break;
		case 7 : area[iup[i]][jdown[j]]+=area[i][j];slopex[iup[i]][jdown[j]]+=slopex[i][j];slopey[iup[i]][jdown[j]]+=slopey[i][j];break;
		case 8 : area[i][jdown[j]]+=area[i][j];slopex[i][jdown[j]]+=slopex[i][j];slopey[i][jdown[j]]+=slopey[i][j];break;}
}

void computecontributingarea()
{   unsigned long i,j,m;
    float slopexl,slopeyl,slopel;
	
	maxarea=0;
	for (j=1;j<=lattice_size_y;j++)
     for (i=1;i<=lattice_size_x;i++)
	  {calculatedrainagedirections(i,j);
       topovec[(j-1)*lattice_size_x+i]=topo[i][j];
	   area[i][j]=deltaxinkm*deltaxinkm;
	   mask[i][j]=0;}
	indexx(lattice_size_x*lattice_size_y,topovec,topovecind);
	m=lattice_size_x*lattice_size_y+1;
	while (m>1)
	 {m--;
      i=(topovecind[m])%lattice_size_x;
      if (i==0) i=lattice_size_x;
      j=(topovecind[m])/lattice_size_x+1;
      if (i==lattice_size_x) j--;
	  d8route(i,j);
	  if (area[i][j]>maxarea) maxarea=area[i][j];}
	for (j=1;j<=lattice_size_y;j++)
     for (i=1;i<=lattice_size_x;i++)
	  {slopex[i][j]/=area[i][j]/(deltaxinkm*deltaxinkm);
	   slopey[i][j]/=area[i][j]/(deltaxinkm*deltaxinkm);
	   slope[i][j]=sqrt(slopex[i][j]*slopex[i][j]+slopey[i][j]*slopey[i][j]);
	   if (slopex[i][j]!=0) angle[i][j]=-180/PI*atan2(slopey[i][j],slopex[i][j]); else {if (slopey[i][j]>0) angle[i][j]=-90; else angle[i][j]=90;}
       slopexl=0.5*oneoverdeltax*(topo[idown[i]][j]-topo[iup[i]][j]);
	   slopeyl=0.5*oneoverdeltax*(topo[i][jdown[j]]-topo[i][jup[j]]);
	   slopel=sqrt(slopexl*slopexl+slopeyl*slopeyl); 
	   if ((area[i][j]<thresharea)||(intern[i][j]==1)) area[i][j]=0;}
}

void identifyjunctions()
{   unsigned long i,j,il,jl,i1,j1,i2,j2;
    float max,max1;

	for (j=1;j<=lattice_size_y;j++)
     for (i=1;i<=lattice_size_x;i++)
      if (area[i][j]>0)
	   {//identify pixel in direction of upstream trib 1 (largest trib)
        max=0;
        if ((area[iup[i]][j]>max)&&(area[i][j]>area[iup[i]][j])) {max=area[iup[i]][j];i1=iup[i];j1=j;}
        if ((area[idown[i]][j]>max)&&(area[i][j]>area[idown[i]][j])) {max=area[idown[i]][j];i1=idown[i];j1=j;}
		if ((area[i][jup[j]]>max)&&(area[i][j]>area[i][jup[j]])) {max=area[i][jup[j]];i1=i;j1=jup[j];}
		if ((area[i][jdown[j]]>max)&&(area[i][j]>area[i][jdown[j]])) {max=area[i][jdown[j]];i1=i;j1=jdown[j];}
		//identify pixel in direction of upstream trib 2
		max1=max;max=0;
		if ((area[iup[i]][j]>max)&&(area[iup[i]][j]<max1)) {max=area[iup[i]][j];i2=iup[i];j2=j;}
		if ((area[idown[i]][j]>max)&&(area[idown[i]][j]<max1)) {max=area[idown[i]][j];i2=idown[i];j2=j;}
		if ((area[i][jup[j]]>max)&&(area[i][jup[j]]<max1)) {max=area[i][jup[j]];i2=i;j2=jup[j];}
		if ((area[i][jdown[j]]>max)&&(area[i][jdown[j]]<max1)) {max=area[i][jdown[j]];i2=i;j2=jdown[j];}
	    //if area of upstream trib 2 is above a threshold, define a junction
	    if ((max>0)&&(fabs((area[i1][j1]+area[i2][j2])/area[i][j])-1)<threshhillslopecontrib) mask[i][j]=1;}
    for (j=1;j<=lattice_size_y;j++)
     for (i=1;i<=lattice_size_x;i++)
      if ((area[i][j]>0)&&(mask[i][j]==0)&&(mask[iup[i]][j]==0)&&(mask[idown[i]][j]==0)&&(mask[i][jup[j]]==0)&&(mask[i][jdown[j]]==0)&&(mask[iup[i]][jup[j]]==0)&&(mask[idown[i]][jup[j]]==0)&&(mask[iup[i]][jdown[j]]==0)&&(mask[idown[i]][jdown[j]]==0))
	   {max=0;
        if ((area[iup[i]][j]>max)&&(area[i][j]>area[iup[i]][j])) {max=area[iup[i]][j];i1=iup[i];j1=j;}
        if ((area[idown[i]][j]>max)&&(area[i][j]>area[idown[i]][j])) {max=area[idown[i]][j];i1=idown[i];j1=j;}
		if ((area[i][jup[j]]>max)&&(area[i][j]>area[i][jup[j]])) {max=area[i][jup[j]];i1=i;j1=jup[j];}
		if ((area[i][jdown[j]]>max)&&(area[i][j]>area[i][jdown[j]])) {max=area[i][jdown[j]];i1=i;j1=jdown[j];}
		if ((area[iup[i]][jup[j]]>max)&&(area[i][j]>area[iup[i]][jup[j]])) {max=area[iup[i]][jup[j]];i1=iup[i];j1=jup[j];}
		if ((area[iup[i]][jdown[j]]>max)&&(area[i][j]>area[iup[i]][jdown[j]])) {max=area[iup[i]][jdown[j]];i1=iup[i];j1=jdown[j];}
		if ((area[idown[i]][jup[j]]>max)&&(area[i][j]>area[idown[i]][jup[j]])) {max=area[idown[i]][jup[j]];i1=idown[i];j1=jup[j];}
        if ((area[idown[i]][jdown[j]]>max)&&(area[i][j]>area[idown[i]][jdown[j]])) {max=area[idown[i]][jdown[j]];i1=idown[i];j1=jdown[j];}
		//identify pixel in direction of upstream trib 2
		max1=max;max=0;
		if ((area[iup[i]][j]>max)&&(area[iup[i]][j]<max1)) {max=area[iup[i]][j];i2=iup[i];j2=j;}
		if ((area[idown[i]][j]>max)&&(area[idown[i]][j]<max1)) {max=area[idown[i]][j];i2=idown[i];j2=j;}
		if ((area[i][jup[j]]>max)&&(area[i][jup[j]]<max1)) {max=area[i][jup[j]];i2=i;j2=jup[j];}
		if ((area[i][jdown[j]]>max)&&(area[i][jdown[j]]<max1)) {max=area[i][jdown[j]];i2=i;j2=jdown[j];}
		if ((area[iup[i]][jup[j]]>max)&&(area[iup[i]][jup[j]]<max1)) {max=area[iup[i]][jup[j]];i2=iup[i];j2=jup[j];}
		if ((area[iup[i]][jdown[j]]>max)&&(area[iup[i]][jdown[j]]<max1)) {max=area[iup[i]][jdown[j]];i2=iup[i];j2=jdown[j];}
		if ((area[idown[i]][jup[j]]>max)&&(area[idown[i]][jup[j]]<max1)) {max=area[idown[i]][jup[j]];i2=idown[i];j2=jup[j];}
        if ((area[idown[i]][jdown[j]]>max)&&(area[idown[i]][jdown[j]]<max1)) {max=area[idown[i]][jdown[j]];i2=idown[i];j2=jdown[j];}
	    //if area of upstream trib 2 is above a threshold, define a junction
		if ((max>0)&&(fabs((area[i1][j1]+area[i2][j2])/area[i][j])-1)<threshhillslopecontrib) mask[i][j]=1;}
}

void analyzejunctiongeometry()
{   unsigned long i,j,n,ilast,jlast,il,jl,i1,j1,i2,j2,i3,j3,i2k,j2k,i1k,j1k;
    float drop,max,max1,slope1,slope2,slope3,slope13,slope23,dropl,dist,diag,angle1,angle2,angle3,angle13,angle23,meanangle1,meanangle2,meanangle3,meanangle13,meanangle23;
	FILE *fp0;
	
	fp0=fopen("./conusjunctions.txt","w");	
	for (j=1;j<=lattice_size_y;j++)
     for (i=1;i<=lattice_size_x;i++)
	  if (mask[i][j]==1)
	   {//identify pixel in direction of upstream trib 2
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
		i3=i;j3=j;
		dist=0;
		max=0;
	    if (area[iup[i3]][j3]>max) {max=area[iup[i3]][j3];il=iup[i3];jl=j3;diag=1;}
	    if (area[idown[i3]][j3]>max) {max=area[idown[i3]][j3];il=idown[i3];jl=j3;diag=1;}
		if (area[i3][jup[j3]]>max) {max=area[i3][jup[j3]];il=i3;jl=jup[j3];diag=1;}
		if (area[i3][jdown[j3]]>max) {max=area[i3][jdown[j3]];il=i3;jl=jdown[j3];diag=1;}
		if (area[iup[i3]][jup[j3]]>max) {max=area[iup[i3]][jup[j3]];il=iup[i3];jl=jup[j3];diag=sqrt2;}
		if (area[iup[i3]][jdown[j3]]>max) {max=area[iup[i3]][jdown[j3]];il=iup[i3];jl=jdown[j3];diag=sqrt2;}
		if (area[idown[i3]][jup[j3]]>max) {max=area[idown[i3]][jup[j3]];il=idown[i3];jl=jup[j3];diag=sqrt2;}
		if (area[idown[i3]][jdown[j3]]>max) {max=area[idown[i3]][jdown[j3]];il=idown[i3];jl=jdown[j3];diag=sqrt2;}
		drop=topo[i][j]-topo[i3][j3];
		dist+=diag*deltax;
		i3=il;j3=jl;
		count=0;
		while ((drop<threshdrop)&&(i3>1)&&(i3<lattice_size_x)&&(j3>1)&&(j3<lattice_size_y)&&(count<10000))
		 {max=0;
	      count++;
	      if (area[iup[i3]][j3]>max) {max=area[iup[i3]][j3];il=iup[i3];jl=j3;diag=1;}
	      if (area[idown[i3]][j3]>max) {max=area[idown[i3]][j3];il=idown[i3];jl=j3;diag=1;}
		  if (area[i3][jup[j3]]>max) {max=area[i3][jup[j3]];il=i3;jl=jup[j3];diag=1;}
		  if (area[i3][jdown[j3]]>max) {max=area[i3][jdown[j3]];il=i3;jl=jdown[j3];diag=1;}
		  if (area[iup[i3]][jup[j3]]>max) {max=area[iup[i3]][jup[j3]];il=iup[i3];jl=jup[j3];diag=sqrt2;}
		  if (area[iup[i3]][jdown[j3]]>max) {max=area[iup[i3]][jdown[j3]];il=iup[i3];jl=jdown[j3];diag=sqrt2;}
		  if (area[idown[i3]][jup[j3]]>max) {max=area[idown[i3]][jup[j3]];il=idown[i3];jl=jup[j3];diag=sqrt2;}
		  if (area[idown[i3]][jdown[j3]]>max) {max=area[idown[i3]][jdown[j3]];il=idown[i3];jl=jdown[j3];diag=sqrt2;}
		  i3=il;j3=jl;
		  drop=topo[i][j]-topo[i3][j3];
		  dist+=diag*deltax;}
	    if (drop>=threshdrop) slope3=(topo[i][j]-topo[i3][j3])/dist; else slope3=0;
	    //angles are defined with positive values increasing counterclockwise from E, -180 to 180 degrees
	    if (i3!=i) angle3=-180/PI*atan2(1.0*j3-j,1.0*i3-i); else {if (j3>j) angle3=-90; else angle3=90;}
		//search upslope along tributary 2 to compute slope and direction
		dist=0;
	    i2=i2k;j2=j2k;
		drop=topo[i2][j2]-topo[i][j];
	    ilast=i2;jlast=j2;
		max=large;
		count=0;
		while ((max>0)&&(drop<threshdrop)&&(i2>1)&&(i2<lattice_size_x)&&(j2>1)&&(j2<lattice_size_y)&&(area[i2][j2]>0)&&(count<10000))
		 {max=0;
	      count++;
	      if ((area[iup[i2]][j2]>max)&&(area[iup[i2]][j2]<area[ilast][jlast])) {max=area[iup[i2]][j2];il=iup[i2];jl=j2;diag=1;}
	      if ((area[idown[i2]][j2]>max)&&(area[idown[i2]][j2]<area[ilast][jlast])) {max=area[idown[i2]][j2];il=idown[i2];jl=j2;diag=1;}
	      if ((area[i2][jup[j2]]>max)&&(area[i2][jup[j2]]<area[ilast][jlast])) {max=area[i2][jup[j2]];il=i2;jl=jup[j2];diag=1;}
		  if ((area[i2][jdown[j2]]>max)&&(area[i2][jdown[j2]]<area[ilast][jlast])) {max=area[i2][jdown[j2]];il=i2;jl=jdown[j2];diag=1;}
		  if ((area[iup[i2]][jup[j2]]>max)&&(area[iup[i2]][jup[j2]]<area[ilast][jlast])) {max=area[iup[i2]][jup[j2]];il=iup[i2];jl=jup[j2];diag=sqrt2;}
		  if ((area[iup[i2]][jdown[j2]]>max)&&(area[iup[i2]][jdown[j2]]<area[ilast][jlast])) {max=area[iup[i2]][jdown[j2]];il=iup[i2];jl=jdown[j2];diag=sqrt2;}
		  if ((area[idown[i2]][jup[j2]]>max)&&(area[idown[i2]][jup[j2]]<area[ilast][jlast])) {max=area[idown[i2]][jup[j2]];il=idown[i2];jl=jup[j2];diag=sqrt2;}
          if ((area[idown[i2]][jdown[j2]]>max)&&(area[idown[i2]][jdown[j2]]<area[ilast][jlast])) {max=area[idown[i2]][jdown[j2]];il=idown[i2];jl=jdown[j2];diag=sqrt2;}
		  if (max>small) {ilast=i2;jlast=j2;i2=il;j2=jl;}
		  dist+=diag*deltax;
		  drop=topo[i2][j2]-topo[i][j];}
	    if (drop>=threshdrop) slope2=(topo[i2][j2]-topo[i][j])/dist; else slope2=0; 
		if (i2!=i) angle2=-180/PI*atan2(1.0*j-j2,1.0*i-i2); else {if (j2>j) angle2=90; else angle2=-90;}
		dist=0;
	    i1=i1k;j1=j1k;
		drop=topo[i1][j1]-topo[i][j];
	    ilast=i1;jlast=j1;
	    max=large;
		count=0;
		while ((max>0)&&(drop<threshdrop)&&(i1>1)&&(i1<lattice_size_x)&&(j1>1)&&(j1<lattice_size_y)&&(area[i1][j1]>0)&&(count<10000))
	 	 {max=0;
	      count++;
	      if ((area[iup[i1]][j1]>max)&&(area[iup[i1]][j1]<area[ilast][jlast])) {max=area[iup[i1]][j1];il=iup[i1];jl=j1;diag=1;}
		  if ((area[idown[i1]][j1]>max)&&(area[idown[i1]][j1]<area[ilast][jlast])) {max=area[idown[i1]][j1];il=idown[i1];jl=j1;diag=1;}
		  if ((area[i1][jup[j1]]>max)&&(area[i1][jup[j1]]<area[ilast][jlast])) {max=area[i1][jup[j1]];il=i1;jl=jup[j1];diag=1;}
		  if ((area[i1][jdown[j1]]>max)&&(area[i1][jdown[j1]]<area[ilast][jlast])) {max=area[i1][jdown[j1]];il=i1;jl=jdown[j1];diag=1;}
		  if ((area[iup[i1]][jup[j1]]>max)&&(area[iup[i1]][jup[j1]]<area[ilast][jlast])) {max=area[iup[i1]][jup[j1]];il=iup[i1];jl=jup[j1];diag=sqrt2;}
		  if ((area[iup[i1]][jdown[j1]]>max)&&(area[iup[i1]][jdown[j1]]<area[ilast][jlast])) {max=area[iup[i1]][jdown[j1]];il=iup[i1];jl=jdown[j1];diag=sqrt2;}
		  if ((area[idown[i1]][jup[j1]]>max)&&(area[idown[i1]][jup[j1]]<area[ilast][jlast])) {max=area[idown[i1]][jup[j1]];il=idown[i1];jl=jup[j1];diag=sqrt2;}
          if ((area[idown[i1]][jdown[j1]]>max)&&(area[idown[i1]][jdown[j1]]<area[ilast][jlast])) {max=area[idown[i1]][jdown[j1]];il=idown[i1];jl=jdown[j1];diag=sqrt2;}
		  if (max>small) {ilast=i1;jlast=j1;i1=il;j1=jl;}
		  dist+=diag*deltax;
		  drop=topo[i1][j1]-topo[i][j];}
	    if (drop>=threshdrop) slope1=(topo[i1][j1]-topo[i][j])/dist; else slope1=0; 
		if (i1!=i) angle1=-180/PI*atan2(1.0*j-j1,1.0*i-i1); else {if (j1>j) angle1=90; else angle1=-90;}
	    angle13=fabs(angle1-angle3);
	    angle23=fabs(angle2-angle3);
	    if (angle13>180) angle13=fabs(360-angle13);
	    if (angle23>180) angle23=fabs(360-angle23);	
		meanangle13=fabs(angle[i1k][j1k]-angle[i][j]);
		meanangle23=fabs(angle[i2k][j2k]-angle[i][j]);
		if (meanangle13>180) meanangle13=fabs(360-meanangle13);
	    if (meanangle23>180) meanangle23=fabs(360-meanangle23);	
		//next two lines output location of each junction, along-channel slope ratio, along-channel junction angle, basin slope ratio, basin junction angle
        if ((slope3>0)&&(slope1>0)&&(slope3<large)&&(slope1<large)&&(intern[i][j]==0)) fprintf(fp0,"%d %d %f %f %f %f %f\n",i,lattice_size_y-j+1,area[i][j],slope3/slope1,angle13,slope[i][j]/slope[i1k][j1k],meanangle13);
	    if ((slope3>0)&&(slope2>0)&&(slope3<large)&&(slope2<large)&&(intern[i][j]==0)) fprintf(fp0,"%d %d %f %f %f %f %f\n",i,lattice_size_y-j+1,area[i][j],slope3/slope2,angle23,slope[i][j]/slope[i2k][j2k],meanangle23);
		//next four lines compute mean along-channel junction angle versus along-channel slope ratio
		if (intern[i][j]==0)
		 {ind=(int)(numbins*slope3/slope1)+1;
		  if ((ind>=1)&&(ind<=numbins)&&(slope3>0)&&(slope1>0)&&(slope3<large)&&(slope1<large)) {avcount[ind]++;avangle[ind]+=cos(PI/180*angle13);}
		  ind=(int)(numbins*slope3/slope2)+1;
		  if ((ind>=1)&&(ind<=numbins)&&(slope3>0)&&(slope2>0)&&(slope3<large)&&(slope2<large)) {avcount[ind]++;avangle[ind]+=cos(PI/180*angle23);}
		  //next four lines compute mean basin junction angle versus basin slope ratio
		  ind=(int)(numbins*slope[i][j]/slope[i1k][j1k])+1;
		  if ((ind>=1)&&(ind<=numbins)&&(slope[i][j]>0)) {avcount2[ind]++;avangle2[ind]+=cos(PI/180*meanangle13);}
		  ind=(int)(numbins*slope[i][j]/slope[i2k][j2k])+1;
	      if ((ind>=1)&&(ind<=numbins)&&(slope[i][j]>0)) {avcount2[ind]++;avangle2[ind]+=cos(PI/180*meanangle23);}}}
	fclose(fp0);
}

int main()
{
     FILE *fp0;
	 unsigned long i,j;
     int dum;
	  
	 lattice_size_x=92401;
	 lattice_size_y=58059;
	 deltaxinkm=0.05;               // km
	 thresharea=0.1;                // km^2
	 deltax=deltaxinkm*1000;
	 oneoverdeltax=1.0/deltax;
	 threshhillslopecontrib=0.01;   // fraction
	 threshdrop=10;                 // m
	 threshfilled=10;               // m
	 numbins=30;
	 setupgridneighbors();
	 setupgrids();
	 fp0=fopen("./conussurfgeo.txt","r");
	 for (j=1;j<=lattice_size_y;j++)
	  for (i=1;i<=lattice_size_x;i++)
	   {fscanf(fp0,"%d",&dum);
        if ((dum==0)||(dum==9)) intern[i][j]=1; else intern[i][j]=0;}
     fclose(fp0);
     fp0=fopen("./conus.txt","r");
	 for (j=1;j<=lattice_size_y;j++)
	  for (i=1;i<=lattice_size_x;i++)
	   {fscanf(fp0,"%f",&topo[i][j]);
        if ((topo[i][j]>-1*large)&&(topo[i][j]<large)); else topo[i][j]=0;
	    topoold[i][j]=topo[i][j];
		mask[i][j]=0;}
     fclose(fp0);		
     for (j=1;j<=numbins;j++)
	  {avcount[j]=0;
       avcount2[j]=0;
       avangle[j]=0;
	   avangle2[j]=0;}	
	 hydrologiccorrection();  
	 computecontributingarea();
	 identifyjunctions();
	 analyzejunctiongeometry();
	 fp0=fopen("./conusmeanjunctions.txt","w");
	 //plots mean along-channel junction angle binned by along-channel slope ratio 
	 for (j=1;j<=numbins;j++) if (avcount[j]>0) fprintf(fp0,"%f %f %f\n",(j-0.5)/numbins,180/PI*acos(avangle[j]/avcount[j]),180/PI*acos((j-0.5)/numbins)); else printf("0.0 0.0\n");	 
	 //plots mean basin junction angle binned by basin slope ratio 
	 for (j=1;j<=numbins;j++) if (avcount2[j]>0) fprintf(fp0,"%f %f %f\n",(j-0.5)/numbins,180/PI*acos(avangle2[j]/avcount2[j]),180/PI*acos((j-0.5)/numbins)); else printf("0.0 0.0\n");	 
     //outputs grid of contributing area with junctions assigned the maximum contributing area 
     fclose(fp0);
	 fp0=fopen("./conusnetwork.pbm","w");	
	 fprintf(fp0,"P1\n92401 58059\n");
	 for (j=1;j<=lattice_size_y;j++)
	  {for (i=1;i<=lattice_size_x;i++)
	    {if ((area[i][j]>thresharea)&&(intern[i][j]==0)) fprintf(fp0,"1 "); else fprintf(fp0,"0 ");}
	   fprintf(fp0,"\n");}
     fclose(fp0);
	 return 0;
}		 