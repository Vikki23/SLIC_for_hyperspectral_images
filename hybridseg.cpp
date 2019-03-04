/**
 * mex -largeArrayDims hybridseg.cpp
 */

#include<mex.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

// ���弸�ֲ�ͬ�Ĺ��������Բ�ȷ�ʽ

//1. ŷʽ����
// ���룺��������ʽ��ʾ�ĸ߹�����Ԫ����������Ϊbands
double calEuclidean(double* v1, double* v2, int bands)
{
    int i;
    double dist = 0.0;
    for(i=0; i<bands; i++)
    {
        dist += (v1[i]-v2[i])*(v1[i]-v2[i]);
    }
    return sqrt(dist);
}

//2. SAD ���׽Ǿ��룬��Χ [0---PI]����
double calSAD(double* v1, double* v2, int bands)
{
    int i;
    double nor_v1, nor_v2, inner, theta;
    nor_v1 = nor_v2 = 0;
    inner = 0;
    for(i=0; i<bands; i++)
    {
        nor_v1 += v1[i]*v1[i];
        nor_v2 += v2[i]*v2[i];
        inner += v1[i]*v2[i];
    }
    inner = inner/(sqrt(nor_v1)*sqrt(nor_v2));
    theta = acos(inner);
    
    return theta; //����
}

//3. SID spectral information diffusion ������Ϣɢ�� ����Χ���[0---1]
double calSID(double* v1, double* v2, int bands)
{
    int i;
    double sum_v1, sum_v2;
    double* p1, *p2;
    double sid;
    
    p1 = (double*)mxMalloc(sizeof(double) * bands);
    p2 = (double*)mxMalloc(sizeof(double) * bands);
    
    sum_v1 = sum_v2 = 0;
    for(i=0; i<bands; i++)
    {
        sum_v1 += v1[i];
        sum_v2 += v2[i];
    }
    
    for(i=0; i<bands; i++)
    {
        p1[i] = v1[i]/sum_v1;
        p2[i] = v2[i]/sum_v2;
    }
    
    sid = 0;
    for(i=0; i<bands; i++)
    {
        sid += (p1[i]-p2[i])*log(p1[i]/p2[i]); 
    }
       
    mxFree(p1);
    mxFree(p2);
    
    return sid;
}

//4. SID * tan(SAD) ����SID��SAD
double cal_SAD_SID(double* v1, double* v2, int bands)
{
    double sad, sid;
    sad = calSAD(v1, v2, bands);
    sid = calSID(v1, v2, bands);
	return sid*tan(sad);  // �޸����� 2016-10-16 ��ʽ�޸�
    //return sad*tan(sid);
}

////////////////////////////////////////////////////////////////////////////////////

/*
�õ���ʼ��������ӣ��������Լ����ӵ���Ŀ
seedIndices����СΪnumseeds,���ÿ�����ӵ�����ֵ
*/
void getSeeds(int STEP, int width, int height, int* seedIndices, int* numseeds)
{
    const bool hexgrid = false;
	int n;
    int xstrips, ystrips;
    int xerr, yerr;
    double xerrperstrip,yerrperstrip;
    int xoff,yoff;
    int x,y;
    int xe,ye;
    int seedx,seedy;
    int i;

	xstrips = (int)(0.5+(double)(width)/(double)(STEP));
	ystrips = (int)(0.5+(double)(height)/(double)(STEP));
    
    xerr = width  - STEP*xstrips;if(xerr < 0){xstrips--;xerr = width - STEP*xstrips;}//�����δ�������Ԫ
    yerr = height - STEP*ystrips;if(yerr < 0){ystrips--;yerr = height- STEP*ystrips;}
    
	xerrperstrip = (double)(xerr)/(double)(xstrips); //�����ÿһ��strip�����
	yerrperstrip = (double)(yerr)/(double)(ystrips);
    
	xoff = STEP/2;
	yoff = STEP/2;
    
    n = 0;
	for( y = 0; y < ystrips; y++ )
	{
		ye = y*yerrperstrip;
		for( x = 0; x < xstrips; x++ )
		{
			xe = x*xerrperstrip;
            seedx = (x*STEP+xoff+xe);
            if(hexgrid){ seedx = x*STEP+(xoff<<(y&0x1))+xe; if(seedx >= width)seedx = width-1; }//for hex grid sampling
            seedy = (y*STEP+yoff+ye);
            i = seedy*width + seedx;
			seedIndices[n] = i;
			n++;
		}
	}
    *numseeds = n;
}


/*��SLIC���ڸ߹�������ʱ����ͳ��5D������Ҫ��չ��(bands+2)ά���ݣ�bands��ʾ����ά�Ĳ�����*/
void PerformSuperpixelSLIC(
						   double* spectral_vec, double* kseeds, int bands, //����������						   
						   double* kseedsx, double* kseedsy,                  //�ռ�������
						   int width, int height, 
						   int numseeds, int* klabels, int STEP, double compactness, int dist_type)
{
    int x1, y1, x2, y2; //������������߽磬��2*Step��С�������ռ�
	
	double dist;
	double distxy;
    int itr; // ��������default=10
    int n, x, y;
    int i;
    int ind;
    int r,c;
    int k;
    int b;
    
    int sz = width*height;     //sz��ʾ��Ԫ����
	const int numk = numseeds; //numk��ʾ������
	int offset = STEP;         //offset��ʾ�ָ��߳�
    
    
    double* distvec = (double*)mxMalloc(sizeof(double) * sz);
// 	double invwt = 1.0/((STEP/compactness)*(STEP/compactness));
    double invwt = compactness/STEP;

	double* clustersize = (double*)mxMalloc(sizeof(double) * numk);
    double* inv         = (double*)mxMalloc(sizeof(double) * numk);

    double* sigmabs     = (double*)mxMalloc(sizeof(double) * numk * bands);
    double* sigmax      = (double*)mxMalloc(sizeof(double) * numk);
    double* sigmay      = (double*)mxMalloc(sizeof(double) * numk);
    
	for( itr = 0; itr < 10; itr++ )  //����10��
	{
		for(i = 0; i < sz; i++){distvec[i] = DBL_MAX;}
     
		for( n = 0; n < numk; n++ ) //��ÿ������Ԫ������
		{
            //2 * STEP ��С�������ռ�
            x1 = kseedsx[n]-offset; if(x1 < 0) x1 = 0;
            y1 = kseedsy[n]-offset; if(y1 < 0) y1 = 0;
            x2 = kseedsx[n]+offset; if(x2 > width)  x2 = width;
            y2 = kseedsy[n]+offset; if(y2 > height) y2 = height;
            
			//����2*Step��С�������ռ���о�����
			for( y = y1; y < y2; y++ )
			{
				for( x = x1; x < x2; x++ )
				{
					i = y * width + x; //��������ص�
                    
                    //������׾���
                    // spectral_vec[i*bands]----spectral_vec[(i+1)*bands-1] ������Ԫ�Ĺ�������
                    // kseeds[n*bands]----kseeds[(n+1)*bands-1] ��ʾ����������Ԫ�Ĺ�������
                    switch(dist_type)
                    {
                        case 1:   //ŷʽ����
                            dist = calEuclidean(&spectral_vec[i*bands], &kseeds[n*bands], bands);
                            break;
                        case 2:   //SAD
                            dist = calSAD(&spectral_vec[i*bands], &kseeds[n*bands], bands);
                            break;
                        case 3:   //SID
                            dist = calSID(&spectral_vec[i*bands], &kseeds[n*bands], bands);
                            break;
                        case 4:   // SAD-SID
                            dist = cal_SAD_SID(&spectral_vec[i*bands], &kseeds[n*bands], bands);
                            break;
                    }
                    
                    // ����ռ����
					distxy = sqrt((x - kseedsx[n])*(x - kseedsx[n]) + (y - kseedsy[n])*(y - kseedsy[n]));
                    
					dist = dist*(1-compactness) + distxy*invwt;  //SLIC���빫ʽ
                    
					if(dist < distvec[i])
					{
						distvec[i] = dist;
						klabels[i] = n; //���ڵ� n �ž���
					}
				}
			}
		}

		//-----------------------------------------------------------------
		// Recalculate the centroid and store in the seed values
		//-----------------------------------------------------------------
		
        //���³�ʼ����������
        for(k = 0; k < numk; k++)
        {
			for(b = 0; b < bands ; b++)
				sigmabs[k*bands+b] = 0;
            sigmax[k] = 0;
            sigmay[k] = 0;
            clustersize[k] = 0;
        }
        
		ind = 0;
        for( r = 0; r < height; r++ )
        {
            for( c = 0; c < width; c++ )
            {
                if(klabels[ind] >= 0) //�ѻ��ֺõľ���
                {
					for(b=0; b<bands; b++)
						sigmabs[klabels[ind]*bands+b] += spectral_vec[ind*bands+b];
                    sigmax[klabels[ind]] += c;
                    sigmay[klabels[ind]] += r;
                    clustersize[klabels[ind]] += 1.0;
                }
                ind++;
            }
        }
        
		for( k = 0; k < numk; k++ )
		{
			if( clustersize[k] <= 0 ) clustersize[k] = 1; //���ܴ��ڵ�δ�ɹ�����ĵ�
			inv[k] = 1.0/clustersize[k];//computing inverse now to multiply, than divide later
		}
		//�����µľ������ģ��˴������ʣ�����ľ����������겻����ʵ���ص�Ĺ���ֵ
		for( k = 0; k < numk; k++ )
		{
			for(b=0; b<bands; b++)
				kseeds[k*bands+b] = sigmabs[k*bands+b]*inv[k];
			kseedsx[k] = sigmax[k]*inv[k];
			kseedsy[k] = sigmay[k]*inv[k];
		}
	}
    mxFree(sigmabs);
    mxFree(sigmax);
    mxFree(sigmay);
    mxFree(clustersize);
    mxFree(inv);
    mxFree(distvec);
}

/*��һЩδ������Ԫ���¾���*/
void EnforceSuperpixelConnectivity(int* labels, int width, int height, int numSuperpixels,int* nlabels, int* finalNumberOfLabels)
{
    int i,j,k;
    int n,c,count;
    int x,y;
    int ind;
    int oindex, adjlabel;
    int label;
	
	// 4 ����    
    const int dx4[4] = {-1,  0,  1,  0};
	const int dy4[4] = { 0, -1,  0,  1};
    
	const int sz = width*height;
    const int SUPSZ = sz/numSuperpixels;
    
	int* xvec = (int*)mxMalloc(sizeof(int)*SUPSZ*10);
	int* yvec = (int*)mxMalloc(sizeof(int)*SUPSZ*10);

	for( i = 0; i < sz; i++ ) nlabels[i] = -1;

    oindex = 0;
    adjlabel = 0;//adjacent label
    label = 0;

	for( j = 0; j < height; j++ )
	{
		for( k = 0; k < width; k++ )
		{
			if( 0 > nlabels[oindex] )
			{
				nlabels[oindex] = label;
				//--------------------
				// Start a new segment
				//--------------------
				xvec[0] = k;
				yvec[0] = j;
				//-------------------------------------------------------
				// Quickly find an adjacent label for use later if needed
				//-------------------------------------------------------
				for( n = 0; n < 4; n++ )
				{
					int x = xvec[0] + dx4[n];
					int y = yvec[0] + dy4[n];
					if( (x >= 0 && x < width) && (y >= 0 && y < height) )
					{
						int nindex = y*width + x;
						if(nlabels[nindex] >= 0) adjlabel = nlabels[nindex];
					}
				}
                
				count = 1;
				for( c = 0; c < count; c++ )
				{
					for( n = 0; n < 4; n++ )
					{
						x = xvec[c] + dx4[n];
						y = yvec[c] + dy4[n];
                        
						if( (x >= 0 && x < width) && (y >= 0 && y < height) )
						{
							int nindex = y*width + x;
                            
							if( 0 > nlabels[nindex] && labels[oindex] == labels[nindex] )
							{
								xvec[count] = x;
								yvec[count] = y;
								nlabels[nindex] = label;
								count++;
							}
						}
                        
					}
				}
				//-------------------------------------------------------
				// If segment size is less then a limit, assign an
				// adjacent label found before, and decrement label count.
				//-------------------------------------------------------
				if(count <= SUPSZ >> 2)
				{
					for( c = 0; c < count; c++ )
					{
                        ind = yvec[c]*width+xvec[c];
						nlabels[ind] = adjlabel;
					}
					label--;
				}
				label++;
			}
			oindex++;
		}
	}
	*finalNumberOfLabels = label;
    
	mxFree(xvec);
	mxFree(yvec);
}

/*
 * ���������8������1-4����ң�����ݣ�һά��������ʽ����
 * ��5���ǳ����ظ�������6���ǳ���mֵ��
 * ��7���ǲ��õ����ƶȲ�������,
 * ��8���ǿ��������ʽ���Ƿ���δ��������أ�
 * 
 * ���������4��, 
 * ��1���Ǿ����ǩ����2����ʵ�ʾ�������
 * ��3��4���Ǿ���������Ԫ����
 **/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])      
{
    int width, height, bands, sz;
	int numSuperpixels = 200;//default value
    double compactness = 0.1;//default value	
    int dist_type = 1;  //���������ͣ�1��ED��ȱʡֵ����2��SAD, 3:SID, 4:SAD-SID
    int handle_all = 1; //�Ƿ���������Ԫ��1������ȱʡֵ���� 2�������������ص�labels�д���δ�������Ԫ
    
	int i, ii, x, y, k ,t;
    
	int* klabels;
    int* clabels;
    
	int step;
    int* seedIndices;
    int numseeds;
    
	double* kseedsx;
    double* kseedsy;
    double* kseeds;  //���ÿ��������Ԫ�Ĺ�������
    
    int* outputNumSuperpixels;
    int* outlabels;
    int finalNumberOfLabels;
    double* outseedx;
    double* outseedy;
    
    int numelements;
    int* dims;
	double* spectral_vec; //���ԭʼͼ������    
    

    if (nrhs < 4) {
        mexErrMsgTxt("At least 4 argument is required.") ;
    } else if(nrhs > 8) {
        mexErrMsgTxt("Too many input arguments.");
    }
    if(nlhs!=4) {
        mexErrMsgIdAndTxt("SLIC:nlhs","Four outputs required, a labels, the number of labels and the seed positions.");
    }

    //��ȡ�������
    numelements = mxGetNumberOfElements(prhs[0]) ; //ң����������ά���ݣ���˴��ݽ�������Ҫ�ؽ�
    dims  = (int*)mxGetDimensions(prhs[0]) ;
    spectral_vec  = (double*)mxGetData(prhs[0]) ;//ң�����ݱ�����double���͵ģ��˴�ע����άתΪһάʱ��˳��
    //---------------------------
    width  = mxGetScalar(prhs[2]); //����
    height  = mxGetScalar(prhs[1]);//����
    bands  = mxGetScalar(prhs[3]); //������
    sz = width*height; //���ظ���
    numSuperpixels  = mxGetScalar(prhs[4]);
    compactness     = mxGetScalar(prhs[5]);
    dist_type = mxGetScalar(prhs[6]);
    handle_all = mxGetScalar(prhs[7]);
	    
    //---------------------------
    // Allocate memory
    //---------------------------    
    klabels = (int*)mxMalloc( sizeof(int)         * sz );//original k-means labels
    clabels = (int*)mxMalloc( sizeof(int)         * sz );//corrected labels after enforcing connectivity
    seedIndices = (int*)mxMalloc( sizeof(int)     * sz ); //������ӵ�����ֵ����ʵ����Ҫ��ô��
    
    //---------------------------
    // Find seeds
    //---------------------------
    step = sqrt((double)(sz)/(double)(numSuperpixels))+0.5; //��������ı߳�����������������
    getSeeds(step, width, height, seedIndices, &numseeds);	
    
    kseedsx    = (double*)mxMalloc( sizeof(double) * numseeds ) ;
    kseedsy    = (double*)mxMalloc( sizeof(double) * numseeds ) ;
    kseeds     = (double*)mxMalloc( sizeof(double) * numseeds * bands ) ; //ÿ�д���һ��������Ԫ
    for(k = 0; k < numseeds; k++)
    {
        kseedsx[k] = seedIndices[k]%width; //�к�
        kseedsy[k] = seedIndices[k]/width; //�к�
        for(t=0; t<bands; t++)
        {
            kseeds[k*bands+t] = spectral_vec[seedIndices[k]*bands+t];
        }
    }
    //---------------------------
    // Compute superpixels
    //---------------------------
    for(i=0; i<sz; i++)
    {
        klabels[i] = clabels[i] = -1;
    }
    
    PerformSuperpixelSLIC(spectral_vec, kseeds,bands, kseedsx,kseedsy,width,height,numseeds,klabels,step,compactness,dist_type);
    
    //---------------------------
    // Enforce connectivity
    //---------------------------
    if(handle_all == 1)
    {
        EnforceSuperpixelConnectivity(klabels,width,height,numseeds,clabels,&finalNumberOfLabels);
    }
    else 
    {
        for(i=0; i<sz; i++)
        {
            clabels[i] = klabels[i];
            finalNumberOfLabels = numseeds;
        }
    }
        
    //---------------------------
    // Assign output labels
    //---------------------------
    plhs[0] = mxCreateNumericMatrix(1, sz, mxINT32_CLASS,mxREAL);
    outlabels = (int*)mxGetData(plhs[0]);
    for(x = 0; x < sz; x++)
    {
       outlabels[x] = clabels[x]+1;  // ���ص��������ȸ�ʽ
    }
    //---------------------------
    // Assign number of labels/seeds
    //---------------------------
    plhs[1] = mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
    outputNumSuperpixels = (int*)mxGetData(plhs[1]);//gives a void*, cast it to int*
    *outputNumSuperpixels = finalNumberOfLabels;

	plhs[2] = mxCreateDoubleMatrix(1, finalNumberOfLabels, mxREAL);
    outseedx = (double*)mxGetData(plhs[2]);
    for(x = 0; x < finalNumberOfLabels; x++)//copying data from row-major C matrix to column-major MATLAB matrix (i.e. perform transpose)
    {
		outseedx[x] = kseedsx[x]+1; //add 1 because matlab is different with C for array index���к�        
    }

	plhs[3] = mxCreateDoubleMatrix(1, finalNumberOfLabels, mxREAL);
    outseedy = (double*)mxGetData(plhs[3]);
    for(x = 0; x < finalNumberOfLabels; x++)//copying data from row-major C matrix to column-major MATLAB matrix (i.e. perform transpose)
    {
		outseedy[x] = kseedsy[x]+1; //add 1 because matlab is different with C for array index���к�       
    }


    //---------------------------
    // Deallocate memory
    //---------------------------     
    mxFree(klabels);
    mxFree(clabels);
    mxFree(seedIndices);
    mxFree(kseedsx);
    mxFree(kseedsy);
    mxFree(kseeds);
}
