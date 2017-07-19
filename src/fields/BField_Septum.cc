// ********************************************************************
// $Id: BField_Septum.cc,v 3.0, 2011/1/19  G2P Exp $
// Implementation of the BField_Septum class.
//
//////////////////////////////////////////////////////////////////////

#include "BField_Septum.hh"
#include "G4UImanager.hh"

#include "g4hrsUsageManager.hh"

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include  <math.h>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"


//#define CREATE_MAP_NTUPLE 1
//#define BFIELD_SEPTUM_DEBUG 1



using namespace std;

BField_Septum* BField_Septum::fInstance=0;
BField_Septum* BField_Septum::GetInstance()
{ 
	if(!fInstance)  
	{
		//new BField_Septum();
		cout<<"BField_Septum is not initialized yet...exit...\n";
		exit(1);
	}
	return fInstance; 
}



//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

BField_Septum::BField_Septum(double pMomentum, const char *mapfile)
{
#ifdef BFIELD_SEPTUM_DEBUG
	if(BFIELD_SEPTUM_DEBUG > Global_Debug_Level)
		SetGlobalDebugLevel("BField_Septum::BField_Septum()", (int)BFIELD_SEPTUM_DEBUG);
#endif

	fInstance=this;

	this->ReadMap(mapfile);

	//update the Current Ratio if necessary
	if(fabs(pMomentum)>1.0E-5) 
	  SetMomentum(pMomentum);
	//SetMomentum(mDefaultMomentumL,mDefaultMomentumR);
	
}

BField_Septum::~BField_Septum()
{
	int i,j,k;
	for(i=0;i<mXNum;i++)
	{
		for (j=0;j<mYNum;j++)
		{
			for (k=0;k<mZNum;k++)
			{
				delete [] mBField[i][j][k];
			}
			delete [] mBField[i][j];
		}
		delete [] mBField[i];
	}
	delete mRotL2F;
	delete mRotF2L;
}

/////////////////////////////////////////////////////////////////////
bool BField_Septum::ReadMap(const char *filename)
{

	////////////////////
	//START INI BLOCK //
	////////////////////
	// This block was added to replace the BField_Septum.ini file		       
	// Initialization parameters are now read in from a header in the map itself   
	mNPara = 6;
	mInterpolateOutOfRange = 0;

	//FIXME: Set in macro!
	mDefaultMomentum = 1.063;

	//FIXME: Are these necessary?  Hard-coding for now to test, fix/delete later.
	mUseUniformB = 0;
	mUniformB[0] = 0.; mUniformB[1] = 0.; mUniformB[2] = 0.; 
	mFieldUnit = 1.0;	
	mRotAxis[0] = 0.; mRotAxis[1] = 0.; mRotAxis[2] = 0.;
	mRotAngle[0] = 0.; mRotAngle[1] = 0.; mRotAngle[2] = 0.;

	char strLog[1024];
	sprintf(strLog,"BField_Septum::ReadMap() is loading field map %s......\n",filename);
		printf(strLog);
	ifstream ins;
	int indexX=0,indexY=0,indexZ=0,col=0;
	double tempLine[10];
	ins.open(filename);
	if (ins.fail())
	{
		sprintf(strLog,"***ERROR! Can not open field map %s...exit!***\n",filename);
		exit(1);
		return false;
	}

	if(ins >> mXmin >> mXmax >> mStepX) {
		G4cout << "xmin, xmax, xstep = " << mXmin << ", " << mXmax << ", " << mStepX << endl;
	} else {
		G4cerr << "Error " << __FILE__ << " line " << __LINE__ 
		<< ": File " << filename << " contains unreadable header.  Aborting" << G4endl;
	    	exit(1);
	}
	
	if(ins >> mYmin >> mYmax >> mStepY) {
		G4cout << "ymin, ymax, ystep = " << mYmin << ", " << mYmax << ", " << mStepY << endl;
	} else {
		G4cerr << "Error " << __FILE__ << " line " << __LINE__ 
		<< ": File " << filename << " contains unreadable header.  Aborting" << G4endl;
	    	exit(1);
	}

	if(ins >> mZmin >> mZmax >> mStepZ) {
		G4cout << "zmin, zmax, zstep = " << mZmin << ", " << mZmax << ", " << mStepZ << endl;
	} else {
		G4cerr << "Error " << __FILE__ << " line " << __LINE__ 
		<< ": File " << filename << " contains unreadable header.  Aborting" << G4endl;
	    	exit(1);
	}

	if(ins >> mOrigin[0] >> mOrigin[1] >> mOrigin[2]) {
		G4cout << "septumX, septumY, septumZ = " << mOrigin[0] << ", " << mOrigin[1] << ", " << mOrigin[2] << endl;
	} else {
		G4cerr << "Error " << __FILE__ << " line " << __LINE__ 
		<< ": File " << filename << " contains unreadable header.  Aborting" << G4endl;
	    	exit(1);
	}

	if(ins >> mTarget[0] >> mTarget[1] >> mTarget[2]) {
		G4cout << "septumX, septumY, septumZ = " << mTarget[0] << ", " << mTarget[1] << ", " << mTarget[2] << endl;
	} else {
		G4cerr << "Error " << __FILE__ << " line " << __LINE__ 
		<< ": File " << filename << " contains unreadable header.  Aborting" << G4endl;
	    	exit(1);
	}
	int LineNum = 6;
	char tempname[256];

	//Burn empty line, column headers
	ins.getline(tempname,256);
	ins.getline(tempname,256);
	
	//////////////////
	//END INI BLOCK //
	//////////////////

	////////////////////////////
	//START CONSTRUCTOR BLOCK //
	////////////////////////////
	// This block is a cut/paste of commands that were previously executed in the constructor
	// These commands were run AFTER the BField_Septum.ini was read, but BEFORE the map was read
	// Therefore they are placed after the initialization header is read, but before the field map is read

	if(fabs(mDefaultMomentum)<1.0E-5 ) 
	{
		cerr<<"\n##DefaultMomentumL is invaid ("<<mDefaultMomentum
			<<") in BField_Septum.ini. I quit... \n";
		exit(1);
	}

	mDoShift=(mOrigin[0]*mOrigin[0]+mOrigin[1]*mOrigin[1]+mOrigin[2]*mOrigin[2]>0)?true:false;
	//mDoRotation=(fabs(mRotAngle[0])+fabs(mRotAngle[1])+fabs(mRotAngle[2])>0)?true:false;

	mDoRotation=false;
	// Default constructor. Gives a unit matrix
	CLHEP::HepRotation pRot[3];

	for(int i=0;i<3;i++)
	{
		if(mRotAxis[i]==0 || fabs(mRotAngle[i])<1.0E-10) continue;
		mDoRotation=true;
		if(mRotAxis[i]==1) pRot[i].rotateX(mRotAngle[i]);
		else if(mRotAxis[i]==2) pRot[i].rotateY(mRotAngle[i]);
		else if(mRotAxis[i]==3) pRot[i].rotateZ(mRotAngle[i]);
	}

	mRotL2F=new CLHEP::HepRotation();
	mRotF2L=new CLHEP::HepRotation();
	*mRotL2F=pRot[2]*pRot[1]*pRot[0];
	*mRotF2L=mRotL2F->inverse(); 
#ifdef BFIELD_SEPTUM_DEBUG
	double rad2deg=180./acos(-1.0);
	if(Global_Debug_Level>=1)
	{
		cout<<"\nCLHEP Lab2Field EulerAngles: phi="<<mRotL2F->getPhi()*rad2deg
			<<"  theta="<<mRotL2F->getTheta()*rad2deg
			<<"  psi="<<mRotL2F->getPsi()*rad2deg<<endl;
	}
#endif

	mXNum=int((mXmax-mXmin)/mStepX)+1;
	if (mXNum<2)
	{
		mXNum=3; //set the minimum to 3
		mXmax=mXmin+3*mStepX;
	}
	mYNum=int((mYmax-mYmin)/mStepY)+1;
	if (mYNum<2)
	{
		mYNum=3; //set the minimum to 3
		mYmax=mYmin+3*mStepY;
	}
	mZNum=int((mZmax-mZmin)/mStepZ)+1;
	if (mZNum<2)
	{
		mZNum=3; //set the minimum to 3
		mZmax=mZmin+3*mStepZ;
	}

	///////////allocate start//////////
	int i,j,k,l;
	mBField=new double ***[mXNum];
	for(i=0;i<mXNum;i++)
	{
		mBField[i]=new double **[mYNum];
		for (j=0;j<mYNum;j++)
		{
			mBField[i][j]=new double *[mZNum];
			for	(k=0;k<mZNum;k++)
			{
				mBField[i][j][k]=new double [mNPara];
				//initial the array
				for (l=0;l<mNPara;l++)    mBField[i][j][k][l]=0.0;
			}
		}
	}
	/////////////allocate end////////////


	//////////////////////////
	//END CONSTRUCTOR BLOCK //
	//////////////////////////

	while (!ins.eof())
	{
		ins.getline(tempname,256);	
		LineNum++;
		//check if it is an empty line	
		if(strlen(tempname) < size_t(2*mNPara-1)) continue;

		istringstream s1(tempname);
		for(col=0;col<mNPara;col++)    s1>>tempLine[col];
		
		//check for X, Y and Z
		//cout << "Check: " << mNPara << " " << tempLine[0] << " " << mXmin << " " << mXmax << " " << tempLine[1] << " " << mYmin << " " << mYmax << " " << tempLine[2] << " " << mZmin << " " << mZmax << endl;

		if (tempLine[0]>=mXmin && tempLine[0]<=mXmax && tempLine[1]>=mYmin && tempLine[1]<=mYmax &&
		tempLine[2]>=mZmin && tempLine[2]<=mZmax)
		  //if(1)
		{//store the value

			//in case there is an empty line, r=Btot=0.0
			if( sqrt(tempLine[0]*tempLine[0]+tempLine[1]*tempLine[1]+tempLine[2]*tempLine[2]) < 1.0E-8 && 
				sqrt(tempLine[3]*tempLine[3]+tempLine[4]*tempLine[4]+tempLine[5]*tempLine[5]) < 1.0E-8 ) 
			{
				cout<<"***Warning: ZERO field at  x="<<tempLine[0]<<"  y="<<tempLine[1]
				<<"  z="<<tempLine[2]<<endl;
				cout<<"***Line "<<LineNum<<" could be an empty line in the map "<<filename<<endl;
			}
			else 
			{
			  indexX=int((tempLine[0]-mXmin)/mStepX);
			  indexY=int((tempLine[1]-mYmin)/mStepY);
			  indexZ=int((tempLine[2]-mZmin)/mStepZ);
				for(col=0;col<mNPara;col++)
				{
					mBField[indexX][indexY][indexZ][col]=tempLine[col];
					//cout << indexX << " " << indexY << " " << indexZ << " " << mBField[indexX][indexY][indexZ][col] << endl;
					if(col>=3 && col<=5) mBField[indexX][indexY][indexZ][col]*=mFieldUnit; //change unit to tesla
				}
			}
		}
	}
	ins.close();

#ifdef CREATE_MAP_NTUPLE
	char rootfile[255];
	double x=0.0,y=0.0,z=0.0,Bx=0.0,By=0.0,Bz=0.0,r=0.0,Br=0.0,Btot=0.0;
	sprintf(rootfile,"%s.root",filename);
	TFile *file=new TFile(rootfile,"RECREATE");
	TTree *field=new TTree("field","field map");
	field->Branch("x",&x,"x/D");
	field->Branch("y",&y,"y/D");
	field->Branch("r",&r,"r/D");
	field->Branch("z",&z,"z/D");
	field->Branch("Bx",&Bx,"Bx/D");
	field->Branch("By",&By,"By/D");
	field->Branch("Br",&Br,"Br/D");
	field->Branch("Bz",&Bz,"Bz/D");
	field->Branch("Btot",&Btot,"Btot/D");


#ifdef BFIELD_SEPTUM_DEBUG
	//if(Global_Debug_Level>=4)
	if(1)
	{
		printf("The Magnetic field Map is:\n");		
		printf("       x         y        z       Bx       By       Bz        r        Br        Btot\n");
	}
#endif
	for (indexX=0;indexX<mXNum;indexX++)
	{
		for (indexY=0;indexY<mYNum;indexY++)
		{
			for (indexZ=0;indexZ<mZNum;indexZ++)
			{
				x=mBField[indexX][indexY][indexZ][0];
				y=mBField[indexX][indexY][indexZ][1];
				r=sqrt(x*x+y*y);
				z=mBField[indexX][indexY][indexZ][2];
				Bx=mBField[indexX][indexY][indexZ][3];
				By=mBField[indexX][indexY][indexZ][4];  
				Bz=mBField[indexX][indexY][indexZ][5];				
				Br=sqrt(Bx*Bx+By*By);
				Btot=sqrt(Bx*Bx+By*By+Bz*Bz);

				if(fabs(r)<1.0E-8 && fabs(z)<1.0E-08 && fabs(Btot)<1.0E-8) 
				{
					cout<<"***Warning: ZERO field at  x="<<x<<"  y="<<y<<"  z="<<z<<endl;
					cout<<"***Warning: There is an empty line in the map buffer..."<<endl;
				}
				else field->Fill();

#ifdef BFIELD_SEPTUM_DEBUG
				if(Global_Debug_Level>=4) 
				{					
					printf("%8.3f %8.3f %8.3f %8.6f %8.6f %8.6f %8.3f %8.6f %8.6f\n",
						x,y,z,Bx,By,Bz,r,Br,Btot);

				}
#endif
				//reset 
				x=y=z=Bx=By=Bz=r=Br=Btot=0.0;
			}
		}
	}
	file->Write();
	file->Close();
	//cout << "Close root file " << rootfile <<endl;  
	file->Delete();
#endif

#ifdef BFIELD_SEPTUM_DEBUG
	//if(Global_Debug_Level>=4)
	if(1)
	{
		printf("The Magnetic field Map is:\n");
		printf("      x         y        z        Bx        By        Bz\n");
		for (indexX=0;indexX<mXNum;indexX++)
		{
			for (indexY=0;indexY<mYNum;indexY++)
			{
				for (indexZ=0;indexZ<mZNum;indexZ++)
				{					
					for(col=0;col<3;col++)  printf(" %8.2f ",mBField[indexX][indexY][indexZ][col]);
					for(col=3;col<mNPara;col++)  printf(" %8.6f ",mBField[indexX][indexY][indexZ][col]);
					printf("\n");
				}
			}
		}
	}
#endif

	return true;
}

/////////////////////////////////////////////////////////////////////
bool BField_Septum::Interpolation(double Pos[3],double B[3],int n)
{/*////////////////////////////////////////
    //function: calculate the nth order interpolation
    //1)found out B[X0-Xn][Y0-Yn][Z0-Zn], (n+1)x(n+1)x(n+1) matrix, 
    //2)then interpolate by X, found 3x(n+1)x(n+1) matrix
    //3)then interpolate by Y, found 3x(n+1) matrix 
    //4)Finally interpolate by Z to get Bx,By,Bz 
    //Input:
    // Pos[3]: the position in x,y,z coordinate in cm,origin variable
    // n : calculate the Nth order
    //Output:
    // B[3]: the magnetic field in B_x,B_y,B_z in tesla
    *//////////////////////////////////////////

#ifdef BFIELD_SEPTUM_DEBU
    //I do not want to check these in release version, just make sure you provide the right parameters
    if(n<1) n=1; 
    if(n>4) n=4; 

    if (n>=mXNum || n>=mYNum|| n>=mZNum)
    {
        printf("\n**Error! Too few points to finish %dth order interpolation!**",n);
        printf("**Please change the config file to load more data points or use lower order Number!**\n");
        return false;
    }
#endif

    //just flip by z
    double x=Pos[0],y=Pos[1],z=Pos[2];

    int StartIndexX=0,StartIndexY=0,StartIndexZ=0;  //the first point to do the interpolation
    StartIndexX=int((x-mBField[0][0][0][0])/mStepX);
    StartIndexY=int((y-mBField[0][0][0][1])/mStepY);
    StartIndexZ=int((z-mBField[0][0][0][2])/mStepZ);


#ifdef BFIELD_SEPTUM_DEBUG
    if (StartIndexX<0 || StartIndexX>=mXNum)
    {
        if(BFIELD_SEPTUM_DEBUG>=1)
            printf("\n**Warning!!! X=%.4f is out of range [%.4f,%.4f) !!!**\n",
                    x,mBField[0][0][0][0],mBField[mXNum-1][0][0][0]);
        B[0]=B[1]=B[2]=0.0;return false;
    }
    if (StartIndexY<0 || StartIndexY>=mYNum)
    {
        if(BFIELD_SEPTUM_DEBUG>=1)
            printf("\n**Warning!!! Y=%.4f is out of range [%.4f,%.4f) !!!**\n",
                    y,mBField[0][0][0][1],mBField[0][mYNum-1][0][1]);
        B[0]=B[1]=B[2]=0.0;return false;
    }
    if (StartIndexZ<0 || StartIndexZ>=mZNum)
    {
        if(BFIELD_SEPTUM_DEBUG>=1)
            printf("\n**Warning!!! Z=%.4f is out of range [%.4f,%.4f) !!!**\n",
                    z,mBField[0][0][0][2],mBField[0][0][mZNum-1][2]);
        B[0]=B[1]=B[2]=0.0;return false;
    }
#endif

    if(mInterpolateOutOfRange==0)
    {
        if ((StartIndexX<0 || StartIndexX>=mXNum) || 
                (StartIndexY<0 || StartIndexY>=mYNum) || 
                (StartIndexZ<0 || StartIndexZ>=mZNum) )
        {
            B[0]=B[1]=B[2]=0.0;return false;
        }
    }

    //choose a best position to interpolate, considering the lower and upper limits
    if (StartIndexX>mXNum-1-n) StartIndexX=mXNum-1-n;
    else if (StartIndexX<0) StartIndexX=0;
    if (StartIndexY>mYNum-1-n) StartIndexY=mYNum-1-n;
    else if (StartIndexY<0) StartIndexY=0;
    if (StartIndexZ>mZNum-1-n) StartIndexZ=mZNum-1-n;
    else if (StartIndexZ<0) StartIndexZ=0;

    // interpolate by X first
    //put the n+1=5 to avoid new and delete, this can save a lot of time
    double temp,tempBX[3][5][5],tempBY[3][5];

    int iX,iY,iZ,ii,jj,kk,ll;
    for (ll=0;ll<=n;ll++)
    {
        iZ=StartIndexZ+ll;
        for (ii=0;ii<=n;ii++)
        {
            //initial
            tempBX[0][ii][ll]=0.;
            tempBX[1][ii][ll]=0.;
            tempBX[2][ii][ll]=0.;
            iY=StartIndexY+ii;		
            for (kk=0;kk<=n;kk++)
            {
                temp=1.;
                for (jj=0;jj<=n;jj++)
                {
                    iX=StartIndexX+jj;
                    if(jj!=kk)
                    {
                        temp*=(x-mBField[iX][iY][iZ][0])/(mBField[StartIndexX+kk][iY][iZ][0]-mBField[iX][iY][iZ][0]);
                    }
                }
                tempBX[0][ii][ll]+=temp*mBField[StartIndexX+kk][iY][iZ][3];
                tempBX[1][ii][ll]+=temp*mBField[StartIndexX+kk][iY][iZ][4];
                tempBX[2][ii][ll]+=temp*mBField[StartIndexX+kk][iY][iZ][5];
            }
        }
    }
    // interpolate by Y 
    iX=StartIndexX;
    for (ii=0;ii<=n;ii++)
    {
        //initial
        for (int i=0;i<3;i++) tempBY[i][ii]=0.;
        iZ=StartIndexZ+ii;
        for (kk=0;kk<=n;kk++)
        {
            temp=1.;
            for (jj=0;jj<=n;jj++)
            {
                iY=StartIndexY+jj;
                if(jj!=kk)
                {
                    temp*=(y-mBField[iX][iY][iZ][1])/(mBField[iX][StartIndexY+kk][iZ][1]-mBField[iX][iY][iZ][1]);
                }
            }
            tempBY[0][ii]+=temp*tempBX[0][kk][ii];			
            tempBY[1][ii]+=temp*tempBX[1][kk][ii];
            tempBY[2][ii]+=temp*tempBX[2][kk][ii];
        }
    }
    // interpolate by Z 
    //initial
    iX=StartIndexX;
    iY=StartIndexY;
    for (int i=0;i<3;i++) B[i]=0.;
    for (kk=0;kk<=n;kk++)
    {
        for(int i=0;i<3;i++ )temp=1.;
        for (jj=0;jj<=n;jj++)
        {
            iZ=StartIndexZ+jj;
            if(jj!=kk)
            {
                temp*=(z-mBField[iX][iY][iZ][2])/(mBField[iX][iY][StartIndexZ+kk][2]-mBField[iX][iY][iZ][2]);
            }
        }
        B[0]+=temp*tempBY[0][kk];			
        B[1]+=temp*tempBY[1][kk];
        B[2]+=temp*tempBY[2][kk];
    }
    return  true;
}

void BField_Septum::Rotate_Lab2Field(const double LabP[3],double FieldP[3])
{
	Hep3Vector pFieldP(LabP[0],LabP[1],LabP[2]);
	pFieldP.transform(*mRotL2F);
	FieldP[0]=pFieldP.x();
	FieldP[1]=pFieldP.y();
	FieldP[2]=pFieldP.z();
}

void BField_Septum::Transform_Lab2Field(const double LabP[3],double FieldP[3])
{
	for(int i=0;i<3;i++) FieldP[i]=LabP[i]-mOrigin[i]-mTarget[i];
	if(mDoRotation) Rotate_Lab2Field(FieldP,FieldP);
}


void BField_Septum::Rotate_Field2Lab(const double FieldP[3],double LabP[3])
{
	//the following 3 lines do the same thing, but the 3rd line will also change pFieldP
	//pLabP=(*mRotF2L)(pFieldP);
	//pLabP=mRotF2L->operator ()(pFieldP);
	//pLabP=pFieldP.transform(*mRotF2L);
	Hep3Vector pLabP(FieldP[0],FieldP[1],FieldP[2]);
	pLabP.transform(*mRotF2L);
	LabP[0]=pLabP.x();
	LabP[1]=pLabP.y();
	LabP[2]=pLabP.z();
}

void BField_Septum::Transform_Field2Lab(const double FieldP[3],double LabP[3])
{
	Rotate_Field2Lab(FieldP,LabP);
	LabP[0]+=mOrigin[0];
	LabP[1]+=mOrigin[1];
	LabP[2]+=mOrigin[2];
}


/////////////////////////////////////////////////////////////////////
bool BField_Septum::GetBField(double Pos[],double B[]){

  //input x,y,z in centimeter, return B field in Tesla
  int i;
  //G4cout << "fancy get" << G4endl;
  /*
    if(mUseUniformB==1)
    {
    if( abs( Pos[0] ) > 8.4 && abs( Pos[0] < 38.8 ) && abs( Pos[1] < 12.2 ) && abs( Pos[2] < 37.0 ) ){// Make sure you are in the septum - quite ideal...
    for (i=0;i<3;i++) B[i]=mUniformB[i];
    return true;
    }else{
    for (i=0;i<3;i++) B[i] = 0.0;
    return true;
    }
    }
  */
  
  double pPos[3],pB[3]={0,0,0},flag[3]={1.0,1.0,1.0};
  //shift and rotate the origin to the field coordinate
  if(mDoShift || mDoRotation) Transform_Lab2Field(Pos,pPos);
  else{
    for (i=0;i<3;i++) pPos[i]=Pos[i];

    if(mUseUniformB==1){
      if( abs( Pos[0] ) > 8.4 && abs( Pos[0] ) < 38.8 && abs( Pos[1] ) < 12.2 && abs( Pos[2] ) < 37.0 ){// Make sure you are in the septum - quite ideal...
	for (i=0;i<3;i++) B[i]=mUniformB[i];
	return true;
      }else{
	for (i=0;i<3;i++) B[i] = 0.0;
	return true;
      }
    }

  }
  
  
  //the map provide fields for the whole range of x and y, but only half of z
  //field is y direction. Need to flip for -z 
  //if z<0, flip Bz. Bx does not need to flip //Not true for Nickie's version! comment out
  //if (pPos[2]<0) {flag[2]*=-1.0;}
  
  //Note that the map is only covers z>0, I have to take the absolute values // This is no longer true for me... Nickie - 24 Aug, 2015
  //double pPosAtMap[]={pPos[0],pPos[1],fabs(pPos[2])};
  double pPosAtMap[]={pPos[0],pPos[1],pPos[2]};
  
  if(!Interpolation(pPosAtMap,pB,1)) {B[0]=B[1]=B[2]=0.0; return false;}
  
#ifdef BFIELD_SEPTUM_DEBUG
  if(Global_Debug_Level>=3){
    printf("Input position(x,y,z) in cm=(%f, %f, %f):==>\n",Pos[0],Pos[1],Pos[2]);
    printf("Map position(x,y,z) in cm=(%f, %f, %f);\n",pPos[0],pPos[1],pPos[2]);
    printf("The raw magnetic field in Tesla without apply Septum_CurrentRatio:\n");
    printf("(B_x=%7.5f, B_y=%7.5f, B_z=%7.5f);\n",
	   pB[0]*=flag[0],pB[1]*=flag[1],pB[2]*=flag[2]);
  }
#endif
  
  //apply real current ratio and the flag
  /*
  G4cout << "Nickie is double checking: " << G4endl;
  G4cout << mCurrentRatioL << " " << mCurrentRatioR << G4endl;
  G4cout << B[0] << " " << B[1] << " " << B[2] << G4endl;
  */
  for (i=0;i<3;i++) {
    //G4cout << flag[i] << G4endl;
    B[i]=pB[i]*(fMomentum/mDefaultMomentum)*flag[i]*1e3;
  }

#ifdef BFIELD_SEPTUM_DEBUG
  if(Global_Debug_Level>=2){
    printf("The magnetic field in Tesla after apply Septum_CurrentRatio: \n");
    printf("(B_x=%7.5f, B_y=%7.5f, B_z=%7.5f);\n",B[0],B[1],B[2]);
  }
#endif
  
  //rotate back to the Hall coordinate
  if(mDoRotation){
    Rotate_Field2Lab(B,B);
    
#ifdef BFIELD_SEPTUM_DEBUG
    if(Global_Debug_Level>=2){
      printf("The magnetic field in Tesla after apply Rotation: \n");
      printf("(B_x=%7.5f, B_y=%7.5f, B_z=%7.5f);\n",B[0],B[1],B[2]);
    }
#endif
  }
  
  return true;
}

/////////////////////////////////////////////////////////////////////
bool BField_Septum::GetBField(float fPos[],float fB[])
{//input x,y,z in centimeter
  G4cout << "simple get" << G4endl;
	bool status=false;
	double dPos[3],dB[3];
	int i;
	for ( i=0;i<3;i++)
	{
		dPos[i]=(double) fPos[i];
	}
	status=GetBField(dPos,dB);

	for ( i=0;i<3;i++)
	{
		fB[i]=(float) dB[i];
	}
	return status;
}

/////////////////////////////////////////////////////////////////////
//
//void CreateNtuple(const char *rootfile="septumfield.root",int verbose=0,
//				  double xmin=-43.0,double xmax=43.0, double xstep=1.0,
//				  double ymin=-12.0,double ymax=12.0, double ystep=1.0,
//				  double zmin=-100.0,double zmax=100.0, double zstep=1.0);
//
//Create root ntuple: if verbose != 0 will print steps out 
void BField_Septum::CreateNtuple(const char *rootfile,int verbose,double xmin,double xmax,double xstep,
								 double ymin,double ymax,double ystep,double zmin,double zmax,double zstep)
{
	double x,y,r,z,Bx,By,Br,Bz,Btot;

	TFile *file=new TFile(rootfile,"RECREATE");
	TTree *field=new TTree("field","field map");
	field->Branch("x",&x,"x/D");
	field->Branch("y",&y,"y/D");
	field->Branch("r",&r,"r/D");
	field->Branch("z",&z,"z/D");
	field->Branch("Bx",&Bx,"Bx/D");
	field->Branch("By",&By,"By/D");
	field->Branch("Br",&Br,"Br/D");
	field->Branch("Bz",&Bz,"Bz/D");
	field->Branch("Btot",&Btot,"Btot/D");

	//if(verbose)
	if(1)
	{
		printf("The Magnetic field Map is:\n");
		printf("       x         y        z       Bx       By       Bz        r        Br        Btot\n");
	}

	int pXNum=(int)ceil((xmax-xmin)/xstep)+1;
	int pYNum=(int)ceil((ymax-ymin)/ystep)+1;
	int pZNum=(int)ceil((zmax-zmin)/zstep)+1;

	double pPos[3],pB[3];

	for (int ix=0;ix<pXNum;ix++)
	{
		pPos[0]=x=xmin+xstep*ix;
		for (int iy=0;iy<pYNum;iy++)
		{
			pPos[1]=y=ymin+ystep*iy;
			r=sqrt(x*x+y*y);
			for (int iz=0;iz<pZNum;iz++)
			{
				pPos[2]=z=zmin+zstep*iz;

				//Get The Field value
				GetBField(pPos,pB);
				Bx=pB[0]; By=pB[1]; Bz=pB[2];  //in Tesla

				Br=sqrt(Bx*Bx+By*By);
				Btot=sqrt(Bx*Bx+By*By+Bz*Bz);

				if(fabs(Btot)<1.0E-13) 
				{
					cout<<"***Warning: Zero field at  x="<<x<<"  y="<<y<<"  z="<<z<<endl;
				}

				field->Fill();

				if(verbose)
				{
					printf("%8.3f %8.3f %8.3f %8.6f %8.6f %8.6f %8.3f %8.6f %8.6f\n",
						x,y,z,Bx,By,Bz,r,Br,Btot);
				}
			}
		}
	}

	file->Write();
	file->Close();
	file->Delete();

}
/////////////////////////////////////////////////////////////////////

