/*  ObjCryst++ Object-Oriented Crystallographic Library
    (c) 2000-2002 Vincent Favre-Nicolin vincefn@users.sourceforge.net
        2000-2001 University of Geneva (Switzerland)

    This program is free software; you can redistibute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#include "ObjCryst/Quirks/VFNStreamFormat.h"

#include "ObjCryst/RefinableObj/LSQNumObj.h"

#ifdef __WX__CRYST__
   #include "ObjCryst/wxCryst/wxLSQ.h"
#endif

#include "newmat/newmatap.h" //for SVD decomposition
#include "newmat/newmatio.h"

#ifdef use_namespace
using namespace NEWMAT;
#endif
using namespace std;

#include <iomanip>
//#include <cstring>
#include <stdexcept>
#include <fstream>

// ----------------------- GSVD ---

/* Generalised SVD decomposition of matrices A and B.
 *
 *   A = U * Sigma1 * [0, R] * Q'
 *   B = V * Sigma2 * [0, R] * Q'
 *
 *   eR = [0, R] ;
 *   r = rank([A; B]) ; k = r-l ; l = rank(B) ;
 */
void GSVD(const Matrix &A, const Matrix &B,
	  Matrix &Sigma1, Matrix &Sigma2,
	  Matrix &U, Matrix &V, Matrix &eR, Matrix &Q,
	  int &K, int &L, vector<int> &ind);

// ----------------------- rest of the file ---

namespace ObjCryst
{

LSQNumObj::LSQNumObj(string objName)
#ifdef __WX__CRYST__
:mpWXCrystObj(0)
#endif
{
   mDampingFactor=1.;
   mSaveReportOnEachCycle=false;
   mName=objName;
   mSaveFileName="LSQrefinement.save";
   mR=0;
   mRw=0;
   mChiSq=0;
   mChiSqReg=0; // Zdenek
   mRex=0; // Zdenek
   mStopAfterCycle=false;
}

LSQNumObj::~LSQNumObj()
{
   #ifdef __WX__CRYST__
   this->WXDelete();
   #endif
}

const string& LSQNumObj::GetClassName()const // Zdenek
{
	static const string className = "LSQNumObj";
	return className;
}

const string& LSQNumObj::GetName()const // Zdenek
{
	return mName;
}

void LSQNumObj::SetParIsFixed(const string& parName,const bool fix)
{
   if(mRefParList.GetNbPar()==0) this->PrepareRefParList();
   mRefParList.SetParIsFixed(parName,fix);
}
void LSQNumObj::SetParIsFixed(const RefParType *type,const bool fix)
{
   if(mRefParList.GetNbPar()==0) this->PrepareRefParList();
   mRefParList.SetParIsFixed(type,fix);
}

void LSQNumObj::SetParIsFixed(RefinablePar &par,const bool fix)
{
   if(mRefParList.GetNbPar()==0) this->PrepareRefParList();
   mRefParList.GetPar(par.GetPointer()).SetIsFixed(fix);
}

void LSQNumObj::SetParIsFixed(RefinableObj &obj,const bool fix)

{
   if(mRefParList.GetNbPar()==0) this->PrepareRefParList();
   for(unsigned int i=0;i<obj.GetNbPar();++i)
      this->SetParIsFixed(obj.GetPar(i),fix);
}

void LSQNumObj::UnFixAllPar()
{
   if(mRefParList.GetNbPar()==0) this->PrepareRefParList();
   mRefParList.UnFixAllPar();
}
   
void LSQNumObj::SetParIsUsed(const string& parName,const bool use)
{
   if(mRefParList.GetNbPar()==0) this->PrepareRefParList();
   mRefParList.SetParIsUsed(parName,use);
}
void LSQNumObj::SetParIsUsed(const RefParType *type,const bool use)
{
   if(mRefParList.GetNbPar()==0) this->PrepareRefParList();
   mRefParList.SetParIsUsed(type,use);
}

void LSQNumObj::Refine (int nbCycle,bool useLevenbergMarquardt,
                        const bool silent, const bool callBeginEndOptimization,
                        const float minChi2var)
{
  TAU_PROFILE("LSQNumObj::Refine()","void ()",TAU_DEFAULT);
  if(callBeginEndOptimization) this->BeginOptimization();
  mObs=this->GetLSQObs();
  mWeight=this->GetLSQWeight();
   
  bool terminateOnDeltaChi2=false;
  if(nbCycle<0)
  {
    nbCycle=-nbCycle;
    terminateOnDeltaChi2=true;
  }

  if(!silent) cout << "LSQNumObj::Refine():Beginning "<<endl;
 //Prepare for refinement (get non-fixed parameters)
  if(mRefParList.GetNbPar()==0) this->PrepareRefParList();
  mRefParList.PrepareForRefinement();
  if(!silent) mRefParList.Print();
  if(mRefParList.GetNbPar()==0) throw ObjCrystException("LSQNumObj::Refine():no parameter to refine !");

 //variables
  long nbVar=mRefParList.GetNbParNotFixed();
  const long nbObs=mObs.numElements();
  CrystVector_REAL calc,calc0,calc1,tmpV1,tmpV2;
  CrystMatrix_REAL M(nbVar,nbVar);
  CrystMatrix_REAL N(nbVar,nbVar);
  CrystMatrix_REAL P(nbVar,nbVar); // Zdenek ... solution projector matrix for constraints
  CrystVector_REAL B(nbVar);
  CrystMatrix_REAL designMatrix(nbVar,nbObs);
  CrystVector_REAL deltaVar(nbVar);
  long i,j,k;
  REAL R_ini,Rw_ini;
  REAL *pTmp1,*pTmp2;
  int nbConvReached = 0; // How many times convergence condition was satisfied ? // Zdenek
  bool bMarquardtConverged = false; // Flag of Levenberg-Marquardt convergence
  
  REAL marquardt=1e-2;
  const REAL marquardtMult=4.;
 //initial Chi^2, needed for Levenberg-Marquardt
  {
    calc=this->GetLSQCalc();
    tmpV1 = mObs;
    tmpV1 -= calc;
    tmpV1 *= tmpV1;
    tmpV1 *= mWeight;
    mChiSq=tmpV1.sum();
    mChiSqReg=this->GetRegularizationChiSq();
    mChiSq += mChiSqReg;
  }
 //initial value of expected R-factor // Zdenek
  {
    tmpV1 = mObs;
    tmpV1 *= tmpV1;
    tmpV1 *= mWeight; // TODO:: Include constraints.
    mRex = sqrt( (nbObs-nbVar) / tmpV1.sum() );
  } 
 //store old values
  mIndexValuesSetInitial=mRefParList.CreateParamSet("LSQ Refinement-Initial Values");
  mIndexValuesSetLast=mRefParList.CreateParamSet("LSQ Refinement-Last Cycle Values");

  mStopAfterCycle = false; // Zdenek
	
 //refine
  int cycle=1; // Zdenek
  for( ; cycle <=nbCycle;cycle++)
  {
    const REAL ChisSqPreviousCycle=mChiSq;
    mRefParList.SaveParamSet(mIndexValuesSetLast);// end of last cycle
    if(!silent) cout << "LSQNumObj::Refine():Cycle#"<< cycle <<endl;
    if(!silent) cout << "LSQNumObj::Refine():Computing initial values" <<endl;
		
  LSQNumObj_Refine_Restart_MarquardtUpdate: //Used in case of update solution for Marquardt (LM-factor set to zero)

   //initial value of function
    calc0=this->GetLSQCalc();
   //R
    tmpV1 =  mObs;
    tmpV1 -= calc0;
    tmpV1 *= tmpV1;
    tmpV2 =  mObs;
    tmpV2 *= mObs;
    R_ini=sqrt(tmpV1.sum()/tmpV2.sum());
   //Rw
    tmpV1 *= mWeight;
    tmpV2 *= mWeight;
    Rw_ini=sqrt(tmpV1.sum()/tmpV2.sum());
   //derivatives
    if(!silent) cout << "LSQNumObj::Refine():Computing derivatives" <<endl; // Zdenek
    designMatrix=0.;
    pTmp2=designMatrix.data();
    //cout <<"obs:"<<FormatHorizVector<REAL>(calc0,10,8);
    //cout <<"calc:"<<FormatHorizVector<REAL>(mObs,10,8);
    //cout <<"weight:"<<FormatHorizVector<REAL>(mWeight,10,8);
    for(i=0;i<nbVar;i++)
    {
      //:NOTE: Real design matrix is the transposed of the one computed here
      //if(!silent) cout << "........." << mRefParList.GetParNotFixed(i).GetName() <<endl;
      
      tmpV1=this->GetLSQDeriv(mRefParList.GetParNotFixed(i));
      pTmp1=tmpV1.data();
      //cout <<"deriv#"<<i<<":"<<FormatHorizVector<REAL>(tmpV1,10,8);
      for(j=0;j<nbObs;j++) *pTmp2++ = *pTmp1++;
    }
    //cout << designMatrix;

  LSQNumObj_Refine_Restart: //Used in case of singular matrix or for Marquardt
      
    if(!silent) cout << "LSQNumObj::Refine():Computing M and B Matrices" <<endl; // Zdenek
   //Calculate M and B matrices
    tmpV1.resize(nbObs);
    tmpV2.resize(nbObs);
    for(i=0;i<nbVar;i++) 
    {
      for(k=0;k<nbObs;k++) tmpV1(k)=designMatrix(i,k);
      tmpV1 *= mWeight;
      for(j=0;j<nbVar;j++)
      {
	for(k=0;k<nbObs;k++) tmpV2(k)=designMatrix(j,k);
	tmpV2 *= tmpV1;
	M(j,i)= tmpV2.sum();
	//M(i,j)=total(designMatrix.row(i)*designMatrix.row(j)*weight);
      }
            
      tmpV2=mObs;
      tmpV2 -= calc0;
      tmpV2 *= tmpV1;
      B(i)=tmpV2.sum();//total((obs-calc0)*weight*designMatrix.row(i))
    }
		
    // Zdenek (begin)
   // Check for singular values
    /* // Zdenek - Done later, but I think it can be done here
       for(i=0;i<nbVar;i++)
       {
       if( M(i,i) < 1e-20) //:TODO: Check what value to use as a limit
       {  
       if(!silent) cout << "LSQNumObj::Refine() Singular parameter !";
				if(!silent) cout << "(null derivate in all points) : ";
				if(!silent) cout << mRefParList.GetParNotFixed(i).GetName() << endl;
				if(!silent) cout << "LSQNumObj::Refine(): Automatically fixing parameter";
				if(!silent) cout << " and re-start cycle..";
				if(!silent) cout << endl;
				mRefParList.GetParNotFixed(i).SetIsFixed(true);
				mRefParList.PrepareForRefinement();
				nbVar=mRefParList.GetNbParNotFixed();
				if(nbVar<=1) throw ObjCrystException("LSQNumObj::Refine(): not enough (1) parameters after fixing one...");
				N.resize(nbVar,nbVar);
				deltaVar.resize(nbVar);

				//Just remove the ith line in the design matrix
				REAL *p1=designMatrix.data();
				const REAL *p2=designMatrix.data();
				p1 += i*nbObs;
				p2 += i*nbObs+nbObs;
				for(long j=i*nbObs+nbObs;j<nbObs*nbVar;j++) *p1++ = *p2++;
				designMatrix.resizeAndPreserve(nbVar,nbObs);
              
        M.resize(nbVar,nbVar);
        B.resize(nbVar);
   			
   			// recalculate expected R-factor
   			mRex *= sqrt((nbObs-nbVar)/(nbObs-(nbVar+1.))); // TODO:: Include constraints
        
        goto LSQNumObj_Refine_Restart;
			}
			}
		*/

   //Used in case of change of LevenBerg-Marquardt factor when solution was not updated (LM-factor was increased)
  LSQNumObj_Refine_Restart_MarquardtNotUpdate:

    // Zdenek (end)
      
   //Apply LevenBerg-Marquardt factor
    if(true==useLevenbergMarquardt)
      {
	const REAL lmfact=1.+marquardt;
	for(i=0;i<nbVar;i++) M(i,i) *= lmfact;
      }
   // Check for singular values
    for(i=0;i<nbVar;i++)
    {
      if( (M(i,i) < 1e-20)||(ISNAN_OR_INF(M(i,i)))) //:TODO: Check what value to use as a limit
      {  
	if(!silent) cout << "LSQNumObj::Refine() Singular parameter !";
	if(!silent) cout << "(null derivate in all points) : "<<M(i,i)<<":";
	if(!silent) cout << mRefParList.GetParNotFixed(i).GetName() << endl;
	/*
	  if(!silent)
               {
                  for(i=0;i<nbVar;i++)
                  {
                     tmpV1=this->GetLSQDeriv(mRefParList.GetParNotFixed(i));
                     cout <<"deriv#"<<i<<":"<<FormatHorizVector<REAL>(tmpV1,10,8);
                  }
               }
	*/
	if(!silent) cout << "LSQNumObj::Refine(): Automatically fixing parameter";
	if(!silent) cout << " and re-start cycle..";
	if(!silent) cout << endl;
	mRefParList.GetParNotFixed(i).SetIsFixed(true);
	mRefParList.PrepareForRefinement();
	nbVar=mRefParList.GetNbParNotFixed();
	if(nbVar<=1)
	{
	  mRefParList.RestoreParamSet(mIndexValuesSetInitial);
	  if(callBeginEndOptimization) this->EndOptimization();
	  if(!silent) mRefParList.Print();
	  throw ObjCrystException("LSQNumObj::Refine(): not enough (1) parameters after fixing one...");
	}
	N.resize(nbVar,nbVar);
	deltaVar.resize(nbVar);

       //Just remove the ith line in the design matrix
	REAL *p1=designMatrix.data();
	const REAL *p2=designMatrix.data();
	p1 += i*nbObs;
	p2 += i*nbObs+nbObs;
	for(long j=i*nbObs+nbObs;j<nbObs*nbVar;j++) *p1++ = *p2++;
	designMatrix.resizeAndPreserve(nbVar,nbObs);
               
               //:TODO: Make this work...
               /* 
               //Remove ith line &Column in M & ith element in B
                  p1=M.data();
                  p2=M.data();
                  for(long j=0;j<=nbVar;j++)
                  {
                     if( (j>=i) && (j<nbVar) ) B(j)=B(j+1);
                     for(long k=0;k<=nbVar;k++)
                     {
                        if((j==i) || (k==i)) p2++;
                        else *p1++ = *p2++;
                     }
                  }
               M.resizeAndPreserve(nbVar,nbVar);
               B.resizeAndPreserve(nbVar);
               */
	M.resize(nbVar,nbVar);
	B.resize(nbVar);

       // recalculate expected R-factor // Zdenek
	mRex *= sqrt((nbObs-nbVar)/(nbObs-(nbVar+1.))); // TODO:: Include constraints
      
	goto LSQNumObj_Refine_Restart;
      }
    } // Check for singular values

    // Zdenek (begin)
   // Get matrix of constraints
    CrystMatrix_REAL Con = GetConstraintsMatrix();
    long nbCon = Con.rows(); // Zdenek
   // Get Regularization Operator matrix
    const CrystMatrix_REAL HReg = GetGlobalWeightedRegularizationMatrix();
   // check if it is zero
    const bool bReg = (fabs(HReg.max())>1.e-7 || fabs(HReg.min())>1.e-7) && !(nbConvReached > 1);
    //cout << "Regularization: " << bReg << endl;
    
    // Zdenek (end)

   //Perform "Eigenvalue Filtering" on normal matrix (using newmat library)
    if (nbCon==0) { // Zdenek
      if(!silent) cout << "LSQNumObj::Refine():Eigenvalue Filtering..." <<endl; // Zdenek
      CrystMatrix_REAL V(nbVar,nbVar);
      CrystVector_REAL invW(nbVar);//diagonal matrix, in fact
      if(!silent) cout << "LSQNumObj::Refine():Eigenvalue Filtering...1" <<endl; // Zdenek
      {
	SymmetricMatrix newmatA(nbVar);
	Matrix   newmatV(nbVar,nbVar),
	         newmatN(nbVar,nbVar);
	DiagonalMatrix newmatW(nbVar);
	ColumnVector newmatB(nbVar);
       //Regularization // Zdenek
	newmatA = 0.;
	newmatB = 0.;
	if(bReg) {
	  SymmetricMatrix newmatH(nbVar);
	  ColumnVector newmatPar(nbVar);
	  for(long i=0; i<nbVar; i++) {
	    newmatPar(i+1) = mRefParList.GetParNotFixed(i).GetValue();
	    for(long j=0; j<=i; j++) newmatH(i+1,j+1)=HReg(i,j);
	  }
	  newmatA += newmatH;
	  //newmatA += newmatH.AsDiagonal();
	  for(int i=1; i<=nbVar; i++)
	    newmatA(i,i) += newmatH(i,i);
	  newmatB +=-newmatH*newmatPar;
	}
       //'Derivative scaling' matrix
	DiagonalMatrix newmatDscale(nbVar);
	for(long i=0;i<nbVar;i++)
	  newmatDscale(i+1,i+1) = 1./sqrt(newmatA(i+1,i+1)+M(i,i)); // Zdenek
	for(long i=0;i<nbVar;i++)
	{
	  newmatB(i+1)+=B(i);//NOT scaled // += Zdenek
	  for(long j=0;j<=i;j++) // for(long j=0;j<nbVar;j++)
	    newmatA(i+1,j+1) = (newmatA(i+1,j+1)+M(i,j)) * newmatDscale(i+1,i+1) * newmatDscale(j+1,j+1); // Zdenek
	}
	if(!silent) cout << "LSQNumObj::Refine():Eigenvalue Filtering...2" <<endl; // Zdenek

	//cout << newmatA.SubMatrix(1,3,1,3) <<endl;
	//if(!silent) cout << "LSQNumObj::Refine():Eigenvalue Filtering...3" <<endl;
            
       //Jacobi(newmatA,newmatW,newmatV);
	try
	{
	  EigenValues(newmatA,newmatW,newmatV);
	}
	catch(const exception & e)
	{
	  cout<<"Caught a Newmat exception :"<<e.what()<<endl;
	  cout<<"A:"<<endl<<newmatA<<endl<<"W:"<<endl<<newmatW<<endl<<"V:"<<endl<<newmatV<<endl;
	  exit(0);
	  //throw ObjCrystException("LSQNumObj::Refine():caught a newmat exception during Eigenvalues computing !");
	}
	ColumnVector newmatDelta(nbVar);
	DiagonalMatrix newmatInvW(nbVar);
        if(!silent) cout << "LSQNumObj::Refine():Eigenvalue Filtering...3" <<endl; // Zdenek
       //Avoid singular values
	{
	  REAL max=newmatW.MaximumAbsoluteValue();
	  REAL minAllowedValue=1e-5*max;// :TODO: Check if reasonable !
	  for(long i=0;i<nbVar;i++) 				
	    if(abs(newmatW(i+1,i+1)) > minAllowedValue) // Zdenek - abs 
	      newmatInvW(i+1,i+1)= 1./newmatW(i+1,i+1);
	    else 
	    {
	      if(!silent) cout << "LSQNumObj::Refine():fixing ill-cond EigenValue "<< i <<endl;
	      newmatInvW(i+1,i+1) = 0.;
	    }
	}
	newmatN=newmatV * newmatInvW * newmatV.t();
            
       //Back 'Derivative Scaling' :TODO:
	if(!silent) cout << "Back-scaling" << endl; // Zdenek
	newmatN = newmatDscale * newmatN * newmatDscale;
	newmatDelta = newmatN * newmatB;
            
	for(long i=0;i<nbVar;i++)
	{
	  invW(i)=newmatInvW(i+1,i+1);   // Zdenek :TODO: Is W-matrix scaled properly?
	  deltaVar(i)=newmatDelta(i+1);
	  for(long j=0;j<nbVar;j++) 
	  {
	    N(i,j) = newmatN(i+1,j+1);
	    V(i,j) = newmatV(i+1,j+1); // Zdenek:TODO: Is V-matrix scaled properly?
	  }
	}
      }
    }//End EigenValue filtering
    // Zdenek (begin)
    else { // We have some constraints
     // Using generalised SVD instead of EigenValue filtering
      if(!silent) cout << "LSQNumObj::Refine():Using constraints and GSVD..." <<endl;
      //CrystMatrix_REAL V(nbVar,nbVar);
      //CrystVector_REAL invW(nbVar);//diagonal matrix, in fact
      if(!silent) cout << "LSQNumObj::Refine():GSVD...1" <<endl;
      {
       // GSVD: A = U * Sig1 * [0,R] * Q'
       //       B = V * Sig2 * [0,R] * Q'
	SymmetricMatrix newmatA(nbVar);
	Matrix   newmatB(nbCon,nbVar),
	  newmatU, newmatV, newmatSig1, newmatSig2, newmatzR, newmatQ, // zR=[0,R]
	  newmatN, newmatW; // W = inv([0,R]*Q'), N = inv(A) = W*inv(Sig1)*U', 
	                    // inv(Sig1) = 1/Sig1(i,i) ... i=1...k, = 0 ... i=k+1...k+l
	Matrix newmatX, newmatP(nbVar,nbVar); // X' = [0,R]*Q', P ... solution projector matrix
	ColumnVector newmatb(nbVar);
	int k, l;
	vector<int> ind;
	newmatA = 0.;
	newmatb = 0;
       //Regularization // Zdenek
	if(bReg) {
	  SymmetricMatrix newmatH(nbVar);
	  ColumnVector newmatPar(nbVar);
	  for(long i=0; i<nbVar; i++) {
	    newmatPar(i+1) = mRefParList.GetParNotFixed(i).GetValue();
	    for(long j=0; j<=i; j++) newmatH(i+1,j+1)=HReg(i,j);
	  }
	  newmatA += newmatH;
	  for(int i=1; i<=nbVar; i++)
	    newmatA(i,i) += newmatH(i,i);
	  newmatb +=-newmatH*newmatPar;
	}
       //'Derivative scaling' matrix
	DiagonalMatrix newmatDscale(nbVar);
	for(long i=0;i<nbVar;i++)
	  newmatDscale(i+1,i+1) = 1./sqrt(newmatA(i+1,i+1)+M(i,i)); // Zdenek
	for(long i=0;i<nbVar;i++)
	{
	  newmatb(i+1)+=B(i);//NOT scaled // += zdenek
	  for(long j=0;j<=i;j++) // TODO:: j=i? // for(long j=0;j<nbVar;j++)
	    newmatA(i+1,j+1) = (M(i,j)+newmatA(i+1,j+1)) * newmatDscale(i+1,i+1) * newmatDscale(j+1,j+1);
	}
       // Constraints
	//memcpy( newmatB.Store(), Con.data(), nbCon*nbVar*sizeof(REAL) );
	for(int i=0; i<nbCon; i++)
	  for(int j=0; j<nbVar; j++) newmatB(i+1,j+1)=Con(i,j)*newmatDscale(j+1,j+1);
       // Right hand side is zero identically
	if(!silent) cout << "LSQNumObj::Refine():GSVD...2" <<endl;

	try
	{
	  GSVD(newmatA,newmatB,newmatSig1,newmatSig2,newmatU,newmatV,newmatzR,newmatQ,k,l,ind);
  
	  if( (k+l)!=newmatA.Ncols() )
	    throw std::logic_error("LSQNumObj::Refine():GSVD: Math error: rank([A; B]) neq nb.ofColms A.");

	}
	catch(const std::exception & e)
	{
	  cout<<"Caught a GSVD exception :"<<e.what()<<endl;
	  cout<<" --- GSVD --- "<<"\n";
	  cout<<" k: "<<k<<", l: "<<l<<"\n";
	  cout<<" --- A ---"<<"\n"<<newmatA;
	  cout<<" --- B ---"<<"\n"<<newmatB;
	  cout<<" --- Sig1 ---"<<"\n"<<newmatSig1;
	  cout<<" --- Sig2 ---"<<"\n"<<newmatSig2;
	  cout<<" --- U ---"<<"\n"<<newmatU;
	  cout<<" --- V ---"<<"\n"<<newmatV;
	  cout<<" --- zR ---"<<"\n"<<newmatzR;
	  cout<<" --- Q ---"<<"\n"<<newmatQ;
	  cout<<" --- ind ---"<<"\n";
	  for(int i=0; i<ind.size(); i++) cout<<" "<<ind[i];
	  cout<<"\n";
	  cout<<" --- Dscale ---"<<"\n"<<newmatDscale;
	  cout<<"This can happen when there is a strong correlation. Check following sets of parameters:\n";
	  for(int j=0; j<(newmatA.Ncols()-k-l); j++) {
	    for(int i=0; i<newmatQ.Nrows(); i++)
	      if (abs(newmatQ(i+1,j+1))>1.e-6) cout << "  " << mRefParList.GetParNotFixed(i).GetName();
	    cout <<"\n";
	  }
	  cout<<flush;
	  exit(0);
	  //throw ObjCrystException("LSQNumObj::Refine():caught a newmat exception when performing GSVD !");
	}
	ColumnVector newmatDelta(nbVar);
	DiagonalMatrix newmatInvSig1(nbVar);
        if(!silent) cout << "LSQNumObj::Refine():GSVD...3" <<endl;
       //Avoid singular values
	{
	  REAL max=newmatSig1.MaximumAbsoluteValue(); // This is probably =1
	  REAL minAllowedValue=1e-5*max;// :TODO: Check if reasonable !
	  for(long i=0;i<k;i++) 				
	    if(abs(newmatSig1(i+1,i+1)) > minAllowedValue) 
	      newmatInvSig1(i+1,i+1)= 1./newmatSig1(i+1,i+1);
	    else 
	    {
	      if(!silent) cout << "LSQNumObj::Refine():fixing ill-cond EigenValue "<< i <<endl;
	      newmatInvSig1(i+1,i+1) = 0.;
	    }
	  for(long i=k;i<nbVar;i++)
	    newmatInvSig1(i+1,i+1)= 0.;
	}
       // Calculate W = inv([0,R]*Q') = Q*inv([0,R])
       // Then N = inv(M) = W*inv(Sig1)*U'
	newmatW=newmatQ*newmatzR.i();
	newmatN=newmatW * newmatInvSig1 * newmatU.t();
       // Projector matrix: P(u,v) = sum(i=1..k) W(u,i) X(v,i), where X = Q*[0;R]
	newmatX=newmatQ*newmatzR.t();
	for(long u=1;u<=nbVar;u++)
	  for(long v=1;v<=nbVar;v++) {
	    REAL t = 0;
	    for(long i=1;i<=k;i++)
	      t+=newmatW(u,i)*newmatX(v,i);
	    newmatP(u,v) = t;
	  }
   
       //Back 'Derivative Scaling' :TODO:
	if(!silent) cout << "Back-scaling" << endl;
	newmatN = newmatDscale * newmatN * newmatDscale;
	newmatDelta = newmatN * newmatb;
	newmatB = newmatB * newmatDscale.i();
	newmatP = newmatDscale * newmatP * newmatDscale.i();

	/*{
	  ColumnVector newmatpar(nbVar);
	  for(int i=0; i<nbVar; i++) newmatpar(i+1) = mRefParList.GetParNotFixed(i).GetValue();
	  cout << " --- B --- \n";
	  cout << newmatB;
	  cout << "--- B*delta --- \n";
	  ColumnVector t = newmatB * newmatDelta;
	  cout << t.t();
	  cout << "--- B*par --- \n";
	  t = newmatB * newmatpar;
	  cout << t.t();
	  cout << "\tmaxAbs: " << t.MaximumAbsoluteValue() << "\n";
	  }*/

	for(long i=0;i<nbVar;i++)
	{
	  //invW(i)=newmatInvSig1(i+1,i+1);   // TODO:: ??? Are we really using this invW?
	  deltaVar(i)=newmatDelta(i+1);
	  for(long j=0;j<nbVar;j++) 
	  {
	    N(i,j) = newmatN(i+1,j+1);
	    P(i,j) = newmatP(i+1,j+1);
	    //V(i,j) = newmatU(i+1,j+1); // TODO:: ??? Are we really using this invW?
	  }
	}
	
      } // GSVD block
    }
    // Zdenek (end)

    if(!bMarquardtConverged) { // Zdenek
			
      /// Applying new computed values :TODO: & Check if a limit has been hit
      if(!silent) cout << "LSQNumObj::Refine():Computing new values for variables" <<endl; // Zdenek
      
      /// Without constraints the non-fixed parameters can be mutated
      if(nbCon==0) {
	for(i=0;i<nbVar;i++)
	{
	  mRefParList.GetParNotFixed(i).Mutate(deltaVar(i));
	}
      } else {
	/// With constraints possible parameters mutations are find first, projected and finally applied
	for(i=0;i<nbVar;i++) {
	  const REAL oldvalue=mRefParList.GetParNotFixed(i).GetValue();
	  const REAL expected=oldvalue+deltaVar(i);
	  mRefParList.GetParNotFixed(i).Mutate(deltaVar(i));
	  const REAL newvalue=mRefParList.GetParNotFixed(i).GetValue();
	  /*if(abs(expected-newvalue)>1.e-4) {
	    cout << "Warning: Parameter i="<<i<<" ("<<mRefParList.GetParNotFixed(i).GetName()<<") ";
	    cout << "not mutated as expected.\n";
	    }*/
	  mRefParList.GetParNotFixed(i).Mutate(-(newvalue-oldvalue));
	  deltaVar(i) = newvalue-oldvalue;
	}
	
	CrystVector_REAL deltaVarP(nbVar); // projected deltaVar
	deltaVarP = 0.;
	for(long i=0; i<nbVar; i++)
	  for(long j=0; j<nbVar; j++) 
	    deltaVarP(i) += P(i,j)*deltaVar(j);
	
	for(i=0;i<nbVar;i++)
	{
	  mRefParList.GetParNotFixed(i).Mutate(deltaVarP(i));
	}
	
	deltaVar = deltaVarP;
      } // else(nbCon==0)

      if(!silent) cout << "LSQNumObj::Refine():Computing statistics for last cycle..." <<endl; // Zdenek
     //for statistics...
      //mRefParList.Print();
      calc=this->GetLSQCalc();
     //Chi^2
      if(!silent) cout << "LSQNumObj::Refine():Computing Chi^2 for last cycle..." <<endl; // Zdenek
      // Zdenek	         
      {
	REAL oldChiSq=mChiSq;
	REAL oldChiSqReg=mChiSqReg;
	tmpV1 = mObs;
	tmpV1 -= calc;
	tmpV1 *= tmpV1;
	tmpV1 *= mWeight;
	mChiSq=tmpV1.sum();
	if(bReg)
	  mChiSqReg = this->GetRegularizationChiSq();
	mChiSq += mChiSqReg;
	if(true==useLevenbergMarquardt)
	{
	  if(mChiSq > (oldChiSq*1.0001))
	  {
	    mRefParList.RestoreParamSet(mIndexValuesSetLast);
	    marquardt *= marquardtMult;
	    if(!silent)
	    {
	      cout << "LSQNumObj::Refine(Chi^2="<<oldChiSq<<"->"<<mChiSq
		   <<")=>new Levenberg-Marquardt factor :"
		   << FormatFloat(marquardt,18,14) <<endl;
	    }
	    mChiSq=oldChiSq;
	    mChiSqReg=oldChiSqReg;
	    if(marquardt>1e8)
	    {
	      mRefParList.RestoreParamSet(mIndexValuesSetInitial);
	      if(callBeginEndOptimization) this->EndOptimization();
	      if(!silent) mRefParList.Print();
	      throw ObjCrystException("LSQNumObj::Refine():Levenberg-Marquardt diverging !");
	    }
	    goto LSQNumObj_Refine_Restart_MarquardtNotUpdate; // Zdenek
	  }
	  else { // Zdenek missing { here - bug?
	    marquardt /= marquardtMult;
	    // Zdenek
	    if(marquardt<1e-2) marquardt=1e-2;
	    if(!silent) cout << "LSQNumObj::Refine():new Levenberg-Marquardt factor :" ;
	    if(!silent) cout << FormatFloat(marquardt,18,14) <<endl;
	  }
	} // if(true==useLevenbergMarquardt) // Zdenek missing } here - bug?
	
	 // convergence condition
	  //if ( ( (oldChiSq-mChiSq)>=0 ) && ( (oldChiSq-mChiSq)<0.01 || (oldChiSq-mChiSq)<0.001*oldChiSq ) ) nbConvReached++;
	if( terminateOnDeltaChi2 && (mChiSq<=1.0001*ChisSqPreviousCycle)
	    && (minChi2var>( (ChisSqPreviousCycle-mChiSq)/abs(ChisSqPreviousCycle+1e-6) ) ) )
	  nbConvReached++; // Zdenek abs?
      } // Chi^2
	// Zdenek (begin)
				
      //} // remove
	     
     // in case of convergence of LevenbergMarquart force matrix recalcualtion with zero factor
      if( true==useLevenbergMarquardt && nbConvReached>1 ) { // Zdenek
	marquardt = 0.;
	if(!silent) cout << "LSQNumObj::Refine():Convergence reached. Force recalculation with Levenberg-Marquardt factor :" ;
	if(!silent) cout << FormatFloat(marquardt,18,14) <<endl;
	bMarquardtConverged=true;
	goto LSQNumObj_Refine_Restart_MarquardtUpdate;
      }
			
    } // if(!bMarquardtConverged)
    // Zdenek (end)

   //Sigmas
    if(nbObs==nbVar) 
      for(i=0;i<nbVar;i++) mRefParList.GetParNotFixed(i).SetSigma(0);
    else
      for(i=0;i<nbVar;i++) mRefParList.GetParNotFixed(i).SetSigma(sqrt(N(i,i)*mChiSq/(nbObs-nbVar)));
   //Correlations :TODO: re-compute M and N if using Levenberg-Marquardt // Zdenek - done at the end
    mCorrelMatrix.resize(nbVar,nbVar);
    mvVarCovar.clear();
  
    for(i=0;i<nbVar;i++)
    {
      RefinablePar *pi=&(mRefParList.GetParNotFixed(i));
      for(j=0;j<nbVar;j++)
      {
	RefinablePar *pj=&(mRefParList.GetParNotFixed(j));
	mCorrelMatrix(i,j)=sqrt(N(i,j)*N(i,j)/N(i,i)/N(j,j));
	if(nbObs!=nbVar) 
	  mvVarCovar[make_pair(pi,pj)]=N(i,j)*mChiSq/(nbObs-nbVar);
      }
    }
   //R-factor
    tmpV1 = mObs;
    tmpV1 -= calc;
    tmpV1 *= tmpV1;
    tmpV2 = mObs;
    tmpV2 *= tmpV2;
    mR=sqrt(tmpV1.sum()/tmpV2.sum());
   //Rw-factor
    tmpV1 *= mWeight;
    tmpV2 *= mWeight;
    mRw=sqrt(tmpV1.sum()/tmpV2.sum());
   //OK, finished
    if(!silent) cout << "finished cycle #"<<cycle <<"/"<<nbCycle <<". Rw="<<Rw_ini<<"->"<<mRw<<",    Chi^2="<<mChiSq<<endl;
    if (mSaveReportOnEachCycle) this->WriteReportToFile();
      
    if(!silent) this->PrintRefResults();
    //if( terminateOnDeltaChi2 && (minChi2var>( (ChisSqPreviousCycle-mChiSq)/abs(ChisSqPreviousCycle+1e-6) ) ) ) break; // Zdenek
   // Zdenek (begin)
    if ( nbConvReached > 1 ) { // Zdenek >=1
      if(!silent) cout << "LSQNumObj::Refine():Convergence reached." <<endl;
      break;
    }
		
    if (mStopAfterCycle) {
      if(!silent) cout << "LSQNumObj::Refine():Refinement stopped." <<endl;
      break;
    }
    
  } // for( ; cycle < ...)

  if(true==useLevenbergMarquardt && (!bMarquardtConverged ) && !silent)
  {
    cout << "LSQNumObj::Refine():Warning: Levenberg-Marquardt has not converged yet.\n";
    cout << "                             Check program output if parameters deviations\n";
    cout << "                             were calculated properly!" <<endl;
	// Zdenek (end)
  }
  if(callBeginEndOptimization) this->EndOptimization();
}

void LSQNumObj::StopAfterCycle () // Zdenek
{
	cout << "LSQNumObj: name='" << this->GetName() << "' requested to stop after current cycle." << endl;
	mStopAfterCycle = true;
}

CrystMatrix_REAL LSQNumObj::CorrelMatrix()const{return mCorrelMatrix;};

REAL LSQNumObj::Rfactor()const{return mR;};

REAL LSQNumObj::RwFactor()const{return mRw;};

REAL LSQNumObj::ChiSquare()const{return mChiSq;};


//void RecursiveMapFunc(RefinableObj &obj,map<RefinableObj*,unsigned int> &themap, const unsigned int value)
void RecursiveMapFunc(RefinableObj &obj,vector< pair<RefinableObj*,unsigned int> > &themap, const unsigned int value) // Zdenek
{
   //themap[&obj]=value;
   themap.push_back( make_pair(&obj,value) );
   ObjRegistry<RefinableObj> *pObjReg=&(obj.GetSubObjRegistry());
   for(int i=0;i<pObjReg->GetNb();i++)
     RecursiveMapFunc(pObjReg->GetObj(i),themap,value);
   return;
}

void LSQNumObj::SetRefinedObj(RefinableObj &obj, const unsigned int LSQFuncIndex, const bool init, const bool recursive)

{
   if(init)
   {
      mvRefinedObjMap.clear();
   }
   if(recursive) RecursiveMapFunc(obj,mvRefinedObjMap,LSQFuncIndex);
   //else mvRefinedObjMap[&obj]=LSQFuncIndex; // Zdenek
   // Zdenek (begin)
   else {
     // try to find the obj
     int ind = -1;
     for(int i=0; i<mvRefinedObjMap.size(); i++)
       if(mvRefinedObjMap[i].first==&obj) { ind = i; break; }
     // set the value
     if(ind==-1) mvRefinedObjMap.push_back( make_pair(&obj,LSQFuncIndex) );
     else mvRefinedObjMap[ind]=make_pair(&obj,LSQFuncIndex);
   }
   // Zdenek (end)
}

//ObjRegistry<RefinableObj> &LSQNumObj::GetRefinedObjList(){return mRecursiveRefinedObjList;}

//const map<RefinableObj*,unsigned int>& LSQNumObj::GetRefinedObjMap() const
const vector< pair<RefinableObj*,unsigned int> >& LSQNumObj::GetRefinedObjMap() const // Zdenek
{
   return mvRefinedObjMap;
}

//map<RefinableObj*,unsigned int>& LSQNumObj::GetRefinedObjMap()
vector< pair<RefinableObj*,unsigned int> >& LSQNumObj::GetRefinedObjMap() // Zdenek
{
   return mvRefinedObjMap;
}


RefinableObj& LSQNumObj::GetCompiledRefinedObj(){return mRefParList;}

const RefinableObj& LSQNumObj::GetCompiledRefinedObj()const{return mRefParList;}

void LSQNumObj::SetUseSaveFileOnEachCycle(bool yesOrNo)
{
   mSaveReportOnEachCycle=yesOrNo;
}

void LSQNumObj::SetSaveFile(string fileName)
{
   mSaveFileName=fileName;
}

void LSQNumObj::PrintRefResults() const
{
   //:TODO:
   //this->PrepareRefParList(); activate this when PrepareRefParList() will be more savy
   cout << "Results after last refinement :(" ;
   cout << mRefParList.GetNbParNotFixed()<< " non-fixed parameters)"<<endl;
   #ifdef __ZDENEK__
   cout << "R-factor  : " << mR<<endl;
   cout << "Rw-factor : " << mRw<<endl;
   cout << "Chi-Square: " << mChiSq-mChiSqReg<<endl;
   cout << "GoF: " << mRw/mRex<<endl; // Zdenek
   #endif
   cout << "Variable information : Initial, last cycle , current values and sigma  + dp = (Initial - current value), (dp / initial), (dp / sigma), extended parameter name"<<endl;
   for (int i=0;i<mRefParList.GetNbPar();i++)
   {
      if( (true==mRefParList.GetPar(i).IsFixed()) 
            || (false == mRefParList.GetPar(i).IsUsed()) ) continue;
      cout << FormatString(mRefParList.GetPar(i).GetName(),30) << "  " ;
      cout << FormatFloat((mRefParList.GetParamSet(mIndexValuesSetInitial))(i)*mRefParList.GetPar(i).GetHumanScale(),15,8) << "  " ;
      cout << FormatFloat((mRefParList.GetParamSet(mIndexValuesSetLast))(i)*mRefParList.GetPar(i).GetHumanScale(),15,8) << "  " ;
      cout << FormatFloat(mRefParList.GetPar(i).GetHumanValue(),15,8) << "  " ;
      cout << FormatFloat(mRefParList.GetPar(i).GetHumanSigma(),15,8) << "  " ;
      double initial = (mRefParList.GetParamSet(mIndexValuesSetInitial))(i) * mRefParList.GetPar(i).GetHumanScale();
      double fitted = mRefParList.GetPar(i).GetHumanValue();
      double sigma = mRefParList.GetPar(i).GetHumanSigma();
      double change_initial = initial - fitted;
      double change_initial_rel = fabs(change_initial) / initial * 100.;
      if (initial < 1e-30) change_initial_rel = 0.0;
      double change_sigma_rel = fabs(change_initial) / sigma * 100.;
      if (sigma < 1e-30) change_sigma_rel = 0.0;
      cout << "    " << FormatFloat(change_initial, 16, 10) << "  " ;
      cout << FormatFloat(change_initial_rel, 7, 2) << "%  " ;
      cout << FormatFloat(change_sigma_rel, 7, 2) << "%  " ;
      cout << FormatString(mRefParList.GetPar(i).GetExtName(),50);
      //
      //cout << FormatFloat(mRefParList(i).DerivStep(),16,12) << "  ";
      //cout << varNames[i] << "  " << var0(i) << "  " << varLast(i) 
      //cout<< "  " << varCurrent(i)<< "  " << sigmaValues(i)<<endl;
      cout << endl;
   }
   #ifndef __ZDENEK__
   cout << "R-factor  : " << mR<<endl;
   cout << "Rw-factor : " << mRw<<endl;
   cout << "Chi-Square: " << mChiSq<<endl;
   //cout << "GoF: " << mChiSq/this->GetLSQWeight().numElements()<<endl; // Zdenek strange ? not missing sqrt ?
   cout << "GoF: " << mRw/mRex<<endl; // Zdenek
   #endif
   cout <<endl;
}

void LSQNumObj::SetDampingFactor(const REAL newDampFact)
{
   mDampingFactor=newDampFact;
}

void LSQNumObj::PurgeSaveFile()
{
   //:TODO:
}

void LSQNumObj::WriteReportToFile() const
{
   //:TODO:
}

void LSQNumObj::OptimizeDerivativeSteps()
{
   //:TODO:
}

const std::map<pair<const RefinablePar*,const RefinablePar*>,REAL > & LSQNumObj::GetVarianceCovarianceMap()const
{ return mvVarCovar;}

void LSQNumObj::PrepareRefParList(const bool copy_param)
{
  /*cout << "map<RefinableObj*,unsigned int>"<<endl;
  int i=0;
  for(vector< pair<RefinableObj*,unsigned int> >::iterator pos=mvRefinedObjMap.begin();pos!=mvRefinedObjMap.end();++pos) 
    {
     cout<<setw(4)<<i; i++;
     cout<<setw(45)<<(*pos).first->GetClassName();
     cout<<setw(40)<<(*pos).first->GetName();
     cout<<" with params.: "<<setw(4)<<(*pos).first->GetNbPar();
     cout<<endl;
   }
   cout << "-------------------------------"<<endl;*/
   mRefParList.ResetParList();
   //for(map<RefinableObj*,unsigned int>::iterator pos=mvRefinedObjMap.begin();pos!=mvRefinedObjMap.end();++pos)
   for(vector< pair<RefinableObj*,unsigned int> >::iterator pos=mvRefinedObjMap.begin();pos!=mvRefinedObjMap.end();++pos) // Zdenek
   {
      VFN_DEBUG_MESSAGE("LSQNumObj::PrepareRefParList():"<<pos->first->GetName(),4);
      //mRecursiveRefinedObjList.GetObj(i).Print();
      //mRefParList.AddPar(*(pos->first),copy_param);
      mRefParList.AddPar(*((*pos).first),copy_param); // Zdenek
   }
   //mRefParList.Print();
   if(copy_param) mRefParList.SetDeleteRefParInDestructor(true);
   else mRefParList.SetDeleteRefParInDestructor(false);

   // prepare list of refined sub-objects with constraints
   // prepare list of refined sub-objects with regularization operators
   mConstrainedRefinedObjList.DeRegisterAll(); // Zdenek
   mRegularizedRefinedObjList.DeRegisterAll(); // Zdenek
   /*for(int i=0;i<mRecursiveRefinedObjList.GetNb();i++)
   {
   		if(mRecursiveRefinedObjList.GetObj(i).GetNbLSQConstraints()>0)
   			mConstrainedRefinedObjList.Register(mRecursiveRefinedObjList.GetObj(i));
			}*/
   for(vector< pair<RefinableObj*,unsigned int> >::iterator pos=mvRefinedObjMap.begin();
       pos!=mvRefinedObjMap.end(); ++pos) // Zdenek
   {
      VFN_DEBUG_MESSAGE("LSQNumObj::PrepareRefParList():"<<pos->first->GetName(),4);
      if( pos->first->GetNbLSQConstraints()>0 )
	mConstrainedRefinedObjList.Register( *(pos->first) ); // Zdenek
      if( pos->first->GetNbLSQRegularizationOperator(pos->second)>0 )
	mRegularizedRefinedObjList.Register( *(pos->first) ); // Zdenek
   }
   if( mRefParList.GetNbLSQConstraints()>0 ) { // Zdenek
     mConstrainedRefinedObjList.Register( mRefParList );
   }
}

const CrystVector_REAL& LSQNumObj::GetLSQCalc() const
{
   const CrystVector_REAL *pV;
   unsigned long nb=0;
   //for(map<RefinableObj*,unsigned int>::const_iterator pos=mvRefinedObjMap.begin();pos!=mvRefinedObjMap.end();++pos)
   for(vector< pair<RefinableObj*,unsigned int> >::const_iterator pos=mvRefinedObjMap.begin();pos!=mvRefinedObjMap.end();++pos) // Zdenek
   {
     //if(pos->first->GetNbLSQFunction()==0) continue;
     //pV=&(pos->first->GetLSQCalc(pos->second));
      if((*pos).first->GetNbLSQFunction()==0) continue;
      pV=&((*pos).first->GetLSQCalc(pos->second));
      const unsigned long n2=pV->numElements();
      if((nb+n2)>mLSQCalc.numElements()) mLSQCalc.resizeAndPreserve(nb+pV->numElements());
      const REAL *p1=pV->data();
      REAL *p2=mLSQCalc.data()+nb;
      for(unsigned long j=0;j<n2;++j) *p2++ = *p1++;
      nb+=n2;
   }
   if(mLSQCalc.numElements()>nb) mLSQCalc.resizeAndPreserve(nb);
   return mLSQCalc;
}

const CrystVector_REAL& LSQNumObj::GetLSQObs() const
{
   const CrystVector_REAL *pV;
   unsigned long nb=0;
   //for(map<RefinableObj*,unsigned int>::const_iterator pos=mvRefinedObjMap.begin();pos!=mvRefinedObjMap.end();++pos)
   for(vector< pair<RefinableObj*,unsigned int> >::const_iterator pos=mvRefinedObjMap.begin();pos!=mvRefinedObjMap.end();++pos) // Zdenek
   {
     //if(pos->first->GetNbLSQFunction()==0) continue;
     //pV=&(pos->first->GetLSQObs(pos->second));
      if((*pos).first->GetNbLSQFunction()==0) continue;
      pV=&((*pos).first->GetLSQObs(pos->second));
      const unsigned long n2=pV->numElements();
      if((nb+n2)>mLSQObs.numElements()) mLSQObs.resizeAndPreserve(nb+pV->numElements());
      const REAL *p1=pV->data();
      REAL *p2=mLSQObs.data()+nb;
      for(unsigned long j=0;j<n2;++j) *p2++ = *p1++;
      nb+=n2;
   }
   if(mLSQObs.numElements()>nb) mLSQObs.resizeAndPreserve(nb);
   return mLSQObs;
}

const CrystVector_REAL& LSQNumObj::GetLSQWeight() const
{
   const CrystVector_REAL *pV;
   unsigned long nb=0;
   //for(map<RefinableObj*,unsigned int>::const_iterator pos=mvRefinedObjMap.begin();pos!=mvRefinedObjMap.end();++pos)
   for(vector< pair<RefinableObj*,unsigned int> >::const_iterator pos=mvRefinedObjMap.begin();pos!=mvRefinedObjMap.end();++pos) // Zdenek
   {
     //if(pos->first->GetNbLSQFunction()==0) continue;
     //pV=&(pos->first->GetLSQWeight(pos->second));
      if((*pos).first->GetNbLSQFunction()==0) continue;
      pV=&((*pos).first->GetLSQWeight(pos->second));
      const unsigned long n2=pV->numElements();
      if((nb+n2)>mLSQWeight.numElements()) mLSQWeight.resizeAndPreserve(nb+pV->numElements());
      const REAL *p1=pV->data();
      REAL *p2=mLSQWeight.data()+nb;
      for(unsigned long j=0;j<n2;++j) *p2++ = *p1++;
      nb+=n2;
   }
   if(mLSQWeight.numElements()>nb) mLSQWeight.resizeAndPreserve(nb);
   return mLSQWeight;
}

const CrystVector_REAL& LSQNumObj::GetLSQDeriv(RefinablePar&par)
{
   const CrystVector_REAL *pV;
   unsigned long nb=0;
   //for(map<RefinableObj*,unsigned int>::iterator pos=mvRefinedObjMap.begin();pos!=mvRefinedObjMap.end();++pos)
   for(vector< pair<RefinableObj*,unsigned int> >::const_iterator pos=mvRefinedObjMap.begin();pos!=mvRefinedObjMap.end();++pos) // Zdenek
   {
     //if(pos->first->GetNbLSQFunction()==0) continue;
     //pV=&(pos->first->GetLSQDeriv(pos->second,par));
      if((*pos).first->GetNbLSQFunction()==0) continue;
      pV=&((*pos).first->GetLSQDeriv(pos->second,par));
      const unsigned long n2=pV->numElements();
      if((nb+n2)>mLSQDeriv.numElements()) mLSQDeriv.resizeAndPreserve(nb+pV->numElements());
      const REAL *p1=pV->data();
      REAL *p2=mLSQDeriv.data()+nb;
      for(unsigned long j=0;j<n2;++j) *p2++ = *p1++;
      nb+=n2;
   }
   if(mLSQDeriv.numElements()>nb) mLSQDeriv.resizeAndPreserve(nb);
   return mLSQDeriv;
}

void LSQNumObj::BeginOptimization(const bool allowApproximations, const bool enableRestraints)
{
  //for(map<RefinableObj*,unsigned int>::iterator pos=mvRefinedObjMap.begin();pos!=mvRefinedObjMap.end();++pos)
    //pos->first->BeginOptimization(allowApproximations, enableRestraints);
  for(vector< pair<RefinableObj*,unsigned int> >::iterator pos=mvRefinedObjMap.begin();pos!=mvRefinedObjMap.end();++pos)
    (*pos).first->BeginOptimization(allowApproximations, enableRestraints);
}

void LSQNumObj::EndOptimization()
{
  //for(map<RefinableObj*,unsigned int>::iterator pos=mvRefinedObjMap.begin();pos!=mvRefinedObjMap.end();++pos)
    //pos->first->EndOptimization();
  for(vector< pair<RefinableObj*,unsigned int> >::iterator pos=mvRefinedObjMap.begin();pos!=mvRefinedObjMap.end();++pos)
    (*pos).first->EndOptimization();
}

CrystMatrix_REAL LSQNumObj::GetConstraintsMatrix() const // Zdenek
{
  unsigned int nc = 0;
  for(int i=0;i<mConstrainedRefinedObjList.GetNb();i++)
    nc += mConstrainedRefinedObjList.GetObj(i).GetNbLSQConstraints();

  unsigned int np = mRefParList.GetNbParNotFixed();
	
  CrystMatrix_REAL C(nc,np);
  C = 0.;
	
  unsigned int n = 0;
  // for all refined sub-objects with some constraints
  for(unsigned int i=0;i<mConstrainedRefinedObjList.GetNb();i++)
    // for all constraints of such an object
    for(unsigned int j=0;j<mConstrainedRefinedObjList.GetObj(i).GetNbLSQConstraints();j++)
    {
      // get constraint
      std::vector< const RefinablePar* > parList;
      CrystVector_REAL coef;
			
      mConstrainedRefinedObjList.GetObj(i).GetLSQConstraint(j,parList,coef);
			
      // build the matrix of constraints
      for(unsigned int k=0;k<np;k++) {
	// identify parameters from their adresses
	const RefinablePar* ppar = &(mRefParList.GetParNotFixed(k));
	//cout<< "ParNotFixed Name:" << ppar->GetName() << ", address:"<<ppar<<"\n";
	for(unsigned int m=0;m<parList.size();m++) {
	  //cout<< "par from parList Name:" << parList[m]->GetName() << ", address:"<<parList[m]<<"\n";
	  if(parList[m]==ppar)
	  {
	    // add constraint coefficient
	    C(n,k) = coef(m);
	    //cout << "IDENTIFIED - coefficient: " << coef(m) << "\n";
	  }
	}
      } // k - 0..np
      n++;
    } // j -  0..mConstrainedRefinedObjList.GetObj(i).GetNbLSQConstraints()

  // reduce empty equations
  CrystVector_int bb(C.rows());
  for(int i=0; i<C.rows(); i++) {
    REAL t = 0.;
    for(int j=0; j<C.cols(); j++) t += pow(C(i,j),2);
    bb(i) = (t>1e-4) ? 1 : 0;
  }
  
  CrystMatrix_REAL CC(bb.sum(),C.cols());
  n = 0;
  for(int i=0; i<C.rows(); i++)
    if(bb(i)>0) {
      for(int j=0; j<C.cols(); j++) CC(n,j) = C(i,j);
      n++;
    }
  
  return CC;
}

CrystMatrix_REAL LSQNumObj::GetGlobalWeightedRegularizationMatrix(const int LSQfunc) const // Zdenek
{
  
  unsigned int nbPar = mRefParList.GetNbParNotFixed();
	
  CrystMatrix_REAL H(nbPar,nbPar);
  H = 0.;
	
  // for all refined sub-objects with regularization operators
  for(unsigned int iObj=0; iObj<mRegularizedRefinedObjList.GetNb(); iObj++)
    // for all operators of such an object
    for(unsigned int iOp=0;
	iOp<mRegularizedRefinedObjList.GetObj(iObj).GetNbLSQRegularizationOperator(LSQfunc);iOp++)
    {
      // get operator
      const LSQRegularizationOperator & Op = mRegularizedRefinedObjList.GetObj(iObj).
	                                         GetLSQRegularizationOperator(iOp,LSQfunc);
      
      // list of parameters of the given operator
      const std::vector< const RefinablePar* > & OpParList = Op.GetParamList();

      // create a map linking indexes of the local operator matrix with the global matrix
      std::map< int, int > mapLoc2Glob;
      
      for(unsigned int iPar=0; iPar<nbPar; iPar++) {

	// identify parameters from their data adresses
	const RefinablePar & par = mRefParList.GetParNotFixed(iPar);
	//cout<< "ParNotFixed Name:" << par.GetName() << ", data adress:"<<par.GetPointer()<<"\n";

	for(unsigned int iOpPar=0; iOpPar<OpParList.size(); iOpPar++)
	  if( OpParList[iOpPar]->GetPointer() == par.GetPointer() ) {
	    //cout << "IDENTIFIED: " << iOpPar << " -> " << iPar << "\n";
	    mapLoc2Glob[iOpPar] = iPar;
	    break;
	  } // iOpPar

      } // iPar
	
      // build the regularization matrix

      const CrystMatrix_REAL & H0 = Op.GetRegularizationOperatorMatrix();
      const REAL Lambda = Op.GetRegularizationOperatorWeight();
      
      for(int i=0; i<OpParList.size(); i++) {
	const int I = mapLoc2Glob[i];
	for(int j=0; j<=i; j++) {
	  const int J = mapLoc2Glob[j];
	  H(I,J) += Lambda * H0(i,j);
	  H(J,I) += Lambda * H0(i,j);
	} // j
      } // i

    } // iOp
  
  return H;
}

REAL LSQNumObj::GetRegularizationChiSq(const int LSQfunc)
{
  REAL ChiSq = 0.;

  for(int iObj=0; iObj<mRegularizedRefinedObjList.GetNb(); iObj++) {
    const RefinableObj & Obj = mRegularizedRefinedObjList.GetObj(iObj);
    const int nbOp = Obj.GetNbLSQRegularizationOperator(LSQfunc);
    for(int iOp=0; iOp<nbOp; iOp++) {
      const LSQRegularizationOperator & Op = Obj.GetLSQRegularizationOperator(iOp,LSQfunc);
      ChiSq += Op.GetRegularizationOperatorWeight()*Op.GetValue();
    } // iOp
  } // iObj

  return ChiSq;
}

#ifdef __WX__CRYST__
WXCrystObjBasic* LSQNumObj::WXCreate(wxWindow* parent)
{
   this->WXDelete();
   mpWXCrystObj=new WXLSQ(parent,this);
   return mpWXCrystObj;
}
WXCrystObjBasic* LSQNumObj::WXGet()
{
   return mpWXCrystObj;
}
void LSQNumObj::WXDelete()
{
   if(0!=mpWXCrystObj)
   {
      VFN_DEBUG_MESSAGE("LSQNumObj::WXDelete()",5)
      delete mpWXCrystObj;
   }
   mpWXCrystObj=0;
}
void LSQNumObj::WXNotifyDelete()
{
   VFN_DEBUG_MESSAGE("LSQNumObj::WXNotifyDelete():"<<mName,5)
   mpWXCrystObj=0;
}
#endif

// Zdenek (begin)

//######################################################################
//
//      LSQRegularizationOperator
//
//######################################################################

ObjRegistry<LSQRegularizationOperator> gLSQRegularizationOperatorRegistry("Global LSQRegularizationOperator registry");

LSQRegularizationOperator::LSQRegularizationOperator (string objName)
  :mLambda(0.)
{
  mName = objName;
  mParamList.clear();
  mRegOpMatrix.resize(0,0);
  gLSQRegularizationOperatorRegistry.Register(*this);
}

LSQRegularizationOperator::LSQRegularizationOperator (const LSQRegularizationOperator & old)
  :mLambda(old.mLambda), mName(old.mName), mParamList(old.mParamList), mRegOpMatrix(old.mRegOpMatrix)
{
  gLSQRegularizationOperatorRegistry.Register(*this);
}

LSQRegularizationOperator::~LSQRegularizationOperator ()
{
  mParamList.clear();
  mRegOpMatrix.resize(0,0);
  gLSQRegularizationOperatorRegistry.DeRegister(*this);
}

const string & LSQRegularizationOperator::GetClassName () const
{
  const static string className="LSQRegularizationOperator";
  return className;
}

const string & LSQRegularizationOperator::GetName () const
{
  return mName;
}

const CrystMatrix_REAL & LSQRegularizationOperator::GetRegularizationOperatorMatrix () const
{
  if (mParamList.size()!=mRegOpMatrix.rows())
    throw ObjCrystException("LSQRegularizationOperator (name="+this->GetName()+"): Number of parametrs "
			    +"not same as the size of the regularization matrix!");
  
  return mRegOpMatrix;
}
  
const std::vector< const RefinablePar * > & LSQRegularizationOperator::GetParamList () const
{
   if (mParamList.size()!=mRegOpMatrix.rows())
     throw ObjCrystException("LSQRegularizationOperator (name="+this->GetName()+"): Regularization matrix size "
			     +"not same as the number of parametrs!");
   
   return mParamList;
}

REAL LSQRegularizationOperator::GetRegularizationOperatorWeight () const
{
  return mLambda;
}

void  LSQRegularizationOperator::SetRegularizationOperatorMatrix (const CrystMatrix_REAL & matrix)
{
  if (matrix.rows()!=matrix.cols())
    throw ObjCrystException("LSQRegularizationOperator (name="+this->GetName()+"): Trying to set a non-rectangular "
			    +"matrix. Regularization operator matrix has to be symmetric!");

  mRegOpMatrix = matrix;
}

void  LSQRegularizationOperator::SetParamList (const std::vector< const RefinablePar * > & paramList)
{
  mParamList.clear();
  mParamList = paramList;
}

void LSQRegularizationOperator::SetRegularizationOperatorWeight (const REAL Lambda)
{
  mLambda = Lambda;
}

REAL LSQRegularizationOperator::GetValue () const
{
  // TODO:: band matrix, save result, speed?

  REAL ChiSq = 0.;

  CrystVector_REAL u(mParamList.size());
  for(int i=0; i<mParamList.size(); i++)
    u(i) = mParamList[i]->GetValue();

  // diagonal terms
  for(int i=0; i<mParamList.size(); i++)
    ChiSq += mRegOpMatrix(i,i)*u(i)*u(i);

  // non-diagonal terms
  for(int i=0; i<mParamList.size(); i++)
    for(int j=0; j<i; j++)
       ChiSq += 2*mRegOpMatrix(i,j)*u(i)*u(j);

  return ChiSq;
}

// Definition of an Empty LSQ Regularization Operator Object
LSQRegularizationOperator EmptyLSQRegularizationOperatorObj("An Empty Global LSQ Regularization Operator");
const LSQRegularizationOperator EmptyLSQRegularizationOperatorObj_const("An Empty Global LSQ Regularization Operator");

// Zdenek (end)

}//namespace

// Zdenek (begin

//***********EXPLICIT INSTANTIATION*******************//
#include <algorithm> // for "find" algorithm in ObjRegistry
#include "ObjCryst/RefinableObj/ObjRegistry_T.cpp"
template class ObjCryst::ObjRegistry<ObjCryst::LSQRegularizationOperator>;

// Zdenek (end)
