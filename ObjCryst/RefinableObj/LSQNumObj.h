/*  ObjCryst++ Object-Oriented Crystallographic Library
    (c) 2000-2002 Vincent Favre-Nicolin vincefn@users.sourceforge.net
        2000-2001 University of Geneva (Switzerland)

    This program is free software; you can redistribute it and/or modify
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
/*   LSQNumObj.h
*  header file for Least-Squares refinement witn Numerical derivatives Objects
*
*/
#ifndef _LSQOBJNUM_H
#define _LSQOBJNUM_H

#include "ObjCryst/CrystVector/CrystVector.h"
#include "ObjCryst/RefinableObj/RefinableObj.h"
#include <string>
#include <map>

namespace ObjCryst
{

class LSQRegularizationOperator; // Zdenek

/** \brief (Quick & dirty) Least-Squares Refinement Object with Numerical derivatives
*
* This is still highly experimental !
*/
class LSQNumObj
{
   public:
      LSQNumObj(std::string objName="Unnamed LSQ object");
      ~LSQNumObj();
      /// Name for this class ("RefinableObj", "Crystal",...)
      const string & GetClassName () const; // Zdenek
      /// Name of the object.
      const string & GetName () const; // Zdenek
      /// Fix one parameter.
      ///
      /// LSQNumObj::PrepareRefParList() must be called first!
      void SetParIsFixed(const std::string& parName,const bool fix);
      /// Fix one family of parameters
      ///
      /// LSQNumObj::PrepareRefParList() must be called first!
      void SetParIsFixed(const RefParType *type,const bool fix);
      /** Fix one parameter
      *
      * Note that this will fix the copied parameter, not the one
      * in the original object. The supplied RefinablePar
      * may be either the copied one or the original.
      *
      * LSQNumObj::PrepareRefParList() must be called first!
      **/ 
      void SetParIsFixed(RefinablePar &par,const bool fix);
      /** Fix all parameters within an object
      *
      * Note that this will fix the copied parameters, not the one
      * in the original objects.
      *
      * LSQNumObj::PrepareRefParList() must be called first!
      **/ 
      void SetParIsFixed(RefinableObj &obj,const bool fix);
      /// UnFix All parameters
      ///
      /// LSQNumObj::PrepareRefParList() must be called first!
      void UnFixAllPar();
      /// Set a parameter to be used
      ///
      /// LSQNumObj::PrepareRefParList() must be called first!
      void SetParIsUsed(const std::string& parName,const bool use);
      /// Set a family of parameters to be used
      ///
      /// LSQNumObj::PrepareRefParList() must be called first!
      void SetParIsUsed(const RefParType *type,const bool use);
      /** Do the refinement
      *
      * \param nbCycle: number of LSQ cycles - if negative, the algorithm will continue
      * until it reaches (-nbcycle) or until the relative variation in Chi2 is less
      * than minChi2var
      * \param useLevenbergMarquardt: enable Levenberg-Marquardt algorithm
      * to ensure that a decrease of Chi^2 will be obtained (actually a 1%
      * increase is allowed)
      * \param callBeginEndOptimization: if true, will call RefinableObj::BeginOptimization(true,...)
      * and RefinableObj::EndOptimization(). You may not want this if the LSQ is done during another
      * (e.g. monte-carlo) optimization - but then the calling function \b must ensure
      * that approximations are disabled (using RefinableObj::SetApproximationFlag)
      * for objects where that would render derivative calculations imprecise.
      * \param minChi2var: used for termination of the refinement if the relative variation
      * of Chi2 between two successive cyles is less than minChi2var
      */
      void Refine (int nbCycle=1,bool useLevenbergMarquardt=false,
                   const bool silent=false, const bool callBeginEndOptimization=true,
                   const float minChi2var=0.01);
      /// Stop after the current cycle. Used for refinement interruption. 
      void StopAfterCycle(); // Zdenek
      CrystVector_REAL Sigma()const;
      CrystMatrix_REAL CorrelMatrix()const;
      REAL Rfactor()const;
      REAL RwFactor()const;
      REAL ChiSquare()const;   //uses the weight if specified
      /** Choose the object to refine. The minimization will be done
      * against its LSQ function and its parameters, as well as the LSQ functions
      * and parameters of its sub-objects (if recursive==true)
      *
      * \param LSQFuncIndex: one object can have a choice of several LSQ
      * functions to minimize- this allows to choose which one to minimize.
      *
      * \param init: if true, the list of refined objects is first cleared. otherwise
      * the new object (and its sub-objects) is just added to the list.
      *
      * \param recursive: if false, only the supplied object is added, and not its sub-objects
      */
      void SetRefinedObj(RefinableObj &obj, const unsigned int LSQFuncIndex=0, const bool init=true, const bool recursive=false);
      // Access to the full list of refined objects. The list is initially built
      // recursively from one object. This function allows to modify the list
      // of sub-objects before refinement (such as for removing certain types
      // of objects).
      //ObjRegistry<RefinableObj> &GetRefinedObjList();
      /** Get the map of refined objects - this is a recursive list of all the objects
      * that are taken into account for the refinement. The key is a pointer to the
      * object and the value is the LSQ function index for that object.
      */
      //const std::map<RefinableObj*,unsigned int>& GetRefinedObjMap() const;
      const std::vector< std::pair<RefinableObj*,unsigned int> >& GetRefinedObjMap() const; // Zdenek
      /** Get the map of refined objects - this is a recursive list of all the objects
      * that are taken into account for the refinement. The key is a pointer to the
      * object and the value is the LSQ function index for that object.
      */
      //std::map<RefinableObj*,unsigned int>& GetRefinedObjMap();
      std::vector< std::pair<RefinableObj*,unsigned int> >& GetRefinedObjMap(); // Zdenek
      /** Access to the RefinableObj which is the compilation of all parameters
      * from the object supplied for optimization and its sub-objects.
      *
      * Since this compilation is only updated from the suplied refinableobj and
      * its sub-objects when SetRefinedObj() and PrepareRefParList() are called,
      * it is possible to alter the fixed/limited status of parameters
      * here without affecting the parameters in the refined objects.
      */
      RefinableObj& GetCompiledRefinedObj();
      /** Access to the RefinableObj which is the compilation of all parameters
      * from the object supplied for optimization and its sub-objects.
      *
      * Since this compilation is only updated from the suplied refinableobj and
      * its sub-objects when SetRefinedObj() and PrepareRefParList() are called,
      * it is possible to alter the fixed/limited status of parameters
      * here without affecting the parameters in the refined objects.
      */
      const RefinableObj& GetCompiledRefinedObj()const;
      void SetUseSaveFileOnEachCycle(bool yesOrNo=true);
      void SetSaveFile(std::string fileName="refine.save");
      void PrintRefResults()const;
      void SetDampingFactor(const REAL newDampFact);
      void PurgeSaveFile();
      void WriteReportToFile()const;
      
      void OptimizeDerivativeSteps();
      const std::map<pair<const RefinablePar*,const RefinablePar*>,REAL > &GetVarianceCovarianceMap()const;
      /** Prepare the full parameter list for the refinement
      * \param copy_param: if false (the default), then the lsq algorithm will work directly
      * on the parameters of the refined object and sub-object. So that any modification
      * to the fixed/used/limited status applies permanently to the parameters.
      * if true, then the parameters are copied and therefore only the value of the
      * parameter is changed (and the clocks are ticked).
      *
      * \note: if copy_param==true, then any modification to the parameters (fixed, limited, used
      * status) only affects the copy and not the original. Also, calling again PrepareRefParList
      * cancels any such modification.
      *
      * \note This will be called automatically before starting the refinement only if
      * the parameter list is empty. Otherwise it should be called before refinement.
      */
      void PrepareRefParList(const bool copy_param=false);
      
      /// Get the LSQ calc vector (using either only the top or the hierarchy of object)
      const CrystVector_REAL& GetLSQCalc() const;
      /// Get the LSQ obs vector (using either only the top or the hierarchy of object)
      const CrystVector_REAL& GetLSQObs() const;
      /// Get the LSQ weight vector (using either only the top or the hierarchy of object)
      const CrystVector_REAL& GetLSQWeight() const;
      /// Get the LSQ deriv vector (using either only the top or the hierarchy of object)
      const CrystVector_REAL& GetLSQDeriv(RefinablePar&par);
      /** Tell all refined object that the refinement is beginning
      */
      void BeginOptimization(const bool allowApproximations=false, const bool enableRestraints=false);
      /** Tell all refined object that the refinement is finished
      */
      void EndOptimization();
   protected:
   private:
#ifdef __ZDENEK__
      /// Get the matrix of constraints for refined parameters
      CrystMatrix_REAL GetConstraintsMatrix() const; // Zdenek
      /** Get the sum of all weighted regularization operators matrices.
       *
       * In addition to the ChiSq value for the given LSQ function \param LSQfunc
       * a regularization term is minimised during refinemet. Minimised quantity is:
       *
       *          ChiSq(a0+da) + Lambda * (a0+da)*H*(a0+da) ,
       *
       * where a0 ... states for a vector of model parameters, da ... is the mutation
       * of model parameters refined by LSQs, Lambda ... is the weighting factor and
       * H ... is the symmetric regularization operator. Multiple such regularisation
       * operators can be possibly applied in general. The matrix returned is the sum
       * of all such operators required from sub-refinable objects. The weighting
       * factors (Lambdas) are already included in the final matrix.
      */
      CrystMatrix_REAL GetGlobalWeightedRegularizationMatrix(const int LSQfunc=0) const; // Zdenek
#endif // __ZDENEK__
      // Refined object
         /// The recursive list of all refined sub-objects
         ObjRegistry<RefinableObj> mRecursiveRefinedObjList;
#ifdef __ZDENEK__
      /// The list of all refined sub-objects with LSQ Constraints
      ObjRegistry<RefinableObj> mConstrainedRefinedObjList;
      /// The list of all refined sub-objects with LSQ Regularization Operators
      ObjRegistry<RefinableObj> mRegularizedRefinedObjList;
      REAL mChiSqReg;
      /// Get the ChiSq value associated with the regularization part
      REAL GetRegularizationChiSq(const int LSQfunc=0);
#endif // __ZDENEK__
      /** The refinable par list used during refinement. It is only a compilation
      * of the parameters in RefinableObj and its sub-objects
      *
      * This list is only updated from the suplied refinableobj and
      * its sub-objects when SetRefinedObj() and PrepareRefParList() are called,
      * so it is possible to alter the fixed/limited status of parameters
      * here without affecting the parameters in the refined objects.
      */
      mutable RefinableObj mRefParList;
      /// Damping factor for the refinement (unused yet...)
      REAL mDampingFactor;
      ///Save result to file after each cycle ?
      bool mSaveReportOnEachCycle;   
      /// Name of the refined object
      std::string mName;
      /// File name where refinement info is saved
      std::string mSaveFileName;
      REAL mR,mRw,mChiSq,mRex;
      /// Correlation matrix between all refined parameters.
      CrystMatrix_REAL mCorrelMatrix;
      ///Variance-Covariance matrix, as a std::map
      std::map<pair<const RefinablePar*,const RefinablePar*>,REAL > mvVarCovar;
      /// Observed values.
      CrystVector_REAL mObs;
      /// Weight corresponding to all observed values.
      CrystVector_REAL mWeight;
      /// Index of the set of saved values for all refinable parameters, before refinement
      /// and before the last cycle.
      int mIndexValuesSetInitial, mIndexValuesSetLast;
      /// If true, then stop at the end of the cycle. Used in multi-threading environment
      bool mStopAfterCycle;
      // The optimized object
      //RefinableObj *mpRefinedObj;
      // The index of the LSQ function in the refined object (if there are several...)
      //unsigned int mLSQFuncIndex;
      
      /** Map of the recursive list of the objects to be refined. The key is the pointer
      * to the object and the value the LSQ function index
      *
      * Individual LSQ functions can be changed using GetRefinedObjMap().
      */
      //std::map<RefinableObj*,unsigned int> mvRefinedObjMap;
      std::vector< std::pair<RefinableObj*,unsigned int> > mvRefinedObjMap; // Zdenek
      /// If true, then parameters to be refined will be copied instead of referenced.
      /// Therefore only their values and the parameter's clocks are affected when
      /// working on the copy.
      bool mCopyRefPar;
      /// Temporary arrays for LSQ functions evaluation - used when
      /// using recursive LSQ function
      mutable CrystVector_REAL mLSQObs,mLSQCalc,mLSQWeight,mLSQDeriv;
#ifdef __WX__CRYST__
   public:
      virtual WXCrystObjBasic* WXCreate(wxWindow* parent);
      WXCrystObjBasic* WXGet();
      void WXDelete();
      void WXNotifyDelete();
   protected:
      WXCrystObjBasic *mpWXCrystObj;
#endif
};

// Zdenekk(begin)

/** \brief Simple Regularization Operator for the Least-Squares Refinement.
*
* During the common LSQ refinement a ChiSq value representing the difference of
* observed and simulated data is minimized. The matter of this object is to
* provide an operator that is used to regularize refinement solution in some
* way. For example it is assumed that the solution should be smooth. This means
* a that each model parameter from a given group of parametres, representing e.g.
* a continous function, should not differ much from neighbouring parameters
* in the group. Such a condition can be expressed as an additional operator 
* that should be optimised. Hence insteed of minimizing only the ChiSq value
* the quantity:
*                      ChiSq(a0+da) + Lambda * (a0+da)*H*(a0+da) 
*
* is minimised now. Here a0 ... states for a vector of model parameters,
* da ... is the mutation of model parameters refined by LSQs,
* Lambda ... is the weighting factor and H ... is the symmetric regularization
* operator.
* 
*/
class LSQRegularizationOperator
{
 public:
  /// Constructor
  LSQRegularizationOperator (std::string objName="Unnamed LSQ Regularization Obj");
  /// Copy constructor
  LSQRegularizationOperator (const LSQRegularizationOperator & );
  /// Descructor
  ~LSQRegularizationOperator();
  /// Name for this class ("LSQRegularizationOperator",...)
  const string & GetClassName () const;
  /// Name of the object.
  const string & GetName () const;

  /** Get the matrix reprezenting the regularization operator.
   *
   * The elements of the matrix descibe relations between parameters returned
   * by GetParamsList(). Hence this list defines the order of matrix rows and columns.
   * 
   * The matrix is symmetric. The lower triangular part of the matrix should be used
   * preferably (becouse of compatibility with Newmat).
   */
  const CrystMatrix_REAL & GetRegularizationOperatorMatrix () const;
  /// Get the list of parameters. This list links the matrix elements with the appropriate parameters.
  const std::vector< const RefinablePar * > & GetParamList () const;
  /// Get the weighting factor of the regularization operator in the minimized ChiSq quantity.
  REAL GetRegularizationOperatorWeight () const;
  
  /** Set the matrix reprezenting the regularization operator.
   *
   * The elements of the matrix descibe relations between parameters returned
   * by GetParamsList(). Hence this list defines the order of matrix rows and columns.
   * 
   * The matrix must be symmetric. The lower triangular part of the matrix should be used
   * preferably (becouse of compatibility with Newmat).
   */
  void SetRegularizationOperatorMatrix (const CrystMatrix_REAL & matrix);
  /// Set the list of parameters on which the operator is working. This list links the matrix elements with the appropriate parameters.
  void SetParamList (const std::vector< const RefinablePar * > & paramList);
  /// Set the weighting factor of the regularization operator in the minimized ChiSq quantity.
  void SetRegularizationOperatorWeight (const REAL Lambda);

  /// Get (non-weighted) value of Regularization Operator applied to the current parameters set.
  REAL GetValue () const;

 private:
  /// Name of the object
  std::string mName;

  /// List of parameters on which the regularization operator is working.
  std::vector< const RefinablePar * > mParamList;
  /// Regularization operator matrix.
  CrystMatrix_REAL mRegOpMatrix;
  /// Weight of the regularization operator.
  REAL mLambda;
};

/// An empty LSQ Regularization Operator Object
extern LSQRegularizationOperator EmptyLSQRegularizationOperatorObj;
/// Constant empty LSQ Regularization Operator Object
extern const LSQRegularizationOperator EmptyLSQRegularizationOperatorObj_const;
/// Global Registry for all LSQRegularizationOperators
extern ObjRegistry<LSQRegularizationOperator> gLSQRegularizationOperatorRegistry;
// Zdenek (end)

}//namespace
#endif //_LSQOBJNUM_H
