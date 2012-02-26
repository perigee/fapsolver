#ifndef COLALGO_HPP
#define COLALGO_HPP

#define NDEBUG
#include <cassert>


// boost library
#include <boost/function_output_iterator.hpp>
#include <boost/iterator/filter_iterator.hpp>
#include <boost/bind.hpp>
#include <boost/lambda/core.hpp>
#include <boost/ref.hpp>
#include <boost/function.hpp>
#include <boost/numeric/conversion/cast.hpp>


#include "Base.hpp"





using namespace std;




// ------------------------------------------- Pre-declaration of Constraint -----
template <typename T> struct VariableBase;
typedef VariableBase<AlgoRelated<dummyT> > Variable;
typedef boost::shared_ptr<Variable> VarPtr;



struct AbsConstraint;
typedef boost::shared_ptr<AbsConstraint> CtrPtr;
// -------------------------------------------------------------------------------



//================================================= Value class =================
// ==============================================================================

template <typename T>
        struct ValueBase : public T{


    ValueBase(int v):v_(v),p_(Neu){}
    ~ValueBase(){}

    void setPolar(Polar const p){p_ =p;}

    bool operator ==(const ValueBase& val) const{return v_ == val.v_ && p_ == val.p_;}
    bool operator < (ValueBase& val) const{return v_ < val.v_ ;}
    bool operator > (ValueBase& val) const{return v_ > val.v_;}

  void print()
  {
    if(p_ == V)
      cout << "+" ;
    
    if(p_ == H)
      cout << "-" ;
    
    if(p_ == Neu)
      cout << "o" ;
    
    cout << v_ << "\t";
  }
    int v_;
    Polar p_; // polarization of frequency
};

// Value class Declaration -------------------------------------------------------

template <typename T>  struct ValueBase;
typedef ValueBase<AlgoRelated<dummyT> >  Value;
typedef boost::shared_ptr<Value> ValPtr;

typedef std::vector<ValPtr> Domain;
// -------------------------------------------------------------------------------



/**
  *
  * general interface of constraint
  */
struct AbsConstraint{

    virtual std::string getID() = 0;
    //virtual void resign()=0;
    virtual void print() = 0;
    virtual void printViz() = 0;
    //virtual void recover() = 0;
    virtual void setCtrLevel(int level = 0) = 0;

    //
    virtual bool isAssigned() const = 0;
    virtual bool hasAssigned() const = 0;
    virtual void    resign() = 0;
    virtual bool satisfy() const = 0;
    virtual bool isViolated() const = 0;
    // check if given variable/value is satisfied
    virtual bool isViolated(VarPtr const v, ValPtr const a) const = 0;
    virtual bool findSolution() const = 0;



    virtual bool isIn(VarPtr var) const = 0;
    virtual bool isInStr(std::string& id) const = 0;
    virtual VarPtr getOther(VarPtr const v) = 0;



    // Gamma Arc
    virtual void initVarsGamma()=0;
    virtual void updateVarsGamma(VarPtr const var, ValPtr const a, ValPtr const b)=0;


    //virtual void setGamma(int i)=0;


    // Propagation
    virtual bool fc(VarPtr  var, bool& deadend) const = 0;
    virtual bool ac(VarPtr  var, bool& deadend) const = 0;



    // static gamma
    virtual Level getLevel() const= 0;
    virtual void  setLevel(Level l) = 0;
    virtual void setGammaS(int g) = 0;
    virtual int  getGammaS() const = 0;
    virtual void resetGammaS() = 0;
    virtual void setGammaD(int g) = 0;
    virtual int  getGammaD() const= 0;
    virtual void resetGamma() = 0;

    // static weight for learning
    virtual void setWeight(int w) = 0;
    virtual void addWeight() = 0;
    virtual int  getWeight() const = 0;
    virtual void resetWeight() = 0;

    // tabu tenure
    virtual void setTenure(int d) = 0;
    virtual int  getTenure() const= 0;
    virtual void resetTenure() = 0;
    virtual void reduceTenure() = 0;
    virtual bool inTenure() const = 0;

    virtual void resetAtts() = 0;
    virtual void clearAtts() = 0;


    // variable related operation
    virtual bool  sameVarsLevel() const = 0;
    virtual bool  varsLevelgreaterThan(Level const l) const = 0;
    virtual bool  varLevelgreaterThan(Level const l) const = 0;
    virtual Level getVarsLevel() const = 0;
    virtual void  setVarsLevel(Level l) = 0;
};



//================================================= Variable class ==============
// ==============================================================================


/*!
  Variable class: implemented basic functions
*/
template <typename T>
        struct VariableBase : public T{


    VariableBase(std::string& id):id_(id),
    assigned_(false),
    hasSolution_(false){}


    ValPtr  current(){assert(assigned_); return current_;}
    void    assign(ValPtr val){current_ = val; assigned_ = true;}
    void    resign(){if(assigned_) assigned_ = false;}
    bool    isAssigned() const{if(assigned_) assert(current_);return assigned_;}



    // only should be used in
    bool    isViolated() const{

        // is variable assigned
        assert(assigned_);

        // if the constraint level is normal
        BOOST_FOREACH(CtrPtr const c, ctrs_){
            if(c->getLevel() > disable){
                if(c->isViolated())
                    return true;
            }
        }

        return false;
    }


    int degree() const{return vars_.size();}
    int degree(Level l) const {
        //cout << "in degree" << endl;
        int cnt = 0;
        BOOST_FOREACH(VarPtr const v, vars_){
            if(v->getLevel() == l)
                ++cnt;
        }

        return cnt;
    }

  int domSize() const{return values_.size();}

    void    resetDomain(){

        std::for_each(values_.begin(),values_.end(),
                      boost::bind(&Value::resetAtts,_1));

        this->clearSolution();
    }

    void    clearDomain(){
        std::for_each(values_.begin(),values_.end(),
                      boost::bind(&Value::clearAtts,_1));

        this->clearSolution();
    }




    // clear domain and resign the variable
    void clear(){resetDomain();resign();}



//    void clearValuesGamma(){
//        std::for_each(values_.begin(), values_.end(),
//                      boost::bind(&Value::resetGamma,_1));

//    }


    bool linkToLevelVar(Level const l){
        BOOST_FOREACH(VarPtr const v, vars_){
            if(v->getLevel() >= l)
                return true;
        }

        return false;
    }


    bool underLevelCtr(Level const l){
        BOOST_FOREACH(CtrPtr const c, ctrs_){
            if(c->getLevel() == l)
                return true;
        }

        return false;
    }


    // solution related operations ==
    void setSolution(){if(isAssigned()) sol_ = current();if(!hasSolution_) hasSolution_ = true;}
    void getSolution(){if(hasSolution_) assign(sol_); hasSolution_ = false;}
    void clearSolution(){if(hasSolution_) hasSolution_ = false;}
    bool hasSolution() const{return hasSolution_;}



    void print(){
        cout << "VAR[" << id_ << "]=(" ;
        BOOST_FOREACH(ValPtr const a, values_){
            //if(a->getGamma() < 1)
            if(a->p_ == H)
                cout << "+" << a->v_ << " ";
            else
                cout << "-" << a->v_ << " ";
        }
        cout << ")" << endl;
    }

    friend std::ostream& operator<<(std::ostream& out, const Variable& v) {
        out << "VAR[" << v.id_ << "]" << v.values_.size();

        return out;
    }



    std::vector<VarPtr> vars_; ///< neighbor variables
    std::vector<CtrPtr> ctrs_; ///< neighbor constraints
    Domain values_; ///< values in domain
    std::string id_; ///< unique id of variable


private:
    bool assigned_; ///< see if variable is assigned
    bool hasSolution_; ///< store the solution during search

    ValPtr current_; ///< current assignment of variable
    ValPtr sol_;     ///< solution stored during tabu search


};






//================================================= Constraint Class=============
// ==============================================================================


/**
 * Coloring problem arc class which is implemented with the
 * gamma structure
 *
 */
template <class Policy>
        struct ColorComparator : public Policy{


    ColorComparator(){}// constructor for arc constraint



    //======================================== Compare function ===  BGN
    // only compare with two given values;
    bool satisfy(ValPtr const a, ValPtr const b) const{
        //cout << "call satisfy" << endl;
        return a->v_ != b->v_;
    }


  void print(string& a, string& b){}

    void setCtrLevel(int level = 0){}

    //======================================== Compare function ===  END

};


/**
 * FAPP
 *
 */
template <class Policy>
        struct FAPP_CI : public Policy{


    FAPP_CI():freq_(true),equal_(true),interval_(0){}// constructor for arc constraint

    void init(bool freq, bool equal, int interval = 0){
        freq_ = freq;
        equal_ = equal;
        interval_ = interval;
    }


    void setCtrLevel(int level = 0){}
    //======================================== Compare function ===  BGN
    // only compare with two given values;
    bool satisfy(ValPtr const a, ValPtr const b) const{

        if(freq_){
            if(equal_){
                return std::abs(a->v_-b->v_) == interval_;
            }

            return std::abs(a->v_-b->v_) != interval_;

        }else{
            if(equal_)
                return a->p_ == b->p_;

            return a->p_ != b->p_;
        }

    }



  void print(string& a, string& b){
    cout << "CI\t" << a << "\t" << b; 
        if(freq_)
            cout << "F\t";
        else
            cout << "P\t";

        if(equal_)
            cout << "E\t";
        else
            cout << "I\t";

        cout  << interval_ << endl;

    }

    //======================================== Compare function ===  END

private:
    bool freq_;
    bool equal_;
    int interval_;

};



template <class Policy>
        struct FAPP_CP : public Policy{


    FAPP_CP():level_(0){}// constructor for arc constraint

    void init(int level = 0){
        level_ = level;
    }


    void setCtrLevel(int level = 0){level_ = level;}



    void push_back(bool equ, int interval){
        if(equ)
            equVec_.push_back(interval);
        else
            inqVec_.push_back(interval);
    }

    //======================================== Compare function ===  BGN
    // only compare with two given values;
    bool satisfy(ValPtr const a, ValPtr const b) const{

        assert(level_ < inqVec_.size());

        // in case polarization inequality
        if(a->p_ != b->p_)
                return std::abs(a->v_-b->v_) >= inqVec_[level_];

        return std::abs(a->v_-b->v_) >= equVec_[level_];
    }

    //======================================== Compare function ===  END
    
  
  void print(string& a, string& b){

    //      cout << "level: " << level_ << endl;

      cout << "CE\t" << a << "\t" << b; 
        BOOST_FOREACH(int const v, equVec_){
	  cout << "\t" << v;
        }


        cout << endl;

      cout << "CD\t" << a << "\t" << b; 
        BOOST_FOREACH(int const v, inqVec_){
	  cout << "\t" << v;
        }

        cout << endl;

    }

private:
    int level_;
    std::vector<int> equVec_;
    std::vector<int> inqVec_;
};





/**
 * Coloring problem arc class which is implemented with the
 * gamma structure
 *
 */
template <class Policy>
        struct LRFAPComparator : public Policy{


    LRFAPComparator(){}// constructor for arc constraint

    void setCtrLevel(int level = 0){}


    //======================================== Compare function ===  BGN
    // only compare with two given values;
    bool satisfy(ValPtr const a, ValPtr const b) const{
        if(opt_ == '=')
            return std::abs(a->v_-b->v_) == diff_;

        if(opt_ == '>')
            return std::abs(a->v_-b->v_) > diff_;


        cout << "error" << endl;
        exit(1);
        return false;
    }

    void setParams( char opt, int diff){
        diff_ = diff;
        opt_ = opt;
    }

    void print(string& a, string& b){}



private:
    int diff_;
    char opt_; // operator


};





/**
 * Binary constraint base class, contains only necessary functions
 *
 */
template <class Policy>
        struct ConstraintBase2 :
        public Policy{


    ConstraintBase2(std::string id):id_(id){}

    template <typename InputIterator>
            void setVars(InputIterator first,InputIterator last,
                         string id1, string id2, CtrPtr ctr){
        bool hasSet = false;
        for(;first!=last;++first){
            if(id1.compare((*first)->id_) == 0){
                //if((*first)->id_ == id1){
                //cout << "fnd 1" << endl;
                v1_ = *first;
                if(hasSet)
                    break;
                hasSet = true;
            }

            if(id2.compare((*first)->id_) == 0){
                //if((*first)->id_ == id2){
                //cout << "fnd 2" << endl;
                v2_ = *first;
                if(hasSet)
                    break;
                hasSet = true;
            }
        }

        if(!hasSet){
            cout << "not set" << endl;
            exit(1);
        }

        //cout << "after set" <<endl;

        v1_->ctrs_.push_back(ctr);
        v1_->vars_.push_back(v2_);
        v2_->ctrs_.push_back(ctr);
        v2_->vars_.push_back(v1_);
    }


    std::string getID(){return id_;}


    bool satisfy() const{
        if(isAssigned())
            return Policy::satisfy(v1_->current(),v2_->current());
        return true;
    }

    void setCtrLevel(int level = 0){Policy::setCtrLevel(level);}

    bool isViolated() const{
        assert(this->isAssigned());
        return !satisfy();
    }

    bool isViolated(VarPtr const v, ValPtr const a) const{
        if(!this->isAssigned())
            return false;

        if(v == v1_){
            if(v2_->isAssigned())
                return ! Policy::satisfy(a,v2_->current());
        }else{
            if(v1_->isAssigned())
                return ! Policy::satisfy(v1_->current(),a);

        }

        return false;
    }


    bool findSolution() const{
        bool consistent = false;
        for(std::vector<ValPtr>::iterator air = v1_->values_.begin();
        air != v1_->values_.end();++air){
            for(std::vector<ValPtr>::iterator bir = v2_->values_.begin();
            bir != v2_->values_.end();++bir){
                if(Policy::satisfy(*air,*bir)){
                    consistent = true;
                    v1_->assign(*air);
                    v2_->assign(*bir);
                }
            }
        }

        return consistent;
    }



    // is given variable is under this contraint
    bool isIn(VarPtr const  var) const{
        return var == v1_ || var == v2_;
    }

    bool isInStr(std::string& id) const{
        if(id == v1_->id_)
            return true;
        if(id == v2_->id_)
            return true;
        return false;
    }

    // get opposite variable, try not to adopt this function
    VarPtr getOther(VarPtr const v){if(v==v1_) return v2_; return v1_;}


    bool    isAssigned() const{return v1_->isAssigned() && v2_->isAssigned();}
    bool    hasAssigned() const{return v1_->isAssigned() || v2_->isAssigned();}
    void    resign() {v1_->resign();v2_->resign();}




    // initialize the gamma system
    void initVarsGamma(){
        //cout << "f:  initVarsgamma" << endl;
        assert(v1_->isAssigned() && v2_->isAssigned());
        for(std::vector<ValPtr>::iterator a = v1_->values_.begin();
        a != v1_->values_.end(); ++a){


            if((*a)->getLevel() > disable && !Policy::satisfy(*a, v2_->current()))
                (*a)->setGammaS((*a)->getGammaS() + this->getGammaS());
        }

        for(std::vector<ValPtr>::iterator b = v2_->values_.begin();
        b != v2_->values_.end(); ++b){
            if((*b)->getLevel() > disable && !Policy::satisfy(v1_->current(), *b))
                (*b)->setGammaS((*b)->getGammaS() + this->getGammaS());
        }


    }

    void updateVarsGamma(VarPtr const var, ValPtr const a, ValPtr const b){
        assert(v1_->isAssigned() && v2_->isAssigned());
        // update opposite variables gamma system

        //cout << "f:  updateVarsGamma" << endl;

        bool order = true;
        VarPtr opposite;
        if(var == v1_){
            opposite = v2_;
        }else{
            order = false;
            opposite = v1_;
        }

        for(std::vector<ValPtr>::iterator ir = opposite->values_.begin();
        ir != opposite->values_.end(); ++ir){
            if((*ir)->getLevel() > disable){
                if(order){
                    if(!Policy::satisfy(a, *ir))
                        (*ir)->setGammaS((*ir)->getGammaS() - this->getGammaS());

                    if(!Policy::satisfy(b, *ir))
                        (*ir)->setGammaS((*ir)->getGammaS() + this->getGammaS());
                }else{
                    if(!Policy::satisfy(*ir, a))
                        (*ir)->setGammaS((*ir)->getGammaS() - this->getGammaS());

                    if(!Policy::satisfy(*ir, b))
                        (*ir)->setGammaS((*ir)->getGammaS() + this->getGammaS());
                }
            }
        }


    }



    // interface for foreign classes
    //@return true: no deadend, false: has deadend
    bool fc(VarPtr var, bool& deadend) const{

        if(this->isAssigned())
            return satisfy();

        assert(isIn(var));

        VarPtr opposite;
        if(var == v1_)
            opposite = v2_;
        else
            opposite = v1_;



        if(var->isAssigned()){
            return fcGamma(var, opposite, deadend);
        }else
            return fcGamma(opposite, var, deadend);
    }


    //@return true: no deadend, false: has deadend
    bool ac(VarPtr var, bool& hasDeleted) const{
        assert(isIn(var));
		

        VarPtr opposite;
        if(var == v1_)
            opposite = v2_;
        else
            opposite = v1_;

        return acGamma(var, opposite, hasDeleted);
    }

    
    bool  sameVarsLevel() const{return v1_->getLevel() == v2_->getLevel();}
    bool  varsLevelgreaterThan(Level const l) const{return (v1_->getLevel() > l && v2_->getLevel() > l);}
    bool  varLevelgreaterThan(Level const l) const{return (v1_->getLevel() > l || v2_->getLevel() > l);}
    Level getVarsLevel() const{return v1_->getLevel();}
    void  setVarsLevel(Level l){v1_->setLevel(l);v2_->setLevel(l);}

    void print(){
        //verify(v1_->current(), v2_->current());
       

        Policy::print(v1_->id_,v2_->id_);
        //Policy::printComparator();
        //cout << endl;
        //cout << v1_->id_ << " = " << (v1_->current())->v_ << " -- " << v2_->id_ << " = " <<  (v2_->current())->v_ << endl;
    }

    void printViz(){


//        if(Policy::getLevel() != critical)
//            cout << v1_->id_ <<  " -- " << v2_->id_ << "[color = grey];"<< endl;
//        else
//            cout << v1_->id_ <<  " -- " << v2_->id_ << ";"<< endl;

        cout << v1_->id_ <<  " -- " << v2_->id_ ;
    }


    friend std::ostream& operator<<(std::ostream& out, const ConstraintBase2& v) {
        out << "CTR[" << v.id_ << "]:" << (v.v1_)->id_ << " -- " << (v.v2_)->id_ << endl;
        return out;
    }



private:


    // @return true: consistent, false: deadend
    // @param a: the assigned variable
    // @param b: opposite variable
    bool fcGamma(VarPtr a, VarPtr b,  bool& hasDeleted) const{
        assert(a->isAssigned());

        hasDeleted = false;
        bool hasValid = false; // if there is consistent value
        for(std::vector<ValPtr>::iterator ir = b->values_.begin();
        ir != b->values_.end(); ++ir){

            // disable or gammaed value
            if((*ir)->getLevel() < normal || (*ir)->getGammaS() > 0)
                continue;

            bool consistent = true;

            // correct the order of variables in constraint
            if(a == v1_){
                consistent = Policy::satisfy(a->current(), *ir);
            }else{
                consistent = Policy::satisfy(*ir,a->current());
            }

            if(!consistent){
                if(!hasDeleted)
                    hasDeleted = true;

                // set static gamma on values
                (*ir)->setGammaS(this->getGammaS());
            }else{
                if(!hasValid)
                    hasValid = true;
            }

        }

        return hasValid;

    }


    // @param a: the variable need to be filtered
      bool acGamma(VarPtr a, VarPtr b, bool& hasDeleted) const{
          assert(!a->isAssigned());
          assert(!b->isAssigned());

          hasDeleted = false;
          bool hasValid = false; // if there is consistent value

          for(std::vector<ValPtr>::iterator xir = a->values_.begin();
          xir != a->values_.end(); ++xir){

              // disable or gammaed value
              if((*xir)->getLevel() < normal || (*xir)->getGammaS() > 0)
                  continue;


              bool hasConsistent = false;
              for(std::vector<ValPtr>::iterator ir = b->values_.begin();
              ir != b->values_.end(); ++ir){

                  // disable or gammaed value
                  if((*ir)->getLevel() < normal || (*ir)->getGammaS() > 0)
                      continue;

                  bool consistent = false;

                  // correct the order of variables in constraint
                  if(a == v1_){
                      consistent = Policy::satisfy(*xir, *ir);
                  }else{
                      consistent = Policy::satisfy(*ir,*xir);
                  }

                  if(consistent){
                      hasConsistent = true;
                      //(*xir)->print();
		      //(*ir)->print();
		      //cout << std::abs((*xir)->v_ - (*ir)->v_)
		      // << "\t ok"<< endl;
                      break;
                  }else{
                    //(*xir)->print();
                    //(*ir)->print();
                    //cout << " inconsistent"<< endl;
                    //cout << std::abs((*xir)->v_ - (*ir)->v_) << "\t INC"<< endl;
                  }

              }

              if(!hasConsistent){
                  (*xir)->setGammaS(this->getGammaS());
                  if(!hasDeleted)
                      hasDeleted = true;
              }else{
                  if(!hasValid)
                      hasValid = true;
              }



          }

          return hasValid; // return false if there were deadend
      }



    std::string id_; /// identification of constraint
    VarPtr v1_; ///< first variable
    VarPtr v2_; ///< second variable

};


// -------------------------------------------- Binary Constraint Class -----------
//typedef ConstraintBase2<WeightablePolicy<MarkablePolicy<Propagation2Gamma<GammaArc<Gamma<ColorComparator<EnablePolicyT<AbsConstraint> > > > > > > >
//        Constraint2;


typedef ConstraintBase2<AlgoRelated<ColorComparator<AbsConstraint> > >
        Constraint2;

typedef boost::shared_ptr<Constraint2> Constraint2Ptr;


/// FAPP
typedef ConstraintBase2<AlgoRelated<FAPP_CI<AbsConstraint> > >
        Constraint2CI;


typedef boost::shared_ptr<Constraint2CI> Constraint2CIPtr;


typedef ConstraintBase2<AlgoRelated<FAPP_CP<AbsConstraint> > >
        Constraint2CP;


typedef boost::shared_ptr<Constraint2CP> Constraint2CPPtr;



//typedef ConstraintBase2<WeightablePolicy<MarkablePolicy<Propagation2Gamma<GammaArc<Gamma<LRFAPComparator<EnablePolicyT<AbsConstraint> > > > > > > >
//        Constraint2FAP;


typedef ConstraintBase2<AlgoRelated<LRFAPComparator<AbsConstraint> > >
        Constraint2FAP;

typedef boost::shared_ptr<Constraint2FAP> Constraint2FAPPtr;
// --------------------------------------------------------------------------------



//================================================= Random Search ===============
// ==============================================================================

int getRandom(int min, int max){

    return min + int((double(rand())/RAND_MAX)*max);
}

/*! Randomly generate an initial solution
  */
template <typename InputIterator>
        void randomSolver(InputIterator first, InputIterator last){

    //srand((unsigned)time(0));

    for(;first!=last;++first){
        int randomNb = 0;
        do{
            randomNb = getRandom(0, ((*first)->values_).size()-1);
        }while((*first)->values_[randomNb]->getLevel() < normal);

        (*first)->assign((*first)->values_[randomNb]);

        // test
        //cout << (*first)->id_ << " = " << (*first)->current()->v_ << endl;
    }

}



/// ========================================================================
/// Statistic class
/// ========================================================================


struct Statistic{

    Statistic(){
        timeVec.reserve(CateEND);
        nbVec.reserve(CateEND);
        init();
    }

    void init(){
        total.restart();reset();
        for(int i = 0; i != CateEND;++i){
            timeVec.push_back(0.0);
            nbVec.push_back(0);}
    }
    //
    void restart(){t.restart();}

    void setTime(Cate const cate){
        timeVec[cate] += t.elapsed();
        nbVec[cate] += 1;
    }


    double getTotal(){return total.elapsed();}

    void print(){


        cout << "t:\t";
        BOOST_FOREACH(double const d,timeVec)
                cout << d << "\t";

        BOOST_FOREACH(int const d,nbVec)
                cout << d << "\t";

        cout << total.elapsed() << endl;
    }

private:

    void reset(){timeVec.clear();nbVec.clear();}

    std::vector<double> timeVec;
    std::vector<int> nbVec;
    boost::timer t;
    boost::timer total;
};

#endif // COLALGO_HPP
