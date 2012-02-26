#ifndef HEURISTIC_HPP
#define HEURISTIC_HPP


#include <functional>
#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <ctime>            // std::time
#include <deque>


//#include <boost/accumulators/accumulators.hpp>
//#include <boost/accumulators/statistics/stats.hpp>
//#include <boost/accumulators/statistics/mean.hpp>
//#include <boost/accumulators/statistics/moment.hpp>
//using namespace boost::accumulators;




#include "colAlgo.hpp"



///=============================== ONLY PRINT THE INFORMATION ================ BGN




template <typename Itr, typename Pred>
void printGraph(Itr first, Itr last, Pred pred, std::string& mode){
    cout << "=========== " << mode << endl;
    cout << "graph G { node [shape=point];" << endl;
    for(;first!=last;first++){
        if(pred(*first)){
            (*first)->printViz();
            cout << ";" << endl;

        }else{
            (*first)->printViz();
            cout << "[color = grey];" << endl;
        }
    }
    cout << "}" << endl;
}


///=============================== ONLY PRINT THE INFORMATION ================ END

struct Acfc{


    template <typename InputIterator>
    bool solve(InputIterator first, InputIterator last){
        bool deadend = false;

        bool hasDeleted = false;

        init(first,last);
        //std::sort(qlist.begin(),qlist.end(),compareVarsDegrees); // sort
        //std::sort(qlist.begin(),qlist.end(),
        //          boost::bind(std::greater<int>(),
        //			  boost::bind(&Variable::degree, _1),
        //			  boost::bind(&Variable::degree, _2))); // sort


        while(!qlist.empty()){

            //cout << qlist.size() << endl;
            // get last variable and pop it
            VarPtr current = qlist.back();
            current->resetTenure();
            qlist.pop_back();

            //verify if there is any value removed.
            hasDeleted = false;
            if(revise(current,hasDeleted)){

                // if no deleted value, then continue
                if(!hasDeleted)
                    continue;

                //add all its valid neighbors which is not in Qlist
                BOOST_FOREACH(VarPtr const var, current->vars_){
                    if(!var->inTenure()){
                        var->setTenure(10);
                        qlist.push_back(var);
                    }
                }

            }else{
                //if the variable becomes a deadend
                deadend = true;
                break;
            }

        }

        // reset all effects
        //std::for_each(first,last, boost::bind(&Variable::resetDomain,_1));
        std::for_each(first,last, boost::bind(&Variable::resetAtts,_1));
        qlist.clear();
        return !deadend;

    }



private:

    bool revise(VarPtr const var, bool& hasDeleted){
        hasDeleted = false;

        bool consistent = false;

        //traverse all valid values
        BOOST_FOREACH(CtrPtr const ctr, var->ctrs_){
            bool hasD = false;
            if(ctr->ac(var,hasD)){
                if(hasD && !hasDeleted)
                    hasDeleted = true;
            }else
                return false;

        }

        return true;
    }


    template <typename InputIterator>
    void init(InputIterator first, InputIterator last){
        if(!qlist.empty())
            qlist.clear();

        for(;first!=last;++first){
            (*first)->setTenure(10);
            qlist.push_back(*first);

        }
    }


    std::deque<VarPtr>  qlist;

};





template <typename Tptr>
struct Marked{

    bool operator () (Tptr const t){
        return t->getLevel() > normal;
    }
};


template <typename Tptr>
struct NonMarked{

    bool operator () (Tptr const t){
        return t->getLevel() != marked;
    }
};


template <typename Tptr>
struct Normal{

    bool operator () (Tptr const t){
        return t->getLevel() > disable;
    }
};


template <typename Tptr>
struct Critical{

    bool operator () (Tptr const t){
        return t->getLevel() == critical;
    }
};

template <typename Tptr>
struct NonCritical{

    bool operator () (Tptr const t){
        return t->getLevel() != critical;
    }
};


template <typename Tptr>
struct All{

    bool operator () (Tptr const t){
        return true;
    }
};


template <typename Element>
struct Counter_Normal
{

    int operator ()(int result, Element const e){
        if(e->getLevel() > disable)
            ++result;

        return result;
    }
};

template <typename Element>
struct Counter_NormalNoGamma
{

    int operator ()(int result, Element const e){
        if(e->getLevel() > disable && e->getGammaS() < 1)
            ++result;

        return result;
    }
};


template <typename Element>
struct Counter_Marked
{

    int operator ()(int result, Element const e){
        if(e->getLevel() == marked)
            ++result;

        return result;
    }
};


template <typename Element>
struct Counter_Critical
{

    int operator ()(int result, Element const e){
        if(e->getLevel() == critical)
            ++result;

        return result;
    }
};



struct Counter_DegreeMean
{

    int operator ()(int result, VarPtr const e){
        int d = e->vars_.size();
        if(result < 1)
            return d;
        else
            return (d+result)/2;

    }
};


/// ===============================================================================================
/// TabuCol
/// Constructive method
/// ===============================================================================================

// verify the minimal gamma value which is not in tabu
inline bool lessValue(ValPtr const a, ValPtr const b){
    if(a->getLevel() < normal)
        return false;

    if(b->getLevel() < normal)
        return true;


    if(!a->inTenure() && !b->inTenure()){
        if(a->getGammaS() == b->getGammaS())
            return rand()/static_cast<double>(RAND_MAX) < 0.5;

        return a->getGammaS() < b->getGammaS();
    }



    if(a->inTenure() && b->inTenure())
        return rand()/static_cast<double>(RAND_MAX) < 0.5;

    return !a->inTenure();
}

// compare two delta values of var,val pairs
inline bool  min_delta_value(VarPtr const va, ValPtr const a,
                             VarPtr const vb, ValPtr const b){

    // in case both values in tabu is false
    if(!a->inTenure() || !b->inTenure()){
        if(a->inTenure())
            return false;

        if(b->inTenure())
            return true;

        if((a->getGammaS()-va->current()->getGammaS())
                == (b->getGammaS() - vb->current()->getGammaS()))
            return rand()/static_cast<double>(RAND_MAX) < 0.5;

        return (a->getGammaS()-va->current()->getGammaS())
                < (b->getGammaS() - vb->current()->getGammaS());
    }

    // in case both values in tabu
    return rand()/static_cast<double>(RAND_MAX) < 0.5;

}

// get minimal gamma variable/value
template <bool weightable = false>
struct MinGamma{



    MinGamma():set_(false),
        hasViolated_(false),
        nbViolated_(0),
        domSize_(0),
        freq_(0),
        tenure_(0){}

    void find(VarPtr v){
        //cout << "in miniGamma" << endl;
        // skip the low level variable
        if(v->getLevel() <= disable)
            return;

        // find
        if(v->isViolated()){
            if(!hasViolated_)
                hasViolated_ = true;
        }else
            return;

        ++nbViolated_;

        if(weightable)
            v->addWeight();

        // only to avoid current value
        v->current()->setTenure(7);

        if(!set_){
            domSize_ = v->values_.size();
        }

        ValPtr a = *std::min_element(v->values_.begin(),
                                     v->values_.end(), lessValue);

        //        if(a == v->current())
        //            cout << "find current" << endl;

        // reduce tenure
        std::for_each(v->values_.begin(), v->values_.end(),
                      boost::bind(&Value::reduceTenure,_1));



        if(!set_ || min_delta_value(v,a,var_,val_)){
            val_ = a;
            var_ = v;
            if(!set_)
                set_ = true;
        }

    }



    bool isSet() const{return set_;}
    int  getDelta() const{return val_->getGammaS() - var_->current()->getGammaS();}
    int  getNbViolated() const{return nbViolated_;}

    bool hasViolated(){return hasViolated_;}
    void reset(){set_ =false; hasViolated_=false;nbViolated_ = 0;}
    void allreset(){reset();freq_ = 0;tenure_=0;}



    // return the gamma difference of two values
    template <typename Pred>
    void update(Pred pred){


        if(set_){

            //            if(val_->inTenure())
            //                cout << "tabu value"<< endl;

            // update gamma
            for(std::vector<CtrPtr>::iterator ir = var_->ctrs_.begin();
                ir != var_->ctrs_.end(); ++ir){
                if(pred(*ir)){
                    // update opposite variable's value gamma == under construction
                    (*ir)->updateVarsGamma(var_, var_->current(), val_);
                }
            }

            ValPtr tmp = var_->current();
            val_->resetTenure();
            var_->assign(val_);
            val_ = tmp;
        }
    }


    // should be invoked after update
    void setTenure(bool improved){val_->setTenure(getTenure(improved));}



private:


    int getTenure(bool improved){

        /// strategy 1: coloring L + lambda(ViolatedNb), L is randomly chosen between 0,9
        return getRandom(0,9) + 0.6 * nbViolated_;

        /// strategy 2: for fap
        if(!improved){



            ++freq_;
            if(tenure_ == 0 || tenure_/nbViolated_ > domSize_)
                tenure_ = domSize_;

            if(freq_ > nbViolated_ ){
                if((tenure_/domSize_) < nbViolated_ ){
                    tenure_ += domSize_;
                    //++tenure_;
                    freq_ = 0;
                }
            }


        }else{
            freq_ = 0;
            if(tenure_ == 0)
                tenure_ = domSize_;
        }

        //return nbViolated_;

        return tenure_;


    }

    bool set_;  // flag
    bool hasViolated_;
    int  nbViolated_;
    int  domSize_;
    int  freq_;
    int  tenure_;
    VarPtr var_;
    ValPtr val_;
};



// tabucol implementation,
// before entering the function, all gamma on constraints
// should be set.
template <typename VarItr, typename CtrItr,
          typename Move, typename PredVar, typename PredCtr>
bool tabuCol(VarItr vfirst, VarItr vlast,
             CtrItr cfirst, CtrItr clast,
             Move& move,PredVar predvar, PredCtr predctr){


    // filter ===========
    typedef boost::filter_iterator<PredCtr, CtrItr> CtrItrFilter;
    CtrItrFilter filter_cfirst(predctr, cfirst, clast);
    CtrItrFilter filter_clast(predctr, clast,clast);

    typedef boost::filter_iterator<PredVar, VarItr> VarItrFilter;
    VarItrFilter filter_vfirst(predvar, vfirst, vlast);
    VarItrFilter filter_vlast(predvar, vlast,vlast);

    // initialize gamma system
    std::for_each(filter_cfirst,filter_clast,
                  boost::bind(&AbsConstraint::initVarsGamma,_1));

    //cout << "after initial" << endl;

    // set of parameters
    register int MAXIteration = 1000;
    register int iteration = 0;
    register int cost = 0;
    register int currentCost = 0;
    // iterations
    bool firstUp = false;
    while(iteration < MAXIteration){
        std::for_each(filter_vfirst, filter_vlast,
                      boost::bind(&Move::find,&move,_1));

        // calculate the difference of gamma

        if(!move.hasViolated()){
            // reserve current assignment and clear the previous solution
            std::for_each(filter_vfirst, filter_vlast,
                          boost::bind(&Variable::clearSolution,_1));
            break;
        }

        int delta = move.getDelta();
        currentCost += delta;
        move.update(predctr);


        // update  solution if necessary

        if(delta < 0){


            if(!firstUp || cost > currentCost){
                std::for_each(filter_vfirst,filter_vlast,
                              boost::bind(&Variable::setSolution,_1));
                if(firstUp){
                    cost = currentCost;
                    move.setTenure(true);
                }



                iteration = 0;

            }else{
                move.setTenure(false);
            }



        }else{
            if(!firstUp){
                firstUp = true;
                cost = 100;
                currentCost = 100;
            }

            move.setTenure(false) ;
            ++iteration;
        }

	//cout << cost << "\t "
	//              << currentCost << endl;

        move.reset();


    }

    //bool feasible = ! move.hasViolated();
    move.allreset();

    // clear gamma
    std::for_each(filter_vfirst, filter_vlast, boost::bind(&Variable::clearDomain,_1));
    std::for_each(filter_vfirst, filter_vlast, boost::bind(&Variable::getSolution,_1));
    // verify the solution
    for(;filter_cfirst != filter_clast;++filter_cfirst){
        if((*filter_cfirst)->isViolated())
            return false;
    }

    return true;

}




template <typename VarItr, typename CtrItr,
          typename Move, typename PredVar, typename PredCtr>
bool localsearch(VarItr vfirst, VarItr vlast,
                 CtrItr cfirst, CtrItr clast,
                 Move& move,PredVar predvar, PredCtr predctr){


    // filter ===========
    typedef boost::filter_iterator<PredCtr, CtrItr> CtrItrFilter;
    CtrItrFilter filter_cfirst(predctr, cfirst, clast);
    CtrItrFilter filter_clast(predctr, clast,clast);

    typedef boost::filter_iterator<PredVar, VarItr> VarItrFilter;
    VarItrFilter filter_vfirst(predvar, vfirst, vlast);
    VarItrFilter filter_vlast(predvar, vlast,vlast);

    // initialize gamma system
    std::for_each(filter_cfirst,filter_clast,
                  boost::bind(&AbsConstraint::initVarsGamma,_1));


    //cout << "size: " << std::distance(cfirst,clast) << "\t" << std::distance(filter_cfirst,filter_clast) << endl;
    //cout << "size: " << std::distance(vfirst,vlast) << "\t" << std::distance(filter_vfirst,filter_vlast) << endl;

    if(std::distance(cfirst,clast) != std::distance(filter_cfirst,filter_clast)){
        for(;cfirst != clast; ++cfirst){
            cout << "level: " << (*cfirst)->getLevel() << endl;
        }
        exit(0);
    }
    
    // set of parameters
    register int MAXIteration = 1000;
    register int iteration = 0;
    register int cost = 0;
    register int currentCost = 0;
    // iterations
    bool firstUp = false;
    while(iteration < MAXIteration){
        std::for_each(filter_vfirst, filter_vlast,
                      boost::bind(&Move::find,&move,_1));

        // calculate the difference of gamma

        if(!move.hasViolated()){
            cout << "local search found no violated" << endl;
            // reserve current assignment and clear the previous solution
            std::for_each(vfirst, vlast,
                          boost::bind(&Variable::clearSolution,_1));
            break;
        }

        int delta = move.getDelta();
        currentCost += delta;
        move.update(predctr);


        // update  solution if necessary

        if(delta < 0){


            if(!firstUp || cost > currentCost){
                std::for_each(vfirst,vlast,
                              boost::bind(&Variable::setSolution,_1));
                if(firstUp){
                    cost = currentCost;
                    move.setTenure(true);
                }



                iteration = 0;

            }else{
                move.setTenure(false);
            }



        }else{
            if(!firstUp){
                firstUp = true;
                cost = 100;
                currentCost = 100;
                break; // local search
            }

            move.setTenure(false) ;
            ++iteration;
        }

        //        cout << cost << "\t "
        //                << currentCost << endl;

        move.reset();


    }

    //bool feasible = ! move.hasViolated();
    move.allreset();

    // clear gamma
    std::for_each(vfirst, vlast, boost::bind(&Variable::clearDomain,_1));
    std::for_each(vfirst, vlast, boost::bind(&Variable::getSolution,_1));
    // verify the solution
    for(;cfirst != clast;++cfirst){
        //cout << "gamma: " << (*cfirst)->getGammaS() << endl;
        if((*cfirst)->isViolated())
            return false;
    }

    return true;

}



inline bool hasValidNeighbor(VarPtr const v){
    for(std::vector<VarPtr>::iterator ir = v->vars_.begin();
        ir !=  v->vars_.end(); ++ir){
        if((*ir)->getLevel() > disable)
            return true;
    }

    return false;
}


// function object for counting the available values

struct Counter_validValue{
    int operator () (int result, ValPtr const t){
        if(t->getGammaS() > 0 || t->getLevel() < normal)
            return result;

        return ++result;
    }
};
// make for constructive Algorithm,
// for choosing a variable according to MRV
// only verify non-assigned variable
// only choose the variable adjacent with assigned one


int validDomainSize(VarPtr const v){



    // counting only the valid and non-gamma value
    return std::accumulate(v->values_.begin(), v->values_.end(),
                           0, Counter_NormalNoGamma<ValPtr>());


}


/// ===============================================================================================
/// breakout algorithm
/// ===============================================================================================

inline void increaseViolatedGamma(CtrPtr c, int g = 1){
    if(c->isViolated()){
        if(g < 2)
            c->setGammaS(c->getGammaS() +1);
        else{
            //cout << "setGamma " << g << endl;
            c->setGammaS(g);
        }
    }
}

template <typename VarItr, typename CtrItr, typename GammaSys, typename PredVar, typename PredCtr>
bool   breakout(VarItr vfirst, VarItr vlast, CtrItr cfirst, CtrItr clast,
                GammaSys& minGamma,PredVar varPred, PredCtr ctrPred){

    //cout << "in breakout " << endl;
    // filter ===========
    typedef boost::filter_iterator<PredVar, VarItr> VarItrFilter;
    VarItrFilter filter_vfirst(varPred, vfirst, vlast);
    VarItrFilter filter_vlast(varPred, vlast,vlast);

    typedef boost::filter_iterator<PredCtr, CtrItr> CtrItrFilter;
    CtrItrFilter filter_cfirst(ctrPred, cfirst, clast);
    CtrItrFilter filter_clast(ctrPred, clast,clast);


    // find maximal degree element
    VarItrFilter itr = std::max_element(filter_vfirst, filter_vlast, 
					boost::bind(std::less<int>(),
						    boost::bind(&Variable::degree,_1),
                                                    boost::bind(&Variable::degree,_2)));




    register int breakoutItr = (*itr)->vars_.size();
    //cout << "maxDegree: "  << breakoutItr;
    //breakoutItr = std::accumulate(filter_vfirst,filter_vlast,0,Counter_DegreeMean());
    // cout << "\tMeanDegree: " << breakoutItr << endl;
    // exit(0);
    //int l = breakoutItr;
    bool feasible = false;
    for(;breakoutItr > 0; --breakoutItr){

        //        cout << "iteration: " << breakoutItr
        //                << "\t" << std::accumulate(filter_cfirst, filter_clast, 0,Counter_violated<CtrPtr>())
        //                << endl;


        /// strategy 1: using tabucol as core
        //        feasible = tabuCol(filter_vfirst,filter_vlast, filter_cfirst, filter_clast,
        //                           minGamma,varPred, ctrPred);

        /// strategy 2: using tabucol as core
        //cout << "before local search in breakout" << endl;
        feasible = localsearch(vfirst,vlast, cfirst, clast,
                               minGamma,varPred, ctrPred);

        if(feasible)
            return true;


        // increment the gamma on violated constraints
        std::for_each(filter_cfirst,filter_clast,
                      boost::bind(increaseViolatedGamma,_1,1));


    }

    return false;
}


template <typename VarItr, typename CtrItr, typename GammaSys, typename PredVar, typename PredCtr>
bool   breakoutViolated(VarItr vfirst, VarItr vlast, CtrItr cfirst, CtrItr clast,
                        GammaSys& minGamma,PredVar varPred, PredCtr ctrPred){


    // filter ===========
    //cout << "in breakoutViolated" << endl;
    typedef boost::filter_iterator<PredVar, VarItr> VarItrFilter;
    VarItrFilter filter_vfirst(varPred, vfirst, vlast);
    VarItrFilter filter_vlast(varPred, vlast,vlast);

    typedef boost::filter_iterator<PredCtr, CtrItr> CtrItrFilter;
    CtrItrFilter filter_cfirst(ctrPred, cfirst, clast);
    CtrItrFilter filter_clast(ctrPred, clast,clast);

    // apply gamma on critical elements
    //bool hasMarked = false;
    //    for(CtrItrFilter current = filter_cfirst; current != filter_clast;++current){
    //      CtrPtr c = *current;
    //      if(c->getLevel() > normal){
    //	      c->setGammaS(100);
    //	      c->setLevel(marked);
    //	      c->setVarsLevel(marked);
    //      }

    //      if(c->getLevel() > normal){
    //	if(!hasMarked)
    //	  hasMarked = true;
    //      }
    //    }


    bool feasible = false;
    while(!feasible){

        //cout << "before tabu" << endl;
        feasible = tabuCol(vfirst,vlast,cfirst,clast,
                           minGamma,varPred, ctrPred);

        if(feasible)
            return true;
        else{
            // evaluate if there is gammaed values violated
            for(CtrItrFilter current = filter_cfirst; current != filter_clast;++current){
                if((*current)->getLevel() == marked && (*current)->isViolated())
                    return false;
            }
        }



        for(CtrItrFilter ir = filter_cfirst; ir != filter_clast; ++ir){
            if((*ir)->getLevel() != marked && (*ir)->isViolated()){
                (*ir)->setGammaS(100);
                (*ir)->setLevel(marked);
                (*ir)->setVarsLevel(marked);
            }
        }
	
        //	if(hasMarked){
        //            for(CtrItrFilter ir = filter_cfirst; ir != filter_clast; ++ir){
        //                if((*ir)->getLevel() != marked && (*ir)->isViolated() && (*ir)->varLevelgreaterThan(normal)){
        //                    (*ir)->setGammaS(100);
        //                    (*ir)->setLevel(marked);
        //                    (*ir)->setVarsLevel(marked);
        //                }
        //            }
        //	}else{
        //            CtrPtr chosenCtr ;
        //            bool isSetCtr = false;
        //            for(CtrItrFilter ir = filter_cfirst; ir != filter_clast; ++ir){
        //                if((*ir)->getLevel() != marked && (*ir)->isViolated()){

        //                    if(!isSetCtr || rand()/static_cast<double>(RAND_MAX) < 0.5)
        //                        chosenCtr = *ir;
        //                }
        //            }

        //            chosenCtr->setGammaS(100);
        //            chosenCtr->setLevel(marked);
        //            chosenCtr->setVarsLevel(marked);
        //            hasMarked = true;
        //	}

        /// applying the saturation

        for(CtrItrFilter current = filter_cfirst; current != filter_clast; ++current){
            CtrPtr c = *current;
	    
	    if(c->getLevel() > normal)
                continue;

            if(c->sameVarsLevel()){
                if(c->getVarsLevel() > normal){
                    c->setGammaS(100);
                    c->setLevel(marked);
                }
            }
        }

    }


    return true;
}





template <typename VarItr, typename CtrItr, typename GammaSys, typename PredVar, typename PredCtr>
bool   breakoutViolatedConstraint(VarItr vfirst, VarItr vlast, CtrItr cfirst, CtrItr clast,
                                  GammaSys& minGamma,PredVar varPred, PredCtr ctrPred){


    // filter ===========
    //cout << "in breakoutViolated" << endl;
    typedef boost::filter_iterator<PredVar, VarItr> VarItrFilter;
    VarItrFilter filter_vfirst(varPred, vfirst, vlast);
    VarItrFilter filter_vlast(varPred, vlast,vlast);

    typedef boost::filter_iterator<PredCtr, CtrItr> CtrItrFilter;
    CtrItrFilter filter_cfirst(ctrPred, cfirst, clast);
    CtrItrFilter filter_clast(ctrPred, clast,clast);


    for(CtrItrFilter current = filter_cfirst; current != filter_clast;++current){
        if((*current)->getLevel() == marked && (*current)->isViolated()){
            cout << "before into mark c" << endl;
            return false;
        }
    }


    bool feasible = false;
    while(!feasible){

        //cout << "before tabu" << endl;
        feasible = tabuCol(vfirst,vlast,cfirst,clast,
                           minGamma,varPred, ctrPred);

        if(feasible){
            //cout << "constraint true" << endl;
            return true;
        }else{
            // evaluate if there is gammaed values violated
            for(CtrItrFilter current = filter_cfirst; current != filter_clast;++current){
                if((*current)->getLevel() == marked && (*current)->isViolated()){
                    //cout << "mark criteria" << endl;
                    return false;
                }
            }
        }


        //        int count = 0;
        for(CtrItrFilter ir = filter_cfirst; ir != filter_clast; ++ir){
            if((*ir)->getLevel() != marked && (*ir)->isViolated()){
                (*ir)->setGammaS(100);
                (*ir)->setLevel(marked);
                (*ir)->setVarsLevel(marked);
            }

            //            if((*ir)->getLevel() == marked){
            //                ++count;
            //            }
        }

        //        cout << "total marked: " << count << endl;


    }


    return true;
}





template <typename VarItr, typename CtrItr, typename GammaSys, typename PredVar, typename PredCtr>
bool   prefiltering(VarItr vfirst, VarItr vlast, CtrItr cfirst, CtrItr clast,
                    GammaSys& minGamma,PredVar varPred, PredCtr ctrPred){


    randomSolver(vfirst,vlast);

    // filter ===========
    //cout << "in breakoutViolated" << endl;
    typedef boost::filter_iterator<PredVar, VarItr> VarItrFilter;
    VarItrFilter filter_vfirst(varPred, vfirst, vlast);
    VarItrFilter filter_vlast(varPred, vlast,vlast);

    typedef boost::filter_iterator<PredCtr, CtrItr> CtrItrFilter;
    CtrItrFilter filter_cfirst(ctrPred, cfirst, clast);
    CtrItrFilter filter_clast(ctrPred, clast,clast);


    for(CtrItrFilter current = filter_cfirst; current != filter_clast;++current){
        if((*current)->getLevel() == marked && (*current)->isViolated()){
            cout << "before into mark c" << endl;
            return false;
        }
    }


    bool feasible = false;
    //int countC = 0;
    while(!feasible){

        cout << "before tabu" << endl;
//        feasible = localsearch(vfirst,vlast,cfirst,clast,
//                               minGamma,varPred, ctrPred);

        feasible = tabuCol(vfirst,vlast,cfirst,clast,
                               minGamma,varPred, ctrPred);

        if(feasible){
            //cout << "constraint true" << endl;
            return true;
        }else{
            // evaluate if there is gammaed values violated
            for(CtrItrFilter current = filter_cfirst; current != filter_clast;++current){
                if((*current)->getLevel() == marked && (*current)->isViolated()){
                    cout << "mark criteria" << endl;
                    return false;
                }
            }
        }


        // mark the violated constraints
        for(CtrItrFilter ir = filter_cfirst; ir != filter_clast; ++ir){
            if((*ir)->getLevel() != marked && (*ir)->isViolated()){
                (*ir)->setGammaS(100);
                (*ir)->setLevel(marked);
                (*ir)->setVarsLevel(marked);
                //cout << countC++ << endl;
            }

        }


        // saturation
        for(CtrItrFilter ir = filter_cfirst; ir != filter_clast; ++ir){
            if((*ir)->getLevel() != marked){
                if((*ir)->sameVarsLevel() && (*ir)->getVarsLevel() == marked){
                    (*ir)->setGammaS(100);
                    (*ir)->setLevel(marked);
                }
            }

        }

    }


    return true;
}



inline bool validValue(ValPtr const a){
    return a->getLevel() > disable && a->getGammaS() < 1;
}

/// before enteing this function, one or two variables should be fixed
/// the rest variables should be disable
template <typename VarItr, typename CtrItr>
bool constructMRV(VarItr vfirst, VarItr vlast, CtrItr cfirst, CtrItr clast){

    // resign all disable variavble
    for(VarItr current = vfirst; current != vlast; ++current){
        if((*current)->getLevel() < normal){
            if((*current)->isAssigned()){
                (*current)->setSolution();
                (*current)->resign();
            }
        }
    }

    // fc
    fc_consistency(vfirst,vlast,Normal<VarPtr>(),All<CtrPtr>());
    
    bool consistentSolution = false;
    while(true){
        
        int validSize = 1;
        bool isSet = false;
        VarPtr current;
        for (VarItr xcurrent = vfirst;xcurrent != vlast;++xcurrent ){

	    if ((*xcurrent)->getLevel() > disable  || !(*xcurrent)->linkToLevelVar(normal))
                continue;

            int size = std::accumulate((*xcurrent)->values_.begin(), (*xcurrent)->values_.end(),
                                       0, Counter_NormalNoGamma<ValPtr>());

	    
            if (!isSet || validSize > size){
                if(!isSet)
                    isSet = true;

                current = *xcurrent;
                validSize = size;
            }
        }

        if(!isSet){
            //            cout << "no set " << endl;
            //            for (CtrItr current = cfirst; current != clast; ++current){
            //                if ((*current)->getLevel() < normal && (*current)->varsLevelgreaterThan(disable)){
            //                    (*current)->setLevel(normal);
            //                }
            //            }

            //            consistentSolution = true;
            //            int cnt = 0;
            //            for (CtrItr current = cfirst; current != clast; ++current){
            //                if ((*current)->getLevel() > disable && (*current)->isViolated()){
            //                    consistentSolution = false;
            //                    break;
            //                }else{
            //                    if ((*current)->getLevel() > disable )
            //                        ++cnt;
            //                }
            //            }

            //            int cntVar = 0;
            //            for(VarItr ir = vfirst; ir != vlast; ++ir){
            //                if((*ir)->getLevel() > disable)
            //                    ++cntVar;
            //            }


            //            if(consistentSolution){
            //                cout << "feasible by MRV var: " << cntVar << "\t ctr: " << cnt <<  endl;
            //            }else{
            //                cout << "nonfeasible by MRV" << endl;
            //            }

            consistentSolution = true;
            break;
        }


        // assigne a value on chosen variable
        if((current)->getLevel() == disable){
            (current)->setLevel(normal);

            // if find consistent value
            std::vector<ValPtr>::iterator aItr = std::find_if((current)->values_.begin(),
                                                              (current)->values_.end(),
                                                              validValue);



            if(aItr != (current)->values_.end()){
                //cout << "find value: " << (*aItr)->v_ << "\t" << (*aItr)->getGammaS() << endl;
                (current)->assign(*aItr);
                if(current->hasSolution())
                    current->clearSolution();

                //cout << "sign " << (*aItr)->v_ << " on " << (current)->id_ << endl;
                // fc consistency
                bool delted = false;
                BOOST_FOREACH(CtrPtr const c, (current)->ctrs_){
                    assert(c->getGammaS() == 1);
                    c->fc(current,delted);
                }

            }else{
                if(current->hasSolution())
                    current->getSolution();
                else
                    current->assign(current->values_.front());
                break;
            }
        }
    }


    for(VarItr current = vfirst; current != vlast; ++current){

        if((*current)->getLevel() < normal){
            if((*current)->hasSolution()){
                (*current)->getSolution();
            }
        }

        (*current)->resetDomain();
    }

    return consistentSolution;
}






/// ===============================================================================================
/// consistency verificator
/// MAC algorithm
/// ===============================================================================================
// entering the assigned variables
// forward checking
template <typename Iter, typename PredVar, typename PredCtr>
bool fc_consistency(Iter first, Iter last,
                    PredVar varPred, PredCtr ctrPred){

    //cout << "in fc " << endl;

    bool hasDeleted = false;
    for(;first != last; ++ first){

        if(!varPred(*first) || !(*first)->isAssigned())
            continue;

        for(std::vector<CtrPtr>::iterator ir = (*first)->ctrs_.begin();
            ir != (*first)->ctrs_.end(); ++ir){
            if(!ctrPred(*ir))
                continue;

            if(!(*ir)->fc(*first, hasDeleted))
                return false;
        }
    }

    return true;

}


// @return true: no deadend, false: has deadend
template <typename PredVar, typename PredCtr>
bool revise(VarPtr var, PredVar varPred, PredCtr ctrPred){


    //cout << "on: " << var->id_ << endl;
    //    if(!varPred(var))
    //        return true;

    bool hasDeleted = false;

    // in revise procedure
    std::vector<CtrPtr>::iterator first = var->ctrs_.begin();
    std::vector<CtrPtr>::iterator last = var->ctrs_.end();

    // verify all valid constraints, see if there is value erased

    for(;first != last; ++first){

        if(ctrPred(*first)){
            // only verify the unsigned constraint
            if(!(*first)->hasAssigned()){
                bool hasD = false;
                //(*first)->print();
                if(!(*first)->ac(var,hasD)){
                    return false;		    
                }else{
                    if(hasD && !hasDeleted){
                        hasDeleted = true;
                    }
                }

            }

        }
    }

    // if there is pruning, then revise
    if(hasDeleted){
        BOOST_FOREACH(VarPtr const v, var->vars_){
            if(varPred(v)){
                if(!v->isAssigned()){
                    if(!revise(v,varPred, ctrPred))
                        return false;
                }
            }
        }
    }

    return true;
}

/// disable  gammaed value and reset its gamma
void removeValue(ValPtr a){if(a->getGammaS() > 0){a->setLevel(disable); a->resetGammaS();}}

// entering non assigned variables, dealing with marked variable only
// @return true: arc-consistent; false: deadend
template <typename Iter, typename PredVar, typename PredCtr>
bool arc_consistency(Iter first, Iter last,
                     PredVar varPred, PredCtr ctrPred){

    //cout << "in ac " << endl;

    // only deal with the gamma free value and test only gamma constraints
    for(Iter current = first; current != last; ++current){

        if(varPred(*current)){

            assert(!(*current)->isAssigned());

            if(!revise(*current, varPred, ctrPred))
                return false;
        }
    }


    return true;

}


template <typename VarItr, typename PredVar, typename PredCtr>
bool mac_solver(VarItr first, VarItr last,
                PredVar varPred, PredCtr ctrPred,bool& depass){



    boost::timer timer;

    // filter the variables
    typedef boost::filter_iterator<PredVar, VarItr> VarItrFilter;
    VarItrFilter filter_vfirst(varPred, first, last);
    VarItrFilter filter_vlast(varPred, last, last);




    int size = std::distance(filter_vfirst, filter_vlast);
    if( size == 0) return true;

    std::for_each(filter_vfirst, filter_vlast,
                  boost::bind(&Variable::clear,_1));

    std::vector<VarPtr> future;
    future.reserve(size);
    std::vector<VarPtr> past;
    past.reserve(size);


    std::copy(filter_vfirst, filter_vlast, std::back_inserter(future));
    //cout << "mac: vars " << future.size() << endl;









    /// weight heuristic to order the variables ================================= BGN
    //    std::sort(future.begin(), future.end(),
    //              boost::bind(std::greater<int>(),
    //                          boost::bind(&Variable::getWeight,_1),
    //                          boost::bind(&Variable::getWeight,_2)));

    //    std::for_each(future.begin(), future.end(),
    //                  boost::bind(&Variable::resetWeight,_1));


    // only for coloring
    std::sort(future.begin(), future.end(),
              boost::bind(std::greater<int>(),
                          boost::bind(&Variable::degree,_1,critical),
                          boost::bind(&Variable::degree,_2,critical)));

    /// weight heuristic to order the variables ================================= END

    VarPtr current = future.back();

    while(!future.empty()){

        //cout << "variable id: " << current->id_ << endl;

        if(timer.elapsed() > 60.0){
            depass = true;
            return true;

            cout << "time out by mac ";
            exit(0);
        }



        for(std::vector<ValPtr>::iterator ir = current->values_.begin();
            ir != current->values_.end(); ++ir){

            if((*ir)->getGammaS() > 0)
                continue;

            current->assign(*ir);
            // for fc and arc consistency reason
            future.pop_back();
            past.push_back(current);

            if(future.empty()){
                std::for_each(first, last, boost::bind(&Variable::resetDomain,_1));
                return true;
            }

            if(!fc_consistency(past.begin(), past.end(),varPred, ctrPred) ||
                    !arc_consistency(future.begin(), future.end(), varPred, ctrPred)){

                (*ir)->setGammaS(100);
                current->resign();

                // if current value is not consistent, then recover all
                // arc effects - a little bit overload, but simple to implement
                // since we only need to recover the last arc effect, but
                // the last arc effect is difficult to be verified
                std::for_each(future.begin(), future.end(),
                              boost::bind(&Variable::resetDomain, _1));

                future.push_back(current);
                past.pop_back();
            }else{
                break;
            }
        }

        if(!current->isAssigned()){

            /// the problem is not feasible
            if(past.empty()){
                std::for_each(filter_vfirst, filter_vlast, boost::bind(&Variable::resetDomain,_1));
                return false;
            }


            std::for_each(future.begin(), future.end(),
                          boost::bind(&Variable::resetDomain, _1));



            current = past.back();
            past.pop_back();
            current->current()->setGammaS(100);
            current->resign();
            future.push_back(current);

        }else{
            current = future.back();
        }

    }


    std::for_each(filter_vfirst, filter_vlast, boost::bind(&Variable::resetDomain,_1));
    return true;
}


/// ===============================================================================================
/// Extract routine
/// ===============================================================================================

void markedToCriticalCtr(CtrPtr v){
    if(v->getLevel() == normal){
        v->setLevel(disable);
    }


    if(v->getLevel() == marked){
        v->setLevel(critical);
        v->setVarsLevel(critical);
        v->resign();
    }

    v->setGammaS(1);
}

void criticalToMarked(CtrPtr c){
    if(c->getLevel() == critical){
        c->setLevel(marked);
	c->setVarsLevel(marked);
	c->setGammaS(100);
    }else{
        c->setGammaS(1);
    }

    
}

void normalToDisableVar(VarPtr v){
    if(v->getLevel() == normal)
        v->setLevel(disable);
}


void NormalToCriticalCtr(CtrPtr v){
    if(v->getLevel() == normal){
        v->setLevel(critical);
        v->setVarsLevel(critical);
        v->resign();
    }

    v->setGammaS(1);
}

void criticalToNormalCtr(CtrPtr v){
    if(v->getLevel() > normal){
        v->setGammaS(1);
        v->setLevel(normal);
        v->setVarsLevel(normal);
    }
}


void normalToDisableCtr(CtrPtr v){
    if(v->getLevel() == normal){
        v->setGammaS(1);
        v->setLevel(disable);
        v->setVarsLevel(disable);
    }
}



void disableCtr(CtrPtr v){
    v->setGammaS(1);
    v->setLevel(disable);
    v->setVarsLevel(disable);
}





void resetVar(VarPtr v){

    v->resetAtts();
    v->setGammaS(1);
    v->resetDomain();
}

void resetCtr(CtrPtr v){
    v->resetAtts();
    v->setGammaS(1);
}



void clearVar(VarPtr v){

    v->clearAtts();
    v->setGammaS(1);
    v->resetDomain();

    if(v->getLevel() > normal)
        v->setLevel(normal);
}

void clearCtr(CtrPtr c){

    c->clearAtts();
    c->setGammaS(1);

    if(c->getLevel() > normal)
        c->setLevel(normal);
}

inline bool emptyDomain(VarPtr const v){
    BOOST_FOREACH(ValPtr const a, v->values_){
        if(a->getLevel() > disable && a->getGammaS() < 1)
            return false;
    }

    return true;
}

template <typename VarItr, typename CtrItr >
bool extract(VarItr vfirst, VarItr vlast, CtrItr cfirst, CtrItr clast){

    Statistic stat;
    stat.init();

    // filtering procedure
    // step 0: arc-consistency ===== BGN
    //cout << "arc-consitency" << endl;

    CtrPtr c; // first located constraint

    bool isLocated = false;


        stat.restart();
        bool feasibleArc = arc_consistency(vfirst,vlast,Normal<VarPtr>(), Normal<CtrPtr>());
        stat.setTime(CateFilter);


    // disable gamma value



    if(!feasibleArc){
        isLocated = true;
        //cout << "arc-consitency infeasible" << endl;
        // step 1.1: locate first constraint in case of infeasbility of arc-consistency
        // find the critical constraint by the means of deadend
        VarItr vItr = std::find_if(vfirst, vlast, emptyDomain);

        VarPtr v = *vItr;

        VarPtr opposite = *std::min_element(v->vars_.begin(), v->vars_.end(),
                                            boost::bind(std::less<int>(),
                                                        boost::bind(validDomainSize,_1),
                                                        boost::bind(validDomainSize,_2)));



        // set normal variables
        BOOST_FOREACH(CtrPtr const cx, v->ctrs_){
            if(cx->isIn(opposite)){
                c = cx;
                break;
            }
        }

        //
        std::for_each(vfirst,vlast,boost::bind(&Variable::setLevel,_1,disable));
        std::for_each(vfirst,vlast,boost::bind(&Variable::resetDomain,_1));
        std::for_each(cfirst,clast,boost::bind(disableCtr,_1));

        c->setLevel(normal);
        c->setVarsLevel(normal);


        if(!c->findSolution()){
            c->setLevel(critical);
            c->setVarsLevel(critical);
            stat.print();
            return false;
        }
    }else{
        //cout << "arc-consitency feasible" << endl;

        /// strategy 1: ===== Effect the arc-consistency
        for(VarItr current = vfirst; current != vlast; ++current){
            std::for_each((*current)->values_.begin(),(*current)->values_.end(),
                          boost::bind(removeValue, _1));
        }

        /// strategy 2: ===== Not Effect the arc-consistency
        //        for(VarItr current = vfirst; current != vlast; ++current){
        //            std::for_each((*current)->values_.begin(),(*current)->values_.end(),
        //                          boost::bind(&Value::resetAtts, _1));
        //        }
        // give a random solution

    }
    // step 0 ====================== END




    // correctness routine
    bool feasible = true;
    //boost::progress_timer t;
    MinGamma<false> minGamma;

    while(feasible){


        if(stat.getTotal() > 2000.0){
            stat.print();
            cout << "EXTTIMEOUT" ;
            return false;
        }


        // step 1.2 : locate first constraint in case of feasbility of arc-consistency
        if(!isLocated){
            //isLocated = true;

            //cout << "Never located !!!!!!!!!!!!!" << endl;
            std::for_each(vfirst,vlast,boost::bind(&Variable::setLevel,_1,normal));
            std::for_each(cfirst,clast,boost::bind(&AbsConstraint::setLevel,_1,normal));

            randomSolver(vfirst, vlast);
            //cout << "breakout locating" << endl;

            // breakout to find the critical constraint
            stat.restart();
            feasible = breakout(vfirst,vlast,cfirst,clast,minGamma,
                                Normal<VarPtr>(),Normal<CtrPtr>());
            stat.setTime(CateLocate);

            // only locate the normal violated constraint;

            if(feasible){
                cout << "breakout feasible" << endl;
                stat.print();
                return true;
            }

            bool isSet = false;
            for(CtrItr current = cfirst; current != clast; ++current){
                // skip non-effective and satisfied constraints
                if((*current)->getLevel() < normal || (*current)->satisfy())
                    continue;

                if(!isSet || c->getGammaS() < (*current)->getGammaS()){
                    c = *current;
                    if(!isSet)
                        isSet = true;
                }

            }


            assert(isSet);

            //cout << "constrat gamma: " << c->getGammaS() << endl;

            // recover all diable values and reset all
            std::for_each(cfirst,clast,boost::bind(disableCtr,_1));
            std::for_each(vfirst,vlast,boost::bind(&Variable::resetDomain,_1));

            c->setLevel(normal);
            c->setVarsLevel(normal);

            if(!c->findSolution()){
                c->setLevel(critical);
                c->setVarsLevel(critical);
                stat.print();
                return false;
            }
        }

        feasible = true;
        // step 2: construct the core arround the constraint

        //bool breakoutViolatedFeasible = false;
        while(feasible){
            //cout << "constructing" << endl;
            stat.restart();
	    

            feasible = constructMRV(vfirst,vlast, cfirst,clast);

            //// print graph
            //            std::string ext = "extension";
            //            printGraph(cfirst,clast,Normal<CtrPtr>(),ext);

            stat.setTime(CateConstruct);

	    if(feasible){
                //cout << "Go To locate by MRV" << endl;

                isLocated = false;
                break;
	    }	    

            //cout << "constructing end" << endl;
            if(!feasible){
                //cout << "constructing not feasible" << endl;

	        for(CtrItr current = cfirst; current != clast; ++current){
                    if ((*current)->getLevel() < normal && (*current)->varsLevelgreaterThan(disable)){
                        (*current)->setLevel(normal);
                    }
		}

                //for(CtrItr current = cfirst; current != clast; ++current){
                //    if((*current)->sameVarsLevel()){
                //        if((*current)->getVarsLevel() > disable){
                //            if((*current)->getLevel()< normal)
                //                (*current)->setLevel(normal);
                //        }
                //    }
                //}

                std::for_each(vfirst,vlast,boost::bind(&Variable::resetDomain,_1));

                // mark routine: one by one
                //cout << "breakout marking" << endl;

		// convert critical to marked
		std::for_each(cfirst,clast,boost::bind(criticalToMarked,_1));

                stat.restart();
                feasible = breakoutViolated(vfirst,vlast,cfirst,clast,minGamma,
                                            Normal<VarPtr>(),Normal<CtrPtr>());
                stat.setTime(CateConstruct);

		
		
                //cout << "after breakout marking" << endl;

            }


	    if(feasible){
                //cout << "Go To locate by BreakoutViolated" << endl;
                std::for_each(vfirst,vlast,boost::bind(&Variable::resetDomain,_1));
                std::for_each(cfirst,clast,boost::bind(&AbsConstraint::setGammaS,_1,1));
                continue;

                // verify all constraint
                //                cout << "constructing feasible" << endl;
                //         bool hasViolated = false;
                //                for(CtrItr current = cfirst; current != clast;++current){
                //    if((*current)->getLevel() == disable && (*current)->isViolated() ){
                //      hasViolated = true;
                //      break;
                //    }
                //                }

                //                if(!hasViolated){
                //                    stat.print();
                //                    return true; // extract find feasible solution
                //                }else
                //  continue;
            }
        }


        //// print graph
        //        std::string ext = "centralization";
        //        printGraph(cfirst,clast,Marked<CtrPtr>(),ext);


        //feasible = mac_solver(vfirst,vlast, Critical<VarPtr>(), Critical<CtrPtr>(),t);
        if(!feasible){
            // gamma system reset before entering the mac_solver
            //cout << "in Mac" << endl;

            // report the reduce level
            //            int normalNb = std::accumulate(vfirst, vlast, 0,
            //                                           Counter_Normal<VarPtr>() );

            //            int markNb = std::accumulate(vfirst, vlast, 0,
            //                                         Counter_Marked<VarPtr>() );

            //            cout << "\n reduce: " << double(normalNb - markNb)/double(normalNb) << endl;

            /// in case the breakoutViolated is infeasible,
            /// only treat the marked elements, and disable all normal elements
            std::for_each(cfirst,clast,boost::bind(markedToCriticalCtr,_1));
            std::for_each(vfirst,vlast,boost::bind(normalToDisableVar,_1));


            //// saturation
            //            for(CtrItr current = cfirst; current != clast; ++current){
            //                if((*current)->sameVarsLevel()){
            //                    if((*current)->getVarsLevel() > normal){
            //                        if((*current)->getLevel()< critical)
            //                            (*current)->setLevel(critical);
            //                    }
            //                }
            //            }



            /// only for coloring problem, check clique ================== BGN
//            int sizeVar = std::accumulate(vfirst,vlast,0,Counter_Critical<VarPtr>());
//            int sizeCtr = std::accumulate(cfirst,clast,0,Counter_Critical<CtrPtr>());

//            //            cout << "sizeVars: " << sizeVar << endl;
//            //            cout << "values  : " << (*vfirst)->values_.size() << endl;



//            if(sizeCtr == (sizeVar-1)*sizeVar /2){

//                if(sizeVar <= (*vfirst)->values_.size()){

//                    feasible =  true;
//                    //isLocated = false;

//                    std::for_each(vfirst,vlast,boost::bind(resetVar,_1));
//                    std::for_each(cfirst,clast,boost::bind(resetCtr,_1));

//                }else{
//                    feasible = false;
//                    break;
//                }
//            }
            /// only for coloring problem, check clique ================== END


            bool depass = false;
            stat.restart();
            if(!feasible){

                return false;

                feasible = mac_solver(vfirst,vlast,Critical<VarPtr>(), Critical<CtrPtr>(),depass);

            }
            stat.setTime(CateExact);

            if(depass){
                stat.print();
                cout << "MACTIMEOUT" ;
                return false;
            }

            if(!feasible){
                stat.print();
                return false;
            }else{
                // for feasible by mac, the harder core is located
                isLocated = true;
            }
        }

        if(isLocated){
            std::for_each(vfirst,vlast,boost::bind(clearVar,_1));
            std::for_each(cfirst,clast,boost::bind(clearCtr,_1));
	}
	
	//else{
        //    std::for_each(vfirst,vlast,boost::bind(resetVar,_1));
        //    std::for_each(cfirst,clast,boost::bind(resetCtr,_1));
        //}

    }

    // setp 3: recover

    stat.print();
    return false;

}


template <typename VarItr, typename CtrItr >
bool extractVars(VarItr vfirst, VarItr vlast, CtrItr cfirst, CtrItr clast){

    Statistic stat;
    stat.init();


    randomSolver(vfirst, vlast);

    // correctness routine
    bool feasible = true;
    //boost::progress_timer t;
    MinGamma<false> minGamma;

    bool isLocated = false;
    while(feasible){


        if(stat.getTotal() > 600.0){
            stat.print();
            cout << "EXTTIMEOUT" ;
            return false;
        }


        stat.restart();
        feasible = breakoutViolated(vfirst,vlast,cfirst,clast,minGamma,
                                    Normal<VarPtr>(),Normal<CtrPtr>());
        stat.setTime(CateConstruct);


        if(feasible){
            bool feas = true;
            for(CtrItr current = cfirst; current != clast;  ++current){
                if((*current)->isViolated()){
                    feas = false;
                    break;
                }
            }


            if(feas)
                return true;
        }


        while(!feasible){





            bool depass = false;
            stat.restart();
            feasible = mac_solver(vfirst,vlast,Marked<VarPtr>(), Marked<CtrPtr>(),depass);
            stat.setTime(CateExact);


            if(depass){
                stat.print();
                cout << "MACTIMEOUT" ;
                return false;
            }

            if(feasible){

                bool feas = true;
                for(CtrItr current = cfirst; current != clast;  ++current){
                    if((*current)->isViolated()){
                        feas = false;
                        break;
                    }
                }


                if(feas)
                    return true;


                std::for_each(vfirst,vlast,boost::bind(&Variable::resetDomain,_1));
            }else{
                stat.print();
                return false;
            }


        }




    }

    // setp 3: recover

    stat.print();
    return false;

}





template <typename VarItr, typename CtrItr >
bool extractVarsColoring(VarItr vfirst, VarItr vlast, CtrItr cfirst, CtrItr clast){

    Statistic stat;
    stat.init();


    randomSolver(vfirst, vlast);

    // correctness routine
    bool feasible = true;
    //boost::progress_timer t;
    MinGamma<false> minGamma;

    bool isLocated = false;
    while(feasible){


        if(stat.getTotal() > 600.0){
            stat.print();
            cout << "EXTTIMEOUT" ;
            return false;
        }


        stat.restart();
        feasible = breakoutViolated(vfirst,vlast,cfirst,clast,minGamma,
                                    Normal<VarPtr>(),Normal<CtrPtr>());
        stat.setTime(CateConstruct);


        if(feasible){
            bool feas = true;
            for(CtrItr current = cfirst; current != clast;  ++current){
                if((*current)->isViolated()){
                    feas = false;
                    break;
                }
            }


            if(feas)
                return true;
        }


        //        while(!feasible){





        //            bool depass = false;
        //            stat.restart();
        //            feasible = mac_solver(vfirst,vlast,Marked<VarPtr>(), Marked<CtrPtr>(),depass);
        //            stat.setTime(CateExact);


        //            if(depass){
        //                stat.print();
        //                cout << "MACTIMEOUT" ;
        //                return false;
        //            }

        //            if(feasible){

        //                bool feas = true;
        //                for(CtrItr current = cfirst; current != clast;  ++current){
        //                    if((*current)->isViolated()){
        //                        feas = false;
        //                        break;
        //                    }
        //                }


        //                if(feas)
        //                    return true;


        //                std::for_each(vfirst,vlast,boost::bind(&Variable::resetDomain,_1));
        //            }else{
        //                stat.print();
        //                return false;
        //            }


        //        }




    }

    // setp 3: recover

    stat.print();
    return false;

}




template <typename VarItr, typename CtrItr >
bool extractConstraint(VarItr vfirst, VarItr vlast, CtrItr cfirst, CtrItr clast){

    Statistic stat;
    stat.init();


    randomSolver(vfirst, vlast);

    // correctness routine
    bool feasible = true;
    //boost::progress_timer t;
    MinGamma<false> minGamma;

    bool isLocated = false;
    while(feasible){


        if(stat.getTotal() > 600.0){
            stat.print();
            cout << "EXTTIMEOUT" ;
            return false;
        }


        stat.restart();
        feasible = breakoutViolatedConstraint(vfirst,vlast,cfirst,clast,minGamma,
                                              Normal<VarPtr>(),Normal<CtrPtr>());
        stat.setTime(CateConstruct);


        while(!feasible){





            bool depass = false;
            stat.restart();
            feasible = mac_solver(vfirst,vlast,Marked<VarPtr>(), Marked<CtrPtr>(),depass);
            stat.setTime(CateExact);


            if(depass){
                stat.print();
                cout << "MACTIMEOUT" ;
                return false;
            }

            if(feasible){



                for(CtrItr current=cfirst;current!=clast;++current){
                    if((*current)->getLevel() !=marked && (*current)->isViolated()){
                        //cout << "has violated" << endl;
                        (*current)->setLevel(marked);
                        (*current)->setVarsLevel(marked);
                        (*current)->setGammaS(100);
                        feasible = false;
                        break;
                    }

                }

                std::for_each(vfirst,vlast,boost::bind(&Variable::resetDomain,_1));
            }else{
                stat.print();
                return false;
            }


        }




    }

    // setp 3: recover

    stat.print();
    return false;

}


#endif // HEURISTIC_HPP
