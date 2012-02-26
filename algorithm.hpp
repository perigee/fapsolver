#ifndef ALGORITHM_HPP
#define ALGORITHM_HPP


#include <functional>
#include <iostream> 
#include <vector> 
#include <algorithm> 
#include <numeric>
#include <ctime>            // std::time

// boost library
#include <boost/function_output_iterator.hpp>
#include <boost/iterator/filter_iterator.hpp>
#include <boost/bind.hpp>
#include <boost/lambda/core.hpp>
#include <boost/ref.hpp>
#include <boost/function.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

// custom library
#include "colAlgo.hpp"






/// =================================================================================================
/// Exact Method MAC
/// =================================================================================================




/// import, why weighted as neigborhood, since concerned connection
inline bool sort_descent_weight(VarPtr const a, VarPtr const b){
    int aw = std::accumulate(a->vars_.begin(), a->vars_.end(),
                             0,Counter_markedWeight<VarPtr>());

    int bw = std::accumulate(b->vars_.begin(), b->vars_.end(),
                             0,Counter_markedWeight<VarPtr>());


    return (a->getWeight() + aw) > (b->getWeight() + bw) ;
}

inline bool zero_gamma_value(ValPtr const a){
    return a->getGamma() == 0;
}




// @return true: no deadend, false: has deadend
template <typename PredCtr>
bool revise(VarPtr var, PredCtr ctrPred){

    bool hasDeleted = false;

    // in revise procedure
    std::vector<CtrPtr>::iterator first = var->ctrs_.begin();
    std::vector<CtrPtr>::iterator last = var->ctrs_.end();

    // verify all valid constraints, see if there is value erased

    for(;first != last; ++first){

        if((*first)->isMarked()){

            if(!(*first)->getOther(var)->isAssigned()){
                bool hasD = false;
                if(!(*first)->ac(var,hasD)){
                    return false;
                }else{
                    if(hasD){
                        if(!hasDeleted)
                            hasDeleted = true;
                    }
                }

            }

        }
    }

    // if there is pruning, then revise
    if(hasDeleted){
        BOOST_FOREACH(VarPtr const v, var->vars_){

            if(v->isMarked()){
                if(!v->isAssigned()){
                    if(!revise(v))
                        return false;
                }
            }
        }
    }

    return true;
}








// entering non assigned variables, dealing with marked variable only
// @return true: arc-consistent; false: deadend
template <typename Iter>
        bool arc_consistency(Iter first, Iter last){

    //cout << "in ac " << endl;

    // only deal with the gamma free value and test only gamma constraints
    for(Iter current = first; current != last; ++current){

        if((*current)->isMarked()){

            if((*current)->isAssigned()){
                cout << "Error: impossible assigned var " << endl;
                exit(0);
            }

            if(!revise(*current))
                return false;
        }
    }


    return true;

}

inline void resetVarBeforeMAC(VarPtr var){
    var->resetDom();
    var->resign();
}



template <typename Iter>
        bool mac_solver(Iter first, Iter last, boost::progress_timer& timer){

    //cout << "in mac_solver"<< endl;



    int size = std::distance(first, last);

    if( size == 0)
        return true;

    std::for_each(first, last, resetVarBeforeMAC);

    std::vector<VarPtr> future;
    future.reserve(size);
    std::vector<VarPtr> past;
    past.reserve(size);


    std::copy(first, last, std::back_inserter(future));
    std::sort(future.begin(), future.end(), sort_descent_weight);
    std::for_each(future.begin(), future.end(), boost::bind(&Variable::clearWeight,_1));

    VarPtr current = future.back();

    while(!future.empty()){

        //cout << "variable id: " << current->id_ << endl;

        if(timer.elapsed() > 800.0){
            cout << "time out by mac\n" << endl;

            exit(0);
        }



        for(std::vector<ValPtr>::iterator ir = current->values_.begin();
        ir != current->values_.end(); ++ir){

            if((*ir)->getGamma() > 0)
                continue;

            current->assign(*ir);
            // for fc and arc consistency reason
            future.pop_back();
            past.push_back(current);

            if(future.empty()){
                std::for_each(first, last, boost::bind(&Variable::resetDom,_1));

                return true;
            }

            if(!fc_consistency(past.begin(), past.end()) ||
               !arc_consistency(future.begin(), future.end())){

                (*ir)->setGamma(100);
                current->resign();

                // if current value is not consistent, then recover all
                // arc effects - a little bit overload, but simple to implement
                // since we only need to recover the last arc effect, but
                // the last arc effect is difficult to be verified
                std::for_each(future.begin(), future.end(),
                              boost::bind(&Variable::resetDom, _1));

                future.push_back(current);
                past.pop_back();
            }else{
                break;
            }
        }

        if(!current->isAssigned()){

            /// the problem is not feasible
            if(past.empty()){
                std::for_each(first, last, boost::bind(&Variable::resetDom,_1));
                return false;
            }


            std::for_each(future.begin(), future.end(),
                          boost::bind(&Variable::resetDom, _1));



            current = past.back();
            past.pop_back();            
            current->current()->setGamma(100);
            current->resign();
            future.push_back(current);

        }else{
            current = future.back();
        }

    }


    std::for_each(first, last, boost::bind(&Variable::resetDom,_1));
    return true;
}



typedef boost::minstd_rand base_generator_type;

/// Filter for violated gammaed variable only
struct is_gamma_variable{
    bool operator()(VarPtr const v){return v->getGamma() > 0 ;}
};

struct is_marked_variable{
    bool operator()(VarPtr const v){return v->isMarked();}
};

template <typename Element>
        struct not_marked_element{
    bool operator ()(Element const e){
        return !e->isMarked();}
};



/// Find best value with minimal gamma
/// used in tabu search
inline bool  min_gamma_value(ValPtr const a, ValPtr const b){

    // in case both values in tabu is false
    if(!a->isTabu() || !b->isTabu()){
        if(a->isTabu())
            return false;

        if(b->isTabu())
            return true;

        if(a->getGamma() == b->getGamma())
            return rand()/static_cast<double>(RAND_MAX) < 0.5;

        return a->getGamma() < b->getGamma();
    }

    // in case both values in tabu
    return rand()/static_cast<double>(RAND_MAX) < 0.5;

}


inline bool  min_delta_value(VarPtr const va, ValPtr const a, 
			     VarPtr const vb, ValPtr const b){

    // in case both values in tabu is false
    if(!a->isTabu() || !b->isTabu()){
        if(a->isTabu())
            return false;

        if(b->isTabu())
            return true;

        if((a->getGamma()-va->current()->getGamma())
            == (b->getGamma() - vb->current()->getGamma()))
            return rand()/static_cast<double>(RAND_MAX) < 0.5;

        return (a->getGamma()-va->current()->getGamma())
                < (b->getGamma() - vb->current()->getGamma());
    }

    // in case both values in tabu
    return rand()/static_cast<double>(RAND_MAX) < 0.5;

}

/// get the best move among the neighbor solutions
struct Get_best_move{


    Get_best_move():isSet_(false),count_(0){}

    void evaluate(VarPtr var){
        //        cout << "in evaluate var: " << var->id_
        //                << " domain size: " << var->values_.size()
        //                << " found conflict: " << count_ << endl;


        if (!var->isViolated())
            return;


        ++count_;


        /// acumulate the weight on violated variables ===== BGN
        var->addWeight();
        /// acumulate the weight on violated variables ===== END

        std::vector<ValPtr>::iterator ir
                = std::min_element(var->values_.begin(),
                                   var->values_.end(),
                                   min_gamma_value);

        // verify if the chosen one is current one
        // if so, get the neighbor value

        //                if(ir == var->values_.end())
        //                    cout << "ir not found " << endl;

        //cout << "before ir " << endl;
        if(*ir == var->current()){
            //cout << "found end" << endl;
            if(++ir == var->values_.end()){
                --ir; // return to current iterator
                --ir; // get previous iterator
            }
        }

        std::for_each(var->values_.begin(), var->values_.end(),
                      boost::bind(&Value::reduceTabu,_1));

        // initialization if needed
        if(!isSet_){
            //cout << "in no set" << endl;
            isSet_ = true;
            bestVal_ = *ir;
            bestVar_ = var;
        }else{
            // comparing procedure
            //cout << "in set" << endl;
            if(min_delta_value(var, *ir, bestVar_, bestVal_)){
                bestVar_ = var;
                bestVal_ = *ir;
            }

        }
        //cout << "variable: " << bestVar_->id_ << " = " << bestVal_->v_ << endl;

    }

    // get the best move pair
    std::pair<VarPtr,ValPtr> getBestMove(){
        //        cout << "in getBestMove()" << endl;
        //        cout << "variable: " << bestVar_->id_ << " = " << bestVal_->v_ << endl;
        //        cout << "after getBestMove()" << endl;
        std::pair<VarPtr,ValPtr> p;
        p.first = bestVar_;
        p.second = bestVal_;
        return p;
    }


    int getConflictNb(){return count_;}

    void reset(){
        isSet_ = false;
        count_ = 0;
    }

private: 
    bool isSet_;
    int count_;
    VarPtr bestVar_; // best variable found
    ValPtr bestVal_; // best value found in best variable
};


/// get the best move among the neighbor solutions
struct Get_best_move_WithoutWeights{


    Get_best_move_WithoutWeights():isSet_(false),count_(0){}

    void evaluate(VarPtr var){
        //        cout << "in evaluate var: " << var->id_
        //                << " domain size: " << var->values_.size()
        //                << " found conflict: " << count_ << endl;


        if (!var->isViolated())
            return;


        ++count_;


        /// acumulate the weight on violated variables ===== BGN
        //var->addWeight();
        /// acumulate the weight on violated variables ===== END

        std::vector<ValPtr>::iterator ir
                = std::min_element(var->values_.begin(),
                                   var->values_.end(),
                                   min_gamma_value);

        // verify if the chosen one is current one
        // if so, get the neighbor value

        //                if(ir == var->values_.end())
        //                    cout << "ir not found " << endl;

        //cout << "before ir " << endl;
        if(*ir == var->current()){
            //cout << "found end" << endl;
            if(++ir == var->values_.end()){
                --ir; // return to current iterator
                --ir; // get previous iterator
            }
        }

        std::for_each(var->values_.begin(), var->values_.end(),
                      boost::bind(&Value::reduceTabu,_1));

        // initialization if needed
        if(!isSet_){
            //cout << "in no set" << endl;
            isSet_ = true;
            bestVal_ = *ir;
            bestVar_ = var;
        }else{
            // comparing procedure
            //cout << "in set" << endl;
            if(min_delta_value(var, *ir, bestVar_, bestVal_)){
                bestVar_ = var;
                bestVal_ = *ir;
            }

        }
        //cout << "variable: " << bestVar_->id_ << " = " << bestVal_->v_ << endl;

    }

    // get the best move pair
    std::pair<VarPtr,ValPtr> getBestMove(){
        //        cout << "in getBestMove()" << endl;
        //        cout << "variable: " << bestVar_->id_ << " = " << bestVal_->v_ << endl;
        //        cout << "after getBestMove()" << endl;
        std::pair<VarPtr,ValPtr> p;
        p.first = bestVar_;
        p.second = bestVal_;
        return p;
    }


    int getConflictNb(){return count_;}

    void reset(){
        isSet_ = false;
        count_ = 0;
    }

private:
    bool isSet_;
    int count_;
    VarPtr bestVar_; // best variable found
    ValPtr bestVal_; // best value found in best variable
};




/// function set the solution
template <typename Iter>
        void setSolution(Iter first, Iter last){
    for(;first!=last;++first)
        (*first)->setSolution((*first)->current());

}

template <typename Iter>
        void getSolution(Iter first, Iter last){
    for(;first!=last;++first){
	if((*first)->hasSolution()){
            (*first)->assign((*first)->getSolution());
            (*first)->clearSolution();
        }
        (*first)->resetDom();
    }

}



template <typename Iter>
        bool isSolutionFeasible(Iter first, Iter last){
    for(;first!=last;++first){
        if((*first)->isViolated())
            return false;
    }

    return true;

}


// initalize the gamma for each value in domain
struct InitGamma{
    void init(VarPtr v){
        //cout << "initValuesGamma" << endl;
        //VarPtr v = VarPtr(this);
        BOOST_FOREACH(ValPtr const a, v->values_){
            int ct = 0;
            //cout << "ctrs size: " << ctrs_.size() << endl;
            BOOST_FOREACH(CtrPtr const c, v->ctrs_){
                //cout << "before initValuesGamma " << ++ct << endl;
                if(c->getGamma() >0){
                    if(c->isViolated(v,a)){
                        //cout << "inconsistent" << endl;
                        a->setGamma(a->getGamma() + c->getGamma());
                    }
                }
            }
        }

        //cout << "initValuesGamma finished " << endl;
    }
};

/// tabu search
template <typename Itr, typename Pred>
        bool tabu(Itr first, Itr last, Pred pred){


    //return true;
    //cout << "in tabu solver" << endl;

    base_generator_type generator(42u);
    generator.seed(static_cast<unsigned int>(std::time(0)));
    boost::uniform_int<> degen_dist(0,9);
    boost::variate_generator<base_generator_type&, boost::uniform_int<> >
            deg(generator, degen_dist);

    boost::uniform_int<> degen_distx(70,90);
    boost::variate_generator<base_generator_type&, boost::uniform_int<> >
            degCELAR(generator, degen_distx);


    register int currentCost = 0;
    register int bestCost = 0;
    bool foundFirst = false;




    /// Initialize the variables filter
    typedef boost::filter_iterator<Pred, Itr> FilterIter;
    FilterIter filter_first(pred, first, last);
    FilterIter filter_last(pred, last,last);

    //register int MaxIter = std::distance(first,last);
    register int MaxIter = 1000;
    register int currentIter = MaxIter;
    register int Freq = 100;
    register int freq = 0;
    int duration = 0;


    // intialize values' gamma
    InitGamma initG;
    std::for_each(filter_first,filter_last,
                  boost::bind(&InitGamma::init,
                              &initG,_1));



    Get_best_move bestMove; //functor consists of finding best move 
    
    
    while(true){

        //cout << "before bestMove" << endl;
        /// Find best move
     	std::for_each(filter_first, filter_last,
		      boost::bind(&Get_best_move::evaluate, 
				  &bestMove,_1));
        //cout << "after bestMove" << endl;



        //        cout << "after getBest" << endl;
        //        cout << "before getConflict" << endl;
        int conflictNb = bestMove.getConflictNb();


        //cout << "best Delta: " << bestDelta << endl;
        //        cout << "conflictNb: " << conflictNb << endl;
        if(conflictNb < 1){
            std::for_each(filter_first, filter_last, boost::bind(&Variable::resetDom,_1));
            std::for_each(filter_first, filter_last, boost::bind(&Variable::clearSolution,_1));
            bestMove.reset();
            //            cout << "found feasible " << endl;
            return true;
        }


        /// Update tabu tenure

        std::pair<VarPtr, ValPtr> pMove = bestMove.getBestMove();
        if((pMove.second)->isTabu())
            (pMove.second)->popTabu();

        int bestDelta = (pMove.second)->getGamma()
                        - (pMove.first)->current()->getGamma();

        bool findBetter = false;
	if(!foundFirst){
            // try to find first move which accumulates the cost
            if(bestDelta >= 0){
                foundFirst = true;
                // update best move
                setSolution(filter_first, filter_last);
                currentCost += bestDelta;
            }
	}else{
            --currentIter;
            currentCost += bestDelta;
            if(currentCost < bestCost){
                // udpate best move
                setSolution(filter_first, filter_last);
                bestCost = currentCost;
                currentIter = MaxIter;
                findBetter = true;
            }

	}



        /// ====================================================
        /// ================================================ BGN

        /// stratgey 1:
	// compute the tabu tenure and current objective funtion
        // duration = deg() + boost:: numeric_cast<int>(0.6*conflictNb);

        /// stratgey 2:
        //duration = degCELAR();


        /// stratgey 3:
        int domSize = (pMove.first)->values_.size();

        if(!foundFirst){
            duration = domSize;
        }else{
            ++freq;

            if(findBetter)
                duration = domSize;

            if(freq > Freq){
                if((duration/domSize) < conflictNb ){
                    duration += domSize;
                    freq = 0;
                }
            }
        }

        //        cout << currentCost
        //                << "\t" << bestCost
        //                << "\t" << conflictNb
        //                << "\t" << duration
        //                << endl;

        /// ================================================ END
        /// ====================================================





        (pMove.first)->current()->pushTabu(duration);  
	
	/// Update best move
        BOOST_FOREACH(CtrPtr const ctr, (pMove.first)->ctrs_){
            if(ctr->getGamma() > 0){
                ctr->updateGamma(pMove.first,(pMove.first)->current(),pMove.second);
            }
        }


	(pMove.first)->assign(pMove.second);

        /// add weights only on one movement variable ================= BGN
        //        if((pMove.second)->getGamma() > 0)
        //            (pMove.first)->addWeight();
        /// add weights only on one movement variable ================= END

        //        if((pMove.first)->varWeight_ > 0)
        //            cout << "w: current: " << (pMove.first)->varWeight_ << endl;

	bestMove.reset();

        /// Stop criteria
        if(currentIter < 0){
            // set the best solution found
            getSolution(filter_first, filter_last);

            return isSolutionFeasible(filter_first, filter_last);
	}



    }

}


/// tabu search with short term of non-improvement
/// also without increasing the weights on variables
template <typename Itr, typename Pred>
        bool tabuShort(Itr first, Itr last, Pred pred){


    //return true;
    //cout << "in tabu solver" << endl;

    base_generator_type generator(42u);
    generator.seed(static_cast<unsigned int>(std::time(0)));
    boost::uniform_int<> degen_dist(0,9);
    boost::variate_generator<base_generator_type&, boost::uniform_int<> >
            deg(generator, degen_dist);

    boost::uniform_int<> degen_distx(70,90);
    boost::variate_generator<base_generator_type&, boost::uniform_int<> >
            degCELAR(generator, degen_distx);


    register int currentCost = 0;
    register int bestCost = 0;
    bool foundFirst = false;




    /// Initialize the variables filter
    typedef boost::filter_iterator<Pred, Itr> FilterIter;
    FilterIter filter_first(pred, first, last);
    FilterIter filter_last(pred, last,last);

    //register int MaxIter = std::distance(first,last);
    register int MaxIter = 100;
    register int currentIter = MaxIter;
    register int Freq = 30;
    register int freq = 0;
    int duration = 0;


    // intialize values' gamma
    InitGamma initG;
    std::for_each(filter_first,filter_last,
                  boost::bind(&InitGamma::init,
                              &initG,_1));



    Get_best_move_WithoutWeights bestMove; //functor consists of finding best move without incrementing the weights on variables


    while(true){

        //cout << "before bestMove" << endl;
        /// Find best move
        std::for_each(filter_first, filter_last,
                      boost::bind(&Get_best_move_WithoutWeights::evaluate,
                                  &bestMove,_1));
        //cout << "after bestMove" << endl;



        //        cout << "after getBest" << endl;
        //        cout << "before getConflict" << endl;
        int conflictNb = bestMove.getConflictNb();


        //cout << "best Delta: " << bestDelta << endl;
        //        cout << "conflictNb: " << conflictNb << endl;
        if(conflictNb < 1){
            std::for_each(filter_first, filter_last, boost::bind(&Variable::resetDom,_1));
            std::for_each(filter_first, filter_last, boost::bind(&Variable::clearSolution,_1));
            bestMove.reset();
            //            cout << "found feasible " << endl;
            return true;
        }


        /// Update tabu tenure

        std::pair<VarPtr, ValPtr> pMove = bestMove.getBestMove();
        if((pMove.second)->isTabu())
            (pMove.second)->popTabu();

        int bestDelta = (pMove.second)->getGamma()
                        - (pMove.first)->current()->getGamma();

        bool findBetter = false;
        if(!foundFirst){
            // try to find first move which accumulates the cost
            if(bestDelta >= 0){
                foundFirst = true;
                // update best move
                setSolution(filter_first, filter_last);
                currentCost += bestDelta;
            }
        }else{
            --currentIter;
            currentCost += bestDelta;
            if(currentCost < bestCost){
                // udpate best move
                setSolution(filter_first, filter_last);
                bestCost = currentCost;
                currentIter = MaxIter;
                findBetter = true;
            }

        }



        /// ====================================================
        /// ================================================ BGN

        /// stratgey 1:
        // compute the tabu tenure and current objective funtion
        // duration = deg() + boost:: numeric_cast<int>(0.6*conflictNb);

        /// stratgey 2:
        //duration = degCELAR();


        /// stratgey 3:
        int domSize = (pMove.first)->values_.size();

        if(!foundFirst){
            duration = domSize;
        }else{
            ++freq;

            if(findBetter)
                duration = domSize;

            if(freq > Freq){
                if((duration/domSize) < conflictNb ){
                    duration += domSize;
                    freq = 0;
                }
            }
        }

        //        cout << currentCost
        //                << "\t" << bestCost
        //                << "\t" << conflictNb
        //                << "\t" << duration
        //                << endl;

        /// ================================================ END
        /// ====================================================





        (pMove.first)->current()->pushTabu(duration);

        /// Update best move
        BOOST_FOREACH(CtrPtr const ctr, (pMove.first)->ctrs_){
            if(ctr->getGamma() > 0){
                ctr->updateGamma(pMove.first,(pMove.first)->current(),pMove.second);
            }
        }


        (pMove.first)->assign(pMove.second);

        /// add weights only on one movement variable ================= BGN
        //        if((pMove.second)->getGamma() > 0)
        //            (pMove.first)->addWeight();
        /// add weights only on one movement variable ================= END

        //        if((pMove.first)->varWeight_ > 0)
        //            cout << "w: current: " << (pMove.first)->varWeight_ << endl;

        bestMove.reset();

        /// Stop criteria
        if(foundFirst || currentIter < 0){
            //if(currentIter < 0){
            //if(currentIter < 0){
            // set the best solution found
            getSolution(filter_first, filter_last);

            return isSolutionFeasible(filter_first, filter_last);
        }



    }

}




struct Breakout_UpdateWeight{

    Breakout_UpdateWeight():weight_(0){}

    void addWeight(VarPtr  v){
        //weight_ = 0;

        if(!v->isViolated())
            return;

        // traveling the weight of constraints
        BOOST_FOREACH(CtrPtr const c, v->ctrs_){

            //c->setGamma(c->getGamma() + 1);
            if(c->isViolated()){
                ++weight_;
                BOOST_FOREACH(ValPtr const a, v->values_){
                    if(c->isViolated(v, a)){
                        a->setGamma(a->getGamma() + 2);
                        //if(a == v->current())
                    }
                }
            }

        }

    }

    void addGamma(VarPtr  v){
        if(!v->isViolated())
            return;

        // traveling the weight of constraints
        BOOST_FOREACH(CtrPtr const c, v->ctrs_){
            if(c->isViolated())
                c->setGamma(c->getGamma() + 1);
        }
    }

    int getWeight(){
        return weight_;
    }

    int resetWeight(){
        weight_ = 0;
    }

private:
    int weight_;
};


/// breakout algorithm -- not working well since the way of implementation
/// take account of incremental weights
template <typename Itr, typename Pred>
        bool breakoutOLD(Itr first, Itr last, Pred pred){


    //return true;
    //cout << "in tabu solver" << endl;

    base_generator_type generator(42u);
    generator.seed(static_cast<unsigned int>(std::time(0)));
    boost::uniform_int<> degen_dist(0,9);
    boost::variate_generator<base_generator_type&, boost::uniform_int<> >
            deg(generator, degen_dist);

    boost::uniform_int<> degen_distx(70,90);
    boost::variate_generator<base_generator_type&, boost::uniform_int<> >
            degCELAR(generator, degen_distx);


    register int currentCost = 0;
    register int bestCost = 0;
    bool foundFirst = false;


    /// Initialize the variables filter
    typedef boost::filter_iterator<Pred, Itr> FilterIter;
    FilterIter filter_first(pred, first, last);
    FilterIter filter_last(pred, last,last);

    //register int MaxIter = std::distance(first,last);
    register int MaxIter = 1000;
    register int currentIter = MaxIter;
    register int Freq = 100;
    register int freq = 0;
    int duration = 0;

    register int breakoutCnt = std::distance(filter_first,filter_last)/2;
    //register int breakoutCnt = 000;


    // intialize values' gamma
    InitGamma initG;
    std::for_each(filter_first,filter_last,
                  boost::bind(&InitGamma::init,
                              &initG,_1));



    Get_best_move bestMove; //functor consists of finding best move
    Breakout_UpdateWeight breakoutWeight;


    while(true){

        //cout << "before bestMove" << endl;
        /// Find best move
        std::for_each(filter_first, filter_last,
                      boost::bind(&Get_best_move::evaluate,
                                  &bestMove,_1));
        //cout << "after bestMove" << endl;



        //        cout << "after getBest" << endl;
        //        cout << "before getConflict" << endl;
        int conflictNb = bestMove.getConflictNb();


        //cout << "best Delta: " << bestDelta << endl;
        //        cout << "conflictNb: " << conflictNb << endl;
        if(conflictNb < 1){
            std::for_each(filter_first, filter_last, boost::bind(&Variable::resetDom,_1));
            std::for_each(filter_first, filter_last, boost::bind(&Variable::clearSolution,_1));
            bestMove.reset();
            //            cout << "found feasible " << endl;
            return true;
        }


        /// Update tabu tenure

        std::pair<VarPtr, ValPtr> pMove = bestMove.getBestMove();
        if((pMove.second)->isTabu())
            (pMove.second)->popTabu();

        int bestDelta = (pMove.second)->getGamma()
                        - (pMove.first)->current()->getGamma();

        bool findBetter = false;
        if(!foundFirst){
            // try to find first move which accumulates the cost
            if(bestDelta >= 0){
                foundFirst = true;
                // update best move
                setSolution(filter_first, filter_last);
                currentCost += bestDelta;
            }
        }else{
            --currentIter;
            currentCost += bestDelta;
            if(currentCost < bestCost){
                //cout << "find better ========================" << endl;
                // udpate best move
                setSolution(filter_first, filter_last);
                bestCost = currentCost;
                currentIter = MaxIter;
                findBetter = true;
            }

        }



        /// ====================================================
        /// ================================================ BGN

        /// stratgey 1:
        // compute the tabu tenure and current objective funtion
        // duration = deg() + boost:: numeric_cast<int>(0.6*conflictNb);

        /// stratgey 2:
        //duration = degCELAR();


        /// stratgey 3:
        int domSize = (pMove.first)->values_.size();

        if(!foundFirst){
            duration = domSize;
        }else{
            ++freq;

            if(findBetter)
                duration = domSize;

            if(freq > Freq){
                if((duration/domSize) < conflictNb ){
                    duration += domSize;
                    freq = 0;
                }
            }
        }

        //        cout << currentCost
        //                << "\t" << bestCost
        //                << "\t" << conflictNb
        //                << "\t" << duration
        //                << endl;

        /// ================================================ END
        /// ====================================================





        (pMove.first)->current()->pushTabu(duration);

        /// Update best move
        BOOST_FOREACH(CtrPtr const ctr, (pMove.first)->ctrs_){
            if(ctr->getGamma() > 0){
                ctr->updateGamma(pMove.first,(pMove.first)->current(),pMove.second);
            }
        }


        (pMove.first)->assign(pMove.second);



        /// combine the weights in gamma system ======================= BGN
        if(foundFirst && !findBetter){
            std::for_each(filter_first,filter_last, boost::bind(&Breakout_UpdateWeight::addWeight, &breakoutWeight, _1));
            std::for_each(filter_first,filter_last, boost::bind(&Breakout_UpdateWeight::addGamma, &breakoutWeight, _1));
            bestCost += breakoutWeight.getWeight();


            currentCost += breakoutWeight.getWeight();
            //cout << "best: " << bestCost << "\tcurrent: " << currentCost << "\tWeight: " << breakoutWeight.getWeight() << endl;
            breakoutWeight.resetWeight();
            -- breakoutCnt;
        }

        /// combine the weights in gamma system ======================= END

        /// add weights only on one movement variable ================= BGN
        //        if((pMove.second)->getGamma() > 0)
        //            (pMove.first)->addWeight();
        /// add weights only on one movement variable ================= END

        //        if((pMove.first)->varWeight_ > 0)
        //            cout << "w: current: " << (pMove.first)->varWeight_ << endl;

        bestMove.reset();

        /// Stop criteria
        if(currentIter < 0 || breakoutCnt < 0){
            // set the best solution found
            getSolution(filter_first, filter_last);

            return isSolutionFeasible(filter_first, filter_last);
        }



    }

}



void addWeightOnVar(VarPtr const v){
    bool hasViolated = false;
    BOOST_FOREACH(CtrPtr const c, v->ctrs_){
        if(c->isViolated()){
            c->setGamma(c->getGamma() + 1);
            if(!hasViolated)
                hasViolated = true;
        }
    }

    if(hasViolated)
        v->addWeight();
}

/// breakout algorithm 2, regarding breakout as a suite of iterative tabu searches
template <typename Itr, typename Pred>
        bool breakout(Itr first, Itr last, Pred pred){

    register int breakoutItr = 10; // define the breakout bound
    for(Itr current = first; current!=last; ++current){
        int d = (*first)->vars_.size();
        if(d>breakoutItr)
            breakoutItr = d;
    }

    ++breakoutItr;

    bool feasible = false;
    while(breakoutItr > 0){
        feasible = tabuShort(first,last,pred);

        if(feasible)
            return feasible;

        // increase the gamma on violated constraints
        //std::for_each(first,last,addWeightOnVar);
        --breakoutItr;
    }


    return feasible;
}



/// =================================================================================================
/// Growing solver with only one gamma system
/// =================================================================================================

template <typename Iter>
        bool hasCore(Iter first, Iter last){

    for(;first != last; ++first){
        if((*first)->getGamma() > 1)
            if((*first)->isViolated())
                return true;
    }

    return false;

}





struct OneViolatedVars{

    OneViolatedVars():isSet_(false){}

    void locate(VarPtr  v){
        if(!v->isViolated()){
            v->setGamma(0);
            return;
        }

        if(isSet_ && maxWeight_->getWeight() >= v->getWeight()){
            v->setGamma(0);
            return;
        }

        if(!isSet_){
            isSet_ = true;
        }else{
            maxWeight_->setGamma(0);
        }

        maxWeight_ = v;
    }

    void set(){
        if(isSet_){
            isSet_ = false;
            bool found = false;
            VarPtr nextVar;
            BOOST_FOREACH(CtrPtr const c, maxWeight_->ctrs_){
                if(c->getGamma() < 1 || !c->isViolated())
                    continue;

                if(found && c->getOther(maxWeight_)->getWeight() <= nextVar->getWeight())
                    continue;

                if(!found)
                    found = true;

                nextVar = c->getOther(maxWeight_);

            }

            if(found){
                nextVar->setGamma(1);
            }

        }


    }

    bool isSet_;
    VarPtr maxWeight_;
};


struct AllViolatedVars{

    AllViolatedVars(){}

    void locate(VarPtr  v){
        if(!v->isViolated()){
            v->setGamma(0);
            return;
        }
    }

    void set(){
        return;


    }


};

// return true if a < b
//inline bool compareByWeightedVar(VarPtr const a, VarPtr const b){
//    if(a->isViolated()){
//        if(b->isViolated())
//            return a->getWeight() < b->getWeight();
//        else
//            return false;
//    }else{
//        return b->isViolated();
//    }
//}
//
bool hasGammaNeighbor(VarPtr const v){
    int great2 = 0;
    BOOST_FOREACH(VarPtr const a, v->vars_){
        if(a->getGamma() > 0){
            return true;
            ++great2;
        }

        if(great2 > 1)
            return true;
    }

    return false;
}

struct LinkToGammaVars{

    LinkToGammaVars(){}

    void extend(VarPtr const var){
        if(var->getGamma() > 0)
            return;

        if(hasGammaNeighbor(var))
            vars_.push_back(var);
    }

    void set(){

        if(vars_.empty())
            return;

        //cout << "find linked: " << vars_.size() << endl;
        BOOST_FOREACH(VarPtr const v, vars_)
                v->setGamma(1);

        vars_.clear();
    }


    std::vector<VarPtr> vars_;

};



struct LinkToWeightVars{

    LinkToWeightVars():isSet_(false){}

    void extend(VarPtr const var){
        if(var->getGamma() > 0 || !hasGammaNeighbor(var))
            return;

        if(isSet_ && maxWeight_->getWeight() >= var->getWeight())
            return;

        if(!isSet_)
            isSet_ = true;


        maxWeight_ = var;
    }

    void set(){

        if(isSet_){
            isSet_ = false;
            maxWeight_->setGamma(1);
        }

    }


    bool isSet_;
    VarPtr maxWeight_;

};



template <typename Element>
        struct PrintMark{

    void operator ()(Element const e){
        if(e->isMarked())
            e->print();
    }
};



void printCore(std::vector<VarPtr>& vars, std::vector<CtrPtr>& ctrs){

    std::for_each(vars.begin(),vars.end(), PrintMark<VarPtr>());
    std::for_each(ctrs.begin(),ctrs.end(), PrintMark<CtrPtr>());

}


template <typename InputIterator>
        void verify(InputIterator first, InputIterator last){
    int cnt = 0;
    int vcnt = 0;
    for(; first != last; ++first){
        if((*first)->isEnable()){
            ++cnt;
            if((*first)->isAssigned()){
                if((*first)->isViolated())
                    ++vcnt;
            }else{
                cout << "Solution not valid" << endl;
                return;
            }
        }

    }
    cout << "total_ctrs: " << cnt << " violated_ctrs: " << vcnt << endl;
    if(cnt == vcnt){
        cout << "feasible" << endl;
        exit(0);
    }
}


template <typename LocalizeVar, typename ConnectedVars>
        bool grow_solver(std::vector<VarPtr>& vars, std::vector<CtrPtr>& ctrs,
                         LocalizeVar localizeVar, ConnectedVars connectedVars){


    // first iteration

    int loopCnt = 1;
    bool feasible = tabu(vars.begin(), vars.end(), is_gamma_variable());
    if(feasible)
        return true;


    //========================= Locator
    /// The locator locates one or many violated constraints as the location of the IIS
    std::for_each(vars.begin(), vars.end(), boost::bind(&LocalizeVar::locate,&localizeVar, _1));
    localizeVar.set();

    //printCore(vars,ctrs,loopCnt);

    do{
        ++loopCnt;
        //=========== Clear all weight
        std::for_each(vars.begin(), vars.end(), boost::bind(&Variable::clearWeight, _1));

        // filter constraints
        std::for_each(ctrs.begin(), ctrs.end(), boost::bind(&AbsConstraint::regardGamma, _1));


        feasible = tabu(vars.begin(), vars.end(), is_gamma_variable());

        //printCore(vars,ctrs,loopCnt);

        if(feasible){
            std::for_each(vars.begin(), vars.end(), boost::bind(&ConnectedVars::extend,&connectedVars, _1));
            connectedVars.set();
        }

    }while(feasible);

    //printCore(vars,ctrs,loopCnt);

    return feasible;

}


/// =================================================================================================
/// Growing solver with mark system
/// marked elements with gamma = 100
/// elements arround with gamma = 1
/// others with gamma = 0
/// =================================================================================================


template <typename Iter>
        bool isMarkedCtrViolated(Iter first, Iter last){
    for(;first!=last; ++first){
        if((*first)->isMarked() && (*first)->isViolated())
            return true;
    }

    return false;
}

///


inline void gammaVar(VarPtr v){

    if(!v->isMarked())
        v->setGamma(1);
}


/// put strong gamma and non-violated constraint first
inline bool compareConstraintsWeights(CtrPtr const a, CtrPtr const b){



    //
    //    if(a->isMarked() && b->isMarked())
    //        return a->getVarsWeights() < b->getVarsWeights();

    if(a->isMarked())
        return true;

    if(b->isMarked())
        return false;


    //    if(a->getGamma() < 1 && b->getGamma() < 1)
    //        return a->getVarsWeights() < b->getVarsWeights();
    //
    //
    //    if(a->getGamma() < 1)
    //        return true;
    //
    //    if(b->getGamma() < 1)
    //        return false;


    if(a->getGamma() == 1 && b->getGamma() == 1){

        if(a->isViolated() && b->isViolated()){
            return a->getVarsWeights() < b->getVarsWeights();
        }

        if(b->isViolated())
            return true;

        if(a->isViolated())
            return false;

        return a->getVarsWeights() < b->getVarsWeights();
    }


    if(a->getGamma() > 1)
        return true;
    else
        return false;

}


inline bool compareVariablesWeights(VarPtr const a, VarPtr const b){
    if(a->isMarked())
        return true;

    if(b->isMarked())
        return false;


    /// add for icai2011 === begin
    //    if(a->getGamma() < 1)
    //        return false;

    //    if(b->getGamma() < 1)
    //        return true;
    /// add for icai2011 === end

    int aw = std::accumulate(a->vars_.begin(), a->vars_.end(),
                             0,Counter_markedWeight<VarPtr>())
    +a->getWeight();


    int bw = std::accumulate(b->vars_.begin(), b->vars_.end(),
                             0,Counter_markedWeight<VarPtr>())
    +b->getWeight();


    if(a->isViolated() && b->isViolated())
        return aw < bw;

    if(a->isViolated())
        return false;


    if(b->isViolated())
        return true;


    return aw < bw;


}

struct OneViolatedVarsWithExtension{

    OneViolatedVarsWithExtension():isSet_(false),weight_(0){}

    void locate(VarPtr  v){
        if(!v->isViolated()){
            v->setGamma(0);
            return;
        }

        int weightCurrent = std::accumulate(v->vars_.begin(),v->vars_.end(),
                                            0,Counter_gammaWeight<VarPtr>())
        + v->getWeight();

        if(isSet_ && weight_ >= weightCurrent){
            v->setGamma(0);
            return;
        }

        if(!isSet_){
            isSet_ = true;
        }else{
            maxWeight_->setGamma(0);
        }

        maxWeight_ = v;
        weight_ = weightCurrent;
    }

    void set(){
        if(isSet_){
            isSet_ = false;


            CtrPtr chosenCtr = *std::max_element(maxWeight_->ctrs_.begin(), maxWeight_->ctrs_.end(),
                                                 compareConstraintsWeights);

            VarPtr nextVar = chosenCtr->getOther(maxWeight_);

            //            if(chosenCtr->isViolated())
            //                cout << "chosen cter is violated" << endl;


            if(nextVar->getGamma() == 0){

                nextVar->setGamma(100);
                nextVar->mark();

                maxWeight_->mark();
                maxWeight_->setGamma(100);

                //                cout << "set(): " << maxWeight_->id_ << endl;
                //                cout << "set(): " << nextVar->id_ << endl;

                // valid surround variable
                std::for_each(nextVar->vars_.begin(), nextVar->vars_.end(),gammaVar);
                std::for_each(maxWeight_->vars_.begin(), maxWeight_->vars_.end(),gammaVar);

            }else{
                cout << "couldn't find constraint: OneViolatedVarsWithExtension::set()" << endl;
                exit(0);
            }

        }


    }

    bool isSet_;
    int weight_;
    VarPtr maxWeight_;

};


struct AllViolatedVarsWithExtension{

    AllViolatedVarsWithExtension(){}

    void locate(VarPtr  v){
        if(!v->isViolated()){
            v->setGamma(0);
            return;
        }

        v->setGamma(100);
        v->mark();
        marked_.push_back(v);
    }

    void set(){
        if(!marked_.empty()){

            BOOST_FOREACH(VarPtr const v, marked_){
                std::for_each(v->vars_.begin(), v->vars_.end(),
                              gammaVar);
            }


            marked_.clear();
        }


    }

    //bool isSet_;
    std::vector<VarPtr> marked_;
    //CtrPtr chosenCtr_;
};



struct OneViolatedVarsWithExtension_2function{

    OneViolatedVarsWithExtension_2function():isSet_(false){}

    void locate(VarPtr  v){
        if(!v->isViolated()){
            v->setGamma(100);
            return;
        }

        if(isSet_ && maxWeight_->getWeight() >= v->getWeight()){
            v->setGamma(100);
            return;
        }

        if(!isSet_){
            isSet_ = true;
        }else{
            maxWeight_->setGamma(100);
        }

        maxWeight_ = v;
    }

    void set(){
        if(isSet_){
            isSet_ = false;


            CtrPtr chosenCtr = *std::max_element(maxWeight_->ctrs_.begin(), maxWeight_->ctrs_.end(),
                                                 compareConstraintsWeights);

            VarPtr nextVar = chosenCtr->getOther(maxWeight_);

            //            if(chosenCtr->isViolated())
            //                cout << "chosen cter is violated" << endl;


            if(nextVar->getGamma() > 1){

                nextVar->mark();

                maxWeight_->mark();
                maxWeight_->setGamma(100);

                //                cout << "set(): " << maxWeight_->id_ << endl;
                //                cout << "set(): " << nextVar->id_ << endl;

                // valid surround variable
                std::for_each(nextVar->vars_.begin(), nextVar->vars_.end(),gammaVar);
                std::for_each(maxWeight_->vars_.begin(), maxWeight_->vars_.end(),gammaVar);

            }else{
                cout << "couldn't find constraint: OneViolatedVarsWithExtension::set()" << endl;
                exit(0);
            }

        }


    }

    bool isSet_;
    VarPtr maxWeight_;

};




inline void resetVar(VarPtr var){
    if(var->isMarked())
        var->demark();

    var->setGamma(1);
    var->resetDom();
    var->clearWeight();
}

inline void resetCtr(CtrPtr ctr){
    if(ctr->isMarked())
        ctr->demark();

    ctr->setGamma(1);
}




/// Locator: which locates the violated constraint/constraints
/// Objective: find the constraint/s inside IIS
template <typename LocalizeVar>
        bool locator(std::vector<VarPtr>& vars, std::vector<CtrPtr>& ctrs,
                     LocalizeVar& localizeVar){


    bool feasible = tabu(vars.begin(), vars.end(), is_gamma_variable());
    if(feasible){
        cout << "locator find feasible"<< endl;
        exit(0);
        return true;
    }
    /// how many constraints violated


    cout << "\ttabu\t"
            << std::accumulate( vars.begin(), vars.end(), 0, Counter_violated<VarPtr>() )
            << "\t"
            <<std::accumulate( ctrs.begin(), ctrs.end(), 0, Counter_violated<CtrPtr>() ) ;//<<  endl;


    //========================= Locator
    /// The locator locates one or many violated constraints as the location of the IIS
    std::for_each(vars.begin(), vars.end(), boost::bind(&LocalizeVar::locate,&localizeVar, _1));
    localizeVar.set();

    return false;

}



struct Extension{


    Extension():isSet_(false){}

    // link to marked
    void extend(VarPtr var){

        // only verify marked variables
        if(var->isMarked()){
            CtrPtr ctr = *std::max_element(var->ctrs_.begin(), var->ctrs_.end(),
                                           compareConstraintsWeights);
            // only verify the constraint whose gamma = 1
            if(ctr->getGamma() == 1){
                if(isSet_ && !compareConstraintsWeights(chosenCtr_, ctr))
                    return;

                if(!isSet_)
                    isSet_ = true;

                chosenCtr_ = ctr;
                maxWeight_ = ctr->getOther(var);

            }
        }

    }


    void extendZone(VarPtr var){

        // only verify marked variables
        if(var->getGamma() == 1){
            std::for_each(var->vars_.begin(), var->vars_.end(), gammaVar);
        }

    }


    bool mark(){
        if(isSet_){
            isSet_ = false;
            maxWeight_->setGamma(100);
            maxWeight_->mark();

            //            cout << "mark(): " << maxWeight_->id_ << endl;

            std::for_each(maxWeight_->vars_.begin(), maxWeight_->vars_.end(),
                          gammaVar);
            return true;

        }else{

            return false;
            // maybe all marked constraints is a component

            cout << "couldn't find constraint: OneViolatedVarsWithExtension::mark()" << endl;
            exit(0);
        }

    }



    void verify(VarPtr const var){
        cout << "variable: " << var->id_ ;
        BOOST_FOREACH(VarPtr const v, var->vars_){
            if(v->getGamma() == 1){
                cout << " has no gamma neighbors !!! " << endl;
                return;
            }
        }
        cout << " has gamma neighbors" << endl;
    }


    void markall(VarPtr v){
        if(v->getGamma() == 1){
            std::for_each(v->vars_.begin(), v->vars_.end(),
                          gammaVar);
        }
    }

    bool isSet_;
    VarPtr maxWeight_;
    CtrPtr chosenCtr_;
};




struct ExtensionVars{


    ExtensionVars():isSet_(false){}

    // link to marked
    void extend(VarPtr var){

        // only verify marked variables
        if(var->isMarked()){
            VarPtr varx = *std::max_element(var->vars_.begin(), var->vars_.end(),
                                            compareVariablesWeights);
            // only verify the constraint whose gamma = 1
            if(var->getGamma() == 1){
                if(isSet_ && !compareVariablesWeights(maxWeight_, varx))
                    return;

                if(!isSet_)
                    isSet_ = true;

                maxWeight_ = varx;

            }
        }

    }


    void extendZone(VarPtr var){

        // only verify marked variables
        if(var->getGamma() == 1){
            std::for_each(var->vars_.begin(), var->vars_.end(), gammaVar);
        }

    }


    bool mark(){
        if(isSet_){
            isSet_ = false;
            maxWeight_->setGamma(100);
            maxWeight_->mark();

            //            cout << "mark(): " << maxWeight_->id_ << endl;

            std::for_each(maxWeight_->vars_.begin(), maxWeight_->vars_.end(),
                          gammaVar);
            return true;

        }else{

            return false;
            // maybe all marked constraints is a component

            cout << "couldn't find constraint: OneViolatedVarsWithExtension::mark()" << endl;
            exit(0);
        }

    }



    void markall(VarPtr v){
        if(v->getGamma() == 1){
            std::for_each(v->vars_.begin(), v->vars_.end(),
                          gammaVar);
        }
    }

    bool isSet_;
    VarPtr maxWeight_;
};


/// Extender: which extends subproblem arround located constraint/constraints
/// Objective:

struct OnlyclearWeightOnNonMarked{

    void operator () (VarPtr var){
        if(var->isMarked())
            return;

        var->clearWeight();

    }

};


//! ========================================================================
//! The
//! ========================================================================
template <typename LocalizeVar, typename ExtendVars>
        bool extender(std::vector<VarPtr>& vars, std::vector<CtrPtr>& ctrs,
                      LocalizeVar& localizeVar,ExtendVars& extension, boost::progress_timer& timer){

    //cout << "in extender" << endl;

    bool feasible = false;
    //int loopCnt = 0;
    int cntGamma = 0;
    while(true){

        if(timer.elapsed() > 800.0){
            cout << "Time out by extender" << endl;
            exit(0);
        }

        //cout << "loop: " << ++loopCnt << endl;


        //=========== Clear all weight
        std::for_each(vars.begin(), vars.end(), OnlyclearWeightOnNonMarked());

        //=========== affect the constraints according to variables
        std::for_each(ctrs.begin(), ctrs.end(), boost::bind(&AbsConstraint::regardMark, _1));


        //        cout << "\nBBBB\t"
        //                << std::accumulate( vars.begin(), vars.end(), 0, Counter_gamma<VarPtr>() )
        //                << "\t"
        //                <<std::accumulate( ctrs.begin(), ctrs.end(), 0, Counter_gamma<CtrPtr>() )<< endl;



        feasible = tabu(vars.begin(), vars.end(), is_gamma_variable());

        // print all the cores information =============================
        //printCore(vars,ctrs,loopCnt);

        // in case the tabu find a feasible subproblem
        if(feasible){
            //cout << "tabu feasible" << endl;
            std::for_each(vars.begin(), vars.end(),
                          boost::bind(&ExtendVars::extendZone,&extension, _1));

            int cnt = std::accumulate(vars.begin(), vars.end(),
                                      0, Counter_gamma<VarPtr>());
            if(cntGamma != cnt){
                cntGamma = cnt;
                continue;
            }else{
                //cout << "in re-localize phase needed" << endl;
                std::for_each(vars.begin(), vars.end(), resetVar);
                std::for_each(ctrs.begin(), ctrs.end(), resetCtr);
                //OneViolatedVarsWithExtension violatedVars;
                feasible = locator(vars, ctrs,localizeVar);
                if(feasible)
                    return true;
            }


            //continue;
        }


        if(!feasible && isMarkedCtrViolated(ctrs.begin(), ctrs.end()))
            break;


        /// Extension: enable surround the variables
        std::for_each(vars.begin(), vars.end(), boost::bind(&ExtendVars::extend,&extension, _1));


        /// Mark another constraint/ constraints
        if(!extension.mark()){
            // reset all gamma
            //cout << "in re-localize phase needed" << endl;
            std::for_each(vars.begin(), vars.end(), resetVar);
            std::for_each(ctrs.begin(), ctrs.end(), resetCtr);
            //OneViolatedVarsWithExtension violatedVars;
            feasible = locator(vars, ctrs,localizeVar);
            if(feasible)
                return true;
            //exit(0);

        }

        //        cout << "\nMMMM\t"
        //                << std::accumulate( vars.begin(), vars.end(), 0, Counter_marked<VarPtr>() )
        //                << "\t"
        //                <<std::accumulate( ctrs.begin(), ctrs.end(), 0, Counter_marked<CtrPtr>() )
        //                <<  endl;


    }// check if there is marked constraints violated

    cout << "\tsubp\t"
            << std::accumulate( vars.begin(), vars.end(), 0, Counter_gamma<VarPtr>() )
            << "\t"
            <<std::accumulate( ctrs.begin(), ctrs.end(), 0, Counter_gamma<CtrPtr>() );//<<  endl;



    //        cout << "\tcore\t"
    //                << std::accumulate( vars.begin(), vars.end(), 0, Counter_marked<VarPtr>() )
    //                << "\t"
    //                <<std::accumulate( ctrs.begin(), ctrs.end(), 0, Counter_marked<CtrPtr>() )
    //                ;//<<  endl;

    return false;

    //exit(0);

}


void enableMarkArround(VarPtr var){
    if(var->isMarked()){
        std::for_each(var->vars_.begin(), var->vars_.end(),gammaVar );
    }
}


bool verificator(std::vector<VarPtr>& vars, std::vector<CtrPtr>& ctrs,
                 boost::progress_timer& timer){

    //cout << "in verificator" << endl;

    std::vector<VarPtr> core;
    if(!core.empty())
        core.clear();


    std::remove_copy_if(vars.begin(), vars.end(),std::back_inserter(core),
                        not_marked_element<VarPtr>());
    std::for_each(ctrs.begin(), ctrs.end(),
                  boost::bind(&AbsConstraint::markByVars,_1));

    cout << "\tcore\t"
            << std::accumulate( vars.begin(), vars.end(), 0, Counter_marked<VarPtr>() )
            << "\t"
            <<std::accumulate( ctrs.begin(), ctrs.end(), 0, Counter_marked<CtrPtr>() )
            ;//<<  endl;

    bool feasible = mac_solver(core.begin(), core.end(),timer);

    return feasible;
}




/// Extension strategy:
/// 1. enable surround variable
/// 2. mark the violated one
template <typename LocalizeVar, typename ExtendVars>
        bool breakout_mark_solver(std::vector<VarPtr>& vars, std::vector<CtrPtr>& ctrs,
                                  LocalizeVar& localizeVar, ExtendVars& extension){

    //cout << "in breakout_mark_solver"<< endl;

    // first iteration

    int loopCnt = 0;
    bool feasible = breakout(vars.begin(), vars.end(), is_gamma_variable());
    std::for_each(ctrs.begin(),ctrs.end(),boost::bind(&AbsConstraint::resetGamma,_1));
    if(feasible)
        return true;

    /// how many constraints violated

    
    cout << "\ttabu\t"
            << std::accumulate( vars.begin(), vars.end(), 0, Counter_violated<VarPtr>() )
            << "\t"
            <<std::accumulate( ctrs.begin(), ctrs.end(), 0, Counter_violated<CtrPtr>() ) ;//<<  endl;


    //========================= Locator
    /// The locator locates one or many violated constraints as the location of the IIS
    std::for_each(vars.begin(), vars.end(), boost::bind(&LocalizeVar::locate,&localizeVar, _1));
    localizeVar.set();

    //printCore(vars,ctrs,loopCnt);

    int cntGamma = 0;
    while(true){

        ++loopCnt;


        //=========== Clear all weight
        std::for_each(vars.begin(), vars.end(), boost::bind(&Variable::clearWeight, _1));

        //=========== affect the constraints according to variables
        std::for_each(ctrs.begin(), ctrs.end(), boost::bind(&AbsConstraint::regardMark, _1));



        feasible = tabu(vars.begin(), vars.end(), is_gamma_variable());

        // print all the cores information =============================
        //printCore(vars,ctrs,loopCnt);

        if(feasible){
            std::for_each(vars.begin(), vars.end(),
                          boost::bind(&ExtendVars::extendZone,&extension, _1));

            int cnt = std::accumulate(vars.begin(), vars.end(),
                                      0, Counter_gamma<VarPtr>());
            if(cntGamma != cnt){
                cntGamma = cnt;
                continue;
            }else{
                //cout << "in re-localize phase needed" << endl;
                std::for_each(vars.begin(), vars.end(), resetVar);
                std::for_each(ctrs.begin(), ctrs.end(), resetCtr);
                //OneViolatedVarsWithExtension violatedVars;
                feasible = breakout(vars.begin(), vars.end(), is_gamma_variable());
                std::for_each(ctrs.begin(),ctrs.end(),boost::bind(&AbsConstraint::resetGamma,_1));
                if(feasible)
                    return true;

                std::for_each(vars.begin(), vars.end(),
                              boost::bind(&LocalizeVar::locate,&localizeVar, _1));
                localizeVar.set();
            }
        }


        if(!feasible && isMarkedCtrViolated(ctrs.begin(), ctrs.end()))
            break;





        /// Extension: enable surround the variables
        std::for_each(vars.begin(), vars.end(),
                      boost::bind(&ExtendVars::extend,&extension, _1));


        /// Mark another constraint/ constraints
        if(!extension.mark()){
            // reset all gamma
            cout << "in re-localize phase" << endl;
            std::for_each(vars.begin(), vars.end(), resetVar);
            std::for_each(ctrs.begin(), ctrs.end(), resetCtr);
            feasible = tabu(vars.begin(), vars.end(), is_gamma_variable());
            if(feasible)
                return true;

            std::for_each(vars.begin(), vars.end(),
                          boost::bind(&LocalizeVar::locate,&localizeVar, _1));
            localizeVar.set();
        }


    }// check if there is marked constraints violated

    /// only for coloring ============ BGN
    int nbVars = std::accumulate(vars.begin(), vars.end(), 0, Counter_marked<VarPtr>());
    int nbCtrs = std::accumulate(ctrs.begin(), ctrs.end(), 0, Counter_marked<CtrPtr>());

    cout << "Mvars: " << nbVars << " Mctrs: " << nbCtrs << endl;
    if(nbCtrs == (nbVars-1)*nbVars/2){
        cout << "find clique" << endl;
    }else{
        cout << "not  clique" << endl;
    }
    /// only for coloring ============ END

    //exit(0);

    return false;

}





//! ========================================================================
//! The algorithm -- in aaai2011 paper
//! ========================================================================
template <typename LocalizeVar, typename ExtendVars>
        bool extension_mark_solver(std::vector<VarPtr>& vars, std::vector<CtrPtr>& ctrs,
                                   LocalizeVar& localizeVar, ExtendVars& extension){

    //cout << "in extension_mark_solver"<< endl;

    // first iteration

    int loopCnt = 0;
    bool feasible = tabu(vars.begin(), vars.end(), is_gamma_variable());

    if(feasible)
        return true;

    /// how many constraints violated


    cout << "\ttabu\t"
            << std::accumulate( vars.begin(), vars.end(), 0, Counter_violated<VarPtr>() )
            << "\t"
            <<std::accumulate( ctrs.begin(), ctrs.end(), 0, Counter_violated<CtrPtr>() ) ;//<<  endl;


    //========================= Locator
    /// The locator locates one or many violated constraints as the location of the IIS
    std::for_each(vars.begin(), vars.end(), boost::bind(&LocalizeVar::locate,&localizeVar, _1));
    localizeVar.set();

    //printCore(vars,ctrs,loopCnt);

    int cntGamma = 0;
    while(true){

        ++loopCnt;


        //=========== Clear all weight
        std::for_each(vars.begin(), vars.end(), boost::bind(&Variable::clearWeight, _1));

        //=========== affect the constraints according to variables
        std::for_each(ctrs.begin(), ctrs.end(), boost::bind(&AbsConstraint::regardMark, _1));



        feasible = tabu(vars.begin(), vars.end(), is_gamma_variable());

        // print all the cores information =============================
        //printCore(vars,ctrs,loopCnt);

        if(feasible){
            std::for_each(vars.begin(), vars.end(),
                          boost::bind(&ExtendVars::extendZone,&extension, _1));

            int cnt = std::accumulate(vars.begin(), vars.end(),
                                      0, Counter_gamma<VarPtr>());
            if(cntGamma != cnt){
                cntGamma = cnt;
                continue;
            }else{
                //cout << "in re-localize phase needed" << endl;
                std::for_each(vars.begin(), vars.end(), resetVar);
                std::for_each(ctrs.begin(), ctrs.end(), resetCtr);
                //OneViolatedVarsWithExtension violatedVars;
                feasible = tabu(vars.begin(), vars.end(), is_gamma_variable());
                if(feasible)
                    return true;

                std::for_each(vars.begin(), vars.end(),
                              boost::bind(&LocalizeVar::locate,&localizeVar, _1));
                localizeVar.set();
            }
        }


        if(!feasible && isMarkedCtrViolated(ctrs.begin(), ctrs.end()))
            return false;





        /// Extension: enable surround the variables
        std::for_each(vars.begin(), vars.end(),
                      boost::bind(&ExtendVars::extend,&extension, _1));


        /// Mark another constraint/ constraints
        if(!extension.mark()){
            // reset all gamma
            cout << "in re-localize phase" << endl;
            std::for_each(vars.begin(), vars.end(), resetVar);
            std::for_each(ctrs.begin(), ctrs.end(), resetCtr);
            feasible = tabu(vars.begin(), vars.end(), is_gamma_variable());
            if(feasible)
                return true;

            std::for_each(vars.begin(), vars.end(),
                          boost::bind(&LocalizeVar::locate,&localizeVar, _1));
            localizeVar.set();
        }


    }// check if there is marked constraints violated

    return false;

}



struct ReachInfeasible{

    ReachInfeasible():isSet_(false),hasFirst_(false){}


    /// in case of feasible, include the neighbors
    void extend(VarPtr var){

        // only verify marked variables
        if(var->getGamma()>0){
            BOOST_FOREACH(VarPtr const v, var->vars_){
                if(v->getGamma() < 1){
                    neighbors_.push_back(v);
                    if(!isSet_)
                        isSet_ = true;
                }
            }
        }

    }


    void validVars(){
        if(isSet_){
            std::for_each(neighbors_.begin(),neighbors_.end(),gammaVar);
            isSet_ = false;
            neighbors_.clear();
        }
    }


    /// methods to extract IS
    void identify(VarPtr var){
        // in case there is no first violated
        if(!hasFirst_){
            if(var->getGamma()>0 && var->isViolated()){


                if(isSet_ && !compareVariablesWeights(maxWeight_, var))
                    return;

                if(!isSet_)
                    isSet_ = true;

                maxWeight_ = var;
            }

        }else{
            if(var->isMarked()){
                BOOST_FOREACH(VarPtr v, var->vars_){
                    if(v->getGamma() == 1 && v->isViolated()){
                        if(isSet_ && !compareVariablesWeights(maxWeight_, v))
                            continue;

                        if(!isSet_)
                            isSet_ = true;

                        maxWeight_ = v;
                    }
                }
            }


        }


    }

    bool mark(){
        if(isSet_){
            isSet_ = false;

        }else{
            return false;
            cout << "couldn't find constraint: ReachInfeasible::mark()" << endl;
            exit(0);
        }




        maxWeight_->setGamma(100);
        maxWeight_->mark();

        //            cout << "mark(): " << maxWeight_->id_ << endl;

        if(!hasFirst_){
            hasFirst_ = true;
            VarPtr varx ;

            bool set = false;
            BOOST_FOREACH(VarPtr v, maxWeight_->vars_){
                if(v->getGamma() > 0 && v->isViolated()){
                    if(set && !compareVariablesWeights(varx, v))
                        continue;

                    if(!set)
                        set = true;

                    varx = v;
                }
            }

            varx->setGamma(100);
            varx->mark();

        }



        return true;
    }


    bool isSet_;
    bool hasFirst_;
    std::vector<VarPtr> neighbors_;// the neighbors of valid variables
    VarPtr maxWeight_;

};



//! ========================================================================
//! The algorithm -- in icai2011 paper
//! 1. identify critical variables by means of breakout
//! 2. tabu until reaches unsolvable subproblem (by means heuristic)
//! 3. force (marked constraints) satisfiability on subproblem
//! ========================================================================
template <typename LocalizeVar, typename ExtendVars>
        bool breakout_heuristic_mark_solver(std::vector<VarPtr>& vars, std::vector<CtrPtr>& ctrs,
                                            LocalizeVar& localizeVar, ExtendVars& extension){

    cout << "in breakout_mark_solver"<< endl;

    // first iteration


    /// applying breakout algorithm to idenify the hard subproblem
    int loopCnt = 0;
    bool feasible = breakout(vars.begin(), vars.end(), is_gamma_variable());
    std::for_each(ctrs.begin(),ctrs.end(),boost::bind(&AbsConstraint::resetGamma,_1));
    if(feasible)
        return true;


    /// how many constraints violated
    cout << "\ttabu\t"
            << std::accumulate( vars.begin(), vars.end(), 0, Counter_violated<VarPtr>() )
            << "\t"
            <<std::accumulate( ctrs.begin(), ctrs.end(), 0, Counter_violated<CtrPtr>() ) ;//<<  endl;


    //========================= Locator
    /// The locator locates one or many violated constraints as the location of the IIS
    std::for_each(vars.begin(), vars.end(), boost::bind(&LocalizeVar::locate,&localizeVar, _1));
    localizeVar.set();

    //printCore(vars,ctrs,loopCnt);

    int cntGamma = 0;
    while(true){

        ++loopCnt;


        //=========== Clear all weight
        std::for_each(vars.begin(), vars.end(), boost::bind(&Variable::clearWeight, _1));

        //=========== affect the constraints according to variables
        std::for_each(ctrs.begin(), ctrs.end(), boost::bind(&AbsConstraint::regardMark, _1));



        feasible = tabu(vars.begin(), vars.end(), is_gamma_variable());

        // print all the cores information =============================
        //printCore(vars,ctrs,loopCnt);

        /// True: in case the tabu find feasible subproblem, extend the subproblem by
        /// adding all its neighbors varaibles
        /// False: in extract routine, extend the IS by marking violated constraints
        /// one by one
        if(feasible){
            std::for_each(vars.begin(), vars.end(),
                          boost::bind(&ExtendVars::extend,&extension, _1));


            extension.validVars();


            int cnt = std::accumulate(vars.begin(), vars.end(),
                                      0, Counter_gamma<VarPtr>());
            if(cntGamma != cnt){
                cntGamma = cnt;
                //cout << "cntGamma " << cntGamma << endl;
                continue;
            }else{
                //cout << "in re-localize phase needed" << endl;
                std::for_each(vars.begin(), vars.end(), resetVar);
                std::for_each(ctrs.begin(), ctrs.end(), resetCtr);
                //OneViolatedVarsWithExtension violatedVars;
                feasible = breakout(vars.begin(), vars.end(), is_gamma_variable());
                std::for_each(ctrs.begin(),ctrs.end(),boost::bind(&AbsConstraint::resetGamma,_1));
                if(feasible)
                    return true;

                std::for_each(vars.begin(), vars.end(),
                              boost::bind(&LocalizeVar::locate,&localizeVar, _1));
                localizeVar.set();
            }
        }

        if(!feasible && isMarkedCtrViolated(ctrs.begin(), ctrs.end()))
            break;





        /// Extension: enable surround the variables
        std::for_each(vars.begin(), vars.end(),
                      boost::bind(&ExtendVars::identify,&extension, _1));


        /// Mark another constraint/ constraints
        if(!extension.mark()){
            // reset all gamma
            //cout << "in re-localize phase heuristic" << endl;
            std::for_each(vars.begin(), vars.end(), resetVar);
            std::for_each(ctrs.begin(), ctrs.end(), resetCtr);
            feasible = tabu(vars.begin(), vars.end(), is_gamma_variable());
            if(feasible)
                return true;

            std::for_each(vars.begin(), vars.end(),
                          boost::bind(&LocalizeVar::locate,&localizeVar, _1));
            localizeVar.set();
        }


    }// check if there is marked constraints violated


    // only for coloring problem
    //    int nbVars = std::accumulate(vars.begin(), vars.end(), 0, Counter_marked<VarPtr>());
    //    int nbCtrs = std::accumulate(ctrs.begin(), ctrs.end(), 0, Counter_marked<CtrPtr>());

    //    cout << "Mvars: " << nbVars << " Mctrs: " << nbCtrs << endl;
    //    if(nbCtrs == (nbVars-1)*nbVars/2){
    //        cout << "find clique" << endl;
    //    }else{
    //        cout << "not  clique" << endl;
    //    }

    //    exit(0);

    return false;

}





#endif // ALGORITHM_HPP
