#ifndef COLOR_HPP
#define COLOR_HPP

#include <set>

// lambda if then need to be combined with lambda::bind
#include <boost/lambda/if.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/lambda/lambda.hpp>


#include "heuristic.hpp"


/// ====================================
/// Locator for coloring problem
/// ====================================

template <typename Iter, typename PredCtr>
bool colorHasViolatedCtr(Iter first,Iter last,PredCtr ctrPred){
    for(;first!=last;++first){
        if(ctrPred(*first) && (*first)->isViolated())
            return true;
    }

    return false;
}

template <typename VarItr, typename CtrItr,typename PredVar, typename PredCtr>
bool colorLocator(VarItr vfirst, VarItr vlast,
                  CtrItr cfirst, CtrItr clast,
                  PredVar varPred, PredCtr ctrPred){


    MinGamma<false> minGamma;

    int gammaValue = 0;
    for(VarItr current = vfirst; current != vlast; ++current){
        int d = (*current)->vars_.size();
        if(d>gammaValue)
            gammaValue = d;
    }

    int cnt = -1;
    while(true){

        cout << "tabu\t" << ++cnt << endl;
        bool feasible = tabuCol(vfirst,vlast,cfirst,clast,minGamma,varPred,ctrPred);
        if(feasible)
            return true;

        // mark violated constraints and check the stop criteria

        if(cnt < 1){
            for(CtrItr current = cfirst; current != clast;++current){
                if((*current)->isViolated()){

                    (*current)->setGammaS(gammaValue);
                    (*current)->setLevel(marked);
                    (*current)->setVarsLevel(marked);

                }
            }

        }else{

            for(CtrItr current = cfirst; current != clast;++current){

                if((*current)->isViolated()){
                    if((*current)->getLevel() == marked)
                        return false;
                    else{
                        // only mark the violated constraints linked to marked
                        if((*current)->varLevelgreaterThan(normal)){
                            (*current)->setGammaS(gammaValue);
                            (*current)->setLevel(marked);
                            (*current)->setVarsLevel(marked);
                        }
                    }
                }

            }
        }

        // saturation
        for(CtrItr current = cfirst; current != clast;++current){
            if((*current)->getLevel() != marked){
                if((*current)->sameVarsLevel()
                        && (*current)->getVarsLevel() == marked){
                    (*current)->setGammaS(gammaValue);
                    (*current)->setLevel(marked);
                }
            }
        }

    }

}


/// ====================================
/// Constructor for coloring problem
/// remove all non marked elements
/// marked one
/// ====================================



bool isDeadEnd(VarPtr const v){
    std::set<int> colors;
    BOOST_FOREACH(VarPtr const a, v->vars_){
        colors.insert(a->current()->v_);
    }

    return colors.size() == v->values_.size();

}


template <typename VarItr, typename CtrItr>
bool colorConstructor(VarItr vfirst, VarItr vlast, CtrItr cfirst, CtrItr clast){

    for(CtrItr current = cfirst; current != clast; ++current){
        (*current)->setGammaS(1);


        (*current)->setLevel(normal);
        (*current)->setVarsLevel(normal);

    }


    for(CtrItr current = cfirst; current != clast; ++current){
        if((*current)->isViolated()){
            (*current)->setLevel(marked);
            (*current)->setVarsLevel(marked);
            break;
        }
    }


    std::for_each(cfirst,clast,boost::bind(&AbsConstraint::initVarsGamma,_1));

    // find deadend in marked constraint
    // only check its marked neighbors
    bool foundDeadEnd = false;

    for(VarItr current = vfirst; current != vlast; ++current){
        if((*current)->getLevel() == marked){
            if(isDeadEnd(*current)){
                foundDeadEnd = true;
                // mark all neighbors
                std::for_each((*current)->vars_.begin(),(*current)->vars_.end(),
                              boost::bind(&Variable::setLevel,_1,marked));
                break;
            }
        }
    }

    // exit when failure occurs
    if(!foundDeadEnd){
        cout << "Constructor: cannot find deadend" << endl;
        exit(0);
    }


    // saturate the contraints
    for(CtrItr current = cfirst; current != clast; ++current){
        if((*current)->sameVarsLevel()){
            if((*current)->getVarsLevel() == marked){
                (*current)->setGammaS(1);
                (*current)->setLevel(marked);
            }
        }
    }


    // resign all non marked variables
    for(VarItr current = vfirst; current!=vlast;++current){
        if((*current)->getLevel() != marked)
            (*current)->resign();
    }


    // extension procedure
    MinGamma<false> minGamma;
    while(tabuCol(vfirst,vlast,cfirst,clast,minGamma,
                  Marked<VarPtr>(),Marked<CtrPtr>())){

        // extend the neighbor using,
        // return true if there is no variable to extend

        // forward checking, prune the values from domains
        fc_consistency(vfirst,vlast,Normal<VarPtr>(),All<CtrPtr>());

        bool hasNonMarked = false;

        // add new variable into core
        int mrv = -1;
        VarPtr candidate;
        for(VarItr current = vfirst; current != vlast; ++current){
            if((*current)->getLevel() != marked){

                bool hasMarked = false;
                BOOST_FOREACH(VarPtr const v, (*current)->vars_){
                    if(v->getLevel() == marked){
                        hasMarked = true;
                        break;
                    }
                }

                if(!hasMarked)
                    continue;

                if(!hasNonMarked) hasNonMarked = true;
                int m = validDomainSize(*current);
                if(mrv < 0 || mrv > m){
                    mrv = m;
                    candidate = *current;
                }
            }
        }

        if(!hasNonMarked){
            cout << "Constructor: cannot find non marked" << endl;
            exit(0);
        }

        // reset all domain
        std::for_each(vfirst,vlast,boost::bind(&Variable::resetDomain,_1));
        candidate->setLevel(marked);
        candidate->assign(candidate->values_.front());

        // saturate the constraints
        for(CtrItr current = cfirst; current != clast; ++current){
            if((*current)->sameVarsLevel()){
                if((*current)->getVarsLevel() == marked){
                    (*current)->setGammaS(1);
                    (*current)->setLevel(marked);
                }
            }
        }

    }

    return false;

}



/// ====================================
/// Verificator for coloring problem
/// ====================================
template <typename VarItr, typename CtrItr, typename PredVar, typename PredCtr>
bool colorVerificator(VarItr vfirst, VarItr vlast,
                      CtrItr cfirst, CtrItr clast,
                      PredVar varPred, PredCtr ctrPred){




}





bool colorExtract(std::vector<VarPtr>& vars, std::vector<CtrPtr>& ctrs){
    //bool colorExtract(VarItr vfirst, VarItr vlast, CtrItr cfirst, CtrItr clast){


    randomSolver(vars.begin(), vars.end());

    Statistic stat;
    stat.init();

    // locator
    stat.restart();
    bool feasible = colorLocator(vars.begin(),vars.end(),ctrs.begin(),ctrs.end(),
                                 Normal<VarPtr>(), Normal<CtrPtr>());
    stat.setTime(CateLocate);

    if(feasible)
        return true;

    // remove all non marked elements
    cout << "problem:\t" << vars.size()
         << "\t" << ctrs.size() << endl;

    vars.erase(std::remove_if(vars.begin(), vars.end(), NonMarked<VarPtr>()),vars.end());
    ctrs.erase(std::remove_if(ctrs.begin(),ctrs.end(), NonMarked<CtrPtr>()),ctrs.end());

    for(std::vector<VarPtr>::iterator current = vars.begin(); current != vars.end();
        ++current){
        ((*current)->ctrs_).erase(
                    remove_if((*current)->ctrs_.begin(), (*current)->ctrs_.end(),
                              NonMarked<CtrPtr>() ),(*current)->ctrs_.end() );

        ((*current)->vars_).erase(
                    remove_if((*current)->vars_.begin(), (*current)->vars_.end(),
                              NonMarked<VarPtr>() ),(*current)->vars_.end() );
    }


    // =========================================== TIME 2 =========================== BGN

//    BOOST_FOREACH(CtrPtr const c, ctrs){
//        c->setGammaS(1);
//        c->setLevel(normal);
//        c->setVarsLevel(normal);
//    }

//    feasible = colorLocator(vars.begin(),vars.end(),ctrs.begin(),ctrs.end(),
//                            Normal<VarPtr>(), Normal<CtrPtr>());


//    vars.erase(std::remove_if(vars.begin(), vars.end(), NonMarked<VarPtr>()),vars.end());
//    ctrs.erase(std::remove_if(ctrs.begin(),ctrs.end(), NonMarked<CtrPtr>()),ctrs.end());

//    for(std::vector<VarPtr>::iterator current = vars.begin(); current != vars.end();
//        ++current){
//        ((*current)->ctrs_).erase(
//                    remove_if((*current)->ctrs_.begin(), (*current)->ctrs_.end(),
//                              NonMarked<CtrPtr>() ),(*current)->ctrs_.end() );

//        ((*current)->vars_).erase(
//                    remove_if((*current)->vars_.begin(), (*current)->vars_.end(),
//                              NonMarked<VarPtr>() ),(*current)->vars_.end() );
//    }


    // =========================================== TIME 2 =========================== END
    std::for_each(ctrs.begin(),ctrs.end(),boost::bind(&AbsConstraint::setGammaS,_1,1));

    cout << "subproblem:\t" << vars.size()
         << "\t" << ctrs.size() << endl;


    stat.restart();
    feasible = colorConstructor(vars.begin(),vars.end(),ctrs.begin(),ctrs.end());
    stat.setTime(CateConstruct);


    if(feasible){
        cout << "feasible" << endl;
        exit(0);
    }


    // =========================================== TIME 3 =========================== BGN
    vars.erase(std::remove_if(vars.begin(), vars.end(), NonMarked<VarPtr>()),vars.end());
    ctrs.erase(std::remove_if(ctrs.begin(),ctrs.end(), NonMarked<CtrPtr>()),ctrs.end());

    for(std::vector<VarPtr>::iterator current = vars.begin(); current != vars.end();
        ++current){
        ((*current)->ctrs_).erase(
                    remove_if((*current)->ctrs_.begin(), (*current)->ctrs_.end(),
                              NonMarked<CtrPtr>() ),(*current)->ctrs_.end() );

        ((*current)->vars_).erase(
                    remove_if((*current)->vars_.begin(), (*current)->vars_.end(),
                              NonMarked<VarPtr>() ),(*current)->vars_.end() );
    }

    BOOST_FOREACH(CtrPtr const c, ctrs){
        c->setGammaS(1);
        c->setLevel(normal);
        c->setVarsLevel(normal);
    }

    stat.restart();
    feasible = colorLocator(vars.begin(),vars.end(),ctrs.begin(),ctrs.end(),
                                 Normal<VarPtr>(), Normal<CtrPtr>());
    stat.setTime(CateLocate);


    if(feasible){
        cout << "feasible" << endl;
        exit(0);
    }

    stat.restart();
    feasible = colorConstructor(vars.begin(),vars.end(),ctrs.begin(),ctrs.end());
    stat.setTime(CateConstruct);


    if(feasible){
        cout << "feasible" << endl;
        exit(0);
    }


    // =========================================== TIME 3 =========================== BGN



    int nbVars = std::accumulate(vars.begin(),vars.end(),0,Counter_Marked<VarPtr>());
    int nbCtrs = std::accumulate(ctrs.begin(),ctrs.end(),0,Counter_Marked<CtrPtr>());


    cout << "subgraph:\t" << nbVars << "\t" << nbCtrs << endl;
    stat.print();
    return true;
}



#endif // COLOR_HPP
