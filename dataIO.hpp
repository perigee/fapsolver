#ifndef DATAIO_HPP
#define DATAIO_HPP

#include <iostream>
#include <fstream>
#include <string>
#include <deque>
#include <map>
#include <exception>
#include <stdexcept>
#include <ostream>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/algorithm/string/trim.hpp>

//#include "statistics.hpp"
//#include "wave.hpp"


//#include "heuristic.hpp"
#include "color.hpp"

using namespace std;


typedef std::pair<int, int> valuePair;

// implement GCP of DIMACS benchmark
template <typename SolverT>
struct FAPP : public SolverT{

    // read benchmark file
    void read(char* filename, int relax = 0){

        ifstream infile;
        infile.open(filename); // open file
        if(infile) {
            string s="";
            //std::deque<std::string> lineDeque;

            int ctrID = 0;

            bool samePreviousState = false;
            bool hasCE = false;
            ValPtr v;
            VarPtr var;
            CtrPtr ctr;
            int ctrId = 0;
            std::string firstStrE;
            std::string secondStrE;
            std::string firstStrI;
            std::string secondStrI;
            Constraint2CPPtr polarCtr;
            while(getline(infile, s)) {


                typedef boost::tokenizer<boost::char_separator<char> >
                        tokenizer;
                boost::char_separator<char> sep(" ");
                tokenizer tok(s, sep);

                //cout << s << endl;
                //for(tokenizer::iterator xbeg = tokx.begin(); xbeg != tokx.end();++xbeg)
                //    cout << *xbeg;

                //cout << endl;

                //exit(0);

                //boost::tokenizer<>  tok(s);

                // create domain list
                if(*(tok.begin()) == "DM"){

                    tokenizer::iterator beg = tok.begin();
                    ++beg;
                    int domainId = std::atoi((*beg).c_str());
                    ++beg;
                    int value = std::atoi((*beg).c_str());

                    //v.reset(new Value(value));
                    std::pair<int, int> p;
                    p.first = value;
                    p.second = domainId;
                    doms_.push_back(p);
                }


                // variables' list
                if(*(tok.begin()) == "TR"){
                    //cout << s << endl;
                    tokenizer::iterator beg = tok.begin();
                    ++beg;
                    string varId = *beg;
                    //cout << *beg << endl;
                    ++beg;
                    int domainId = std::atoi((*beg).c_str());
                    //cout << *beg << endl;
                    ++beg;
                    int polarInt = std::atoi((*beg).c_str());
                    //string polarStr = *beg;
                    //cout << *beg << endl;
                    var.reset(new Variable(varId));


                    BOOST_FOREACH(valuePair const p, doms_){

                        if(p.second != domainId)
                            continue;

                        if(polarInt == 0){
                            //if(polarStr == "0"){
                            v.reset(new Value(p.first));
                            v->setPolar(H);
                            var->values_.push_back(v);
                            v.reset(new Value(p.first));
                            v->setPolar(V);
                            var->values_.push_back(v);
                        }

                        if (polarInt > 0){
                            //if(polarStr == "1"){
                            //cout << "H polarity only: " << var->id_ << endl;
                            v.reset(new Value(p.first));
                            v->setPolar(H);
                            var->values_.push_back(v);
                        }

                        if (polarInt < 0){
                            //if(polarStr == "-1"){
                            //cout << "V polarity only: " << var->id_ << endl;
                            v.reset(new Value(p.first));
                            v->setPolar(V);
                            var->values_.push_back(v);
                        }


                    }

                    vars_.push_back(var);
                }


                // constraints list
                if(*(tok.begin()) == "CI"){
                    tokenizer::iterator beg = tok.begin();
                    std::string idStr = boost::lexical_cast<std::string>(static_cast<int>(ctrID));
                    Constraint2CIPtr ctr2(new Constraint2CI(idStr));
                    ++beg;
                    string firstId = *beg;
                    ++beg;
                    string secondId = *beg;
                    ++beg;
                    std::string freq = *beg;
                    ++beg;
                    std::string equ = *beg;
                    ++beg;
                    int interval = std::atoi((*beg).c_str());

                    if(freq == "F"){
                        if(equ == "E"){
                            ctr2->init(true,true,interval);
                        }else{
                            ctr2->init(true,false,interval);
                        }
                    }

                    if(freq == "P"){
                        if(equ == "E"){
                            ctr2->init(false,true);
                        }else{
                            ctr2->init(false,false);
                        }

                    }

                    ctr2->setVars(vars_.begin(),vars_.end(),
                                  firstId,secondId,ctr2);
                    ctrs_.push_back(ctr2);
                    ++ctrID;
                }


                if(*(tok.begin()) == "CE"){
                    std::string idStr = boost::lexical_cast<std::string>(static_cast<int>(ctrID));
                    polarCtr.reset(new Constraint2CP(idStr));
                    tokenizer::iterator beg = tok.begin();
                    ++beg;
                    firstStrE = *beg;
                    ++beg;
                    secondStrE = *beg;
                    ++beg;

                    for(;beg != tok.end();++beg){
                        polarCtr->push_back(true,std::atoi((*beg).c_str()));
                    }

                    ///ctr.reset();
                }

                if(*(tok.begin()) == "CD"){
                    tokenizer::iterator beg = tok.begin();
                    ++beg;
                    firstStrI = *beg;
                    ++beg;
                    secondStrI = *beg;
                    ++beg;

                    assert(firstStrE== firstStrI);
                    assert(secondStrE== secondStrI);

                    for(;beg != tok.end();++beg){
                        polarCtr->push_back(false,std::atoi((*beg).c_str()));
                    }

                    polarCtr->init(relax);


                    polarCtr->setVars(vars_.begin(),vars_.end(),
                                      firstStrI,secondStrI,polarCtr);
                    ctrs_.push_back(polarCtr);
                    ++ctrID;
                }

            }
        }

        infile.close(); // file close

    }

    void solve(){

        // since gamma is set to 1 by default, just clear all gamma in values


        SolverT::solve(vars_, ctrs_);

        //SolverT::testSimpleWave(vars_,ctrs_);
    }

    template <typename InputIterator>
    void verify(InputIterator first, InputIterator last){

        int totalCtr = 0;
        int violCtr = 0;
        for(;first != last; ++first){
            ++totalCtr;
            if((*first)->isAssigned()){
                if((*first)->isViolated())
                    ++violCtr;
            }else{
                cout << "In problem there is no assigned variables" << endl;
                exit(1);
            }
        }

        cout << "total     ctrs: " << totalCtr << endl;
        cout << "violated  ctrs: " << violCtr << endl;
    }

private:
    std::vector<VarPtr> vars_;
    std::vector<CtrPtr> ctrs_;
    std::vector<valuePair> doms_;

};




// implement GCP of DIMACS benchmark
template <typename SolverT>
struct DIMACS_coloring : public SolverT{

    // read benchmark file
    void read(char* filename, int nbCol){

        ifstream infile;
        infile.open(filename); // open file
        if(infile) {
            string s="";
            std::deque<std::string> lineDeque;

            int ctrID = 0;

            while(getline(infile, s)) {
                //cout<<s<<endl;
                boost::tokenizer<> tok(s);
                if(!lineDeque.empty())
                    lineDeque.clear();

                // create variables list
                if(*(tok.begin()) == "p"){

                    std::copy(tok.begin(),tok.end(),back_inserter(lineDeque));

                    //cout << lineDeque.back() << "\n";
                    ctrs_.reserve(std::atoi(lineDeque.back().c_str()));
                    lineDeque.pop_back();
                    int nbVars = std::atoi(lineDeque.back().c_str());
                    vars_.reserve(nbVars);

                    VarPtr var;

                    //for(int i=0;i< nbVars; ++i){
                    for(int i=1;i< nbVars+1 ; ++i){
                        //std::string id = "V" + boost::lexical_cast<std::string>(static_cast<int>(i));
                        std::string id = boost::lexical_cast<std::string>(static_cast<int>(i));
                        //cout << id << endl;
                        var.reset(new Variable(id));
                        // set domains

                        ValPtr v;
                        for(int idCol=0;idCol<nbCol;++idCol){
                            v.reset(new Value(idCol));
                            (var->values_).push_back(v);
                        }


                        vars_.push_back(var);
                    }

                }

                // create constraints list
                if(*(tok.begin()) == "e"){
                    std::copy(tok.begin(),tok.end(),back_inserter(lineDeque));

                    std::string varb = lineDeque.back();
                    lineDeque.pop_back();
                    std::string vara = lineDeque.back();

                    // verify the duplicated constraints,

                    bool duplicated = false;
                    BOOST_FOREACH(CtrPtr const c, ctrs_){
                        int cnt=0;
                        if(c->isInStr(vara))
                            ++ cnt;
                        if(c->isInStr(varb))
                            ++ cnt;

                        if(cnt ==2){
                            //cout << "find duplicated" << endl;
                            duplicated = true;
                            break;
                        }
                    }

                    if(duplicated)
                        continue;





                    Constraint2Ptr ctr2(new Constraint2(boost::lexical_cast<std::string>(ctrID)));
                    ++ctrID;
                    ctr2->setVars(vars_.begin(),vars_.end(),vara,varb,ctr2);
                    ctrs_.push_back(ctr2);
                }
            }
        }

        infile.close(); // file close

        //cout << "vars: " << vars_.size() << endl;
        //setDomain(vars_.begin(), vars_.end(), nbCol);

        //cout << "ctrs: " << ctrs_.size() << endl;
    }

    void solve(){

        // since gamma is set to 1 by default, just clear all gamma in values


        SolverT::solve(vars_, ctrs_);

        //SolverT::testSimpleWave(vars_,ctrs_);
    }

    template <typename InputIterator>
    void verify(InputIterator first, InputIterator last){

        int totalCtr = 0;
        int violCtr = 0;
        for(;first != last; ++first){
            ++totalCtr;
            if((*first)->isAssigned()){
                if((*first)->isViolated())
                    ++violCtr;
            }else{
                cout << "In problem there is no assigned variables" << endl;
                exit(1);
            }
        }

        cout << "total     ctrs: " << totalCtr << endl;
        cout << "violated  ctrs: " << violCtr << endl;
    }

private:
    std::vector<VarPtr> vars_;
    std::vector<CtrPtr> ctrs_;

};

//====================================================================

typedef std::vector<ValPtr> Domain;

template <typename SolverT>
struct LRFAP_solver : public SolverT{

public:
    LRFAP_solver(){}
    ~LRFAP_solver(){}




    void solve(){
        SolverT::solve(vars_, ctrs_);
        //SolverT::testSimpleWave(vars_,ctrs_);
    }




    template <typename InputIterator>
    void verify(InputIterator first, InputIterator last){
        int total = 0;
        int violated = 0;
        for(;first != last; ++first){
            ++total;
            if(!(*first)->isAssigned())
                throw std::runtime_error("Error: has non assigned variables");
            if((*first)->isViolated()){
                ++violated;
            }

            //(*first)->print();
        }

        // total ctrs + violated ctrs
        cout << " " << total << " " << violated << " ";
    }


    /*
   * Read Domain File
   *
   */
    void readDomFile(char* filename){

        //std::cout << "Reading Domain File..." << std::endl;

        ifstream fin;
        fin.open(filename);


        // File reading error
        if(fin.good()==false){
            cerr << "Error: Domain file could not be opened" << endl;
            exit(1);
        }


        string id; // id of domain



        std::string lineString;
        while(getline(fin,lineString)){
            boost::algorithm::trim(lineString);
            if( lineString == "")
                continue;

            boost::tokenizer<> tok(lineString);

            Domain dom;
            ValPtr val;

            id = *(tok.begin());
            int cnt = 0;
            for(boost::tokenizer<>::iterator beg=tok.begin(); beg!=tok.end();++beg){
                ++cnt;
                if(cnt > 2){
                    val.reset(new Value(boost::lexical_cast<int>(*beg)));
                    dom.push_back(val);
                    //cout << (dom.back())->v_ << "\t";
                }
            }
            //cout << endl;

            doms[id] = dom;
        }

        //exit(0);

        fin.close();
    }


    // Reading variables' file
    void readVarFile(char* filename){
        ifstream fin;
        fin.open(filename);

        //std::cout << "Reading Variable File..." << std::endl;

        // File reading error
        if(fin.good()==false){
            cerr << "Error: Variable file could not be opened" << endl;
            exit(1);
        }

        string id; // id of variable
        string idDom; // id of domain of variable

        VarPtr varPtr;

        std::string lineString;
        while(getline(fin,lineString)){

            boost::algorithm::trim(lineString);
            if( lineString == "")
                continue;

            boost::tokenizer<> tok(lineString);

            id = *(tok.begin());

            varPtr.reset(new Variable(id));

            int cnt = 0;
            for(boost::tokenizer<>::iterator beg=tok.begin(); beg!=tok.end();++beg){
                ++cnt;
                if(cnt == 2){

                    Domain dom = doms[*beg];
                    BOOST_FOREACH(ValPtr const a, dom){
                        ValPtr val(new Value(a->v_));
                        varPtr->values_.push_back(val);
                    }

                    //cout << *varPtr << endl;

                    vars_.push_back(varPtr);
                    break;
                }
            }
            //cout << endl;


        }

        //exit(0);

        fin.close();

    }

    void readCtrFile(char* filename){
        ifstream fin;
        fin.open(filename);

        //std::cout << "Reading Constraint File..." << std::endl;

        // File reading error
        if(fin.good()==false){
            cerr << "Error: Constraint file could not be opened" << endl;
            exit(1);
        }


        string firstId; // id of constraint
        string secondId; // id of domain of variable
        char typeCtr;
        char oprSign;
        double diff;
        //int weight;
        // for max-csp instances

        int id = 0; //  id of constraint


        while(fin >> firstId >> secondId >> typeCtr >> oprSign >> diff ){
            //while(fin >> firstId >> secondId >> typeCtr >> oprSign >> diff >> weight){

            ++id;


            std::string idStr = boost::lexical_cast<std::string>(static_cast<int>(id));

            bool duplicated = false;
            BOOST_FOREACH(CtrPtr const c, ctrs_){
                int cnt=0;
                if(c->isInStr(firstId))
                    ++ cnt;
                if(c->isInStr(secondId))
                    ++ cnt;

                if(cnt ==2){
                    cout << "find duplicated" << endl;
                    duplicated = true;
                    break;
                }
            }


            Constraint2FAPPtr ctr(new Constraint2FAP(idStr));
            ctr->setParams(oprSign, diff);
            ctr->setVars(vars_.begin(),vars_.end(),firstId,secondId,ctr);


            ctrs_.push_back(ctr);

            //(ctrs_.back())->print();


        }

        //exit(0);


        fin.close();
    }



    std::map<string,Domain> doms;
    std::vector<VarPtr> vars_;
    std::vector<CtrPtr> ctrs_;


};




struct Solver{


    //    void test_Mac(std::vector<VarPtr>& vars, std::vector<CtrPtr>& ctrs){


    //        //std::for_each(vars.begin(),vars.end(), boost::bind(&Variable::resetGamma,_1));
    //        std::for_each(vars.begin(),vars.end(), boost::bind(&Variable::mark,_1));

    //        std::for_each(ctrs.begin(),ctrs.end(), boost::bind(&AbsConstraint::markByVars,_1));

    //        boost::progress_timer t;

    //        bool feasible = mac_solver(vars.begin(), vars.end(),t );
    //        if(feasible){
    //            cout << " feasible " ;

    //        }else
    //            cout << " infeasible " ;

    //        cout << endl;

    //        std::for_each(vars.begin(), vars.end(), boost::bind(&Variable::print,_1));

    //    }



    //    void tabuTest(std::vector<VarPtr>& vars_, std::vector<CtrPtr>& ctrs_){
    //        boost::progress_timer t;


    //        //==============================================================================

    //        cout << "problem: " << vars_.size() << " " << ctrs_.size() << endl;

    //        srand((unsigned)time(0));


    //        //cout << "in random solver" << endl;
    //        randomSolver(vars_.begin(), vars_.end());
    //        std::for_each(vars_.begin(), vars_.end(), boost::bind(&Variable::setGamma, _1, 1));
    //        std::for_each(ctrs_.begin(), ctrs_.end(), boost::bind(&AbsConstraint::setGamma, _1, 1));


    //        bool feasible = tabu(vars_.begin(), vars_.end(), is_gamma_variable());

    //        // test breakout
    //        //bool feasible = breakout(vars_.begin(), vars_.end(), is_gamma_variable());
    //        //std::for_each(ctrs_.begin(), ctrs_.end(), boost::bind(&AbsConstraint::resetGamma, _1));

    //        if(feasible){
    //            cout << " feasible " ;
    //        }else{



    //            //cout << " infeasible" ;


    //            cout << "r:\t"
    //                    << std::accumulate( vars_.begin(), vars_.end(), 0, Counter_violated<VarPtr>() )
    //                    << "\t"
    //                    <<std::accumulate( ctrs_.begin(), ctrs_.end(), 0, Counter_violated<CtrPtr>() )
    //                    << "\t";

    //        }


    ////        int w = 0;
    ////        BOOST_FOREACH(VarPtr const a, vars_){
    ////            //cout << a->id_ << " " << a->getWeight() << endl;
    ////            if(a->getWeight() > w)
    ////                w = a->getWeight();
    ////        }


    ////        cout << "maximal weight: " << w << endl;

    //        //cout << endl;

    //    }


    //    void grow(std::vector<VarPtr>& vars_, std::vector<CtrPtr>& ctrs_){
    //        boost::progress_timer t;


    //        //==============================================================================

    //        srand((unsigned)time(0));

    //        randomSolver(vars_.begin(), vars_.end());


    //        OneViolatedVars violatedVars;
    //        //AllViolatedVars violatedVars;
    //        //LinkToGammaVars gammaVars;
    //        LinkToWeightVars gammaVars;

    //        bool feasible = grow_solver(vars_,ctrs_,violatedVars,gammaVars);


    //        if(feasible){

    //            cout << " feasible " ;
    //        }else{
    //            std::vector<VarPtr> core;
    //            BOOST_FOREACH(VarPtr const v, vars_){
    //                if(v->getGamma() >0){
    //                    v->mark();
    //                }
    //            }

    //            std::remove_copy_if(vars_.begin(), vars_.end(),std::back_inserter(core), not_marked_element<VarPtr>());
    //            std::for_each(ctrs_.begin(), ctrs_.end(), boost::bind(&AbsConstraint::regardMark,_1));

    //            cout << "core size: " << core.size() << endl;

    //            feasible = mac_solver(core.begin(), core.end(), t);
    //            if(feasible)
    //                cout << " feasible by MAC" ;
    //            else
    //                cout << " infeasible " ;
    //        }


    //    }



    //    void extend(std::vector<VarPtr>& vars_, std::vector<CtrPtr>& ctrs_){
    //        boost::progress_timer t;


    //        //==============================================================================

    //        srand((unsigned)time(0));

    //        randomSolver(vars_.begin(), vars_.end());


    //        OneViolatedVarsWithExtension violatedVars;
    //        //AllViolatedVarsWithExtension violatedVars;
    //        Extension extension;

    //        bool feasible = extension_mark_solver(vars_,ctrs_,violatedVars,extension);


    //        if(feasible){

    //            cout << " feasible " ;
    //        }else{
    //            std::vector<VarPtr> core;



    //            std::remove_copy_if(vars_.begin(), vars_.end(),std::back_inserter(core), not_marked_element<VarPtr>());
    //            std::for_each(ctrs_.begin(), ctrs_.end(), boost::bind(&AbsConstraint::markByVars,_1));

    //            cout << "core size: " << core.size() << endl;

    //            feasible = mac_solver(core.begin(), core.end(),t);
    //            if(feasible)
    //                cout << " feasible by MAC" ;
    //            else
    //                cout << " infeasible " ;
    //        }


    //    }


    //    // arc-consistency verification
    //    void ac(std::vector<VarPtr>& vars, std::vector<CtrPtr>& ctrs){

    //        std::for_each(vars.begin(),vars.end(), boost::bind(&Variable::mark,_1));

    //        std::for_each(ctrs.begin(),ctrs.end(), boost::bind(&AbsConstraint::markByVars,_1));

    //        //printCore(vars,ctrs);

    //        cout << "r: ";
    //        if(arc_consistency(vars.begin(), vars.end()))
    //            cout << "arcconsistent" << endl;
    //        else
    //            cout << "deadend" << endl;
    //        //printCore(vars,ctrs);
    //    }



    //    void extend_Guarrentee(std::vector<VarPtr>& vars_, std::vector<CtrPtr>& ctrs_){
    //        boost::progress_timer t;
    //        boost::progress_timer tcore;
    //        double totalH = 0.0;
    //        double totalE = 0.0;
    //        int count = 0;



    //        //==============================================================================
    //        cout << "b:\tExecution"  << endl;
    //        cout << "p:\tprob\t"
    //                << vars_.size()
    //                << "\t"
    //                << ctrs_.size()
    //                ;//<<  endl;

    //        srand((unsigned)time(0));

    //        randomSolver(vars_.begin(), vars_.end());


    //        OneViolatedVarsWithExtension violatedVars;
    //        //AllViolatedVarsWithExtension violatedVars;
    //        Extension extension;


    //        bool feasible = true;
    //        std::vector<VarPtr> core;
    //        while(feasible){
    //            ++count;

    //            tcore.restart();
    //            feasible = extension_mark_solver(vars_,ctrs_,violatedVars,extension);
    //            totalH += tcore.elapsed();

    //            if(feasible){
    //                cout << "feasible by tabu" << endl;
    //                return;
    //            }

    //            cout << "\tsubp\t"
    //                    << std::accumulate( vars_.begin(), vars_.end(), 0, Counter_gamma<VarPtr>() )
    //                    << "\t"
    //                    <<std::accumulate( ctrs_.begin(), ctrs_.end(), 0, Counter_gamma<CtrPtr>() )
    //                    ;//<<  endl;



    //            if(!core.empty())
    //                core.clear();


    //            std::remove_copy_if(vars_.begin(), vars_.end(),std::back_inserter(core), not_marked_element<VarPtr>());
    //            std::for_each(ctrs_.begin(), ctrs_.end(), boost::bind(&AbsConstraint::markByVars,_1));

    //            cout << "\tcore\t"
    //                    << std::accumulate( vars_.begin(), vars_.end(), 0, Counter_marked<VarPtr>() )
    //                    << "\t"
    //                    <<std::accumulate( ctrs_.begin(), ctrs_.end(), 0, Counter_marked<CtrPtr>() )
    //                    ;//<<  endl;





    //            tcore.restart();
    //            feasible = mac_solver(core.begin(), core.end(),t);
    //            totalE += tcore.elapsed();

    //            if(feasible){

    //                //                BOOST_FOREACH(VarPtr const vv, core){
    //                //                    if(!vv->isAssigned()){
    //                //                        cout << "Error: mac problem" << endl;
    //                //                        exit(0);
    //                //                    }
    //                //
    //                //                }


    //                std::for_each(ctrs_.begin(), ctrs_.end(), resetCtr);
    //                std::for_each(vars_.begin(), vars_.end(), resetVar);



    //            }else{
    //                cout << "\nr:\tcore\t"
    //                        << std::accumulate( vars_.begin(), vars_.end(), 0, Counter_marked<VarPtr>() )
    //                        << "\t"
    //                        <<std::accumulate( ctrs_.begin(), ctrs_.end(), 0, Counter_marked<CtrPtr>() )
    //                        ;//<<  endl;
    //                cout << "\t"
    //                        << count
    //                        << "\t" << totalH
    //                        << "\t" << totalE
    //                        << "\t" << t.elapsed()
    //                        << " \tinfeasible "  ;

    //                //printCore(vars_, ctrs_);
    //            }
    //        }


    //    }



    //    bool extract(std::vector<VarPtr>& vars_, std::vector<CtrPtr>& ctrs_, boost::progress_timer& t){
    //        //boost::progress_timer t;
    ////        boost::progress_timer tcore;
    ////        double totalH = 0.0;
    ////        double totalE = 0.0;
    ////        int count = 0;



    //        //==============================================================================



    //        randomSolver(vars_.begin(), vars_.end());


    //        OneViolatedVarsWithExtension violatedVars;
    //        //AllViolatedVarsWithExtension violatedVars;
    //        Extension extension;


    //        bool feasible = true;
    //        std::vector<VarPtr> core;
    //        while(feasible){
    //            //++count;

    //            //tcore.restart();
    //            feasible = extension_mark_solver(vars_,ctrs_,violatedVars,extension);
    //            //totalH += tcore.elapsed();

    //            if(feasible){
    //                cout << "feasible by tabu" << endl;
    //                return true;
    //            }

    ////            cout << "\tsubp\t"
    ////                    << std::accumulate( vars_.begin(), vars_.end(), 0, Counter_gamma<VarPtr>() )
    ////                    << "\t"
    ////                    <<std::accumulate( ctrs_.begin(), ctrs_.end(), 0, Counter_gamma<CtrPtr>() )
    ////                    ;//<<  endl;



    //            if(!core.empty())
    //                core.clear();


    //            std::remove_copy_if(vars_.begin(), vars_.end(),std::back_inserter(core), not_marked_element<VarPtr>());
    //            std::for_each(ctrs_.begin(), ctrs_.end(), boost::bind(&AbsConstraint::markByVars,_1));

    ////            cout << "\tcore\t"
    ////                    << std::accumulate( vars_.begin(), vars_.end(), 0, Counter_marked<VarPtr>() )
    ////                    << "\t"
    ////                    <<std::accumulate( ctrs_.begin(), ctrs_.end(), 0, Counter_marked<CtrPtr>() )
    ////                    ;//<<  endl;





    //            //tcore.restart();
    //            feasible = mac_solver(core.begin(), core.end(),t);
    //            //totalE += tcore.elapsed();

    //            if(feasible){

    //                //                BOOST_FOREACH(VarPtr const vv, core){
    //                //                    if(!vv->isAssigned()){
    //                //                        cout << "Error: mac problem" << endl;
    //                //                        exit(0);
    //                //                    }
    //                //
    //                //                }


    //                std::for_each(ctrs_.begin(), ctrs_.end(), resetCtr);
    //                std::for_each(vars_.begin(), vars_.end(), resetVar);



    //            }else{
    //                return false;
    ////                cout << "\nr:\tcore\t"
    ////                        << std::accumulate( vars_.begin(), vars_.end(), 0, Counter_marked<VarPtr>() )
    ////                        << "\t"
    ////                        <<std::accumulate( ctrs_.begin(), ctrs_.end(), 0, Counter_marked<CtrPtr>() )
    ////                        ;//<<  endl;
    ////                cout << "\t"
    ////                        << count
    ////                        << "\t" << totalH
    ////                        << "\t" << totalE
    ////                        << "\t" << t.elapsed()
    ////                        << " \tinfeasible "  ;

    //                //printCore(vars_, ctrs_);
    //            }
    //        }


    //    }






    //    bool breakout_procedure(std::vector<VarPtr>& vars_, std::vector<CtrPtr>& ctrs_, boost::progress_timer& t){
    //        //boost::progress_timer t;
    ////        boost::progress_timer tcore;
    ////        double totalH = 0.0;
    ////        double totalE = 0.0;
    ////        int count = 0;



    //        //==============================================================================



    //        randomSolver(vars_.begin(), vars_.end());


    //        OneViolatedVarsWithExtension violatedVars;
    //        //AllViolatedVarsWithExtension violatedVars;

    //        Extension extension;
    //        //ReachInfeasible extension;


    //        bool feasible = true;
    //        std::vector<VarPtr> core;
    //        while(feasible){
    //            //++count;

    //            //tcore.restart();
    //            feasible = breakout_mark_solver(vars_,ctrs_,violatedVars,extension);
    //            //feasible = breakout_heuristic_mark_solver(vars_,ctrs_,violatedVars,extension);
    //            //totalH += tcore.elapsed();

    //            if(feasible){
    //                cout << "feasible by breakout" << endl;
    //                return true;
    //            }

    ////            cout << "\tsubp\t"
    ////                    << std::accumulate( vars_.begin(), vars_.end(), 0, Counter_gamma<VarPtr>() )
    ////                    << "\t"
    ////                    <<std::accumulate( ctrs_.begin(), ctrs_.end(), 0, Counter_gamma<CtrPtr>() )
    ////                    ;//<<  endl;



    //            if(!core.empty())
    //                core.clear();


    //            std::remove_copy_if(vars_.begin(), vars_.end(),std::back_inserter(core), not_marked_element<VarPtr>());
    //            std::for_each(ctrs_.begin(), ctrs_.end(), boost::bind(&AbsConstraint::markByVars,_1));

    ////            cout << "\tcore\t"
    ////                    << std::accumulate( vars_.begin(), vars_.end(), 0, Counter_marked<VarPtr>() )
    ////                    << "\t"
    ////                    <<std::accumulate( ctrs_.begin(), ctrs_.end(), 0, Counter_marked<CtrPtr>() )
    ////                    ;//<<  endl;





    //            //tcore.restart();
    //            feasible = mac_solver(core.begin(), core.end(),t);
    //            //totalE += tcore.elapsed();

    //            if(feasible){

    //                cout << "exact found feasible" << endl;
    //                exit(0);

    //                //                BOOST_FOREACH(VarPtr const vv, core){
    //                //                    if(!vv->isAssigned()){
    //                //                        cout << "Error: mac problem" << endl;
    //                //                        exit(0);
    //                //                    }
    //                //
    //                //                }


    //                std::for_each(ctrs_.begin(), ctrs_.end(), resetCtr);
    //                std::for_each(vars_.begin(), vars_.end(), resetVar);



    //            }else{
    //                return false;
    ////                cout << "\nr:\tcore\t"
    ////                        << std::accumulate( vars_.begin(), vars_.end(), 0, Counter_marked<VarPtr>() )
    ////                        << "\t"
    ////                        <<std::accumulate( ctrs_.begin(), ctrs_.end(), 0, Counter_marked<CtrPtr>() )
    ////                        ;//<<  endl;
    ////                cout << "\t"
    ////                        << count
    ////                        << "\t" << totalH
    ////                        << "\t" << totalE
    ////                        << "\t" << t.elapsed()
    ////                        << " \tinfeasible "  ;

    //                //printCore(vars_, ctrs_);
    //            }
    //        }


    //    }





    //    bool breakout_heuristicReach_procedure(std::vector<VarPtr>& vars_, std::vector<CtrPtr>& ctrs_, boost::progress_timer& t){
    //        //boost::progress_timer t;
    ////        boost::progress_timer tcore;
    ////        double totalH = 0.0;
    ////        double totalE = 0.0;
    ////        int count = 0;



    //        //==============================================================================



    //        randomSolver(vars_.begin(), vars_.end());


    //        OneViolatedVarsWithExtension violatedVars;
    //        //AllViolatedVarsWithExtension violatedVars;

    //        //Extension extension;
    //        ReachInfeasible extension;


    //        bool feasible = true;
    //        std::vector<VarPtr> core;
    //        while(feasible){
    //            //++count;

    //            //tcore.restart();
    //            //feasible = breakout_mark_solver(vars_,ctrs_,violatedVars,extension);
    //            feasible = breakout_heuristic_mark_solver(vars_,ctrs_,violatedVars,extension);
    //            //totalH += tcore.elapsed();

    //            if(feasible){
    //                cout << "feasible by breakout" << endl;
    //                return true;
    //            }


    //            /// only test for coloring  ============= BGN
    //            return false;
    //            /// only test for coloring  ============= END


    ////            cout << "\tsubp\t"
    ////                    << std::accumulate( vars_.begin(), vars_.end(), 0, Counter_gamma<VarPtr>() )
    ////                    << "\t"
    ////                    <<std::accumulate( ctrs_.begin(), ctrs_.end(), 0, Counter_gamma<CtrPtr>() )
    ////                    ;//<<  endl;



    //            if(!core.empty())
    //                core.clear();


    //            std::remove_copy_if(vars_.begin(), vars_.end(),std::back_inserter(core), not_marked_element<VarPtr>());
    //            std::for_each(ctrs_.begin(), ctrs_.end(), boost::bind(&AbsConstraint::markByVars,_1));

    ////            cout << "\tcore\t"
    ////                    << std::accumulate( vars_.begin(), vars_.end(), 0, Counter_marked<VarPtr>() )
    ////                    << "\t"
    ////                    <<std::accumulate( ctrs_.begin(), ctrs_.end(), 0, Counter_marked<CtrPtr>() )
    ////                    ;//<<  endl;





    //            //tcore.restart();
    //            feasible = mac_solver(core.begin(), core.end(),t);
    //            //totalE += tcore.elapsed();

    //            if(feasible){

    //                //                BOOST_FOREACH(VarPtr const vv, core){
    //                //                    if(!vv->isAssigned()){
    //                //                        cout << "Error: mac problem" << endl;
    //                //                        exit(0);
    //                //                    }
    //                //
    //                //                }


    //                std::for_each(ctrs_.begin(), ctrs_.end(), resetCtr);
    //                std::for_each(vars_.begin(), vars_.end(), resetVar);



    //            }else{
    //                return false;
    ////                cout << "\nr:\tcore\t"
    ////                        << std::accumulate( vars_.begin(), vars_.end(), 0, Counter_marked<VarPtr>() )
    ////                        << "\t"
    ////                        <<std::accumulate( ctrs_.begin(), ctrs_.end(), 0, Counter_marked<CtrPtr>() )
    ////                        ;//<<  endl;
    ////                cout << "\t"
    ////                        << count
    ////                        << "\t" << totalH
    ////                        << "\t" << totalE
    ////                        << "\t" << t.elapsed()
    ////                        << " \tinfeasible "  ;

    //                //printCore(vars_, ctrs_);
    //            }
    //        }


    //    }






    //    void extract_routine(std::vector<VarPtr>& vars_, std::vector<CtrPtr>& ctrs_){

    //        boost::progress_timer t;
    //        srand((unsigned)time(0));

    //        cout << "b:\tExecution"  << endl;
    //        cout << "p:\tprob\t"
    //                << vars_.size()
    //                << "\t"
    //                << ctrs_.size()
    //                <<  endl;



    //        bool feasible = false;

    //        int count = 0;
    //        while(!feasible){
    //            ++count;
    //            feasible = extract(vars_,ctrs_,t);

    //            int size = vars_.size();

    //            ctrs_.erase(remove_if(ctrs_.begin(), ctrs_.end(), Non_Marked<CtrPtr>()), ctrs_.end() );
    //            vars_.erase(remove_if(vars_.begin(), vars_.end(), Non_Marked<VarPtr>()), vars_.end() );

    //                    BOOST_FOREACH(VarPtr const var, vars_){
    //                        (var->ctrs_).erase(remove_if(var->ctrs_.begin(), var->ctrs_.end(), Non_Marked<CtrPtr>()), var->ctrs_.end() );
    //                        (var->vars_).erase(remove_if(var->vars_.begin(), var->vars_.end(), Non_Marked<VarPtr>()), var->vars_.end() );
    //                        //var->resetDom();
    //                        //var->resign();
    //                    }

    //            if(size == vars_.size())
    //                break;
    //            else{

    //                std::for_each(ctrs_.begin(), ctrs_.end(), resetCtr);
    //                std::for_each(vars_.begin(), vars_.end(), resetVar);

    //            }


    //        }

    //        cout << "\nr:\tcore\t"
    //                << std::accumulate( vars_.begin(), vars_.end(), 0, Counter_marked<VarPtr>() )
    //                << "\t"
    //                <<std::accumulate( ctrs_.begin(), ctrs_.end(), 0, Counter_marked<CtrPtr>() )
    //                ;//<<  endl;
    //        cout << "\t"
    //                << count
    ////                << "\t" << totalH
    ////                << "\t" << totalE
    //                << "\t" << t.elapsed()
    //                << " \tinfeasible "  ;

    //    }





    //    void breakout_routine(std::vector<VarPtr>& vars_, std::vector<CtrPtr>& ctrs_){

    //        boost::progress_timer t;
    //        srand((unsigned)time(0));

    //        cout << "b:\tExecution"  << endl;
    //        cout << "p:\tprob\t"
    //                << vars_.size()
    //                << "\t"
    //                << ctrs_.size()
    //                <<  endl;



    //        bool feasible = false;

    //        int count = 0;
    //        while(!feasible){
    //            ++count;
    //            feasible = breakout_procedure(vars_,ctrs_,t);

    //            int size = vars_.size();

    //            ctrs_.erase(remove_if(ctrs_.begin(), ctrs_.end(), Non_Marked<CtrPtr>()), ctrs_.end() );
    //            vars_.erase(remove_if(vars_.begin(), vars_.end(), Non_Marked<VarPtr>()), vars_.end() );

    //                    BOOST_FOREACH(VarPtr const var, vars_){
    //                        (var->ctrs_).erase(remove_if(var->ctrs_.begin(), var->ctrs_.end(), Non_Marked<CtrPtr>()), var->ctrs_.end() );
    //                        (var->vars_).erase(remove_if(var->vars_.begin(), var->vars_.end(), Non_Marked<VarPtr>()), var->vars_.end() );
    //                        //var->resetDom();
    //                        //var->resign();
    //                    }

    //            if(size == vars_.size())
    //                break;
    //            else{

    //                std::for_each(ctrs_.begin(), ctrs_.end(), resetCtr);
    //                std::for_each(vars_.begin(), vars_.end(), resetVar);

    //            }


    //        }

    //        cout << "\nr:\tcore\t"
    //                << std::accumulate( vars_.begin(), vars_.end(), 0, Counter_marked<VarPtr>() )
    //                << "\t"
    //                <<std::accumulate( ctrs_.begin(), ctrs_.end(), 0, Counter_marked<CtrPtr>() )
    //                ;//<<  endl;
    //        cout << "\t"
    //                << count
    ////                << "\t" << totalH
    ////                << "\t" << totalE
    //                << "\t" << t.elapsed()
    //                << " \tinfeasible "  ;

    //    }





    //    void breakout_heuristicReach_routine(std::vector<VarPtr>& vars_, std::vector<CtrPtr>& ctrs_){

    //        boost::progress_timer t;
    //        srand((unsigned)time(0));

    //        cout << "b:\tExecution"  << endl;
    //        cout << "p:\tprob\t"
    //                << vars_.size()
    //                << "\t"
    //                << ctrs_.size()
    //                <<  endl;



    //        bool feasible = false;

    //        int count = 0;
    //        while(!feasible){
    //            ++count;

    //            if(count < 2)
    //                feasible = breakout_heuristicReach_procedure(vars_,ctrs_,t);
    //            else
    //                feasible = breakout_procedure(vars_,ctrs_,t);

    //            int size = vars_.size();

    //            ctrs_.erase(remove_if(ctrs_.begin(), ctrs_.end(), Non_Marked<CtrPtr>()), ctrs_.end() );
    //            vars_.erase(remove_if(vars_.begin(), vars_.end(), Non_Marked<VarPtr>()), vars_.end() );

    //                    BOOST_FOREACH(VarPtr const var, vars_){
    //                        (var->ctrs_).erase(remove_if(var->ctrs_.begin(), var->ctrs_.end(), Non_Marked<CtrPtr>()), var->ctrs_.end() );
    //                        (var->vars_).erase(remove_if(var->vars_.begin(), var->vars_.end(), Non_Marked<VarPtr>()), var->vars_.end() );
    //                        //var->resetDom();
    //                        //var->resign();
    //                    }

    //            if(size == vars_.size())
    //                break;
    //            else{

    //                std::for_each(ctrs_.begin(), ctrs_.end(), resetCtr);
    //                std::for_each(vars_.begin(), vars_.end(), resetVar);

    //            }


    //        }

    //        cout << "\nr:\tcore\t"
    //                << std::accumulate( vars_.begin(), vars_.end(), 0, Counter_marked<VarPtr>() )
    //                << "\t"
    //                <<std::accumulate( ctrs_.begin(), ctrs_.end(), 0, Counter_marked<CtrPtr>() )
    //                ;//<<  endl;
    //        cout << "\t"
    //                << count
    ////                << "\t" << totalH
    ////                << "\t" << totalE
    //                << "\t" << t.elapsed()
    //                << " \tinfeasible "  ;

    //    }





    //    void extend_conserver(std::vector<VarPtr>& vars_, std::vector<CtrPtr>& ctrs_){
    //        boost::progress_timer t;
    //        boost::progress_timer tcore;
    //        double totalLocator = 0.0;
    //        double totalExtender = 0.0;
    //        double totalVerificator = 0.0;




    //        //==============================================================================
    //        cout << "b:\tExecution"  << endl;
    //        cout << "p:\tprob\t"
    //                << vars_.size()
    //                << "\t"
    //                << ctrs_.size()
    //                ;//<<  endl;

    //        srand((unsigned)time(0));

    //        randomSolver(vars_.begin(), vars_.end());



    //        OneViolatedVarsWithExtension violatedVars;
    //        //AllViolatedVarsWithExtension violatedVars;
    //        Extension extension;


    //        //=========================
    //        // locator
    //        tcore.restart();
    //        bool feasible = locator(vars_, ctrs_,
    //                                violatedVars);
    //        totalLocator += tcore.elapsed();

    //        feasible = true;
    //        int counter = 0;
    //        while(feasible){
    //            ++counter;
    //            tcore.restart();
    //            feasible = extender(vars_, ctrs_,
    //                                violatedVars,extension,t);
    //            totalExtender += tcore.elapsed();

    //            if(!feasible){
    //                tcore.restart();
    //                feasible = verificator(vars_, ctrs_, t);
    //                totalVerificator += tcore.elapsed();

    //                if(feasible){
    //                    std::for_each(vars_.begin(),vars_.end(), enableMarkArround);
    //                }
    //            }

    //        }


    //        if(!feasible){
    //            cout << "\nr:\tcore\t"
    //                    << std::accumulate( vars_.begin(), vars_.end(), 0, Counter_marked<VarPtr>() )
    //                    << "\t"
    //                    <<std::accumulate( ctrs_.begin(), ctrs_.end(), 0, Counter_marked<CtrPtr>() )
    //                    ;//<<  endl;
    //            cout << "\t"
    //                    << counter
    //                    << "\t" << totalLocator
    //                    << "\t" << totalExtender
    //                    << "\t" << totalVerificator
    //                    << "\t" << t.elapsed()
    //                    << " \tinfeasible "  ;
    //        }



    //    }


    template <typename Iter>
    int totalPruneValues(Iter first, Iter last){
        int total = 0;
        for(;first != last;++first){
            total += validDomainSize(*first);
        }

        return total;
    }

    void testAC(std::vector<VarPtr>& vars, std::vector<CtrPtr>& ctrs){
        std::for_each(ctrs.begin(), ctrs.end(),
                      boost::bind(&AbsConstraint::setGammaS,_1,1));
        std::for_each(vars.begin(), vars.end(),
                      boost::bind(&Variable::setGammaS,_1,1));

        bool feasible = false;

        boost::timer timer;
        Acfc acSolver;
        bool foundFeasible = false;
        for(int i = 0; i < 11; ++i){
            //for(int i = 0; i < 1; ++i){

            
            if(!feasible){
                // set constraint relax level
                std::for_each(ctrs.begin(), ctrs.end(),
                              boost::bind(&AbsConstraint::setCtrLevel,_1,i));

                timer.restart();
                bool ifeasible = arc_consistency(vars.begin(),vars.end(),
                                                 Normal<VarPtr>(),Normal<CtrPtr>());

                cout << "d:\t" << timer.elapsed() << "\t" << totalPruneValues(vars.begin(),vars.end()) << "\t";

                std::for_each(vars.begin(), vars.end(),
                              boost::bind(&Variable::resetDomain,_1));


                timer.restart();
                feasible = acSolver.solve(vars.begin(),vars.end());

                cout << timer.elapsed() << "\t" << totalPruneValues(vars.begin(),vars.end()) << endl;


                if (ifeasible != feasible){
                    cout << "ac problem" << endl;
                    exit(0);
                }


                std::for_each(vars.begin(), vars.end(),
                              boost::bind(&Variable::resetDomain,_1));

                if(feasible){
                    foundFeasible = true;
                    cout << "r:\t feasible at level: " << i << endl;
                    break;
                }

            }

        }


        if(!foundFeasible){
            cout << "r:\t feasible at level:  11" <<  endl;
        }

    }



    void testACwave(std::vector<VarPtr>& vars, std::vector<CtrPtr>& ctrs){
        std::for_each(ctrs.begin(), ctrs.end(),
                      boost::bind(&AbsConstraint::setGammaS,_1,1));
        std::for_each(vars.begin(), vars.end(),
                      boost::bind(&Variable::setGammaS,_1,1));


        Acfc acSolver;
        bool feasible = false;
        //for(int i = 0; i < 11; ++i){
        for(int i = 0; i < 1; ++i){
            
            if(!feasible){
                //std::for_each(ctrs.begin(), ctrs.end(),
                //boost::bind(&AbsConstraint::setCtrLevel,_1,i));

                feasible = acSolver.solve(vars.begin(),vars.end());
                std::for_each(vars.begin(), vars.end(),
                              boost::bind(&Variable::resetDomain,_1));
            }
	    

            if(feasible){
                cout << "\t OK" ;
            }else{
                cout << "\t INC" ;
            }

            

        }

        cout << endl;
    }

    template <typename Itr>
    void printCtrs(Itr first, Itr last){
        for(; first != last; ++first)
            (*first)->print();
    }


    void testTabu(std::vector<VarPtr>& vars, std::vector<CtrPtr>& ctrs){

        std::for_each(ctrs.begin(), ctrs.end(), boost::bind(&AbsConstraint::setGammaS,_1,1));

        // random solution
        randomSolver(vars.begin(), vars.end());



        MinGamma<false> minGamma;

        if(tabuCol(vars.begin(), vars.end(), ctrs.begin(), ctrs.end(), minGamma,
                   Normal<VarPtr>(), Normal<CtrPtr>())){
            int nbViolated = std::accumulate(ctrs.begin(), ctrs.end(), 0,
                                             Counter_violated<CtrPtr>());

            cout << "feasible: " << nbViolated << " " << ctrs.size() <<endl;
        }else{
            int nbViolated = std::accumulate(ctrs.begin(), ctrs.end(), 0,
                                             Counter_violated<CtrPtr>());

            cout << "violated: " << nbViolated <<  " " << ctrs.size() <<endl;

            BOOST_FOREACH(CtrPtr const c, ctrs){
                if(c->isViolated()){
                    cout << c->getID() << "\t" << c->getGammaS() << endl;
                }
            }
        }
    }


    void testMac(std::vector<VarPtr>& vars, std::vector<CtrPtr>& ctrs){

        boost::progress_timer t;

        std::for_each(ctrs.begin(), ctrs.end(), boost::bind(&AbsConstraint::setGammaS,_1,1));
        std::for_each(ctrs.begin(), ctrs.end(), boost::bind(&AbsConstraint::print,_1));

        ctrs[2]->setLevel(disable);
        ctrs[4]->setLevel(disable);
        ctrs[5]->setLevel(disable);

        //         ctrs[2]->setGammaS(0);
        //         ctrs[4]->setGammaS(0);
        //         ctrs[5]->setGammaS(0);

        vars[3]->setLevel(disable);
        vars[3]->setGammaS(0);

        bool depass = false;
        if(mac_solver(vars.begin(), vars.end(), Normal<VarPtr>(), Normal<CtrPtr>(),depass)){
            int nbViolated = 0;

            cout << "feasible: " << nbViolated << " " << ctrs.size() <<endl;
        }else{


            cout << "violated: "  <<endl;


        }
    }


    void info(std::vector<VarPtr>& vars, std::vector<CtrPtr>& ctrs){

        std::vector<VarPtr>::iterator ir = std::max_element(vars.begin(),
                                                            vars.end(),
                                                            boost::bind(std::less<int>(),
                                                                        boost::bind(&Variable::domSize,_1),
                                                                        boost::bind(&Variable::domSize,_2)));

        int freq = 0;
        BOOST_FOREACH(VarPtr const v, vars){
            int f = (v->values_.back())->v_;
            if(f>freq)
                freq = f;
        }
	

        cout << "info:\t" << vars.size()
             << "\t" << ctrs.size()
                //<< "\t" << (*ir)->domSize() << endl;
             << "\t" << freq ; // << endl;

        std::for_each(ctrs.begin(), ctrs.end(),
                      boost::bind(&AbsConstraint::setGammaS,_1,1));
        std::for_each(vars.begin(), vars.end(),
                      boost::bind(&Variable::setGammaS,_1,1));

        bool feasible = arc_consistency(vars.begin(),vars.end(),
                                        Normal<VarPtr>(),Normal<CtrPtr>());

        if(feasible){
            cout << "\t CON" << endl;
        }else{
            cout << "\t INC" << endl;
        }
    }

    void testExtract(std::vector<VarPtr>& vars, std::vector<CtrPtr>& ctrs){

        //boost::progress_timer t;


        //if(!extractIntensive(vars.begin(), vars.end(),ctrs.begin(), ctrs.end())){


        // extract IIS by means of variables

        boost::timer timer;

        bool  feasible = extract(vars.begin(), vars.end(),ctrs.begin(), ctrs.end());


        // final output
        //        if(feasible){
        //            cout << "r:\t FAILURE" << endl;
        //        }else{
        //            cout << "r:\t found" << endl;
        //        }


    }



    void testExtractCorrecting(std::vector<VarPtr>& vars, std::vector<CtrPtr>& ctrs){

        //boost::progress_timer t;


        //if(!extractIntensive(vars.begin(), vars.end(),ctrs.begin(), ctrs.end())){

        int pervious = 0;
        bool feasible = false;

        // extract IIS by means of variables
        bool isFirst = true;
        boost::timer timer;
        while(!feasible && pervious != vars.size()){

            if(isFirst)
                feasible = extract(vars.begin(), vars.end(),ctrs.begin(), ctrs.end());
            else
                feasible = extractVars(vars.begin(), vars.end(),ctrs.begin(), ctrs.end());


            pervious = vars.size();


            if(feasible){
                BOOST_FOREACH(CtrPtr const c, ctrs){
                    if(c->isViolated()){
                        cout << "found violated constraint" << endl;
                        exit(0);
                    }
                }

                cout << "all constraints consistent" << endl;
            }



            if(isFirst){
                ctrs.erase(remove_if(ctrs.begin(), ctrs.end(), NonCritical<CtrPtr>()), ctrs.end() );
                vars.erase(remove_if(vars.begin(), vars.end(), NonCritical<VarPtr>()), vars.end() );




                if(isFirst){
                    isFirst = false;
                    cout << "varis:\t" << vars.size() << "\t" << ctrs.size() << "\t" << timer.elapsed()<< endl;
                }


                BOOST_FOREACH(VarPtr const var, vars){
                    (var->ctrs_).erase(remove_if(var->ctrs_.begin(), var->ctrs_.end(), NonCritical<CtrPtr>() ),
                                       var->ctrs_.end() );
                    (var->vars_).erase(remove_if(var->vars_.begin(), var->vars_.end(), NonCritical<VarPtr>() ),
                                       var->vars_.end() );
                    var->resetDomain();
                    var->resign();
                    var->resetAtts();

                }
            }else{
                ctrs.erase(remove_if(ctrs.begin(), ctrs.end(), NonMarked<CtrPtr>()), ctrs.end() );
                vars.erase(remove_if(vars.begin(), vars.end(), NonMarked<VarPtr>()), vars.end() );


                BOOST_FOREACH(VarPtr const var, vars){
                    (var->ctrs_).erase(remove_if(var->ctrs_.begin(), var->ctrs_.end(), NonMarked<CtrPtr>() ),
                                       var->ctrs_.end() );
                    (var->vars_).erase(remove_if(var->vars_.begin(), var->vars_.end(), NonMarked<VarPtr>() ),
                                       var->vars_.end() );
                    var->resetDomain();
                    var->resign();
                    var->resetAtts();

                }

            }

            std::for_each(ctrs.begin(),ctrs.end(),boost::bind(&AbsConstraint::resetAtts,_1));
            std::for_each(ctrs.begin(),ctrs.end(),boost::bind(&AbsConstraint::setGammaS,_1,1));
            //cout << "next: " << vars.size() << "\t" << ctrs.size() << endl;
        }

        // print result for
        cout << "variis:\t" << vars.size() << "\t" << ctrs.size() << "\t" << timer.elapsed()<< endl;

        // extract IIS by means of constraints
        int perviousCtr = 0;
        while(!feasible && perviousCtr != ctrs.size()){
            perviousCtr = ctrs.size();

            feasible = extractConstraint(vars.begin(), vars.end(),ctrs.begin(), ctrs.end());


            ctrs.erase(remove_if(ctrs.begin(), ctrs.end(), NonMarked<CtrPtr>()), ctrs.end() );
            vars.erase(remove_if(vars.begin(), vars.end(), NonMarked<VarPtr>()), vars.end() );



            BOOST_FOREACH(VarPtr const var, vars){
                (var->ctrs_).erase(remove_if(var->ctrs_.begin(), var->ctrs_.end(),NonMarked<CtrPtr>()),
                                   var->ctrs_.end() );
                (var->vars_).erase(remove_if(var->vars_.begin(), var->vars_.end(),NonMarked<VarPtr>()),
                                   var->vars_.end() );
                var->resetDomain();
                var->resign();
                var->resetAtts();
            }

            std::for_each(ctrs.begin(),ctrs.end(),boost::bind(&AbsConstraint::resetAtts,_1));
            std::for_each(ctrs.begin(),ctrs.end(),boost::bind(&AbsConstraint::setGammaS,_1,1));

        }

        // print result for IIS constraints
        cout << "ctriis:\t" << vars.size() << "\t" << ctrs.size() << "\t"<< timer.elapsed()<< endl;


        // final output
        if(feasible){
            cout << "r:\t FAILURE" << endl;
        }else{
            cout << "r:\t found" << endl;
        }


    }



    void testExtractColoring(std::vector<VarPtr>& vars, std::vector<CtrPtr>& ctrs){

        //boost::progress_timer t;


        //if(!extractIntensive(vars.begin(), vars.end(),ctrs.begin(), ctrs.end())){

        int pervious = 0;
        bool feasible = false;

        // extract IIS by means of variables
        bool isFirst = true;
        boost::timer timer;


        MinGamma<false> minGamma;

        // pprefilter the problem
        feasible = prefiltering(vars.begin(),vars.end(),
                                ctrs.begin(), ctrs.end(),minGamma,
                                Normal<VarPtr>(),Normal<CtrPtr>());

//        int sizeVarx = std::accumulate(vars.begin(),vars.end(),0,Counter_Marked<VarPtr>());
//        int sizeCtrx = std::accumulate(ctrs.begin(),ctrs.end(),0,Counter_Marked<CtrPtr>());
//        cout << "varmarked:\t" << sizeVarx << "\t" << sizeCtrx << "\t" << timer.elapsed()<< endl;


        ctrs.erase(std::remove_if(ctrs.begin(),ctrs.end(),NonMarked<CtrPtr>()),ctrs.end());
        vars.erase(std::remove_if(vars.begin(),vars.end(),NonMarked<VarPtr>()),vars.end());



        BOOST_FOREACH(VarPtr const v,vars){
            v->vars_.erase(std::remove_if(v->vars_.begin(),v->vars_.end(),
                                          NonMarked<VarPtr>()),v->vars_.end());
            v->ctrs_.erase(std::remove_if(v->ctrs_.begin(),v->ctrs_.end(),
                                          NonMarked<CtrPtr>()),v->ctrs_.end());

            v->setGammaS(1);
            v->setLevel(normal);
            v->resign();
            v->resetDomain();
        }


        BOOST_FOREACH(CtrPtr const c, ctrs){
            c->setGammaS(1);
            c->setLevel(normal);
        }


        feasible = prefiltering(vars.begin(),vars.end(),
                                ctrs.begin(), ctrs.end(),minGamma,
                                Normal<VarPtr>(),Normal<CtrPtr>());


        ctrs.erase(std::remove_if(ctrs.begin(),ctrs.end(),NonMarked<CtrPtr>()),ctrs.end());
        vars.erase(std::remove_if(vars.begin(),vars.end(),NonMarked<VarPtr>()),vars.end());



        BOOST_FOREACH(VarPtr const v,vars){
            v->vars_.erase(std::remove_if(v->vars_.begin(),v->vars_.end(),
                                          NonMarked<VarPtr>()),v->vars_.end());
            v->ctrs_.erase(std::remove_if(v->ctrs_.begin(),v->ctrs_.end(),
                                          NonMarked<CtrPtr>()),v->ctrs_.end());

            v->setGammaS(1);
            v->setLevel(normal);
            v->resign();
            v->resetDomain();
        }


        BOOST_FOREACH(CtrPtr const c, ctrs){
            c->setGammaS(1);
            c->setLevel(normal);
        }


        cout << "fi:\t" << vars.size() << "\t" <<  ctrs.size();
        // in critical subgraph identification
        if(!feasible){
            feasible = extract(vars.begin(), vars.end(),ctrs.begin(), ctrs.end());
        }else{
            cout << "lost" << endl;
        }





        if(feasible){
            BOOST_FOREACH(CtrPtr const c, ctrs){
                if(c->isViolated()){
                    cout << "found violated constraint" << endl;
                    exit(0);
                }
            }

            cout << "all constraints consistent" << endl;
        }


        int sizeVar = std::accumulate(vars.begin(),vars.end(),0,Counter_Critical<VarPtr>());
        int sizeCtr = std::accumulate(ctrs.begin(),ctrs.end(),0,Counter_Critical<CtrPtr>());
        cout << "varis:\t" << sizeVar << "\t" << sizeCtr << "\t" << timer.elapsed()<< endl;



        // final output
        if(feasible){
            cout << "r:\t FAILURE" << endl;
        }else{
            cout << "r:\t found" << endl;
        }


    }



    void testColoring(std::vector<VarPtr>& vars, std::vector<CtrPtr>& ctrs){

        // prefiltering the problem
        //        randomSolver(vars.begin(),vars.end());
        //        MinGamma<false> minGamma;
        //        bool feasible = prefiltering(vars.begin(),vars.end(),ctrs.begin(),
        //                                     ctrs.end(),minGamma, Normal<VarPtr>(),
        //                                     Normal<CtrPtr>());

        //        if(feasible){
        //            cout << "prefiltering found feasible" << endl;
        //            exit(0);
        //        }



        //        // removing the no marked elements and reset all effects
        //        ctrs.erase(remove_if(ctrs.begin(), ctrs.end(), NonMarked<CtrPtr>()), ctrs.end() );
        //        vars.erase(remove_if(vars.begin(), vars.end(), NonMarked<VarPtr>()), vars.end() );



        //        BOOST_FOREACH(VarPtr const var, vars){
        //            (var->ctrs_).erase(remove_if(var->ctrs_.begin(), var->ctrs_.end(),NonMarked<CtrPtr>()),
        //                               var->ctrs_.end() );
        //            (var->vars_).erase(remove_if(var->vars_.begin(), var->vars_.end(),NonMarked<VarPtr>()),
        //                               var->vars_.end() );
        //            var->resetDomain();
        //            var->resign();
        //            var->resetAtts();
        //        }

        //        std::for_each(ctrs.begin(),ctrs.end(),boost::bind(&AbsConstraint::resetAtts,_1));
        //        std::for_each(ctrs.begin(),ctrs.end(),boost::bind(&AbsConstraint::setGammaS,_1,1));




        testExtractColoring(vars,ctrs);


        // identifying the IIS


    }


    template <typename CtrItr>
    void graphviz(CtrItr first, CtrItr last){
        cout << "graph G{" << endl;
        cout << "node [shape=point];" << endl;
        for(; first != last; ++first){
            (*first)->printViz();
            cout << ";" << endl;
        }
        cout << "}" << endl;
    }


    void solve(std::vector<VarPtr>& vars_, std::vector<CtrPtr>& ctrs_){
        std::for_each(ctrs_.begin(), ctrs_.end(), boost::bind(&AbsConstraint::setGammaS,_1,1));

        //std::for_each(vars_.begin(), vars_.end(), boost::bind(&Variable::resetDomain, _1));

        srand((unsigned)time(0));


        //boost::progress_timer t;

        // Test 1:
        //this->tabuTest(vars_,ctrs_);


        // Test 1.1: ac testing
        //this->ac(vars_,ctrs_);


        // Test 2:
        //this->grow(vars_,ctrs_);



        // Test 4: mac solver
        //this->test_Mac(vars_,ctrs_);


        // Test 5: extension mark
        //this->extend(vars_,ctrs_);


        // Test 6: extension mark -- with guarrentee routine
        //extend_Guarrentee(vars_,ctrs_);

        // Test 6.1: extension mark -- with conservation
        //extend_conserver(vars_,ctrs_);


        /// Stable routine
        // Test 6.2: extension mark -- with guarrentee routine - extract
        //extract_routine(vars_,ctrs_);


        // Test 6.3: extension mark -- with guarrentee routine - breakout
        //breakout_routine(vars_,ctrs_);


        // Test 6.4: extension mark -- with guarrentee routine - breakout and heuristic Reach
        //breakout_heuristicReach_routine(vars_,ctrs_);



        /// ================================ New version test ====== BGN
        //testTabu(vars_, ctrs_);


        //testMac(vars_, ctrs_);

        //info(vars_, ctrs_);
        //testExtract(vars_, ctrs_);

        // random shuffle the vars ordering
//        std::random_shuffle(vars_.begin(),vars_.end());
//        testExtractCorrecting(vars_, ctrs_);



        /// testColoring(vars_, ctrs_);
        //colorExtract(vars_,ctrs_);

        //std::for_each(ctrs_.begin(),ctrs_.end(),boost::bind(&AbsConstraint::print,_1));
        //graphviz(ctrs_.begin(),ctrs_.end()); // print the dot file


        //testAC(vars_, ctrs_);
        
        //testACwave(vars_, ctrs_);
        //printCtrs(ctrs_.begin(), ctrs_.end());

        /// coloring problem
        testExtractColoring(vars_,ctrs_);
        //graphviz(ctrs_.begin(), ctrs_.end());

    }




    template <typename InputIterator>
    int validCnt(InputIterator first, InputIterator last){
        int cnt = 0;
        for(; first != last; ++first){
            if((*first)->isEnable())
                ++cnt;
        }
        return cnt;
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
    }

};




#endif // DATAIO_HPP
