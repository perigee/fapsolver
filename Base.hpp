#ifndef BASE_HPP
#define BASE_HPP


#include <vector>
#include <deque>

#include <boost/shared_ptr.hpp>
#include <boost/foreach.hpp>
#include <boost/progress.hpp>


// the elements' level
enum Level {disable, normal, marked, critical};
enum Cate {CateFilter, CateLocate, CateConstruct, CateExact, CateEND};
enum Polar {Neu,V,H};
//enum ProblemType

struct dummyT{
    // nothing inside
};

template <typename T>
        struct AlgoRelated : public T{

    AlgoRelated():l_(normal),
                  gammaD_(0),
                  gammaS_(0),
                  weight_(0){}


    void setLevel(Level l){l_=l;}
    Level getLevel() const{return l_;}


    // dynamic gamma
    void setGammaD(int g){gammaD_ = g;}
    int  getGammaD() const{return gammaD_ ;}
    void resetGammaD(){gammaD_ = 0;}

    // static gamma
    void setGammaS(int g){gammaS_ = g;}
    int  getGammaS() const{return gammaS_ ;}
    void resetGammaS(){gammaS_ = 0;}
    void resetGamma(){gammaS_ = 0;gammaD_=0;}

    // static weight for learning
    void setWeight(int w){weight_ = w;}
    void addWeight(){++weight_;}
    int  getWeight() const{return weight_;}
    void resetWeight(){weight_=0;}

    // tabu tenure
    void setTenure(int d){tenure_ = d;}
    int  getTenure() const{return tenure_;}
    void resetTenure(){tenure_ = 0;}
    void reduceTenure(){if(tenure_>0) --tenure_;}
    bool inTenure() const{return tenure_ > 0;}

    // reset all attributes
    void resetAtts(){l_ = normal;clearAtts();}
    void clearAtts(){gammaD_ = 0; gammaS_=0;weight_=0;tenure_=0;}

    // list of attributes
    Level l_;       // level of variable
    int gammaD_;    // dynamic gamma system
    int gammaS_;    // static gamma system used in
    int weight_;    // weight system
    int tenure_;    // tabu tenure
};



template <typename Element>
        struct Counter_markedWeight
{

    int operator ()(int result, Element const e){
        if(e->isMarked())
            result += e->getWeight();

        return result;
    }
};


template <typename Element>
        struct Counter_gammaWeight
{

    int operator ()(int result, Element const e){
        if(e->getGamma() > 0)
            result += e->getWeight();

        return result;
    }
};


template <typename Element>
        struct Counter_violated
{

    int operator ()(int result, Element const e){
        if(e->isViolated())
            ++result;

        return result;
    }
};


template <typename Element>
        struct Counter_marked
{

    int operator ()(int result, Element const e){
        if(e->isMarked())
            ++result;

        return result;
    }
};


template <typename Element>
        struct Counter_gamma
{

    int operator ()(int result, Element const e){
        if(e->getGamma() > 0)
            ++result;

        return result;
    }
};


template <typename Element>
        struct Non_Marked
{

    int operator ()(Element const e){

        return !e->isMarked();
    }
};


////======================================= Only for testing
template <typename Element>
inline bool sortWeight(Element const a, Element const b){
    return a->getWeight()< b->getWeight();
}

#endif // BASE_HPP
