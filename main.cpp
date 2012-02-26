//#include <QtCore/QCoreApplication>


#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <time.h>
#include <ctime>

#include "dataIO.hpp"


//#include "heuristic.hpp";

using namespace std;


struct A{

    void f(int i, int j){
        cout << i << "\t" << j << endl;
    }

    void d(){}

};


int main(int argc, char *argv[])
{


    /// rlfap celar et graph
//    LRFAP_solver<Solver>* pbm = new LRFAP_solver<Solver>();
//    pbm->readDomFile(argv[1]);
//    pbm->readVarFile(argv[2]);
//    pbm->readCtrFile(argv[3]);


    /// dimacs
    DIMACS_coloring<Solver>* pbm = new DIMACS_coloring<Solver>();
    pbm->read(argv[1],
              boost::numeric_cast<int32_t>(boost::lexical_cast<int>(argv[2])));


    /// roadef 2001
    //    FAPP<Solver>* pbm = new FAPP<Solver>();
    //    pbm->read(argv[1],boost::numeric_cast<int32_t>(boost::lexical_cast<int>(argv[2])));




    pbm->solve();

    delete pbm;




    exit(0);



}
