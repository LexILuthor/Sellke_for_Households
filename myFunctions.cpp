//
// Created by popcorn on 06/01/2021.
//


#include <iostream>
#include <fstream>
#include <vector>
#include<random>
#include <time.h>
#include <algorithm>
#include "myFunctions.h"







int in_which(std::vector<size_t> &household_is,int i, int number_of_people_per_household){
    return (int) household_is[i]/number_of_people_per_household;
}


void activeintervels(double from, std::vector<double> &startInfection, std::vector<double> &endInfection,
                     std::vector<size_t> indexEndInfection, int &last) {
    int numberofactiveintervals = 0;
    for (int j = 0; j < indexEndInfection.size(); j++) {
        int i = indexEndInfection[j];
        if (endInfection[i] < from) {
            last = j;
        } else if (startInfection[i] < from && endInfection[i] > from) {
            numberofactiveintervals++;
        }
    }
}


void read_Parameters_From_File(std::string inputpath, int &number_of_households, int &number_of_people_per_household,
                               double &beta, double &betaH, double &ny, double &gamma) {

    std::string line;
    std::ifstream infile(inputpath);
    if (infile.is_open()) {
        getline(infile, line, ':');
        getline(infile, line);
        number_of_households = std::stoi(line);

        getline(infile, line, ':');
        getline(infile, line);
        number_of_people_per_household = std::stoi(line);

        getline(infile, line, ':');
        getline(infile, line);
        beta = std::stod(line);

        getline(infile, line, ':');
        getline(infile, line);
        betaH = std::stod(line);

        getline(infile, line, ':');
        getline(infile, line);
        ny = std::stod(line);

        getline(infile, line, ':');
        getline(infile, line);
        gamma = std::stod(line);

        infile.close();
    } else std::cout << "Unable to open file";

}

void generate_and_decide_the_household(std::vector<double> &Q, std::vector<double> &L, std::vector<double> &I, double ny, double gamma, int N,
                                       std::vector<size_t> & household_is) {

    // for random seed:
    std::default_random_engine generator(time(0));


    std::exponential_distribution<double> exp1_distribution(1.0);
    std::exponential_distribution<double> expNy_distribution(ny);
    std::exponential_distribution<double> expGamma_distribution(gamma);


    // Generate Q, L, I respectively the vectors of the resistance to the infection, expose time, infected time


    for (int i = 0; i < N; i++) {
        Q[i] = exp1_distribution(generator);
        L[i] = expNy_distribution(generator);
        I[i] = expGamma_distribution(generator);

        //function used to guarantee independence (may be not necessary)
        /*
        exp1_distribution.reset();
        expNy_distribution.reset();
        expGamma_distribution.reset();
         */
    }

    household_is = sort_indexes(Q);

    //the first one is already infected
    Q[0] = 0;
    //order Q
    std::sort(Q.begin(), Q.end());
}


void write_the_csv_file(std::string outputpath, std::vector<std::vector<int> > &SEIR, std::vector<double> &temp) {
    std::ofstream outfile(outputpath);
    if (!outfile.is_open()) {
        std::cout << "Unable to open file";
    } else {
        for (int i = 0; i < temp.size(); i++) {

            outfile << SEIR[0][i] << ",\t" << SEIR[1][i] << ",\t" << SEIR[2][i] << ",\t" << SEIR[3][i] << ",\t"
                    << temp[i] << "\n";
        }
        outfile.close();
    }
}



void generate(std::vector<double> &Q, std::vector<double> &L, std::vector<double> &I, double ny, double gamma, int N) {
    // for random seed:
    std::default_random_engine generator(time(0));

    std::exponential_distribution<double> exp1_distribution(1.0);
    std::exponential_distribution<double> expNy_distribution(ny);
    std::exponential_distribution<double> expGamma_distribution(gamma);


    // Generate Q, L, I respectively the vectors of the resistance to the infection, expose time, infected time


    for (int i = 0; i < N; i++) {
        Q[i] = exp1_distribution(generator);
        L[i] = expNy_distribution(generator);
        I[i] = expGamma_distribution(generator);

        //function used to guarantee independence (may be not necessary)
        /*
        exp1_distribution.reset();
        expNy_distribution.reset();
        expGamma_distribution.reset();
         */
    }

    //the first one is already infected
    Q[0] = 0;
    //order Q
    std::sort(Q.begin(), Q.end());
}