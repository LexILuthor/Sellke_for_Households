//
// Created by popcorn on 06/01/2021.
//

#ifndef SELLKE_FOR_HOUSEHOLDS_MYFUNCTIONS_H
#define SELLKE_FOR_HOUSEHOLDS_MYFUNCTIONS_H

#endif //SELLKE_FOR_HOUSEHOLDS_MYFUNCTIONS_H


void generate_and_decide_the_household(std::vector<double> &Q, std::vector<double> &L, std::vector<double> &I, double ny, double gamma, int N,
              std::vector<size_t> &household_is);

void generate(std::vector<double> &Q, std::vector<double> &L, std::vector<double> &I, double ny, double gamma, int N);

int in_which(std::vector<size_t> &household_is, int i, int number_of_people_per_household);


void activeintervels(double from, std::vector<double> &startInfection, std::vector<double> &endInfection,
                     std::vector<size_t> indexEndInfection, int &last);

void read_Parameters_From_File(std::string inputpath, int &number_of_households, int &number_of_people_per_household,
                               double &beta, double &betaH, double &ny, double &gamma);

void write_the_csv_file(std::string outputpath, std::vector<std::vector<int> > &SEIR, std::vector<double> &temp);



//template for the function that returns the ordered index

#include <algorithm>
#include <numeric>      // std::iota

template<typename T>
std::vector<size_t> sort_indexes(const std::vector<T> &v) {

    // initialize original index locations
    std::vector<size_t> idx(v.size());
    iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    // using std::stable_sort instead of std::sort
    // to avoid unnecessary index re-orderings
    // when v contains elements of equal values
    stable_sort(idx.begin(), idx.end(),
                [&v](size_t i1, size_t i2) { return v[i1] < v[i2]; });

    return idx;
}
