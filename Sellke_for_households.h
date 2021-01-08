//
// Created by popcorn on 06/01/2021.
//

#ifndef SELLKE_FOR_HOUSEHOLDS_SELLKE_FOR_HOUSEHOLDS_H
#define SELLKE_FOR_HOUSEHOLDS_SELLKE_FOR_HOUSEHOLDS_H

#endif //SELLKE_FOR_HOUSEHOLDS_SELLKE_FOR_HOUSEHOLDS_H


std::vector<double>
sellke_for_households(int number_of_households, int number_of_people_per_household, double beta, double betaH,
                      double ny, double gamma, std::vector<double> &startInfection,
                      std::vector<double> &endInfection);

std::vector<double> sellke(int N, double beta, double ny, double gamma, std::vector<double> &startInfection,
                           std::vector<double> &endInfection, std::vector<double> &Q, std::vector<double> &L,
                           std::vector<double> &I, std::vector<bool> &has_already_been_exposed);


double integral(double from, int last_infected, std::vector<double> &startInfection,
                std::vector<double> &endInfection, double &integral_value, std::vector<double> &Q, double beta,
                int &numberofactiveintervals, int &last_event_R);

