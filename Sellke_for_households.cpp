//
// Created by popcorn on 06/01/2021.
//


#include <vector>
#include <iostream>
#include "myFunctions.h"
#include "Sellke_for_households.h"


std::vector<double>
sellke_for_households(int number_of_households, int number_of_people_per_household, double beta, double betaH,
                      double ny, double gamma, std::vector<std::vector<int> > &SEIR) {

    int N = number_of_households * number_of_people_per_household;
    // Generate Q, L, I respectively the vectors of the resistance to the infection, expose time, infected time
    std::vector<double> Q(N);
    std::vector<double> L(N);
    std::vector<double> I(N);
    std::vector<size_t> household_is(N);

    std::vector<bool> has_already_been_exposed(N, false);
    has_already_been_exposed[0] = true;

    //generates the vectors and order Q
    generate_and_decide_the_household(Q, L, I, ny, gamma, N, household_is);



    // matrix containiing ar each row (i.e. at each household) the index of the people in that household
    std::vector<std::vector<int>> household_list(number_of_households, std::vector<int>(0, -1));
    //initialize household
    for (int i = 0; i < Q.size(); i++) {
        household_list[in_which(household_is, i, number_of_people_per_household)].push_back(i);
    }


    //vector containing at the i-th position the time of first infectious contact of individual i
    std::vector<double> time_vector(N, -1);
    time_vector[0] = 0;

    //indice della persona infettata più di recente
    int last_infected = 0;



    //valore dell'integrale tra 0 e ts
    double integral_value;
    integral_value = 0;

    //per ogni cella i abbaimo il tempo in cui l'individuo i-esimo comincia ad essere infetto e smette di essere infetto
    //nota che in realtà endInfection=startInfection+I dove I è il vettore degli infected time, mentre
    //startInfection = time_vector+ L dove L è il vettore degli exposed times


    std::vector<double> startInfection;
    std::vector<double> endInfection;

    startInfection.push_back(L[0]);
    endInfection.push_back(L[0] + I[0]);

    // Ultimo tempo t fino al quale l'integrale è conosciuto
    //coinciderà sempre con l'ultimo elemento di "time_vector" tranne che per questo primo step
    //nota: viene inizializzato sotto
    double ts;

    //number of active intervals at time ts
    int number_of_active_intervals = 0;

    //posizione (all'interno del vettore ordinato) dell'ultimo evento del tipo "persona diventa Recovered" che avviene prima del tempo ts
    //NB non è sempre l'ultimo evento, potrebbe essere anche il penultimo
    int last_event_R = 0;


    int infected_household = in_which(household_is, last_infected, number_of_people_per_household);

    std::vector<double> Q_of_infected_household;
    std::vector<double> I_of_infected_household;
    std::vector<double> L_of_infected_household;

    has_already_been_exposed[last_infected] = true;
    Q_of_infected_household.push_back(Q[last_infected]);
    L_of_infected_household.push_back(L[last_infected]);
    I_of_infected_household.push_back(I[last_infected]);

    for (int i = 0; i < number_of_people_per_household; i++) {
        int index_person_in_household = household_list[infected_household][i];
        if (!has_already_been_exposed[index_person_in_household]) {
            Q_of_infected_household.push_back(Q[index_person_in_household]);
            L_of_infected_household.push_back(L[index_person_in_household]);
            I_of_infected_household.push_back(I[index_person_in_household]);
        }
    }

    std::vector<double> startInfection_of_household;
    std::vector<double> endInfection_of_household;
    startInfection_of_household.push_back(L[last_infected] + ts);
    endInfection_of_household.push_back(L[last_infected] + ts + I[last_infected]);

    std::vector<double> time_vector_tmp(Q_of_infected_household.size(), -1);
    time_vector_tmp=sellke(number_of_people_per_household, betaH, ny, gamma, startInfection_of_household,
           endInfection_of_household, Q_of_infected_household, L_of_infected_household,
           I_of_infected_household, has_already_been_exposed);

    for(int i=1;i<startInfection_of_household.size();i++){
        startInfection.push_back(startInfection_of_household[i]);
    }
    for(int i=1;i<endInfection_of_household.size();i++){
        endInfection.push_back(endInfection_of_household[i]);
    }

    int tmp=1;
    for (int i = 0; i < number_of_people_per_household; i++) {
        int index_person_in_household = household_list[infected_household][i];
        if (!has_already_been_exposed[index_person_in_household]) {
            time_vector[index_person_in_household]=time_vector_tmp[tmp];
            tmp++;
        }
    }


    // Ultimo tempo t fino al quale l'integrale è conosciuto
    //coinciderà sempre con l'ultimo elemento di "time_vector" tranne che per questo primo step
    ts = L[0];


    while (last_infected < N - 1) {
        if (!has_already_been_exposed[last_infected + 1]) {
            ts = integral(ts, last_infected, startInfection, endInfection, integral_value, Q, (beta / N),
                          number_of_active_intervals, last_event_R);
            if (ts == -1) {
                //non abbiamo più nuovi contagi
                break;
            } else {
                last_infected++;
                time_vector[last_infected] = ts;
                startInfection.push_back(ts + L[last_infected]);
                endInfection.push_back(ts + L[last_infected] + I[last_infected]);

                //teoricamente il comando seguente è già vero
                integral_value = Q[last_infected];
                has_already_been_exposed[last_infected] = true;

                int infected_household = in_which(household_is, last_infected, number_of_people_per_household);

                std::vector<double> Q_of_infected_household;
                std::vector<double> I_of_infected_household;
                std::vector<double> L_of_infected_household;

                has_already_been_exposed[last_infected] = true;
                Q_of_infected_household.push_back(Q[last_infected]);
                L_of_infected_household.push_back(L[last_infected]);
                I_of_infected_household.push_back(I[last_infected]);

                for (int i = 0; i < number_of_people_per_household; i++) {
                    int index_person_in_household = household_list[infected_household][i];
                    if (!has_already_been_exposed[index_person_in_household]) {
                        Q_of_infected_household.push_back(Q[index_person_in_household]);
                        L_of_infected_household.push_back(L[index_person_in_household]);
                        I_of_infected_household.push_back(I[index_person_in_household]);
                    }
                }

                std::vector<double> startInfection_of_household;
                std::vector<double> endInfection_of_household;
                startInfection_of_household.push_back(L[last_infected] + ts);
                endInfection_of_household.push_back(L[last_infected] + ts + I[last_infected]);

                std::vector<double> time_vector_tmp(Q_of_infected_household.size(), -1);
                time_vector_tmp=sellke(number_of_people_per_household, betaH, ny, gamma, startInfection_of_household,
                                       endInfection_of_household, Q_of_infected_household, L_of_infected_household,
                                       I_of_infected_household, has_already_been_exposed);

                for(int i=1;i<startInfection_of_household.size();i++){
                    startInfection.push_back(startInfection_of_household[i]);
                }
                for(int i=1;i<endInfection_of_household.size();i++){
                    endInfection.push_back(endInfection_of_household[i]);
                }

                int tmp=1;
                for (int i = 0; i < number_of_people_per_household; i++) {
                    int index_person_in_household = household_list[infected_household][i];
                    if (!has_already_been_exposed[index_person_in_household]) {
                        time_vector[index_person_in_household]=time_vector_tmp[tmp];
                        tmp++;
                    }
                }






            }
        } else {
            last_infected++;
        }
    }


    return time_vector;


}


double integral(double from, int last_infected, std::vector<double> &startInfection,
                std::vector<double> &endInfection, double &integral_value, std::vector<double> &Q, double beta,
                int &numberofactiveintervals, int &last_event_R) {
    std::vector<size_t> indexStartInfection = sort_indexes(startInfection);
    std::vector<size_t> indexEndInfection = sort_indexes(endInfection);

    //compute the number of intervals that contains "from", it also sets not_ordered_j to the last individual before from
    //int last_event_R;
    //activeintervels(from, startInfection, endInfection, indexEndInfection, last_event_R);
    int last_event_I = last_event_R;

    double last_from = from;
    bool check = true;
    while (last_event_I < indexStartInfection.size() || last_event_R < indexEndInfection.size()) {
        int i = 0;
        int j = indexEndInfection[last_event_R];
        if (last_event_I >= indexStartInfection.size()) {
            check = false;
        } else {
            i = indexStartInfection[last_event_I];
        }

        if (startInfection[i] < endInfection[j] && check) {
            if (startInfection[i] >= last_from) {
                // nota che "last_from" si potrebbe sostituire con "from"
                if (integral_value + numberofactiveintervals * beta * (startInfection[i] - last_from) >=
                    Q[last_infected + 1]) {
                    // from qui diventa il tempo t al quale il nuovo infetto Q[last_infected+1] subisce il contatto infetto
                    from = (double) (Q[last_infected + 1] - integral_value +
                                     (numberofactiveintervals * beta * last_from)) /
                           (numberofactiveintervals * beta);


                    //comando inutile, è solo di controllo
                    // NB se si toglie controllare che il comando
                    // integral_value = Q[last_infected];
                    // a riga 64 sia attivo
                    //integral_value = integral_value + numberofactiveintervals * beta * (from - last_from);


                    return from;
                } else {
                    from = startInfection[i];
                    integral_value = integral_value + (numberofactiveintervals * beta * (from - last_from));
                    numberofactiveintervals++;
                    last_from = startInfection[i];
                }
            }
            last_event_I++;
        } else {
            if (endInfection[j] > last_from) {

                if (integral_value + numberofactiveintervals * beta * (endInfection[j] - last_from) >=
                    Q[last_infected + 1]) {
                    // from qui diventa il tempo t al quale il nuovo infetto Q[last_infected+1] subisce il contatto infetto
                    from = (double) (Q[last_infected + 1] - integral_value +
                                     numberofactiveintervals * beta * last_from) /
                           (numberofactiveintervals * beta);


                    return from;
                } else {
                    from = endInfection[j];
                    integral_value = integral_value + (numberofactiveintervals * beta * (from - last_from));
                    numberofactiveintervals--;
                    last_from = endInfection[j];
                }
            }
            last_event_R++;
        }
    }
    return -1;
}


std::vector<double> sellke(int N, double beta, double ny, double gamma, std::vector<double> &startInfection,
                           std::vector<double> &endInfection, std::vector<double> &Q, std::vector<double> &L,
                           std::vector<double> &I, std::vector<bool> &has_already_been_exposed) {





    //vector containing at the i-th position the time of first infectious contact of individual i
    std::vector<double> time_vector(N, -1);
    time_vector[0] = 0;

    //indice della persona infettata più di recente
    int last_infected = 0;

    // Ultimo tempo t fino al quale l'integrale è conosciuto
    //coinciderà sempre con l'ultimo elemento di "time_vector" tranne che per questo primo step
    double ts = L[0];

    //valore dell'integrale tra 0 e ts
    double integral_value;
    integral_value = Q[0];

    //per ogni cella i abbaimo il tempo in cui l'individuo i-esimo comincia ad essere infetto e smette di essere infetto
    //nota che in realtà endInfection=startInfection+I dove I è il vettore degli infected time, mentre
    //startInfection = time_vector+ L dove L è il vettore degli exposed times


    startInfection.push_back(L[0]);
    endInfection.push_back(L[0] + I[0]);

    //number of active intervals at time ts
    int number_of_active_intervals = 0;

    //posizione (all'interno del vettore ordinato) dell'ultimo evento del tipo "persona diventa Recovered" che avviene prima del tempo ts
    //NB non è sempre l'ultimo evento, potrebbe essere anche il penultimo
    int last_event_R = 0;
    while (last_infected < N - 1) {
        if (!has_already_been_exposed[last_infected + 1]) {
            ts = integral(ts, last_infected, startInfection, endInfection, integral_value, Q, (beta / N),
                          number_of_active_intervals, last_event_R);
            if (ts == -1) {
                //non abbiamo più nuovi contagi
                break;
            } else {
                last_infected++;
                time_vector[last_infected] = ts;
                startInfection.push_back(ts + L[last_infected]);
                endInfection.push_back(ts + L[last_infected] + I[last_infected]);
                //teoricamente il comando seguente è già vero
                integral_value = Q[last_infected];
                has_already_been_exposed[last_infected] = true;
            }
        } else {
            last_infected++;
        }
    }


    return time_vector;

}

