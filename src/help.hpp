#ifndef HELP
#define HELP

#include "traits.hpp"
#include <string>


/**
 * \brief a struct to import and export data directly from cpp without using R
 * useful to debug and test the program
*/
struct stampante
{
    /**
     * \brief create a txt file containing the matrix/vector in a real and cool way(standard ones) or in a way comfortable for import in r(perimport ones)
     * \param p         a shared pointer to the matrix/vector to print out
     * \param file_id   the name of the file where to print *p
     * \param n         an integer needed to print cool matrices
    */
    static void stampamatrice(const cd::matrixptr &p, const std::string &file_id, const unsigned int n);
    static void stampamatrice(const cd::matrixIptr &p, const std::string &file_id, const unsigned int n);
    static void stampavofvs(const cd::vofvsptr &p, const std::string &file_id, const unsigned int n);
    static void stampamatriceperimport(const cd::matrixptr &p, const std::string &file_id);
    static void stampamatriceperimport(const cd::matrixIptr &p, const std::string &file_id);
    static void stampavettore(const cd::vectorptr &p, const std::string &file_id);

    /**
     * \brief pika pika pikachu
    */
    static void pikapika();

    /**
     * \brief import a matrix or a vector in cpp
     * \param file_id   the name the file to open
    */
    static cd::matrixptr caricamatrice(const std::string &file_id);
    static cd::vectorptr caricavettore(const std::string &file_id);
};


#endif //HELP