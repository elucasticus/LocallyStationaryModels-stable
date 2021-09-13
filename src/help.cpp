#include "help.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace cd;

void stampante::stampamatrice(const cd::matrixptr &p, const std::string &file_id, const unsigned int n)
{
    std::ofstream myfile(file_id);
    for (unsigned int i = 0; i < p->rows(); ++i)
    {
        for (unsigned int j = 0; j < p->cols(); ++j)
        {
            myfile << std::left << std::setw(n) << std::setfill(' ') << p->operator()(i, j);
        }
        myfile << std::endl;
    }
    myfile.close();
}

void stampante::stampamatrice(const cd::matrixIptr &p, const std::string &file_id, const unsigned int n)
{
    std::ofstream myfile(file_id);
    for (unsigned int i = 0; i < p->rows(); ++i)
    {
        for (unsigned int j = 0; j < p->cols(); ++j)
        {
            myfile << std::left << std::setw(n) << std::setfill(' ') << p->operator()(i, j);
        }
        myfile << std::endl;
    }
    myfile.close();
}

void stampante::stampamatriceperimport(const cd::matrixptr &p, const std::string &file_id)
{
    std::ofstream myfile(file_id);
    for (unsigned int i = 0; i < p->rows(); ++i)
    {
        for (unsigned int j = 0; j < p->cols(); ++j)
        {
            myfile << p->operator()(i, j) << " ";
        }
        myfile << std::endl;
    }
    myfile.close();
}

void stampante::stampamatriceperimport(const cd::matrixIptr &p, const std::string &file_id)
{
    std::ofstream myfile(file_id);
    for (unsigned int i = 0; i < p->rows(); ++i)
    { 
        for (unsigned int j = 0; j < p->cols(); ++j)
        {
            myfile << p->operator()(i, j) << " ";
        }
        myfile << std::endl;
    }
    myfile.close();
}


void stampante::stampavettore(const vectorptr &p, const std::string &file_id)
{
    std::ofstream myfile(file_id);
    for (unsigned int i = 0; i < p->rows(); ++i)
    {
        myfile << p->operator[](i) << "\n";
    }
    myfile.close();
}

cd::matrixptr stampante::caricamatrice(const std::string &file_id)
{
    std::ifstream myfile(file_id);
	std::string line;
	matrixptr d = std::make_shared<matrix>();

	std::getline(myfile, line);
	std::stringstream ss(line);

	unsigned int n = 0;
	ss >> n;
	d->resize(n, 2);

	unsigned int i = 0;
	while(std::getline(myfile, line))
	{
		std::stringstream ss(line);
		unsigned int j = 0;
		std::string element;
		while (std::getline(ss, element, ','))
		{
			d->operator()(i, j) = std::stod(element);
			++j;
		}
		++i;
	}

	myfile.close();
	return d;
}

cd::vectorptr stampante::caricavettore(const std::string &file_id)
{
    std::ifstream myfile(file_id);
	std::string line;
	vectorptr d = std::make_shared<vector>();

	std::getline(myfile, line);
	std::stringstream ss(line);

	unsigned int n = 0;
	ss >> n;
	d->resize(n, 1);

	unsigned int i = 0;
	while(std::getline(myfile, line))
	{
		std::stringstream ss(line);
		unsigned int j = 0;
		std::string element;
		while (std::getline(ss, element, ','))
		{
			d->operator()(i, j) = std::stod(element);
			++j;
		}
		++i;
	}

	myfile.close();
	return d;
}

void stampante::stampavofvs(const cd::vofvsptr &p, const std::string &file_id, const unsigned int n)
{
    std::ofstream myfile(file_id);
    for (unsigned int i = 0; i < p->size(); ++i)
    {
        for (unsigned int j = 0; j < p->operator[](i).size(); ++j)
        {
            myfile << std::left << std::setw(n) << std::setfill(' ') << p->operator[](i)[j];
        }
        myfile << std::endl;
    }
    myfile.close();
}


void stampante::pikapika()
{
    std::ifstream myfile("pikachu.csv");
	std::string line;

	while(std::getline(myfile, line))
	{
		std::cout << line << std::endl;
	}

	myfile.close();
}