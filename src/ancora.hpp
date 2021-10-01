#ifndef ANCORA
#define ANCORA

#include "traits.hpp"
#include <vector>
#include <algorithm>

class ancora
{
private:
    cd::matrixptr data;
    double n_cubotti;
    double larghezza = 0;
    double altezza = 0;
    double larghezza_cubo = 0;
    double altezza_cubo = 0;
    cd::vector center{0,0};

    cd::vector cubotti() 
    {
        unsigned int n = data->rows();
        center(0) = (data->col(0)).minCoeff();
        center(1) = (data->col(1)).minCoeff();

        larghezza = (data->col(0)).maxCoeff() - center(0);
        altezza = (data->col(1)).maxCoeff() - center(1);

        larghezza_cubo = larghezza/n_cubotti;
        altezza_cubo = altezza/n_cubotti;



        cd::vector result(n);
        for (unsigned int i=0; i<n; ++i)
        {
            cd::vector coordinates = data->row(i);
            result(i) = ceil((coordinates(0)-center(0))/larghezza_cubo) + n_cubotti*floor((coordinates(1)-center(1))/altezza_cubo);
        }

        return result;
    }

public:
    ancora(const cd::matrixptr &data_, const double h_): data(data_), n_cubotti(h_){};

    const cd::matrix find_anchorpoints()
    {
        std::vector<unsigned int> positions;
        unsigned int n = data->rows();

        cd::vector results = cubotti();

        for (unsigned int i=0; i<n; ++i)
        {
            unsigned int pos = results(i);
            if (std::find(positions.begin(), positions.end(), pos) == positions.end())
                positions.push_back(pos);
        }


        cd::matrix anchorpos(positions.size(), data->cols());

        for (unsigned int i=0; i<anchorpos.rows(); ++i)
        {
            unsigned int I = positions[i];
            anchorpos(i,0) = center(0) + (I - floor(I/n_cubotti)*n_cubotti)*larghezza_cubo - larghezza_cubo/2;
            anchorpos(i,1) = center(1) + ceil(I/n_cubotti)*altezza_cubo - altezza_cubo/2;
        }

        return anchorpos;
    }
};

#endif //ANCORA