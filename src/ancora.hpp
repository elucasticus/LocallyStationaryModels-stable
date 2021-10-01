#ifdef ANCORA
#define ANCORA

#include "traits.hpp"
#include <vector>
#include <algorithm>

class ancora
{
private:
    cd::matrixptr data;
    double n_cubotti;

    cd::vector cubotti() const
    {
        unsigned int n = data->rows();
        double larghezza = (data->col(0)).maxCoeff() - (data->col(0)).minCoeff();
        double altezza = (data->col(1)).maxCoeff() - (data->col(1)).minCoeff();

        double larghezza_cubo = larghezza/n_cubotti;
        double altezza_cubo = altezza/n_cubotti;



        cd::vector result(n);
        for (unsigned int i=0; i<n; ++i)
        {
            cd::vector coordinates = data->row(i);
            result(i) = ceil(coordinates(0)/larghezza_cubo) + 2*n_cubotti*floor(coordinates(1)/altezza_cubo);
        }

        return result;
    }

public:
    ancora(const cd::matrixptr &data_, const double h_): data(data_), n_cubotti(n_cubotti);

    const cd::matrixptr find_anchorpoints() const
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

        double larghezza = (data->col(0)).maxCoeff() - (data->col(0)).minCoeff();
        double altezza = (data->col(1)).maxCoeff() - (data->col(1)).minCoeff();

        double larghezza_cubo = larghezza/n_cubotti;
        double altezza_cubo = altezza/n_cubotti;

        cd::matrixptr anchorpos(positions.size(), data->cols());
        for (unsigned int i=0; i<anchorpos->rows(); ++i)
        {
            unsigned int I = positions[i];
            anchorpos(i,0) = (I - 2*n_cubotti)*larghezza_cubo - larghezza_cubo/2;
            anchorpos(i,1) = floor(I/n_cubotti)*altezza_cubo - altezza_cubo/2;
        }

        return anchorpos;
    }
};

#endif //ANCORA