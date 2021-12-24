//
//  Copyright (c) 2021
//  Erick Suzart Souza
//
//  Permission to use, copy, modify, distribute and sell this software
//  and its documentation for any purpose is hereby granted without fee,
//  provided that the above copyright notice appear in all copies and
//  that both that copyright notice and this permission notice appear
//  in supporting documentation.  The authors make no representations
//  about the suitability of this software for any purpose.
//  It is provided "as is" without express or implied warranty.
//

#include <iostream>
#include <cmath>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "storage_adaptors.hpp"

using boost::numeric::ublas::make_matrix_from_pointer;
using boost::numeric::ublas::matrix;
using boost::numeric::ublas::vector;

/**
 * @brief Retorna o coeficiente de atividade de cada componente
 *
 * @param Alpha Parâmetro não aleatório de cada componente
 * @param EI Parâmetros de energia de interação de cada componente
 * @param T Temperatura em Kelvin
 * @param M Fração molar de cada componente
 * @return vector<double> Coeficiente de ativação para cada componente
 */
vector<double>
nrtl(
    matrix<double> Alpha,
    matrix<double> EI,
    double T,
    vector<double> M)
{
    matrix<double> Tau = EI / T;
    matrix<double> G = prod(-Alpha, Tau);
    // aplicar exponenciação em todos os elementos
    std::transform(G.data().begin(), G.data().end(), G.data().begin(), exp);
    int nElementos = M.size();
    // vetor dos resultados iniciado em zero
    vector<double> resultado(nElementos, 0.0);

    for (int i = 0; i < nElementos; i++)
    {
        // container para somatório
        double soma = 0;

        // inner_prod = Produto interno de dois vetores : O resultado é um número
        // element_prod = Produto elemento-elemento de dois vetores : O resultado é um vetor
        for (int j = 0; j < nElementos; j++)
            soma += M[j] * G(i, j) / inner_prod(column(G, j), M) * (Tau(i, j) - (inner_prod(element_prod(M, column(Tau, j)), column(G, j)) / inner_prod(column(G, j), M)));

        resultado[i] = inner_prod(element_prod(column(Tau, i), column(G, i)), M) / inner_prod(column(G, i), M) + soma;
    }

    // aplicar exponenciação em todos os elementos
    std::transform(resultado.data().begin(), resultado.data().end(), resultado.data().begin(), exp);

    return resultado;
}

static const double Gij_data[3][3] = {
    {0.00, 927.63, 1208.21},
    {-55.35, 0.00, 26.56},
    {141.33, 200.28, 0.00}};

static const double Aij_data[3][3] = {
    {0.0000, 0.2875, 0.2000},
    {0.2875, 0.0000, 0.3700},
    {0.200, 0.3700, 0.0000}};

static const double LLE_organic_data[9][3] = {
    {0.161, 0.000, 0.839},
    {0.200, 0.032, 0.768},
    {0.226, 0.058, 0.716},
    {0.276, 0.089, 0.635},
    {0.330, 0.120, 0.550},
    {0.352, 0.130, 0.518},
    {0.389, 0.143, 0.468},
    {0.490, 0.164, 0.346},
    {0.591, 0.160, 0.249}};

static const double LLE_aqueous_data[9][3] = {
    {0.985, 0.000, 0.015},
    {0.979, 0.006, 0.015},
    {0.976, 0.008, 0.016},
    {0.970, 0.016, 0.014},
    {0.960, 0.023, 0.017},
    {0.955, 0.026, 0.019},
    {0.952, 0.030, 0.018},
    {0.941, 0.039, 0.020},
    {0.917, 0.052, 0.031}};

static const double T = 308.15;

int main()
{
    matrix<double> paramG(3, 3), paramA(3, 3), MolFracOrganic(9, 3), MolFracAqueous(9, 3);

    paramG = make_matrix_from_pointer(3, 3, &(Gij_data[0][0]));
    paramA = make_matrix_from_pointer(3, 3, &(Aij_data[0][0]));
    MolFracOrganic = make_matrix_from_pointer(9, 3, &(LLE_organic_data[0][0]));
    MolFracAqueous = make_matrix_from_pointer(9, 3, &(LLE_aqueous_data[0][0]));

    std::cout << "Energia de ativação em meio aquoso para cada fração molar:" << std::endl;

    for (int i = 0; i < 9; i++)
    {
        std::cout
            << nrtl(paramA, paramG, T, row(MolFracAqueous, i))
            << "\tFracao molar: "
            << row(MolFracAqueous, i)
            << std::endl;
    }
    std::cout << "------------------------------------------------------------" << std::endl;

    std::cout << "Energia de ativação em meio orgânico para cada fração molar:" << std::endl;

    for (int i = 0; i < 9; i++)
    {
        std::cout
            << nrtl(paramA, paramG, T, row(MolFracOrganic, i))
            << "\tFracao molar: "
            << row(MolFracAqueous, i)
            << std::endl;
    }
}
