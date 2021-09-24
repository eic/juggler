/*  Fuzzy K Clustering Algorithms
 *
 *  Author: Chao Peng (ANL)
 *  Date: 10/14/2020
 *
 */

#include "FuzzyKClusters.h"
#include <exception>
#include <iostream>
#include <cmath>


using namespace fkc;
using namespace Eigen;


// =================================================================================================
//  KMeans Algorithm
// =================================================================================================

KMeans::KMeans()
: n_iters(0), variance(0.)
{
    // place holder
}

KMeans::~KMeans()
{
    // place holder
}

MatrixXd KMeans::Fit(const MatrixXd &data, int k, double q, double epsilon, int max_iters)
{
    auto res = Initialize(data, k, q);

    for (n_iters = 0; n_iters < max_iters; ++n_iters) {
        auto old_mems = mems;
        Distances(res, data);
        Memberships(q);
        FormClusters(res, data, q);

        if ((old_mems - mems).cwiseAbs().maxCoeff() < epsilon) {
            break;
        }
    }

    return res;
}

// initialize and guess the clusters
MatrixXd KMeans::Initialize(const MatrixXd &data, int k, double q)
{
    // resize matrices
    dists.resize(k, data.rows());

    // guess the cluster centers
    mems = MatrixXd::Random(k, data.rows());
    for (int j = 0; j < mems.cols(); ++j) {
        auto csum = mems.col(j).sum();
        for (int i = 0; i < mems.rows(); ++i) {
            mems(i, j) = mems(i, j)/csum;
        }
    }

    MatrixXd clusters(k, data.cols());
    FormClusters(clusters, data, q);
    return clusters;
}

// distance matrix (num_clusters, num_data)
void KMeans::Distances(const MatrixXd &centroids, const MatrixXd &data)
{
    for (int i = 0; i < centroids.rows(); ++i) {
        for (int j = 0; j < data.rows(); ++j) {
            dists(i, j) = (centroids.row(i) - data.row(j)).cwiseAbs2().sum();
        }
    }
}

// membership matrix (num_clusters, num_data)
void KMeans::Memberships(double q)
{
    // coeffcient-wise operation
    auto d = dists.array().pow(-1.0/(q - 1.0)).matrix();

    for (int j = 0; j < d.cols(); ++j) {
        auto dsum = d.col(j).sum();
        for (int i = 0; i < d.rows(); ++i) {
            mems(i, j) = d(i, j)/dsum;
        }
    }
}

// rebuild clusters
void KMeans::FormClusters(MatrixXd &clusters, const MatrixXd &data, double q)
{
    auto weights = mems.array().pow(q).matrix();
    for (int i = 0; i < clusters.rows(); ++i) {
        clusters.row(i) *= 0;
        for (int j = 0; j < data.rows(); ++j) {
            clusters.row(i) += data.row(j)*weights(i, j);
        }
        clusters.row(i) /= weights.row(i).sum();
    }
}


// =================================================================================================
//  KRings Algorithm, extended from KMeans
//  Reference:
//      [1] Y. H. Man and I. Gath,
//          "Detection and separation of ring-shaped clusters using fuzzy clustering,"
//          in IEEE Transactions on Pattern Analysis and Machine Intelligence,
//          vol. 16, no. 8, pp. 855-861, Aug. 1994, doi: 10.1109/34.308484.
// =================================================================================================

KRings::KRings()
: KMeans()
{
    // place holder
}

KRings::~KRings()
{
    // place holder
}

MatrixXd KRings::Fit(const MatrixXd &data, int k, double q, double epsilon, int max_iters)
{
    auto res = Initialize(data, k, q);

    for (n_iters = 0; n_iters < max_iters; ++n_iters) {
        auto old_mems = mems;
        Distances(res, data);
        Memberships(q);
        FormRadii(res, q);
        FormClusters(res, data, q);

        if ((old_mems - mems).cwiseAbs().maxCoeff() < epsilon) {
            break;
        }
    }

    return res;
}

// initialize and guess the clusters
MatrixXd KRings::Initialize(const MatrixXd &data, int k, double q)
{
    MatrixXd clusters(k, data.cols() + 1);
    auto centers = clusters.leftCols(data.cols());

    // call KMeans to help initialization
    KMeans fkm;
    centers = fkm.Fit(data, k, q, 1e-4, 5);

    dists.resize(k, data.rows());
    dists_euc = fkm.GetDistances().cwiseSqrt();
    mems = fkm.GetMemberships();
    FormRadii(clusters, q);
    return clusters;
}

// distance matrix (num_clusters, num_data)
void KRings::Distances(const MatrixXd &centroids, const MatrixXd &data)
{
    auto const centers = centroids.leftCols(centroids.cols() - 1);
    auto const radii = centroids.rightCols(1);

    for (int i = 0; i < centroids.rows(); ++i) {
        for (int j = 0; j < data.rows(); ++j) {
            dists_euc(i, j) = std::sqrt((centers.row(i) - data.row(j)).cwiseAbs2().sum());
            dists(i, j) = std::pow(dists_euc(i, j) - radii(i, 0), 2);
        }
    }
}

// rebuild clusters radii
void KRings::FormRadii(MatrixXd &clusters, double q)
{
    auto radii = clusters.rightCols(1);
    auto weights = mems.array().pow(q).matrix();

    for (int i = 0; i < weights.rows(); ++i) {
        radii(i, 0) = 0;
        for (int j = 0; j < weights.cols(); ++j) {
            radii(i, 0) += weights(i, j)*dists_euc(i, j);
        }
        radii(i, 0) /= weights.row(i).sum();
    }
}

// rebuild clusters centers
void KRings::FormClusters(MatrixXd &clusters, const MatrixXd &data, double q)
{
    auto centers = clusters.leftCols(data.cols());
    const auto &radii = clusters.rightCols(1);
    auto weights = mems.array().pow(q).matrix();

    for (int i = 0; i < weights.rows(); ++i) {
        MatrixXd icenter = centers.row(i);
        centers.row(i) *= 0;
        for (int j = 0; j < weights.cols(); ++j) {
            double scale = radii(i, 0)/dists_euc(i, j);
            centers.row(i) += weights(i, j)*(data.row(j) - (data.row(j) - icenter)*scale);
        }
        centers.row(i) /= weights.row(i).sum();
    }
}


// =================================================================================================
//  KEllipses Algorithm, extended from KRings
//  Reference:
//      [1] I. Gath and D. Hoory, Pattern Recognition Letters 16 (1995) 727-741,
//          https://doi.org/10.1016/0167-8655(95)00030-K.
// =================================================================================================
// @TODO

