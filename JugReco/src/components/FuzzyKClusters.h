#pragma once

#include <Eigen/Dense>

namespace fkc {

class KMeans
{
public:
    KMeans();
    virtual ~KMeans();

    virtual Eigen::MatrixXd Fit(const Eigen::MatrixXd &data, int k,
                                double q = 2.0, double epsilon = 1e-4, int max_iters = 1000);

    int NIters() const { return n_iters; }
    double Variance() const { return variance; }
    const Eigen::MatrixXd &GetDistances() const { return dists; }
    const Eigen::MatrixXd &GetMemberships() const { return mems; }

    Eigen::MatrixXd &GetDistances() { return dists; }
    Eigen::MatrixXd &GetMemberships() { return mems; }

protected:
    virtual Eigen::MatrixXd Initialize(const Eigen::MatrixXd &data, int k, double q);
    virtual void Distances(const Eigen::MatrixXd &centroids, const Eigen::MatrixXd &data);
    virtual void Memberships(double q);
    virtual void FormClusters(Eigen::MatrixXd &clusters, const Eigen::MatrixXd &data, double q);

protected:
    int n_iters;
    double variance;
    Eigen::MatrixXd dists, mems;
};

class KRings : public KMeans
{
public:
    KRings();
    ~KRings();

    virtual Eigen::MatrixXd Fit(const Eigen::MatrixXd &data, int k,
                                double q = 2.0, double epsilon = 1e-4, int max_iters = 1000);

protected:
    virtual Eigen::MatrixXd Initialize(const Eigen::MatrixXd &data, int k, double q);
    virtual void Distances(const Eigen::MatrixXd &centroids, const Eigen::MatrixXd &data);
    virtual void FormClusters(Eigen::MatrixXd &clusters, const Eigen::MatrixXd &data, double q);
    virtual void FormRadii(Eigen::MatrixXd &clusters, double g);

protected:
    Eigen::MatrixXd dists_euc;
};

} // namespace fkc

