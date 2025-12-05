#include <iostream>
#include <vector>
#include <functional>
#include <cmath>
#include <cassert>
#include <iomanip>
#include <limits>

// TODO put this in a header and use Eigen Matrices (or maybe just heyokas expressions? combination of both??)
using Vec = std::vector<long double>;

// ---------------- vector utilities ----------------
static Vec add(const Vec &a, const Vec &b){
    assert(a.size() == b.size());
    Vec r(a.size());
    for(size_t i=0;i<a.size();++i) r[i] = a[i] + b[i];
    return r;
}
static Vec sub(const Vec &a, const Vec &b){
    assert(a.size() == b.size());
    Vec r(a.size());
    for(size_t i=0;i<a.size();++i) r[i] = a[i] - b[i];
    return r;
}
static Vec scale(const Vec &a, long double s){
    Vec r(a.size());
    for(size_t i=0;i<a.size();++i) r[i] = a[i]*s;
    return r;
}

static void rk4_step(const std::function<Vec(const Vec&, long double)> &force,
                     std::vector<long double> &y, long double t, long double h){
    int d = (int)(y.size()/2);
    auto eval = [&](const std::vector<long double> &yy, long double tt){
        std::vector<long double> dy(2*d);
        for(int i=0;i<d;i++) dy[i] = yy[d+i];  // x' = v
        Vec x(d);
        for(int i=0;i<d;i++) x[i] = yy[i];
        Vec a = force(x, tt);
        for(int i=0;i<d;i++) dy[d+i] = a[i];  // v' = f(x)
        return dy;
    };

    std::vector<long double> y0 = y;
    auto k1 = eval(y0, t);
    std::vector<long double> ytmp(2*d);
    for(int i=0;i<2*d;i++) ytmp[i] = y0[i] + 0.5L*h*k1[i];
    auto k2 = eval(ytmp, t + 0.5L*h);
    for(int i=0;i<2*d;i++) ytmp[i] = y0[i] + 0.5L*h*k2[i];
    auto k3 = eval(ytmp, t + 0.5L*h);
    for(int i=0;i<2*d;i++) ytmp[i] = y0[i] + h*k3[i];
    auto k4 = eval(ytmp, t + h);
    for(int i=0;i<2*d;i++)
        y[i] = y0[i] + (h/6.0L) * (k1[i] + 2.0L*k2[i] + 2.0L*k3[i] + k4[i]);
}

class StomerCowell {
public:
    StomerCowell(std::function<Vec(const Vec&, long double)> force, int order_p = 4)
        : f(force)
    {
        if(order_p < 2 || (order_p % 2) != 0)
            throw std::invalid_argument("order must be even and >=2");
        p = order_p;
        r = p/2;
        compute_beta();
    }

    std::vector<Vec> integrate(const Vec &x0, const Vec &v0, long double t0, long double t1, long double h){
        int dim = (int)x0.size();
        int steps = (int)std::llround((t1 - t0)/h);
        if(steps < 1) return {x0};

        // --------------- bootstrap ----------------
        std::vector<Vec> pos_buf(r+1); // last r+1 positions
        std::vector<Vec> traj; traj.reserve(steps+1);

        std::vector<long double> y(2*dim); // [x,v]
        for(int i=0;i<dim;i++){ y[i] = x0[i]; y[dim+i] = v0[i]; }

        long double t = t0;
        // x0
        Vec xcur(dim);
        for(int i=0;i<dim;i++) xcur[i] = y[i];
        pos_buf[0] = xcur;
        traj.push_back(xcur);

        // RK4 bootstrap for x1..xr
        for(int k=1;k<=r;k++){
            rk4_step(f, y, t, h);
            t += h;
            Vec xk(dim);
            for(int i=0;i<dim;i++) xk[i] = y[i];
            pos_buf[k] = xk;
            traj.push_back(xk);
        }

        long double h2 = h*h;

        // --------------- main SC loop ----------------
        for(int step = r; step < steps; ++step){
            Vec &x_n   = pos_buf[r];
            Vec &x_nm1 = pos_buf[r-1];

            Vec x_np1 = sub(scale(x_n, 2.0L), x_nm1);

            // compute symmetric force sum
            for(int j=0;j<r;j++){
                Vec f_n_minus_j = f(pos_buf[r - j], t - j*h);
                Vec f_nm1_plus_j;
                if(j==0) f_nm1_plus_j = f(pos_buf[r - 1], t - h); // j=0 special
                else {
                    // use previous positions to avoid out-of-bounds
                    f_nm1_plus_j = f(pos_buf[r - 1 + (j>1 ? 1 : 0)], t + (j-1)*h);
                }
                Vec term = add(f_n_minus_j, f_nm1_plus_j);
                x_np1 = add(x_np1, scale(term, h2*beta[j]));
            }

            // shift buffer left
            for(int k=0;k<r;k++)
                pos_buf[k] = pos_buf[k+1];
            pos_buf[r] = x_np1;

            t += h;
            traj.push_back(x_np1);
        }

        return traj;
    }

    std::vector<long double> get_beta() const { return beta; }

private:
    std::function<Vec(const Vec&, long double)> f;
    int p, r;
    std::vector<long double> beta;

    void compute_beta(){
        std::vector<std::vector<long double>> A(r, std::vector<long double>(r));
        std::vector<long double> b(r);
        for(int k=0;k<r;k++){
            for(int j=0;j<r;j++){
                long double v1 = pow((long double)j, 2LL*k);
                long double v2 = pow((long double)(j+1), 2LL*k);
                A[k][j] = v1 + v2;
            }
            long double denom = 1.0L;
            for(int t=1;t<=2*k+2;t++) denom *= (long double)t;
            b[k] = 1.0L / denom;
        }

        beta.resize(r);
        // naive Gaussian elimination
        for(int k=0;k<r;k++){
            long double pivot = A[k][k];
            for(int j=k;j<r;j++) A[k][j] /= pivot;
            b[k] /= pivot;
            for(int i=k+1;i<r;i++){
                long double factor = A[i][k];
                for(int j=k;j<r;j++) A[i][j] -= factor*A[k][j];
                b[i] -= factor*b[k];
            }
        }
        // back-substitution
        for(int i=r-1;i>=0;i--){
            beta[i] = b[i];
            for(int j=i+1;j<r;j++) beta[i] -= A[i][j]*beta[j];
        }
    }
};

int main(){
    long double omega = 1.5L;
    auto force = [&](const Vec &x, long double t) -> Vec{
        Vec a(1); a[0] = -omega*omega*x[0]; return a;
    };

    Vec x0 = {1.0L};
    Vec v0 = {0.0L};
    long double t0 = 0.0L;
    long double t1 = 10.0L;
    long double h = 0.001L;
    int order_p = 10;

    StomerCowell sc(force, order_p);

    auto beta = sc.get_beta();
    std::cout << "Computed beta coefficients (r=" << order_p/2 << "):\n";
    for(size_t i=0;i<beta.size();i++)
        std::cout << "beta["<<i<<"] = " << (double)beta[i] << "\n";

    auto traj = sc.integrate(x0, v0, t0, t1, h);

    // print every 1000th step
    std::cout << "\nTrajectory sample:\n";
    for(size_t n=0;n<traj.size();n+=1000)
        std::cout << t0 + n*h << " " << (double)traj[n][0] << "\n";

    return 0;
}