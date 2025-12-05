// shampine_wrapper.cpp
// Minimal C++ wrapper around the Shampine-Gordon "ode" solver (C port).
// Requires: ode.c / ode.h compiled & linked into the project (Burkardt/netlib sources).

#include <iostream>
#include <vector>
#include <functional>
#include <cmath>
#include <cstring>

// Forward declare the C API from ode.h (the public driver)
// The real ode.h from Burkardt declares:
// void ode(void f(double t, double y[], double yp[]), int neqn,
//          double y[], double *t, double tout, double relerr,
//          double abserr, int *iflag, double work[], int iwork[] );
extern "C" {
    void ode( void f(double t, double y[], double yp[]),
              int neqn,
              double y[],
              double *t,
              double tout,
              double relerr,
              double abserr,
              int *iflag,
              double work[],
              int iwork[] );
}

// A thin adapter to let C++ call the C-style callback.
struct CallbackData {
    std::function<void(double, const std::vector<double>&, std::vector<double>&)> func;
    int n;
};
static CallbackData *g_cb = nullptr;

// C-style bridge used by ode()
extern "C" void c_f(double t, double y[], double yp[]) {
    if (!g_cb) return;
    std::vector<double> Y(g_cb->n);
    for (int i = 0; i < g_cb->n; ++i) Y[i] = y[i];
    std::vector<double> Yp(g_cb->n);
    g_cb->func(t, Y, Yp);
    for (int i = 0; i < g_cb->n; ++i) yp[i] = Yp[i];
}

// Simple high-level solver interface
class ShampineSolver {
public:
    using RHS = std::function<void(double, const std::vector<double>&, std::vector<double>&)>;

    ShampineSolver(int neq) : neq_(neq) {
        work_.resize(neq_*16 + 200); // heuristic workspace (matches ode's expectations)
        iwork_.resize(200);
    }

    // Integrate from t0 to t1. y is both input and output (size must be neq).
    // relerr/abserr are tolerances (typical: 1e-6).
    bool integrate(RHS f, double &t, double tout, std::vector<double> &y,
                   double relerr = 1e-6, double abserr = 1e-10)
    {
        if ((int)y.size() != neq_) {
            throw std::runtime_error("y size mismatch");
        }
        // Set global callback data
        CallbackData cb;
        cb.func = f;
        cb.n = neq_;
        g_cb = &cb;

        // The C API expects raw arrays
        int iflag = 1; // normal start (1 means start/integrate)
        // convert vectors to C arrays
        // Note: work and iwork must be large enough per ode's documentation.
        // We sized them heuristically; if ode complains, increase sizes.
        // Call the solver
        std::vector<double> yarr = y; // local copy (ode writes directly)
        double tcopy = t;

        ode(c_f, neq_, yarr.data(), &tcopy, tout, relerr, abserr, &iflag,
            work_.data(), iwork_.data());

        // iflag==2 means normal return (integration reached tout)
        bool ok = (iflag == 2 || iflag == -2);
        if (ok) {
            y = yarr;
            t = tcopy;
        } else {
            std::cerr << "ode returned iflag=" << iflag << " (integration issue)\n";
            // Even on failure, ode may return current y and t
            y = yarr;
            t = tcopy;
        }

        g_cb = nullptr;
        return ok;
    }

private:
    int neq_;
    std::vector<double> work_;
    std::vector<int> iwork_;
};

//////////////////////
// Example: simple harmonic oscillator
// y0' = y1
// y1' = -omega^2 * y0
int main() {
    const int neq = 2;
    ShampineSolver solver(neq);

    double t = 0.0;
    double tout = 10.0;
    std::vector<double> y = { 1.0, 0.0 }; // initial: x=1, v=0
    double omega = 1.0;

    auto rhs = [omega](double tt, const std::vector<double>& Y, std::vector<double>& Yp) {
        Yp[0] = Y[1];
        Yp[1] = -omega*omega * Y[0];
    };

    bool ok = solver.integrate(rhs, t, tout, y, 1e-8, 1e-10);
    if (ok) {
        std::cout << "Reached t=" << t << "\n";
        std::cout << "y[0]=" << y[0] << "  y[1]=" << y[1] << "\n";
    } else {
        std::cout << "Integration failed or stopped early. t=" << t << "\n";
    }
    return 0;
}
