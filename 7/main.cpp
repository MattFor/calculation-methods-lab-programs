#include <hprint>

enum class Method
{
    Jacobi,
    GaussSeidel,
    SOR
};

enum class StopCriterion
{
    UpdateNorm,   ///< ||x(k+1) - x(k)||_inf
    ResidualNorm, ///< ||A x(k+1) - b||_inf
    MaxIterations
};

template <std::size_t N>
using Vec = std::array<double, N>;

template <std::size_t N>
using Mat = std::array<std::array<double, N>, N>;


/**
 * Odejmowanie dwóch wektorów o tym samym rozmiarze
 *
 * @tparam N Rozmiar wektora
 * @param a Pierwszy wektor
 * @param b Drugi wektor
 *
 * @return Różnica a - b
 */
template <std::size_t N>
[[nodiscard]] Vec<N> sub(const Vec<N>& a, const Vec<N>& b)
{
    Vec<N> out{};

    for (std::size_t i = 0; i < N; ++i)
    {
        out[i] = a[i] - b[i];
    }

    return out;
}


/**
 * Norma nieskończona wektora
 * Zwraca największą wartość bezwzględną ze wszystkich składowych
 *
 * @tparam N Rozmiar wektora
 * @param v Wektor wejściowy
 *
 * @return Norma nieskończona wektora
 */
template <std::size_t N>
[[nodiscard]] double norm_inf(const Vec<N>& v)
{
    double m = 0.0;
    for (const double x : v)
    {
        m = std::max(m, std::abs(x));
    }

    return m;
}

/**
 * Solver iteracyjny dla układu Ax = b
 * Klasa przechowuje macierz A oraz wektor b
 * Wykonuje iteracyjnie metody Jacobiego Gaussa-Seidla i SOR
 *
 * @tparam N Rozmiar układu równań
 */
template <std::size_t N>
class IterativeSolver
{
    Mat<N> A_{};
    Vec<N> b_{};

    /**
     * Liczy residuum układu dla aktualnego przybliżenia x
     * Residuum ma postać Ax - b i pokazuje jak daleko jest się od rozwiązania
     *
     * @param x Aktualne przybliżenie rozwiązania
     *
     * @return Wektor residuum
     */
    [[nodiscard]] Vec<N> residual(const Vec<N>& x) const
    {
        Vec<N> r{};
        for (std::size_t i = 0; i < N; ++i)
        {
            double s = 0.0;
            for (std::size_t j = 0; j < N; ++j)
            {
                s += A_[i][j] * x[j];
            }

            r[i] = s - b_[i];
        }
        return r;
    }

    /**
     * Wykonuje jeden krok iteracji wybraną metodą
     * Jacobi liczy wszystkie składowe z poprzedniego przybliżenia
     * Gauss Seidel używa od razu nowych wartości gdy są już policzone
     * SOR robi to samo co Gauss-Seidel ale dodatkowo miesza nowe i stare przybliżenia parametrem omega
     *
     * @param method Wybrana metoda iteracyjna
     * @param x_prev Poprzednie przybliżenie
     * @param omega Parametr relaksacji używany tylko w metodzie SOR
     *
     * @return Następne przybliżenie rozwiązania
     */
    [[nodiscard]] Vec<N> step(const Method method, const Vec<N>& x_prev, const double omega) const
    {
        Vec<N> x = x_prev;

        switch (method)
        {
            case Method::Jacobi:
            {
                for (std::size_t i = 0; i < N; ++i)
                {
                    double sum = 0.0;
                    for (std::size_t j = 0; j < N; ++j)
                    {
                        if (j != i)
                        {
                            sum += A_[i][j] * x_prev[j];
                        }
                    }

                    x[i] = ( b_[i] - sum ) / A_[i][i];
                }
                break;
            }

            case Method::GaussSeidel:
            {
                for (std::size_t i = 0; i < N; ++i)
                {
                    double sum_before = 0.0;
                    double sum_after  = 0.0;

                    for (std::size_t j = 0; j < i; ++j)
                    {
                        sum_before += A_[i][j] * x[j];
                    }

                    for (std::size_t j = i + 1; j < N; ++j)
                    {
                        sum_after += A_[i][j] * x_prev[j];
                    }

                    x[i] = ( b_[i] - sum_before - sum_after ) / A_[i][i];
                }
                break;
            }

            case Method::SOR:
            {
                for (std::size_t i = 0; i < N; ++i)
                {
                    double sum_before = 0.0;
                    double sum_after  = 0.0;

                    for (std::size_t j = 0; j < i; ++j)
                    {
                        sum_before += A_[i][j] * x[j];
                    }

                    for (std::size_t j = i + 1; j < N; ++j)
                    {
                        sum_after += A_[i][j] * x_prev[j];
                    }

                    const double gs = ( b_[i] - sum_before - sum_after ) / A_[i][i];
                    x[i]            = ( 1.0 - omega ) * x_prev[i] + omega * gs;
                }

                break;
            }
        }

        return x;
    }

    static const char* method_name(const Method m)
    {
        switch (m)
        {
            case Method::Jacobi:
            {
                return "Jacobi";
            }

            case Method::GaussSeidel:
            {
                return "Gauss-Seidel";
            }

            case Method::SOR:
            {
                return "SOR";
            }
        }

        return "?";
    }

    /**
     *
     * @param s
     * @return
     */
    static const char* stop_name(const StopCriterion s)
    {
        switch (s)
        {
            case StopCriterion::UpdateNorm:
            {
                return "||Δx||∞";
            }

            case StopCriterion::ResidualNorm:
            {
                return "||r||∞";
            }

            case StopCriterion::MaxIterations:
            {
                return "max_iter";
            }
        }

        return "?";
    }

public:
    struct Options
    {
        std::size_t   max_iter = 1000;
        double        tol      = 1e-10;                     ///< Dokładność wymaganej zbieżności
        double        omega    = 0.5;                       ///< Parametr relaksacji dla metody SOR
        StopCriterion stop     = StopCriterion::UpdateNorm; ///< Wybrane kryterium zakończenia iteracji
        bool          verbose  = true;
    };

    struct Result
    {
        Vec<N>      x{};          ///< Ostateczne przybliżenie rozwiązania
        std::size_t iterations{}; ///< Liczba wykonanych iteracji
        bool        converged{};  ///< Czy metoda osiągnęła zbieżność
    };

    /**
     * Konstruuje solver z macierzy A i wektora b
     * Dodatkowo sprawdza, czy na przekątnej nie ma zera
     * bo wtedy metody iteracyjne nie są bezpieczne
     *
     * @param A Macierz współczynników układu
     * @param b Wektor wyrazów wolnych
     */
    IterativeSolver(Mat<N> A, Vec<N> b)
        : A_(std::move(A)), b_(std::move(b))
    {
        for (std::size_t i = 0; i < N; ++i)
        {
            if (std::abs(A_[i][i]) < 1e-15)
            {
                throw std::runtime_error("Zero on diagonal: iterative methods cannot be used safely.");
            }
        }
    }


    /**
     * Rozwiązuje układ równań wybraną metodą iteracyjną
     * W każdej iteracji liczy nowe przybliżenie, residuum i estymator błędu
     * Iteracje kończą się po spełnieniu wybranego kryterium stopu
     *
     * @param method Wybrana metoda iteracyjna
     * @param x0 Przybliżenie początkowe
     * @param opt Parametry obliczeń
     * @param exact Dokładne rozwiązanie używane tylko do porównania, jeśli jest podane
     *
     * @return Struktura zawierająca wynik i liczbę iteracji
     */
    [[nodiscard]] Result solve(const Method method, Vec<N> x0, const Options& opt, std::optional<Vec<N>> exact = std::nullopt) const
    {
        if (opt.verbose)
        {
            hp::print("Iterative solver");
            hp::printf("Method: {}, stop = {}", method_name(method), stop_name(opt.stop));
            hp::printf("tol = {:.3e}, max_iter = {}, omega = {:.6f}", opt.tol, opt.max_iter, opt.omega);
        }

        Vec<N> x_prev = x0;
        Vec<N> x_next = x0;

        for (std::size_t k = 1; k <= opt.max_iter; ++k)
        {
            x_next = step(method, x_prev, opt.omega);

            const double dx = norm_inf(sub(x_next, x_prev));
            const double r  = norm_inf(residual(x_next));

            double     err        = 0.0;
            const bool have_exact = exact.has_value();
            if (have_exact)
            {
                err = norm_inf(sub(x_next, *exact));
            }

            if (opt.verbose)
            {
                std::print("k = {:>3} | x = [", k);
                for (std::size_t i = 0; i < N; ++i)
                {
                    if (i != 0)
                    {
                        std::print(", ");
                    }

                    std::print("{:>15.16e}", x_next[i]);
                }

                std::print("]  |  ||Δx||∞ = {:>12.5e}  ||r||∞ = {:>12.5e}", dx, r);

                if (have_exact)
                {
                    std::print("  ||e||∞ = {:>12.5e}", err);
                }

                std::print("\n");
            }

            bool stop_now = false;
            switch (opt.stop)
            {
                case StopCriterion::UpdateNorm:
                {
                    stop_now = dx <= opt.tol;
                    break;
                }

                case StopCriterion::ResidualNorm:
                {
                    stop_now = r <= opt.tol;
                    break;
                }

                case StopCriterion::MaxIterations:
                {
                    stop_now = k == opt.max_iter;
                    break;
                }

                default:
                {
                    throw std::logic_error("Unknown stop criterion.");
                }
            }

            if (stop_now)
            {
                return Result{ x_next, k, true };
            }

            x_prev = x_next;
        }

        return Result{ x_next, opt.max_iter, false };
    }
};

int main()
{
    constexpr Mat<5> A{ { { { 50, 5, 4, 3, 2 } }, { { 1, 40, 1, 2, 3 } }, { { 4, 5, 30, -5, -4 } }, { { -3, -2, -1, 20, 0 } }, { { 1, 2, 3, 4, 30 } } } };

    constexpr Vec<5> b{ 140, 67, 62, 89, 153 };
    constexpr Vec<5> x0{ 6, 6, 6, 6, 6 };

    constexpr Vec<5> x_exact{ 2, 1, 3, 5, 4 };

    const IterativeSolver solver(A, b);

    IterativeSolver<5>::Options opt{};
    opt.tol      = 1e-16;
    opt.max_iter = 1000;
    opt.omega    = 0.5;
    opt.verbose  = true;

    const auto jacobi_res = solver.solve(Method::Jacobi, x0, opt, x_exact);
    const auto gs_res     = solver.solve(Method::GaussSeidel, x0, opt, x_exact);
    const auto sor_res    = solver.solve(Method::SOR, x0, opt, x_exact);

    std::print("\nSummary:\n");
    std::print("Jacobi        : iter = {}, converged = {}\n", jacobi_res.iterations, jacobi_res.converged);
    std::print("Gauss-Seidel  : iter = {}, converged = {}\n", gs_res.iterations, gs_res.converged);
    std::print("SOR (ω=0.5)   : iter = {}, converged = {}\n", sor_res.iterations, sor_res.converged);
}
