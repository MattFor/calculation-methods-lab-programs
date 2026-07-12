#include <print>
#include <cmath>
#include <numeric>

/**
 * max_iter: maksymalna liczba iteracji
 * verbose: czy wypisywać każdą iterację
 * tol_abs_x: tolerancja bezwzględna na zmianę przybliżenia
 * tol_resid: tolerancja na wartość funkcji w punkcie
 * tol_rel: tolerancja względna na zmianę przybliżenia
 */
struct Options
{
    int  max_iter = 200;
    bool verbose  = true;

    double tol_abs_x = 1e-12;
    double tol_resid = 1e-16;
    double tol_rel   = 1e-14;
};

/**
 * iter: numer iteracji
 * x: aktualne przybliżenie rozwiązania
 * resid: residuum, czyli |f(x)|
 * extra: dodatkowa informacja zależna od metody
 * err_est: estymator błędu, np. |x_{n+1} - x_n|
 *
 * Dla Picarda w extra trafia przybliżona pochodna g'(x)
 * Dla Newtona w extra trafia f'(x)
 * Dla siecznych w extra trafia mianownik f_cur - f_prev
 * Dla bisekcji extra jest nieużywane, więc NaN
 */
struct IterInfo
{
    int    iter{};
    double x{};
    double resid{};
    double extra{};
    double err_est{};
};

static void print_header(std::string_view method)
{
    std::println("=== {} ===\n{:>6}{:>20}{:>18}{:>18}{:>16}", method, "it", "x_n", "err_est", "resid |f(x)|", "extra");
}

static void print_iter(const IterInfo& it)
{
    std::print("{:>6}{:>20.15f}{:>18.10e}{:>18.10e}", it.iter, it.x, it.err_est, it.resid);

    if (std::isfinite(it.extra))
    {
        std::print("{:>18.10e}", it.extra);
    }
    else
    {
        std::print("{:>18}", "X");
    }

    std::println();
}

/**
 * Funkcja musi być w przedziale [a; b] i musi być zmiana znaku f(a) * f(b) <= 0
 * W każdej iteracji przedział [a; b] jest dzielony na pół, a następnie
 * wybierana jest ta połowa, w której nadal występuje zmiana znaku
 *
 * @tparam F Typ funkcji przyjmującej double i zwracającej double
 *
 * @param f Funkcja dla której szukane jest miejsce zerowe
 * @param a Lewy koniec przedziału
 * @param b Prawy koniec przedziału
 * @param opt Options
 *
 * @return std::optional<std::pair<double, int>>
 *
 * @note Estymator błędu to err_est = (b - a) / 2, czyli połowa długości aktualnego przedziału
 *
 * @warning Metoda wymaga, aby funkcja była ciągła na [a, b] oraz aby f(a) i f(b) miały różne znaki (lub jedno z nich było zerem)
 */
template <typename F>
[[nodiscard]] std::optional<std::pair<double, int>>
bisection(F f, double a, double b, const Options& opt)
{
    double fa = f(a);
    double fb = f(b);

    // Funkcja musi być poprawna
    if (!std::isfinite(fa) || !std::isfinite(fb))
    {
        std::println(stderr, "Bisekcja: nieprawidlowe wartosci funkcji na koncach przedzialu!");
        return std::nullopt;
    }

    // Jeśli na końcach przedziału to od razu jest pierwiastek, można zwrócić
    if (fa == 0.0)
    {
        return std::pair{a, 0};
    }

    if (fb == 0.0)
    {
        return std::pair{b, 0};
    }

    // Nie ma zmiany znaku więc bisekcja nie ma pewności działania
    if (fa * fb > 0.0)
    {
        std::println(stderr, "Bisekcja: brak zmiany znaku na [{}, {}]", a, b);
        return std::nullopt;
    }

    double left  = a;
    double right = b;

    print_header("Bisekcja");
    for (int k = 1; k <= opt.max_iter; ++k)
    {
        /**
         * mid: środek przedziału
         * fmid: wartość funkcji w środku
         * err_est: estymator błędu, czyli połowa długości przedziału
         */
        const double mid     = std::midpoint(left, right);
        const double fmid    = f(mid);
        const double err_est = 0.5 * ( right - left );

        // Tabela iteracji
        if (opt.verbose)
        {
            IterInfo it{k, mid, std::abs(fmid), std::numeric_limits<double>::quiet_NaN(), err_est};
            print_iter(it);
        }

        /**
         * 1. Czy mała zmiana |x_{n+1} - x_{n}| < tol_abs_x
         * Czyli czy kolejne iteracje się praktycznie nie zmieniają
         *
         * 2. Małe residuum |f(x_{n+1})| < tol_resid
         * Czyli czy punkt już praktycznie jest miejscem zerowym
         *
         * 3. Mała zmiana względna |x_{n+1} - x_{n}| / max(1, |x_{n+1}) < tol_rel
         * Gdy pierwiastki są duże to sama różnica bezwzgl może nie wystarczyć
         */

        // Jeśli mała zmiana i małe residuum to można zwrócić
        // Następnie osobno sprawdzane
        if (err_est < opt.tol_abs_x && std::abs(fmid) < opt.tol_resid)
        {
            return std::pair{mid, k};
        }

        if (!std::isfinite(fmid))
        {
            std::println(stderr, "Bisekcja: nieprawidlowa wartosc funkcji w srodku przedzialu");
            return std::nullopt;
        }

        // if (std::abs(fmid) < opt.tol_resid)
        // {
        //     return std::pair{mid, k};
        // }
        //
        // if (err_est / std::max(1.0, std::abs(mid)) < opt.tol_rel)
        // {
        //     return std::pair{mid, k};
        // }

        // Wybieranie nowego przedziału gdzie nadal jest zmiana znaku
        if (fa * fmid <= 0.0)
        {
            right = mid;
            fb    = fmid;
        }
        else
        {
            left = mid;
            fa   = fmid;
        }
    }

    std::println(stderr, "Bisekcja: std::max_iter osiągnięte");
    return std::nullopt;
}

/**
 * Najpierw funkcja musi być przekształcona do x = g(x)
 * Następnie trzeba stworzyć ciąg iteracyjny x_{n+1} = g(x_n)
 *   który zbiega do punktu stałego x*, który spełnia g(x*) = x* -> rozwiązanie f(x)=0.
 * Funkcja musi być kontrakycjna czyli musi istnieć stała X < 1 <=> |g'(x)| <= X
 * Dlatego obliczam pochodną g'(x)
 *
 * @tparam G Typ funkcji iteracyjnej g(x)
 * @tparam F Typ funkcji oryginalnej f(x), używanej do residuum |f(x)|
 *
 * @param g Funkcja iteracyjna przekształcona do postaci x = g(x)
 * @param f_resid Oryginalna funkcja f(x) żeby można było sprawdzać dokładność
 * @param x0 Punkt początkowy iteracji
 * @param opt Options
 *
 * @return std::optional<std::pair<double, int>>
 *
 * @note Estymator błędu to |x_{n+1} - x_n|
 *
 * @warning Metoda może być rozbieżna, jeśli |g'(x)| >= 1
 */
template <typename G, typename F>
[[nodiscard]] std::optional<std::pair<double, int>>
picard(G g, F f_resid, const double x0, const Options& opt)
{
    print_header("Picard (punkt-stały)");

    double x = x0;

    for (int k = 1; k <= opt.max_iter; ++k)
    {
        constexpr double h_base  = 1e-8;
        const double     xnext   = g(x);                     // Następne przybliżenie
        const double     err_est = std::abs(xnext - x);      // Zmiana między iteracjami
        const double     resid   = std::abs(f_resid(xnext)); // Residuum oryginalnego równania

        /**
         * Numeryczne wyliczenie tego czy iteracja jest kontrakcyjna
         * No bo jeśli nie zbiega poprawnie to po co?
         * Z wzoru g'(x) ≈ (g(x+h) − g(x−h)) / 2h
         *
         * P.S Wziąłem to z internetu, nie ma szans że to zapamiętam
         */
        const double h          = h_base * std::max(1.0, std::abs(xnext));
        const double gplus      = g(xnext + h);
        const double gminus     = g(xnext - h);
        const double gprime_est = ( gplus - gminus ) / ( 2.0 * h );

        if (opt.verbose)
        {
            IterInfo it{k, xnext, resid, gprime_est, err_est};
            print_iter(it);
        }

        /**
         * 1. Czy mała zmiana |x_{n+1} - x_{n}| < tol_abs_x
         * Czyli czy kolejne iteracje się praktycznie nie zmieniają
         *
         * 2. Małe residuum |f(x_{n+1})| < tol_resid
         * Czyli czy punkt już praktycznie jest miejscem zerowym
         *
         * 3. Mała zmiana względna |x_{n+1} - x_{n}| / max(1, |x_{n+1}) < tol_rel
         * Gdy pierwiastki są duże to sama różnica bezwzgl może nie wystarczyć
         */

        if (err_est < opt.tol_abs_x && resid < opt.tol_resid)
        {
            return std::pair{xnext, k};
        }

        // if (resid < opt.tol_resid)
        // {
        //     return std::pair{xnext, k};
        // }
        //
        // if (err_est / std::max(1.0, std::abs(xnext)) < opt.tol_rel)
        // {
        //     return std::pair{xnext, k};
        // }

        x = xnext;
    }

    std::println(stderr, "picard: std::max_iter osiągnięte");
    return std::nullopt;
}

/**
 * Przybliżanie funkcji linią styczną w punkcie x_n
 * Nowe przybliżenie to punkt przecięcia stycznej z osią OX x_{n+1} = x_n - f(x_n) / f'(x_n)
 *
 * @tparam F  Typ funkcji f(x)
 * @tparam DF Typ funkcji pochodnej f'(x)
 *
 * @param f Funkcja dla której szukane jest miejsce zerowe
 * @param df Pochodna funkcji f(x)
 * @param x0 Punkt początkowy iteracji
 * @param opt Parametry algorytmu (tolerancje, max_iter, verbose)
 *
 * @return std::optional<std::pair<double, int>>
 *
 * @note Estymator błędu to |x_{n+1} - x_n|
 *
 * @warning Jeśli pochodna = 0, czyli styczna jest stała to nie ma jak dalej pójść
 */
template <typename F, typename DF>
[[nodiscard]] std::optional<std::pair<double, int>>
newton(F f, DF df, const double x0, const Options& opt)
{
    print_header("Newton");

    double x = x0;
    for (int k = 1; k <= opt.max_iter; ++k)
    {
        const double fx  = f(x);  // Funckja
        const double dfx = df(x); // Jej pochodna w obecnym punkcie

        if (!std::isfinite(fx) || !std::isfinite(dfx))
        {
            std::println(stderr, "newton: nieprawidlowe wartosci f lub df");
            return std::nullopt;
        }

        // Jeśli pochodna stała, no to nie da się ruszyć (i dzielić przez 0 przy okazji)
        if (dfx == 0.0)
        {
            std::println(stderr, "newton: pochodna zerowa w iteracji {}", k);
            return std::nullopt;
        }

        // Wzór Newtona
        const double xnext   = x - fx / dfx;        // Nowy punkt
        const double err_est = std::abs(xnext - x); // Zmiana
        const double resid   = std::abs(f(xnext));  // Resudiuum

        if (opt.verbose)
        {
            IterInfo it{k, xnext, resid, dfx, err_est};
            print_iter(it);
        }

        /**
         * 1. Czy mała zmiana |x_{n+1} - x_{n}| < tol_abs_x
         * Czyli czy kolejne iteracje się praktycznie nie zmieniają
         *
         * 2. Małe residuum |f(x_{n+1})| < tol_resid
         * Czyli czy punkt już praktycznie jest miejscem zerowym
         *
         * 3. Mała zmiana względna |x_{n+1} - x_{n}| / max(1, |x_{n+1}) < tol_rel
         * Gdy pierwiastki są duże to sama różnica bezwzgl może nie wystarczyć
         */

        if (err_est < opt.tol_abs_x && resid < opt.tol_resid)
        {
            return std::pair{xnext, k};
        }

        // if (resid < opt.tol_resid)
        // {
        //     return std::pair{xnext, k};
        // }
        //
        // if (err_est / std::max(1.0, std::abs(xnext)) < opt.tol_rel)
        // {
        //     return std::pair{xnext, k};
        // }

        x = xnext;
    }

    std::println(stderr, "newton: std::max_iter osiągnięte");
    return std::nullopt;
}

/**
 * Jak Newton ale nie wymaga siecznych
 * Zamiast tego jest przybliżenie nachylenia na podstawie dwóch punktów
 * x_{n+1} = x_n - f(x_n) * (x_n - x_{n-1}) / (f(x_n) - f(x_{n-1}))
 *
 * Nachylenie stycznej jest tu zastąpione nachyleniem siecznej przechodzącej przez punkty (x_{n-1}, f(x_{n-1})) oraz (x_n, f(x_n))
 *
 * @tparam F Typ funkcji f(x)
 *
 * @param f Funkcja dla której szukane jest miejsce zerowe
 * @param x0 Pierwszy punkt początkowy
 * @param x1 Drugi punkt początkowy
 * @param opt Options
 *
 * @return std::optional<std::pair<double, int>>
 *
 * @note Estymator błędu to znowu |x_{n+1} - x_n|
 */
template <typename F>
[[nodiscard]] std::optional<std::pair<double, int>>
Sieczne(F f, const double x0, const double x1, const Options& opt)
{
    print_header("Sieczne");

    // Dwa ostatnie punkty i wartości funkcji w tych punktach
    double x_prev = x0;
    double x      = x1;
    double f_prev = f(x_prev);
    double f_cur  = f(x);

    if (!std::isfinite(f_prev) || !std::isfinite(f_cur))
    {
        std::println(stderr, "Sieczne: nieprawidlowe wartosci funkcji na start");
        return std::nullopt;
    }

    for (int k = 1; k <= opt.max_iter; ++k)
    {
        // Mianownik z wzoru
        const double denom = f_cur - f_prev;
        if (denom == 0.0)
        {
            std::println(stderr, "Sieczne: zerowy mianownik w iteracji {}", k);
            return std::nullopt;
        }

        // Nowe przybliżenie
        const double xnext   = x - f_cur * ( x - x_prev ) / denom;
        const double err_est = std::abs(xnext - x);
        const double resid   = std::abs(f(xnext));

        if (opt.verbose)
        {
            IterInfo it{k, xnext, resid, denom, err_est};
            print_iter(it);
        }

        /**
         * 1. Czy mała zmiana |x_{n+1} - x_{n}| < tol_abs_x
         * Czyli czy kolejne iteracje się praktycznie nie zmieniają
         *
         * 2. Małe residuum |f(x_{n+1})| < tol_resid
         * Czyli czy punkt już praktycznie jest miejscem zerowym
         *
         * 3. Mała zmiana względna |x_{n+1} - x_{n}| / max(1, |x_{n+1}) < tol_rel
         * Gdy pierwiastki są duże to sama różnica bezwzgl może nie wystarczyć
         */

        if (err_est < opt.tol_abs_x && resid < opt.tol_resid)
        {
            return std::pair{xnext, k};
        }

        // if (resid < opt.tol_resid)
        // {
        //     return std::pair{xnext, k};
        // }
        //
        // if (err_est / std::max(1.0, std::abs(xnext)) < opt.tol_rel)
        // {
        //     return std::pair{xnext, k};
        // }

        // Aktualizacja do następnej iteracji
        x_prev = x;
        f_prev = f_cur;
        x      = xnext;
        f_cur  = f(x);
    }

    std::println(stderr, "Sieczne: std::max_iter osiągnięte");
    return std::nullopt;
}

int main()
{
    Options opt;
    opt.max_iter  = 200;
    opt.verbose   = true;
    opt.tol_abs_x = 1e-12;
    opt.tol_resid = 1e-12;
    opt.tol_rel   = 1e-14;

    std::println("\n{:=^128}", " Zadanie 1: ocena zbieznosci metody Picarda ");

    // Funkcja A i jej pochodna
    auto f_a = [](const double x) -> double
    {
        return std::tanh(x) + 2.0 * ( x - 1.0 );
    };

    auto df_a = [](const double x) -> double
    {
        const double t     = std::tanh(x);
        const double sech2 = 1.0 - t * t;
        return sech2 + 2.0;
    };

    // Przekształcenie do postaci x = g(x) dla Picarda
    auto g_a = [](const double x) -> double
    {
        return 1.0 - 0.5 * std::tanh(x);
    };

    std::println("\n{:*^96}", " Przypadek (a): tanh(x) + 2(x-1) = 0 ");
    std::println("Uwaga: dla g(x) = 1 - 0.5*tanh(x) mamy |g'| <= 0.5 -> metoda Picarda powinna zbiegać globalnie.");

    /**
     * Bisekcja na przedziale [0;1]
     * Picard z startem w 0.5
     * Newton z 0.5
     * Sieczne z 0, 1
     */

    if (const auto r = bisection(f_a, 0.0, 1.0, opt))
    {
        std::println("Bisekcja wynik: x = {:.12f} (it={})\n", r->first, r->second);
    }

    if (const auto r = picard(g_a, f_a, 0.5, opt))
    {
        std::println("Picard wynik: x = {:.12f} (it={})\n", r->first, r->second);
    }

    if (const auto r = newton(f_a, df_a, 0.5, opt))
    {
        std::println("Newton wynik: x = {:.12f} (it={})\n", r->first, r->second);
    }

    if (const auto r = Sieczne(f_a, 0.0, 1.0, opt))
    {
        std::println("Sieczne wynik: x = {:.12f} (it={})\n", r->first, r->second);
    }

    // Funkcja B i jej pochodna
    auto f_b = [](const double x) -> double
    {
        return std::sinh(x) + 0.25 * x - 1.0;
    };

    auto df_b = [](const double x) -> double
    {
        return std::cosh(x) + 0.25;
    };

    // Przekształcenie do x = g(x) do Picarda
    auto g_b = [](const double x) -> double
    {
        return std::asinh(1.0 - 0.25 * x);
    };

    std::println("\n********** Przypadek (b): sinh(x) + x/4 - 1 = 0 **********");
    std::println("Uwaga: g(x) = asinh(1 - x/4) ma |g'| <= 0.25 -> metoda Picarda powinna zbiegać.");

    // Bisekcja [0,1] (f(0)=-1; f(1)>0)
    if (const auto r = bisection(f_b, 0.0, 1.0, opt))
    {
        std::println("Bisekcja wynik: x = {:.12f} (it={})\n", r->first, r->second);
    }

    if (const auto r = picard(g_b, f_b, 0.5, opt))
    {
        std::println("Picard wynik: x = {:.12f} (it={})\n", r->first, r->second);
    }

    if (const auto r = newton(f_b, df_b, 0.5, opt))
    {
        std::println("Newton wynik: x = {:.12f} (it={})\n", r->first, r->second);
    }

    if (const auto r = Sieczne(f_b, 0.0, 1.0, opt))
    {
        std::println("Sieczne wynik: x = {:.12f} (it={})\n", r->first, r->second);
    }
}