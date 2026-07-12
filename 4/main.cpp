#include <print>
#include <functional>

using Vec = std::array<double, 3>;
using Mat = std::array<std::array<double, 3>, 3>;

// Układ trzech równań nieliniowych F(x) = 0
struct System3
{
    // Funkcja wektorowa F(x)
    std::function<Vec(const Vec&)> F;

    // Macierz Jacobiego J(x) = dF/dx
    std::function<Mat(const Vec&)> J;
};

// Parametry metody Newtona
struct NewtonConfig
{
    std::size_t maxIter     = 50;    // Maksymalna liczba iteracji
    double      epsResiduum = 1e-12; // Kryterium 1: ||F(x)||_inf
    double      epsStepAbs  = 1e-12; // Kryterium 2: ||dx||_inf
    double      epsStepRel  = 1e-12; // Kryterium 3: ||dx||_inf / max(1, ||x||_inf)
    double      epsPivot    = 1e-14; // Ochrona przed osobliwą macierzą Jacobiego
};

// Wynik działania metody Newtona
struct NewtonResult
{
    Vec         x{};
    std::size_t iterations = 0;
    bool        converged  = false;
    std::string reason;
};

// Norma nieskończona wektora
double norm_inf(const Vec& v)
{
    return std::max({std::abs(v[0]), std::abs(v[1]), std::abs(v[2])});
}

// Norma nieskończona macierzy:
// max z sum bezwzględnych wartości elementów wiersza
double norm_inf(const Mat& A)
{
    double m = 0.0;

    for (const auto& row : A)
    {
        double s = 0.0;
        for (const double a : row)
        {
            s += std::abs(a);
        }

        m = std::max(m, s);
    }

    return m;
}

Vec operator+(const Vec& a, const Vec& b)
{
    return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

Vec operator-(const Vec& a, const Vec& b)
{
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

Vec operator*(const double s, const Vec& v)
{
    return {s * v[0], s * v[1], s * v[2]};
}

// Wypisywanie wektora w czytelnym formacie
void printVec(const Vec& v, const std::string& name)
{
    if (!name.empty())
    {
        std::println("{} = [{:>16.14e}, {:>16.14e}, {:>16.14e}]", name, v[0], v[1], v[2]);
    }
    else
    {
        std::println("[{:>16.14e}, {:>16.14e}, {:>16.14e}]", v[0], v[1], v[2]);
    }
}

// Wypisywanie macierzy 3x3
void printMat(const Mat& A, const std::string& name)
{
    std::println("{} =", name);
    for (const auto& row : A)
    {
        std::println("[{:>16.14e}, {:>16.14e}, {:>16.14e}]", row[0], row[1], row[2]);
    }
}

// Rozwiązanie układu 3x3 metodą eliminacji Gaussa
// z częściowym wyborem elementu głównego
std::optional<Vec> solve3x3(Mat A, Vec b, const double epsPivot)
{
    for (int col = 0; col < 3; ++col)
    {
        // Szukanie pivota w danej kolumnie
        int pivot = col;
        for (int r = col + 1; r < 3; ++r)
        {
            if (std::abs(A[r][col]) > std::abs(A[pivot][col]))
            {
                pivot = r;
            }
        }

        // Zbyt mały pivot oznacza macierz osobliwą lub bliską osobliwej
        if (std::abs(A[pivot][col]) < epsPivot)
        {
            return std::nullopt;
        }

        // Zamiana wierszy
        if (pivot != col)
        {
            std::swap(A[pivot], A[col]);
            std::swap(b[pivot], b[col]);
        }

        // Zerowanie elementów pod pivotem
        for (int r = col + 1; r < 3; ++r)
        {
            const double factor = A[r][col] / A[col][col];

            for (int c = col; c < 3; ++c)
            {
                A[r][c] -= factor * A[col][c];
            }

            b[r] -= factor * b[col];
        }
    }

    // Podstawianie wsteczne
    Vec x{};

    for (int i = 2; i >= 0; --i)
    {
        double sum = b[i];

        for (int j = i + 1; j < 3; ++j)
        {
            sum -= A[i][j] * x[j];
        }

        if (std::abs(A[i][i]) < epsPivot)
        {
            return std::nullopt;
        }

        x[i] = sum / A[i][i];
    }

    return x;
}

// Uogólniona metoda Newtona dla układu trzech równań nieliniowych
// Rozwiązuje układ F(x) = 0 iteracyjnie:
// x_{k+1} = x_k + dx_k, gdzie dx_k spełnia układ liniowy:
// J(x_k) * dx_k = -F(x_k)
NewtonResult generalizedNewton3(const System3& sys, const Vec& x0, const NewtonConfig& cfg)
{
    // Aktualne przybliżenie rozwiązania (start od x0)
    Vec x = x0;

    NewtonResult result;
    result.x = x0;

    // Informacja początkowa
    std::println("START");
    printVec(x, "x0");
    std::println("");

    // Główna pętla iteracyjna metody Newtona
    for (std::size_t k = 0; k < cfg.maxIter; ++k)
    {
        // ========================
        // KROK 1: Obliczenie F(x_k)
        // ========================

        // Wartość układu równań w punkcie x_k
        const Vec Fx = sys.F(x);

        // Residuum: miara "jak bardzo" F(x_k) różni się od zera
        // Używamy normy nieskończonej (max z wartości bezwzględnych składowych)
        const double residuum = norm_inf(Fx);

        std::println("\n\nIteracja {}", k);
        printVec(x, "x_k");
        printVec(Fx, "1. F(x_k)");
        std::println("2. ||F(x_k)||_inf = {:>16.14e}", residuum);

        // Jeśli residuum jest już wystarczająco małe,
        // uznajemy, że znaleźliśmy rozwiązanie
        if (residuum <= cfg.epsResiduum)
        {
            result.x          = x;
            result.iterations = k;
            result.converged  = true;
            result.reason     = "Spełnione kryterium residuum.";
            std::println("STOP: {}", result.reason);
            return result;
        }

        // ========================
        // KROK 2: Obliczenie Jacobiego J(x_k)
        // ========================

        // Macierz Jacobiego zawiera pochodne cząstkowe funkcji F
        // Jest to liniowe przybliżenie układu w punkcie x_k
        const Mat Jx = sys.J(x);
        printMat(Jx, "3. J(x_k)");

        // ========================
        // KROK 3: Rozwiązanie układu liniowego
        // ========================

        // Równanie Newtona:
        // J(x_k) * dx_k = -F(x_k)
        // prawa strona to -F(x_k)
        const Vec rhs{-Fx[0], -Fx[1], -Fx[2]};

        // Rozwiązujemy układ 3x3 metodą eliminacji Gaussa
        auto dxOpt = solve3x3(Jx, rhs, cfg.epsPivot);

        // Jeśli nie udało się rozwiązać układu (np. Jacobian osobliwy)
        if (!dxOpt)
        {
            result.x          = x;
            result.iterations = k;
            result.converged  = false;
            result.reason     = "Macierz Jacobiego osobliwa lub bliska osobliwej.";
            std::println("STOP: {}", result.reason);
            return result;
        }

        // Wektor poprawki Newtona
        const Vec dx = *dxOpt;

        // Nowe przybliżenie rozwiązania
        const Vec xNext = x + dx;

        // ========================
        // KROK 4: Estymacja błędu i residuum
        // ========================

        // Estymator błędu bezwzględnego ~ ||dx||
        // (w metodzie Newtona poprawka jest przybliżeniem błędu)
        const double stepAbs = norm_inf(dx);

        // Estymator błędu względnego
        const double stepRel = stepAbs / std::max(1.0, norm_inf(xNext));

        // Residuum dla nowego punktu x_{k+1}
        const Vec    FNext        = sys.F(xNext);
        const double residuumNext = norm_inf(FNext);

        std::print("4. dx_k = ");
        printVec(dx, "");

        std::println("5. ||dx_k||_inf = {:>16.14e}", stepAbs);

        std::println("6. estymator_błędu_rel = {:>16.14e} = ||dx_k|| / max(1, ||x_{{k+1}}||)", stepRel);

        std::print("7. x_{{k+1}} = ");
        printVec(xNext, "");

        std::print("8. F(x_{{k+1}}) = ");
        printVec(FNext, "");

        std::println("9. ||F(x_{{k+1}})||_inf = {:>16.14e}", residuumNext);

        // ========================
        // KROK 5: Kryteria stopu
        // ========================

        // 3 niezależne kryteria zakończenia iteracji:
        // 1. małe residuum -> układ spełniony
        // 2. mały krok bezwzględny -> brak zmian
        // 3. mały krok względny -> stabilizacja rozwiązania
        const bool stopResiduum = residuumNext <= cfg.epsResiduum;
        const bool stopStepAbs  = stepAbs <= cfg.epsStepAbs;
        const bool stopStepRel  = stepRel <= cfg.epsStepRel;

        std::println("\nKryteria stopu:");
        std::println("- residuum <= epsResiduum\t->  {}", stopResiduum ? "TAK" : "NIE");
        std::println("- krok abs <= epsStepAbs\t->  {}", stopStepAbs ? "TAK" : "NIE");
        std::println("- krok rel <= epsStepRel\t->  {}", stopStepRel ? "TAK" : "NIE");

        // Aktualizacja rozwiązania
        x                 = xNext;
        result.x          = x;
        result.iterations = k + 1;

        // Sprawdzenie warunków zakończenia (priorytet: residuum -> abs -> rel)
        if (stopResiduum)
        {
            result.converged = true;
            result.reason    = "Spełnione kryterium residuum.";
            return result;
        }

        if (stopStepAbs)
        {
            result.converged = true;
            result.reason    = "Spełnione kryterium bezwzględnej poprawki.";
            return result;
        }

        if (stopStepRel)
        {
            result.converged = true;
            result.reason    = "Spełnione kryterium względnej poprawki.";
            return result;
        }
    }

    // Jeśli przekroczono maksymalną liczbę iteracji -> brak zbieżności
    result.converged = false;
    result.reason    = "Osiągnięto maksymalną liczbę iteracji.";
    return result;
}

int main()
{
    System3 sys;

    // Układ z zadania:
    // x^2 + y^2 + z^2 = 4
    // x^2 + y^2/2 = 1
    // x*y = 1/2
    sys.F = [](const Vec& x) -> Vec
    {
        const double X = x[0];
        const double Y = x[1];
        const double Z = x[2];

        return {X * X + Y * Y + Z * Z - 4.0, X * X + 0.5 * Y * Y - 1.0, X * Y - 0.5};
    };

    // Jakobian:
    // [2x  2y  2z]
    // [2x   y   0]
    // [ y   x   0]
    sys.J = [](const Vec& x) -> Mat
    {
        const double X = x[0];
        const double Y = x[1];
        const double Z = x[2];

        return {{{2.0 * X, 2.0 * Y, 2.0 * Z}, {2.0 * X, Y, 0.0}, {Y, X, 0.0}}};
    };

    // Przybliżenie początkowe dobrane tak, aby metoda była zbieżna
    constexpr Vec x0{1.0, 0.5, 1.0};

    NewtonConfig cfg;
    cfg.maxIter     = 30;
    cfg.epsResiduum = 1e-12;
    cfg.epsStepAbs  = 1e-12;
    cfg.epsStepRel  = 1e-12;
    cfg.epsPivot    = 1e-14;

    auto [x, iterations, converged, reason] = generalizedNewton3(sys, x0, cfg);

    std::println("\nWYNIK KONCOWY");
    printVec(x, "x*");
    std::println("Iteracje = {}", iterations);
    std::println("Zbieznosc = {}", converged ? "TAK" : "NIE");
    std::println("Powod = {}", reason);
}