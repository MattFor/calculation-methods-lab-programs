#include <cmath>
#include <hprint>
#include <vector>
#include <chrono>
#include <iostream>

#define timing true

constexpr std::size_t WIDTH = 45;
constexpr double      EPS   = 1e-16;

template<typename Matrix, typename Pivot>
void print_matrix(Pivot &pivot, Matrix& A)
{
    for (int r = 0; r < A.size(); r++) {
        for (int c = 0; c < A[0].size(); c++) {
            std::print(" {} ", A[pivot[r]][c]);
        }
        std::print("\n");
    }
}

template <typename Vec>
void print_vector(const std::string_view name, const Vec& v)
{
    std::print("{:-^{}}\n", std::format(" {} ", name), WIDTH);

    for (std::size_t i = 0; i < v.size(); ++i)
    {
        std::print("{:^{}s}\n", std::format("x[{}] = {:.10f}", i, v[i]), WIDTH);
    }

    std::print("\n");
}

/**
 * Wyznacza rozkład LU z częściowym wyborem elementu głównego
 * Zapisuje macierze L i U w jednej macierzy A oraz permutację w wektorze pivot
 * 1. PA = LU | gdzie P oznacza permutację wierszy, L jest dolnotrójkątna z jedynkami na przekątnej, a U jest górnotrójkątna
 *
 * @tparam Matrix Typ macierzy przechowującej dane
 * @tparam Pivot Typ wektora permutacji wierszy
 * @param A Macierz wejściowa, nadpisywana przez zapis LU
 * @param pivot Wektor przechowujący kolejność wierszy po pivotowaniu
 * @param N Rozmiar układu
 */
template <typename Matrix, typename Pivot>
void lu_decomposition(Matrix& A, Pivot& pivot, const std::size_t N)
{
    // Inicjalizacja wektora permutacji z brakiem zamian wierszy, czyli po prostu 0 1 2 3
    for (std::size_t i = 0; i < N; ++i)
    {
        pivot[i] = i;
    }

    // Przechodzenie kolumna po kolumnie wykonując eliminację Gaussa i budując L oraz U w tej samej macierzy
    for (std::size_t k = 0; k < N; ++k)
    {
        // Wyszukanie w kolumnie k największego elementu od wiersza k w dół, aby zwiększyć stabilność numeryczną, czyli wybór pivotu jako argmax |A[pivot[i]][k]|
        std::size_t pivot_row = k;
        auto        max_value = std::fabs(A[pivot[k]][k]);

        // Szukanie wiersza z największą wartością bezwzględną (pivota)
        for (std::size_t i = k + 1; i < N; ++i)
        {
            // Liczba znormalizowana więc pewnie nie trzeba sprawdzać czy < EPS
            if (auto val = std::fabs(A[pivot[i]][k]); val > max_value)
            {
                pivot_row = i;
                max_value = val;
            }
        }

        /*
        // Jeśli największy element jest bliski zeru to macierz jest osobliwa, czyli brak jednoznacznego rozwiązania
        if (max_value < EPS)
        {
            throw std::runtime_error("Singular matrix!");
        }
        */

        if (pivot_row != k)
        {
            std::swap(pivot[k], pivot[pivot_row]);
        }

        // Zamiana indeksu zamiast fizycznych wierszy co odpowiada mnożeniu przez macierz permutacji P
        // std::swap(pivot[k], pivot[pivot_row]);

        // Eliminacja elementu pod pivotem poprzez utworzenie zera w U oraz zapisując współczynniki w L, czyli factor = A[i][k] / A[k][k]
        for (std::size_t i = k + 1; i < N; ++i)
        {
            double factor = A[pivot[i]][k] / A[pivot[k]][k];

            // W miejscu A[i][k] jest współczynnik L[i][k]
            A[pivot[i]][k] = factor;

            // Aktualizacja pozostałych elementów wiersza, odejmując factor razy wiersz pivotu, czyli A[i][j] = A[i][j] - factor * A[k][j]
            for (std::size_t j = k + 1; j < N; ++j)
            {
                A[pivot[i]][j] -= factor * A[pivot[k]][j];
            }
        }
    }
}

/**
 * Rozwiązuje układ Ax = b na podstawie rozkładu LU:
 * 1. PA = LU
 * 2. Ly = Pb (podstawienie do przodu)
 * 3. x = y (podstawienie do tyłu)
 *
 * @tparam Matrix Typ macierzy z zapisem LU
 * @tparam Pivot Wektor permutacji
 * @tparam VecB Wektor b
 * @tparam VecX Wektor x
 * @param A Macierz LU
 * @param pivot Permutacja wierszy
 * @param b Prawa strona
 * @param x Wynik
 * @param N Rozmiar
 */
template <typename Matrix, typename Pivot, typename VecB, typename VecX>
void lu_solve(const Matrix& A, const Pivot& pivot, const VecB& b, VecX& x, const std::size_t N)
{
    // Wektor pomocniczy y, w którym zapisywane jest rozwiązanie układu L*y = b
    VecX y(N);

    // Rozwiązanie układu L*y = b idąc od góry, ponieważ każdy y[i] zależy tylko od wcześniejszych y[j] czyli y[i] = b[pivot[i]] - suma L[i][j]*y[j]
    for (std::size_t i = 0; i < N; ++i)
    {
        // Uwzględniana permutacja wierszy, więc zamiast b[i] wybierany jest b[pivot[i]]
        y[i] = b[pivot[i]];

        // Odejmowanie wpływu wcześniej policzonych składników wynikający z dolnej macierzy L zapisanej w A
        for (std::size_t j = 0; j < i; ++j)
        {
            y[i] -= A[pivot[i]][j] * y[j];
        }

        // Brak dzielenia przez element diagonalny, ponieważ w L na przekątnej są jedynki
    }

    // Rozwiązanie układu U*x = y idąc od dołu, ponieważ każdy x[i] zależy od już policzonych x[j] dla j > i, czyli x[i] = (y[i] - suma U[i][j]*x[j]) / U[i][i]
    for (std::size_t i = N; i-- > 0;)
    {
        // Początek od y[i], które będzie korygowana
        x[i] = y[i];

        // Odejmowanie wpływu znanych już wartości x[j] wynikających z górnej macierzy U
        for (std::size_t j = i + 1; j < N; ++j)
        {
            x[i] -= A[pivot[i]][j] * x[j];
        }

        // Dzielenie przez element diagonalny U[i][i], aby wyznaczyć końcową wartość x[i]
        x[i] /= A[pivot[i]][i];
    }
}

void run_lu()
{
    std::vector<std::vector<long double>> A{ { 5, 4, 3, 2, 1 }, { 10, 8, 7, 6, 5 }, { -1, 2, -3, 4, -5 }, { 6, 5, -4, 3, -2 }, { 1, 2, 3, 4, 5 } };
    const std::size_t                     N = A.size();

    // ReSharper disable once CppVariableCanBeMadeConstexpr
    const std::vector<long double> b{ 37, 99, -9, 12, 53 };

    std::vector<long double> x(N);
    std::vector<std::size_t> pivot(N);

    hp::print<'-'>(" Task 1: LU Decomposition ", WIDTH);

    lu_decomposition(A, pivot, N);
    lu_solve(A, pivot, b, x, N);

    print_vector("Solution vector (x)", x);
}

// Badanie wpływu bardzo małego parametru e na stabilność numeryczną rozwiązania układu równań liniowych
void run_bonus()
{
    hp::print<'-'>(" Bonus task | Epsilon influence ", WIDTH);

    // Iteracja po malejących wartościach e, w celu zaobserwowania jak zmienia się uwarunkowanie macierzy oraz błąd numeryczny
    for (double e : { 1.e-5, 1.e-6, 1.e-7, 1.e-8, 1.e-9, 1.e-10, 1.e-11, 1.e-12, 1.e-13, 1.e-14, 1.e-15, 1.e-16, 1.e-17, 1.e-18, 1.e-19, 1.e-20 })
    {
        // Macierz A = J + eI, gdzie J to macierz jedynek
        // Dla tej macierzy wartości własne to e (krotność 3) oraz 4 + e,
        // więc wraz ze zmniejszaniem e układ staje się coraz gorzej uwarunkowany
        std::vector<std::vector<double>> A{ { 1. + e, 1., 1., 1. }, { 1., 1. + e, 1., 1. }, { 1., 1., 1. + e, 1. }, { 1., 1., 1., 1. + e } };
        const std::size_t                N = A.size();

        // Wektor b taki, że dokładne rozwiązanie układu Ax = b jest równe x = [1, 2, 2, 1]
        std::vector b{ 6 + e, 6 + 2 * e, 6 + 2 * e, 6 + e };

        // Wektor rozwiązania numerycznego oraz wektor permutacji pivot
        std::vector<double>      x(N);
        std::vector<std::size_t> pivot(N);

        // Rozkład LU z częściowym wyborem elementu głównego, czyli PA = LU gdzie P to macierz permutacji
        lu_decomposition(A, pivot, N);

        // print_matrix(pivot, A);

        // Rozwiązanie układu wykorzystujące rozkład LU, czyli najpierw Ly = Pb, a następnie Ux = y
        lu_solve(A, pivot, b, x, N);

        // Aktualna wartość e. Jaki ma wpływ na wynik?
        std::print("{:=^{}}\n", std::format(" e = {:.1e} ", e), WIDTH);

        // ReSharper disable once CppVariableCanBeMadeConstexpr
        const std::vector<double> exact{ 1, 2, 2, 1 };

        double max_error = 0.0;
        for (std::size_t i = 0; i < N; ++i)
        {
            // Błąd bezwzględny jako |x_num - x_dokładny|
            const double error = std::fabs(x[i] - exact[i]);
            max_error          = std::max(max_error, error);

            std::print("x[{}] = {:>16.32e}    | err = {:>12.32e}\n", i, x[i], error);
        }

        // Dla tej macierzy teoretyczne uwarunkowanie rośnie mniej więcej od (4 + e) / e
        const double cond2 = ( 4.0 + e ) / e;
        std::print("max err = {:>12.32e} | cond2 ≈ {:>12.32e}\n\n", max_error, cond2);
    }
}

int main()
{
    int choice = 1;

    while (choice == 1 || choice == 2)
    {
        std::print("{:+^40}\n1. LU | 2. Bonus | Anything else to end the program.\n: ", " Choose task ");
        std::cin >> choice;

        switch (choice)
        {
#if timing
            case 1:
            {
                const auto start = std::chrono::high_resolution_clock::now();

                run_lu();

                const auto end      = std::chrono::high_resolution_clock::now();
                const auto duration = std::chrono::duration<double, std::milli>(end - start);

                std::print("Execution time: {:.9} ms\n\n", duration.count());
                break;
            }
#else
            case 1:
            {
                run_lu();
                break;
            }
#endif

#if timing
            case 2:
            {
                const auto start = std::chrono::high_resolution_clock::now();

                run_bonus();

                const auto end      = std::chrono::high_resolution_clock::now();
                const auto duration = std::chrono::duration<double, std::milli>(end - start);

                std::print("Execution time: {:.9} ms\n\n", duration.count());
                break;
            }
#else
            case 2:
            {
                run_bonus();
                break;
            }
#endif

            default:
            {
                std::print("Ending program!");
            }
        }
    }
}
