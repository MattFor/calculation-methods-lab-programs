#include <hprint>

#include <cmath>
#include <vector>
#include <future>
#include <fstream>
#include <execution>
#include <filesystem>

using cd  = const double;
using cst = const std::size_t;

namespace fs = std::filesystem;

constexpr double pi = std::numbers::pi_v<double>;

constexpr double explicit_target_lambda = 0.4;
constexpr double implicit_target_lambda = 1.0;
constexpr double jacobi_tolerance       = 1e-14;

constexpr std::size_t snapshot_count        = 20;
constexpr std::size_t jacobi_max_iterations = 200000;

struct grid_config
{
    double length{};
    double t_max{};
    double diffusion{};
    double target_lambda{};

    std::size_t n_intervals{};
    std::size_t n_time_steps{};

    double h{};
    double dt{};
    double actual_lambda{};

    std::vector<double> x_coordinates{};
    std::vector<double> forcing_values{}; // Źródło wymuszenia dobrane do znanego rozwiązania dokładnego
};

struct snapshot_record
{
    std::string method_name;
    std::size_t snapshot_id{};
    std::size_t step_index{};

    double time{};

    std::vector<double> x;
    std::vector<double> numerical;

    std::vector<double> exact;
    std::vector<double> abs_error;
};

struct simulation_result
{
    std::string method_name;

    grid_config grid;

    double runtime_seconds{};
    double max_abs_error_at_t_max{};

    std::vector<double> time_points;
    std::vector<double> max_abs_error_by_time;
    std::vector<double> final_values;

    std::vector<snapshot_record> snapshots;
};

struct convergence_row
{
    std::string method_name;
    std::size_t n_intervals{};

    double h{};
    double dt{};
    double actual_lambda{};

    std::size_t n_time_steps{};

    double runtime_seconds{};
    double max_abs_error_at_t_max{};
};

struct thomas_workspace
{
    std::vector<double> lower;    // Współczynnik podprzekątnej
    std::vector<double> diagonal; // Współczynnik przekątnej głównej
    std::vector<double> upper;    // Współczynnik nadprzekątnej
    std::vector<double> c_prime;  // Zmodyfikowane współczynniki nadprzekątnej
    std::vector<double> denom;    // Zmodyfikowane przekątne po eliminacji w przód
};

/**
 * Oblicza rozwiązanie dokładne testowego problemu dyfuzji
 *
 * @param x Współrzędna przestrzenna
 * @param t Czas obserwacji
 * @param diffusion Współczynnik dyfuzji
 *
 * @return Wartość rozwiązania dokładnego w punkcie (x, t)
 */
double exact_solution(cd x, cd t, cd diffusion)
{
    return ( 1.0 - std::exp(-pi * pi * diffusion * t) ) * std::sin(pi * x); // (1 - e^{-pi^2 D t}) sin(pi x)
}

/**
 * Sprawdza czy podstawowe parametry problemu są dodatnie
 *
 * @param length Długość przedziału przestrzennego
 * @param t_max Czas końcowy symulacji
 * @param diffusion Współczynnik dyfuzji
 * @param target_lambda Docelowa wartość lambda używana do doboru kroku czasu
 *
 * @throws Jeśli którykolwiek parametr nie jest dodatni
 */
void validate_problem_parameters(cd length, cd t_max, cd diffusion, cd target_lambda)
{
    if (!( length > 0.0 ))
    {
        throw std::invalid_argument("Długość musi być dodatnia!");
    }

    if (!( t_max > 0.0 ))
    {
        throw std::invalid_argument("t_max musi być dodatnie!");
    }

    if (!( diffusion > 0.0 ))
    {
        throw std::invalid_argument("Dyfuzja musi być dodatnia!");
    }

    if (!( target_lambda > 0.0 ))
    {
        throw std::invalid_argument("target_lambda musi być dodatnie!");
    }
}

/**
 * Buduje siatkę obliczeniową i wyprowadza z niej kroki czasowe
 *
 * @param n_intervals Liczba przedziałów przestrzennych
 * @param length Długość domeny
 * @param t_max Czas końcowy symulacji
 * @param diffusion Współczynnik dyfuzji
 * @param target_lambda Docelowa wartość lambda używana przy doborze dt
 *
 * @return Wypełniona konfiguracja siatki
 *
 * @throws Jeśli parametry są niepoprawne albo siatka ma za mało przedziałów
 */
grid_config make_grid_config(cst n_intervals, cd length, cd t_max, cd diffusion, cd target_lambda)
{
    validate_problem_parameters(length, t_max, diffusion, target_lambda);

    if (n_intervals < 2)
    {
        throw std::invalid_argument("n_intervals musi być co najmniej 2!");
    }

    grid_config grid;
    grid.length        = length;
    grid.t_max         = t_max;
    grid.diffusion     = diffusion;
    grid.target_lambda = target_lambda;
    grid.n_intervals   = n_intervals;
    grid.h             = grid.length / static_cast<double>(grid.n_intervals); // h = L / N czyli odległość między kolejnymi węzłami

    cd target_dt      = grid.target_lambda * grid.h * grid.h / grid.diffusion; // dt = lambda * h^2 / D a potem zaokrąglamy liczbę kroków do całkowitej wartości
    grid.n_time_steps = static_cast<std::size_t>(std::ceil(grid.t_max / target_dt));
    if (grid.n_time_steps == 0)
    {
        grid.n_time_steps = 1;
    }

    grid.dt            = grid.t_max / static_cast<double>(grid.n_time_steps); // dt po dopasowaniu liczby kroków do całego przedziału czasu
    grid.actual_lambda = grid.diffusion * grid.dt / ( grid.h * grid.h );      // Lambda = D * dt / h^2 dla rzeczywistej siatki

    grid.x_coordinates.resize(grid.n_intervals + 1);
    grid.forcing_values.resize(grid.n_intervals + 1);

    for (std::size_t i = 0; i <= grid.n_intervals; ++i)
    {
        cd x                   = grid.h * static_cast<double>(i);
        grid.x_coordinates[i]  = x;
        grid.forcing_values[i] = grid.diffusion * pi * pi * std::sin(pi * x); // Wymuszenie dobrane tak, aby dokładne rozwiązanie miało postać sinusa w czasie
    }

    return grid;
}

/**
 * Rozkłada migawki możliwie równomiernie po wszystkich krokach czasu
 *
 * @param n_time_steps Liczba kroków czasu w symulacji
 * @param desired_snapshot_count Liczba migawek które chcemy zapisać
 *
 * @return Posortowane i unikalne indeksy kroków dla migawek
 *
 * @throws Jeśli żądana liczba migawek jest mniejsza niż 2
 */
std::vector<std::size_t> make_even_snapshot_steps(cst n_time_steps, cst desired_snapshot_count)
{
    if (desired_snapshot_count < 2)
    {
        throw std::invalid_argument("desired_snapshot_count musi być co najmniej 2!");
    }

    std::vector<std::size_t> steps;
    steps.reserve(desired_snapshot_count);

    for (std::size_t i = 0; i < desired_snapshot_count; ++i)
    {
        cst step = i * n_time_steps / ( desired_snapshot_count - 1 ); // Rzutowanie równego podziału zakresu na całkowite indeksy
        steps.push_back(step);
    }

    std::ranges::sort(steps);
    steps.erase(std::ranges::unique(steps).begin(), steps.end()); // Powtórzenia mogą pojawić się przez dzielenie całkowite

    if (steps.empty() || steps.front() != 0)
    {
        steps.insert(steps.begin(), 0);
    }

    if (steps.back() != n_time_steps)
    {
        steps.push_back(n_time_steps);
    }

    std::ranges::sort(steps);
    steps.erase(std::ranges::unique(steps).begin(), steps.end());

    return steps;
}

/**
 * Liczy maksymalny błąd bezwzględny między rozwiązaniem numerycznym a dokładnym
 *
 * @param numerical Wartości numeryczne na siatce
 * @param grid Konfiguracja siatki z pozycjami węzłów
 * @param time Chwila czasu dla której porównujemy rozwiązania
 *
 * @return Maksymalny błąd bezwzględny w całej siatce
 *
 * @throws Jeśli rozmiar wektora nie pasuje do liczby węzłów
 */
double compute_max_abs_error(const std::vector<double>& numerical, const grid_config& grid, cd time)
{
    if (numerical.size() != grid.x_coordinates.size())
    {
        throw std::runtime_error("Niezgodność rozmiaru rozwiązania!");
    }

    std::atomic max_error = 0.0;

    std::vector<std::size_t> indices(numerical.size());
    std::iota(indices.begin(), indices.end(), 0);

    std::for_each(std::execution::par_unseq, indices.begin(), indices.end(), [&](cst i)
    {
        cd exact = exact_solution(grid.x_coordinates[i], time, grid.diffusion);

        cd err = std::abs(numerical[i] - exact); // |u_num - u_exact| w jednym węźle

        double current = max_error.load(std::memory_order_relaxed);

        while (err > current && !max_error.compare_exchange_weak(current, err, std::memory_order_relaxed))
        {}
    });

    return max_error.load();
}

/**
 * Zapisuje jeden pełny stan symulacji wraz z błędem względem rozwiązania dokładnego
 *
 * @param method_name Nazwa metody numerycznej
 * @param snapshot_id Kolejny numer migawki
 * @param step_index Indeks kroku czasu
 * @param time Czas odpowiadający migawce
 * @param numerical Wartości numeryczne na siatce
 * @param grid Konfiguracja siatki
 *
 * @return Struktura zawierająca współrzędne wartości i błąd bezwzględny
 *
 * @throws Jeśli rozmiar wektora nie pasuje do liczby węzłów
 */
snapshot_record make_snapshot_record(const std::string& method_name, cst snapshot_id, cst step_index, cd time, const std::vector<double>& numerical, const grid_config& grid)
{
    if (numerical.size() != grid.x_coordinates.size())
    {
        throw std::runtime_error("Niezgodność rozmiaru migawki!");
    }

    snapshot_record record;
    record.method_name = method_name;
    record.snapshot_id = snapshot_id;
    record.step_index  = step_index;
    record.time        = time;
    record.x           = grid.x_coordinates;
    record.numerical   = numerical;
    record.exact.resize(numerical.size());
    record.abs_error.resize(numerical.size());

    for (std::size_t i = 0; i < numerical.size(); ++i)
    {
        record.exact[i]     = exact_solution(grid.x_coordinates[i], time, grid.diffusion);
        record.abs_error[i] = std::abs(record.numerical[i] - record.exact[i]);
    }

    return record;
}

/**
 * Przygotowuje przestrzeń roboczą do rozwiązania układu trójdiagonalnego metodą Thomasa
 *
 * @param system_size Rozmiar układu liniowego
 * @param lower_value Wartość podprzekątnej
 * @param diagonal_value Wartość przekątnej głównej
 * @param upper_value Wartość nadprzekątnej
 *
 * @return Gotowa przestrzeń robocza z wyliczonymi współczynnikami pośrednimi
 *
 * @throws Jeśli układ ma rozmiar zero albo prowadzi do dzielenia przez zero
 */
thomas_workspace make_thomas_workspace(cst system_size, cd lower_value, cd diagonal_value, cd upper_value)
{
    if (system_size == 0)
    {
        throw std::invalid_argument("system_size musi być dodatnie!");
    }

    thomas_workspace ws;
    ws.lower.assign(system_size, lower_value);
    ws.diagonal.assign(system_size, diagonal_value);
    ws.upper.assign(system_size, upper_value);
    ws.c_prime.assign(system_size, 0.0);
    ws.denom.assign(system_size, 0.0);

    ws.lower[0]               = 0.0; // Pierwszy węzeł nie ma lewego sąsiada
    ws.upper[system_size - 1] = 0.0; // Ostatni węzeł nie ma prawego sąsiada

    if (std::abs(ws.diagonal[0]) < 1e-30)
    {
        throw std::runtime_error("Nieprawidłowy układ trójdiagonalny: zerowa przekątna!");
    }

    ws.denom[0]   = ws.diagonal[0];
    ws.c_prime[0] = ws.upper[0] / ws.denom[0]; // Start eliminacji w przód dla metody Thomasa

    for (std::size_t i = 1; i < system_size; ++i)
    {
        ws.denom[i] = ws.diagonal[i] - ws.lower[i] * ws.c_prime[i - 1]; // Eliminacja wpływu poprzedniego wiersza
        if (std::abs(ws.denom[i]) < 1e-30)
        {
            throw std::runtime_error("Nieprawidłowy układ trójdiagonalny: zerowy element główny w faktoryzacji Thomasa!");
        }
        ws.c_prime[i] = i + 1 < system_size ? ws.upper[i] / ws.denom[i] : 0.0;
    }

    return ws;
}

/**
 * Rozwiązuje układ trójdiagonalny metodą Thomasa bez iteracji
 *
 * @param ws Wcześniej przygotowana przestrzeń robocza
 * @param rhs Wektor wyrazów wolnych
 * @param solution Wynikowy wektor rozwiązania
 * @param forward_buffer Bufor na przebieg w przód
 *
 * @throws Jeśli wejście ma zły rozmiar albo układ jest osobliwy
 */
void solve_tridiagonal_thomas(const thomas_workspace& ws, const std::vector<double>& rhs, std::vector<double>& solution, std::vector<double>& forward_buffer)
{
    cst n = rhs.size();
    if (n == 0)
    {
        throw std::invalid_argument("rhs nie może być puste!");
    }
    if (ws.denom.size() != n || ws.c_prime.size() != n || ws.lower.size() != n)
    {
        throw std::runtime_error("Niezgodność rozmiaru przestrzeni roboczej Thomasa!");
    }

    forward_buffer.resize(n);
    forward_buffer[0] = rhs[0] / ws.denom[0]; // Pierwszy krok eliminacji w przód

    for (std::size_t i = 1; i < n; ++i)
    {
        forward_buffer[i] = ( rhs[i] - ws.lower[i] * forward_buffer[i - 1] ) / ws.denom[i]; // Recurencja eliminacji w przód dla kolejnych wierszy
    }

    solution.resize(n);
    if (n == 1)
    {
        solution[0] = forward_buffer[0];
        return;
    }

    solution[n - 1] = forward_buffer[n - 1]; // Wzór x_i = d_i - c_i * x_{i+1} zaczyna się od ostatniego węzła
    for (std::size_t i = n - 1; i-- > 0;)
    {
        solution[i] = forward_buffer[i] - ws.c_prime[i] * solution[i + 1]; // Wzór x_i = d_i - c_i * x_{i+1}
    }
}

/**
 * Rozwiązuje układ trójdiagonalny metodą Jacobiego do ustalonej tolerancji
 *
 * @param rhs Wektor wyrazów wolnych
 * @param lower_value Wartość podprzekątnej
 * @param diagonal_value Wartość przekątnej głównej
 * @param upper_value Wartość nadprzekątnej
 * @param initial_guess Startowe przybliżenie rozwiązania
 * @param tolerance Maksymalna dopuszczalna zmiana między iteracjami
 * @param max_iterations Maksymalna liczba iteracji
 *
 * @return Przybliżone rozwiązanie układu
 *
 * @throws Jeśli wejście ma zły rozmiar albo solver nie zbiegnie na czas
 */
std::vector<double> solve_tridiagonal_jacobi(const std::vector<double>& rhs, cd lower_value, cd diagonal_value, cd upper_value, const std::vector<double>& initial_guess, cd tolerance, cst max_iterations)
{
    cst n = rhs.size();
    if (n == 0)
    {
        throw std::invalid_argument("rhs nie może być puste!");
    }

    if (initial_guess.size() != n)
    {
        throw std::invalid_argument("Niezgodność rozmiaru initial_guess!");
    }

    if (std::abs(diagonal_value) < 1e-30)
    {
        throw std::runtime_error("Nieprawidłowy układ trójdiagonalny: zerowa przekątna!");
    }

    std::vector<double> previous = initial_guess;
    std::vector         next(n, 0.0);

    for (std::size_t iter = 0; iter < max_iterations; ++iter)
    {
        std::atomic max_delta = 0.0;

        std::vector<std::size_t> indices(n);
        std::iota(indices.begin(), indices.end(), 0);

        std::for_each(std::execution::par_unseq, indices.begin(), indices.end(), [&](cst i)
        {
            cd left = i == 0 ? 0.0 : lower_value * previous[i - 1];

            cd right = i + 1 == n ? 0.0 : upper_value * previous[i + 1];

            next[i] = ( rhs[i] - left - right ) / diagonal_value; // Nowe przybliżenie x_i^{k+1}

            cd local_delta = std::abs(next[i] - previous[i]); // delta_i = |x_i^{k+1} - x_i^k|

            double current_max = max_delta.load(std::memory_order_relaxed);

            while (local_delta > current_max && !max_delta.compare_exchange_weak(current_max, local_delta, std::memory_order_relaxed))
            {}
        });

        previous.swap(next);

        if (max_delta < tolerance)
        {
            return previous;
        }
    }

    throw std::runtime_error("Solver Jacobiego nie zbiega w limicie iteracji!");
}

/**
 * Wykonuje pełną symulację w czasie i zbiera metryki oraz migawki
 *
 * @tparam Stepper Funkcja wykonująca pojedynczy krok czasowy
 * @param method_name Nazwa metody numerycznej
 * @param grid Konfiguracja siatki
 * @param snapshot_steps Indeksy kroków przeznaczone do zapisywania migawek
 * @param stepper Funkcja aktualizująca stan z kroku n do n + 1
 *
 * @return Wynik symulacji z historią błędu i zarejestrowanymi migawkami
 *
 * @throws Jeśli rozmiar siatki jest niespójny
 */
template <typename Stepper>
simulation_result simulate_solution(const std::string& method_name, const grid_config& grid, const std::vector<std::size_t>& snapshot_steps, Stepper&& stepper)
{
    if (grid.x_coordinates.size() != grid.n_intervals + 1)
    {
        throw std::runtime_error("Niezgodność rozmiaru siatki!");
    }

    simulation_result result;
    result.method_name = method_name;
    result.grid        = grid;
    result.time_points.reserve(grid.n_time_steps + 1);
    result.max_abs_error_by_time.reserve(grid.n_time_steps + 1);

    std::vector current(grid.n_intervals + 1, 0.0);
    std::vector next(grid.n_intervals + 1, 0.0);

    const auto start = std::chrono::steady_clock::now();

    [[maybe_unused]] std::size_t snapshot_cursor = 0;
    auto                         store_state     = [&](cst step_index)
    {
        cd time = static_cast<double>(step_index) * grid.dt;
        result.time_points.push_back(time);
        result.max_abs_error_by_time.push_back(compute_max_abs_error(current, grid, time));
    };

    auto store_snapshots_if_needed = [&](cst step_index)
    {
        while (snapshot_cursor < snapshot_steps.size() && snapshot_steps[snapshot_cursor] == step_index)
        {
            cd time = static_cast<double>(step_index) * grid.dt;
            result.snapshots.push_back(make_snapshot_record(result.method_name, snapshot_cursor, step_index, time, current, grid));
            ++snapshot_cursor;
        }
    };

    store_state(0); // Zapis stanu początkowego przed wykonaniem pierwszego kroku czasowego
    store_snapshots_if_needed(0);

    double time = 0.0;
    for (std::size_t step = 0; step < grid.n_time_steps; ++step)
    {
        stepper(current, next, time); // stepper liczy przejście z u^n do u^{n+1}
        next.front() = 0.0;
        next.back()  = 0.0;
        current.swap(next);

        time = static_cast<double>(step + 1) * grid.dt;
        store_state(step + 1);
        store_snapshots_if_needed(step + 1);
    }

    const auto stop               = std::chrono::steady_clock::now();
    result.runtime_seconds        = std::chrono::duration<double>(stop - start).count();
    result.final_values           = current;
    result.max_abs_error_at_t_max = result.max_abs_error_by_time.back();

    return result;
}

/**
 * Uruchamia jawną metodę FTCS dla równania dyfuzji
 *
 * @param n_intervals Liczba przedziałów przestrzennych
 * @param length Długość domeny
 * @param t_max Czas końcowy symulacji
 * @param diffusion Współczynnik dyfuzji
 * @param snapshot_steps Indeksy kroków przeznaczone do zapisywania migawek
 *
 * @return Wynik symulacji dla metody jawnej
 *
 * @throws Metoda staje się niestabilna
 */
simulation_result run_explicit_method(cst n_intervals, cd length, cd t_max, cd diffusion, const std::vector<std::size_t>& snapshot_steps = {})
{
    const grid_config grid = make_grid_config(n_intervals, length, t_max, diffusion, explicit_target_lambda);

    if (grid.actual_lambda > 0.5)
    {
        throw std::runtime_error("Metoda jawna jest niestabilna: lambda > 0.5!");
    }

    auto stepper = [grid](const std::vector<double>& current, std::vector<double>& next, double)
    {
        next[0]     = 0.0;
        next.back() = 0.0;

        std::vector<std::size_t> indices(grid.n_intervals - 1);

        std::iota(indices.begin(), indices.end(), 1);

        std::for_each(std::execution::par_unseq, indices.begin(), indices.end(), [&](cst i)
        {
            next[i] = current[i] + grid.actual_lambda * ( current[i - 1] - 2.0 * current[i] + current[i + 1] ) + grid.dt * grid.forcing_values[i];
        });
    };

    // Forward Time Center Space czyli czas w przód i przestrzeń w centrum
    // Druga pochodna po przestrzeni jest aproksymowana różnicą centralną
    // Po połączeniu jest jawny schemat dla równania dyfuzji
    return simulate_solution("explicit_ftcs", grid, snapshot_steps, stepper);
}

/**
 * Uruchamia metodę Cranka Nicolsona z solverem Thomasa
 *
 * @param n_intervals Liczba przedziałów przestrzennych
 * @param length Długość domeny
 * @param t_max Czas końcowy symulacji
 * @param diffusion Współczynnik dyfuzji
 * @param snapshot_steps Indeksy kroków przeznaczone do zapisywania migawek
 *
 * @return Wynik symulacji dla schematu Cranka Nicolsona
 */
simulation_result run_crank_nicolson_thomas(cst n_intervals, double length, double t_max, double diffusion, const std::vector<std::size_t>& snapshot_steps = {})
{
    const grid_config grid     = make_grid_config(n_intervals, length, t_max, diffusion, implicit_target_lambda);
    cst               unknowns = grid.n_intervals - 1;

    cd lower_value    = -0.5 * grid.actual_lambda;
    cd diagonal_value = 1.0 + grid.actual_lambda;
    cd upper_value    = -0.5 * grid.actual_lambda; // Współczynniki układu wynikające ze schematu Cranka Nicolsona

    const thomas_workspace ws = make_thomas_workspace(unknowns, lower_value, diagonal_value, upper_value);

    std::vector         rhs(unknowns, 0.0);
    std::vector<double> interior_solution;
    std::vector<double> forward_buffer;

    auto stepper = [&, grid](const std::vector<double>& current, std::vector<double>& next, double) mutable
    {
        for (std::size_t i = 0; i < unknowns; ++i)
        {
            cst node = i + 1;
            rhs[i]   = ( 1.0 - grid.actual_lambda ) * current[node] + 0.5 * grid.actual_lambda * ( current[node - 1] + current[node + 1] ) + grid.dt * grid.forcing_values[node]; // Prawa strona układu liniowego dla węzłów wewnętrznych
        }

        solve_tridiagonal_thomas(ws, rhs, interior_solution, forward_buffer);

        next[0]     = 0.0;
        next.back() = 0.0;
        for (std::size_t i = 0; i < unknowns; ++i)
        {
            next[i + 1] = interior_solution[i];
        }
    };

    return simulate_solution("crank_nicolson_thomas", grid, snapshot_steps, stepper);
}

/**
 * Uruchamia metodę Cranka Nicolsona z solverem Jacobiego
 *
 * @param n_intervals Liczba przedziałów przestrzennych
 * @param length Długość domeny
 * @param t_max Czas końcowy symulacji
 * @param diffusion Współczynnik dyfuzji
 * @param snapshot_steps Indeksy kroków przeznaczone do zapisywania migawek
 *
 * @return Wynik symulacji dla schematu Cranka Nicolsona
 */
simulation_result run_crank_nicolson_jacobi(cst n_intervals, cd length, cd t_max, cd diffusion, const std::vector<std::size_t>& snapshot_steps = {})
{
    const grid_config grid     = make_grid_config(n_intervals, length, t_max, diffusion, implicit_target_lambda);
    cst               unknowns = grid.n_intervals - 1;

    cd lower_value    = -0.5 * grid.actual_lambda;
    cd diagonal_value = 1.0 + grid.actual_lambda;
    cd upper_value    = -0.5 * grid.actual_lambda; // Współczynniki układu wynikające ze schematu Cranka Nicolsona

    std::vector rhs(unknowns, 0.0);
    std::vector initial_guess(unknowns, 0.0);

    auto stepper = [&, grid](const std::vector<double>& current, std::vector<double>& next, double) mutable
    {
        for (std::size_t i = 0; i < unknowns; ++i)
        {
            cst node         = i + 1;
            rhs[i]           = ( 1.0 - grid.actual_lambda ) * current[node] + 0.5 * grid.actual_lambda * ( current[node - 1] + current[node + 1] ) + grid.dt * grid.forcing_values[node]; // Prawa strona układu liniowego dla węzłów wewnętrznych
            initial_guess[i] = current[node];                                                                                                                                             // Start od wartości z poprzedniego kroku czasowego
        }

        const std::vector<double> interior_solution = solve_tridiagonal_jacobi(rhs, lower_value, diagonal_value, upper_value, initial_guess, jacobi_tolerance, jacobi_max_iterations);

        next[0]     = 0.0;
        next.back() = 0.0;
        for (std::size_t i = 0; i < unknowns; ++i)
        {
            next[i + 1] = interior_solution[i];
        }
    };

    return simulate_solution("crank_nicolson_jacobi", grid, snapshot_steps, stepper);
}

/**
 * Upewnia się że katalog wyjściowy istnieje
 *
 * @param output_dir Ścieżka do katalogu na wyniki
 *
 * @throws Jeśli katalogu nie da się utworzyć albo ścieżka nie jest katalogiem
 */
void ensure_output_directory(const fs::path& output_dir)
{
    if (std::error_code ec; !fs::exists(output_dir, ec))
    {
        if (!fs::create_directories(output_dir, ec) || ec)
        {
            throw std::runtime_error("Nie udało się utworzyć katalogu wyjściowego: " + output_dir.string());
        }
    }
    else if (!fs::is_directory(output_dir, ec) || ec)
    {
        throw std::runtime_error("Ścieżka wyjściowa istnieje, ale nie jest katalogiem: " + output_dir.string());
    }
}

/**
 * Zapisuje tabelę z wynikami badania zbieżności do pliku CSV
 *
 * @param file_path Ścieżka do pliku wynikowego
 * @param rows Wiersze z metrykami poszczególnych uruchomień
 *
 * @throws Jeśli pliku nie da się otworzyć do zapisu
 */
void write_convergence_csv(const fs::path& file_path, const std::vector<convergence_row>& rows)
{
    std::ofstream out(file_path);
    out.imbue(std::locale::classic());
    if (!out)
    {
        throw std::runtime_error("Nie udało się otworzyć pliku do zapisu: " + file_path.string());
    }

    out << "method,n_intervals,h,dt,lambda_actual,n_time_steps,runtime_seconds,max_abs_error_at_t_max\n"; // Nagłówek CSV opisuje kolejność kolumn
    for (const auto& [method_name, n_intervals, h, dt, actual_lambda, n_time_steps, runtime_seconds, max_abs_error_at_t_max] : rows)
    {
        out << method_name << ',' << n_intervals << ',' << h << ',' << dt << ',' << actual_lambda << ',' << n_time_steps << ',' << runtime_seconds << ',' << max_abs_error_at_t_max << '\n';
    }
}

/**
 * Zapisuje przebieg błędu w czasie do pliku CSV
 *
 * @param file_path Ścieżka do pliku wynikowego
 * @param results Zestaw rezultatów symulacji
 *
 * @throws Jeśli pliku nie da się otworzyć do zapisu
 */
void write_time_error_csv(const fs::path& file_path, const std::vector<simulation_result>& results)
{
    std::ofstream out(file_path);
    out.imbue(std::locale::classic());
    if (!out)
    {
        throw std::runtime_error("Nie udało się otworzyć pliku do zapisu: " + file_path.string());
    }

    out << "method,step_index,time,max_abs_error\n"; // Każdy wiersz odpowiada jednemu momentowi w czasie
    for (const auto& result : results)
    {
        for (std::size_t i = 0; i < result.time_points.size(); ++i)
        {
            out << result.method_name << ',' << i << ',' << result.time_points[i] << ',' << result.max_abs_error_by_time[i] << '\n';
        }
    }
}

/**
 * Zapisuje pełne migawki rozwiązań do pliku CSV
 *
 * @param file_path Ścieżka do pliku wynikowego
 * @param results Zestaw rezultatów symulacji
 *
 * @throws Jeśli pliku nie da się otworzyć do zapisu
 */
void write_snapshots_csv(const fs::path& file_path, const std::vector<simulation_result>& results)
{
    std::ofstream out(file_path);
    out.imbue(std::locale::classic());
    if (!out)
    {
        throw std::runtime_error("Nie udało się otworzyć pliku do zapisu: " + file_path.string());
    }

    out << "method,snapshot_id,step_index,time,x,numerical,exact,abs_error\n"; // Zapisywanie migawek punkt po punkcie dla każdego węzła
    for (const auto& result : results)
    {
        for (const auto& [method_name, snapshot_id, step_index, time, x, numerical, exact, abs_error] : result.snapshots)
        {
            for (std::size_t i = 0; i < x.size(); ++i)
            {
                out << method_name << ',' << snapshot_id << ',' << step_index << ',' << time << ',' << x[i] << ',' << numerical[i] << ',' << exact[i] << ',' << abs_error[i] << '\n';
            }
        }
    }
}

/**
 * Zamienia wynik symulacji na pojedynczy wiersz tabeli zbieżności
 *
 * @param result Wynik symulacji dla jednej metody i jednej siatki
 *
 * @return Wiersz gotowy do zapisu w CSV
 */
convergence_row to_convergence_row(const simulation_result& result)
{
    convergence_row row;
    row.method_name            = result.method_name;
    row.n_intervals            = result.grid.n_intervals;
    row.h                      = result.grid.h;
    row.dt                     = result.grid.dt;
    row.actual_lambda          = result.grid.actual_lambda;
    row.n_time_steps           = result.grid.n_time_steps;
    row.runtime_seconds        = result.runtime_seconds;
    row.max_abs_error_at_t_max = result.max_abs_error_at_t_max;
    return row;
}

/**
 * Uruchamia skrypt Pythona generujący wykresy z zapisanych wyników
 *
 * @param script_path Ścieżka do skryptu wykonywanego przez interpreter
 *
 * @return Prawda, gdy skrypt zakończy się kodem zero
 */
[[nodiscard]] bool run_python_script(const std::string_view script_path)
{
    const std::string command = "PYTHONWARNINGS=ignore ../.venv/bin/python " + std::string(script_path);
    return std::system(command.c_str()) == 0;
}

int main()
{
    try
    {
        const fs::path output_dir = "../results";
        ensure_output_directory(output_dir);

        constexpr double length    = 1.0;
        constexpr double t_max     = 0.5;
        constexpr double diffusion = 1.0;

        const std::vector<std::size_t> convergence_n_values = { 16, 32, 64, 128, 256, 512/*, 1024, 2048, 4096, 8192*/ };

        std::vector<convergence_row> convergence_rows;
        convergence_rows.reserve(convergence_n_values.size() * 3);

        std::print("Uruchamianie badania zbieżności...\n");

        for (cst n_intervals : convergence_n_values)
        {
            auto explicit_future = std::async(std::launch::async, [=]
            {
                return run_explicit_method(n_intervals, length, t_max, diffusion);
            });

            auto thomas_future = std::async(std::launch::async, [=]
            {
                return run_crank_nicolson_thomas(n_intervals, length, t_max, diffusion);
            });

            auto jacobi_future = std::async(std::launch::async, [=]
            {
                return run_crank_nicolson_jacobi(n_intervals, length, t_max, diffusion);
            });

            const auto explicit_result = explicit_future.get();
            const auto thomas_result   = thomas_future.get();
            const auto jacobi_result   = jacobi_future.get();

            convergence_rows.push_back(to_convergence_row(explicit_result));
            convergence_rows.push_back(to_convergence_row(thomas_result));
            convergence_rows.push_back(to_convergence_row(jacobi_result));

            std::print("{:<24} n={:>4}  h={:.6e}  dt={:.6e}  lambda={:.6f}  błąd={:.6e}  czas={:.3f}s\n", explicit_result.method_name, explicit_result.grid.n_intervals, explicit_result.grid.h, explicit_result.grid.dt, explicit_result.grid.actual_lambda, explicit_result.max_abs_error_at_t_max, explicit_result.runtime_seconds);

            std::print("{:<24} n={:>4}  h={:.6e}  dt={:.6e}  lambda={:.6f}  błąd={:.6e}  czas={:.3f}s\n", thomas_result.method_name, thomas_result.grid.n_intervals, thomas_result.grid.h, thomas_result.grid.dt, thomas_result.grid.actual_lambda, thomas_result.max_abs_error_at_t_max, thomas_result.runtime_seconds);

            std::print("{:<24} n={:>4}  h={:.6e}  dt={:.6e}  lambda={:.6f}  błąd={:.6e}  czas={:.3f}s\n", jacobi_result.method_name, jacobi_result.grid.n_intervals, jacobi_result.grid.h, jacobi_result.grid.dt, jacobi_result.grid.actual_lambda, jacobi_result.max_abs_error_at_t_max, jacobi_result.runtime_seconds);
        }

        constexpr std::size_t plot_n_intervals       = 160;
        const auto            explicit_grid_for_plot = make_grid_config(plot_n_intervals, length, t_max, diffusion, explicit_target_lambda);

        const std::vector<std::size_t> snapshot_steps = make_even_snapshot_steps(explicit_grid_for_plot.n_time_steps, snapshot_count); // Migawek pilnujemy tak aby były rozłożone możliwie równomiernie

        std::print("\nUruchamianie wybranej siatki dla wykresów rozwiązania/czasu...\n");

        auto explicit_plot_future = std::async(std::launch::async, [&]
        {
            return run_explicit_method(plot_n_intervals, length, t_max, diffusion, snapshot_steps);
        });

        auto thomas_plot_future = std::async(std::launch::async, [&]
        {
            return run_crank_nicolson_thomas(plot_n_intervals, length, t_max, diffusion, snapshot_steps);
        });

        auto jacobi_plot_future = std::async(std::launch::async, [&]
        {
            return run_crank_nicolson_jacobi(plot_n_intervals, length, t_max, diffusion, snapshot_steps);
        });

        std::vector<simulation_result> plot_results;
        plot_results.reserve(3);
        plot_results.push_back(explicit_plot_future.get());
        plot_results.push_back(thomas_plot_future.get());
        plot_results.push_back(jacobi_plot_future.get());

        // Zapisywanie trzech plików CSV potrzebnych do analizy i wykresów
        write_convergence_csv(output_dir / "convergence.csv", convergence_rows);
        write_time_error_csv(output_dir / "time_error.csv", plot_results);
        write_snapshots_csv(output_dir / "snapshots.csv", plot_results);

        std::print("\nZapisane pliki:\n");
        std::print("- {}\n", ( output_dir / "convergence.csv" ).string());
        std::print("- {}\n", ( output_dir / "time_error.csv" ).string());
        std::print("- {}\n", ( output_dir / "snapshots.csv" ).string());

        // Skrypt Pythona tworzy wykresy na podstawie wygenerowanych plików CSV
        if (!run_python_script("../plot_results.py"))
        {
            std::print("Nie udało się uruchomić plot_results.py!\n");
        }

        std::print("\nOczekiwany wynik teoretyczny:\n");
        std::print("- Metoda Bezpośrednia  -> O(dt + h^2), a przy stałym lambda dt ~ h^2, więc łącznie ~ O(h^2)\n");
        std::print("- Crank-Nicolson -> O(dt^2 + h^2), a przy stałym lambda dt ~ h^2, więc łącznie ~ O(h^2)\n");
        std::print("- Jacobi zmienia tylko solver liniowy, a nie rząd dyskretyzacji PDE, o ile zbiega dostatecznie dokładnie\n");

        return 0;
    }
    catch (const std::exception& e)
    {
        std::print("błąd: {}\n", e.what());
        return 1;
    }
}
