#include <print>
#include <cmath>
#include <vector>
#include <fstream>

/**
 * Przechowuje węzły interpolacji oraz współczynniki postaci Newtona
 */
struct Interpolant
{
    std::vector<double> nodes;
    std::vector<double> coefficients;
};

/**
 * Zwraca wartość funkcji testowej
 *
 * @param x Argument funkcji
 * @param a Parametr kształtu
 *
 * @return Wartość funkcji w punkcie x
 */
double target_function(const double x, const double a)
{
    return x / ( 1.0 + a * std::pow(x, 4.0) );
}

/**
 * Generuje n węzłów równoodległych na przedziale [-1, 1]
 *
 * @param count Liczba węzłów
 *
 * @return Wektor węzłów równoodległych
 */
std::vector<double> equidistant_nodes(const std::size_t count)
{
    std::vector<double> nodes(count);

    if (count == 1)
    {
        nodes[0] = 0.0;
        return nodes;
    }

    constexpr double left  = -1.0;
    constexpr double right = 1.0;
    const double     step  = ( right - left ) / static_cast<double>(count - 1);

    for (std::size_t i = 0; i < count; ++i)
    {
        nodes[i] = left + static_cast<double>(i) * step;
    }

    return nodes;
}

/**
 * Generuje n węzłów Czebyszewa na przedziale [-1, 1]
 *
 * @param count Liczba węzłów
 *
 * @return Wektor węzłów Czebyszewa
 */
std::vector<double> chebyshev_nodes(const std::size_t count)
{
    std::vector<double> nodes(count);

    for (std::size_t k = 0; k < count; ++k)
    {
        const double angle = ( 2.0 * static_cast<double>(k) + 1.0 ) * std::numbers::pi / ( 2.0 * static_cast<double>(count) );
        nodes[k]           = std::cos(angle);
    }

    std::sort(nodes.begin(), nodes.end());

    return nodes;
}

/**
 * Oblicza ilorazy różnicowe potrzebne do postaci Newtona
 *
 * @param x Węzły interpolacji
 * @param y Wartości funkcji w węzłach
 *
 * @return Współczynniki postaci Newtona
 */
std::vector<double> divided_differences(const std::vector<double>& x, const std::vector<double>& y)
{
    const std::size_t   n      = x.size();
    std::vector<double> coeffs = y;

    for (std::size_t j = 1; j < n; ++j)
    {
        for (std::ptrdiff_t i = static_cast<std::ptrdiff_t>(n) - 1; i >= static_cast<std::ptrdiff_t>(j); --i)
        {
            const auto        current  = static_cast<std::size_t>(i);
            const std::size_t previous = current - j;
            coeffs[current]            = ( coeffs[current] - coeffs[previous] ) / ( x[current] - x[previous] );
        }
    }

    return coeffs;
}

/**
 * Konstruuje wielomian interpolacyjny w postaci Newtona
 *
 * @param count Liczba węzłów
 * @param a Parametr funkcji testowej
 * @param use_chebyshev Informacja czy użyć węzłów Czebyszewa
 *
 * @return Gotowy interpolant z węzłami i współczynnikami
 */
Interpolant build_interpolant(const std::size_t count, const double a, const bool use_chebyshev)
{
    Interpolant result;

    result.nodes = use_chebyshev ? chebyshev_nodes(count) : equidistant_nodes(count);
    result.coefficients.resize(count);

    for (std::size_t i = 0; i < count; ++i)
    {
        result.coefficients[i] = target_function(result.nodes[i], a);
    }

    result.coefficients = divided_differences(result.nodes, result.coefficients);

    return result;
}

/**
 * Oblicza wartość wielomianu Newtona w punkcie x
 *
 * @param x Punkt obliczeń
 * @param nodes Węzły interpolacji
 * @param coefficients Współczynniki Newtona
 *
 * @return Wartość wielomianu interpolacyjnego
 */
double evaluate_newton(double x, const std::vector<double>& nodes, const std::vector<double>& coefficients)
{
    if (coefficients.empty())
    {
        return 0.0;
    }

    double value = coefficients.back();

    for (std::ptrdiff_t i = static_cast<std::ptrdiff_t>(coefficients.size()) - 2; i >= 0; --i)
    {
        const auto index = static_cast<std::size_t>(i);
        value            = value * ( x - nodes[index] ) + coefficients[index];
    }

    return value;
}

/**
 * Wylicza maksymalny błąd bezwzględny na siatce punktów
 *
 * @param exact Wartości dokładne
 * @param approximation Wartości przybliżone
 *
 * @return Maksymalny błąd bezwzględny
 */
double max_abs_error(const std::vector<double>& exact, const std::vector<double>& approximation)
{
    double result = 0.0;

    for (std::size_t i = 0; i < exact.size(); ++i)
    {
        result = std::max(result, std::fabs(exact[i] - approximation[i]));
    }

    return result;
}

/**
 * Zapisuje dane do pliku CSV dla wykresów
 *
 * @param path Ścieżka pliku
 * @param xs Punkty siatki
 * @param exact Wartości funkcji dokładnej
 * @param equidistant Wartości interpolacji na węzłach równoodległych
 * @param chebyshev Wartości interpolacji na węzłach Czebyszewa
 */
void write_plot_csv(const std::string& path, const std::vector<double>& xs, const std::vector<double>& exact, const std::vector<double>& equidistant, const std::vector<double>& chebyshev)
{
    std::ofstream file(path);

    file << "x,exact,equidistant,chebyshev\n";

    for (std::size_t i = 0; i < xs.size(); ++i)
    {
        file << xs[i] << ',' << exact[i] << ',' << equidistant[i] << ',' << chebyshev[i] << '\n';
    }
}

/**
 * Zapisuje węzły i wartości funkcji do pliku CSV
 *
 * @param path Ścieżka pliku
 * @param nodes Węzły interpolacji
 * @param a Parametr funkcji testowej
 */
void write_nodes_csv(const std::string& path, const std::vector<double>& nodes, double a)
{
    std::ofstream file(path);

    file << "x,f\n";

    for (double x : nodes)
    {
        file << x << ',' << target_function(x, a) << '\n';
    }
}

[[nodiscard]] bool run_python_script(const std::string_view script_path)
{
    const std::string command = "PYTHONWARNINGS=ignore ../.venv/bin/python " + std::string(script_path);
    return std::system(command.c_str()) == 0;
}

int main()
{
    constexpr double      a            = 25;
    constexpr std::size_t node_count   = 11;
    constexpr std::size_t sample_count = 2001;

    const auto [nodes, coefficients]   = build_interpolant(node_count, a, false);
    const auto [nodes2, coefficients2] = build_interpolant(node_count, a, true);

    std::vector<double> xs;
    std::vector<double> exact_values;
    std::vector<double> equidistant_values;
    std::vector<double> chebyshev_values;

    xs.reserve(sample_count);
    exact_values.reserve(sample_count);
    equidistant_values.reserve(sample_count);
    chebyshev_values.reserve(sample_count);

    constexpr double left  = -1.0;
    constexpr double right = 1.0;
    constexpr double step  = ( right - left ) / static_cast<double>(sample_count - 1);

    for (std::size_t i = 0; i < sample_count; ++i)
    {
        const double x = left + static_cast<double>(i) * step;

        xs.push_back(x);
        exact_values.push_back(target_function(x, a));
        equidistant_values.push_back(evaluate_newton(x, nodes, coefficients));
        chebyshev_values.push_back(evaluate_newton(x, nodes2, coefficients2));
    }

    const double error_equidistant = max_abs_error(exact_values, equidistant_values);
    const double error_chebyshev   = max_abs_error(exact_values, chebyshev_values);

    write_plot_csv("../output/csv/runge_data.csv", xs, exact_values, equidistant_values, chebyshev_values);
    write_nodes_csv("../output/csv/equidistant_nodes.csv", nodes, a);
    write_nodes_csv("../output/csv/chebyshev_nodes.csv", nodes2, a);

    std::print("Parametr a = {:.6f}\n", a);
    std::print("Liczba węzłów = {}\n", node_count);
    std::print("Maksymalny błąd dla węzłów równoodległych = {:.12f}\n", error_equidistant);
    std::print("Maksymalny błąd dla węzłów Czebyszewa = {:.12f}\n", error_chebyshev);
    std::print("Zapisano pliki runge_data.csv, equidistant_nodes.csv i chebyshev_nodes.csv\n");

    if (!run_python_script("../create_graphs.py"))
    {
        std::print("Failed to run create_graphs.py!\n");
    }
}
