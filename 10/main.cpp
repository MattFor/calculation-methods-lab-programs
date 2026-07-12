#include <print>
#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip>
#include <cstdlib>

// Stała końca przedziału obliczeń
constexpr double T = 1.0;

// Warunek początkowy rozwiązania
constexpr double y_init = 2.0;

/**
 * Współczynnik stojący przy wyrażeniu y(t) - 1 w równaniu różniczkowym
 *
 * @param t Aktualna wartość zmiennej niezależnej
 *
 * @return Wartość współczynnika a(t)
 */
double a(const double t)
{
    return ( 100.0 * t + 10.0 ) / ( t + 1.0 );
}

/**
 * Dokładne rozwiązanie analityczne służące do porównania z rozwiązaniami numerycznymi
 *
 * @param t Aktualna wartość zmiennej niezależnej
 *
 * @return Wartość rozwiązania dokładnego w punkcie t
 */
double exact(const double t)
{
    return 1.0 + std::pow(1.0 + t, 90.0) * std::exp(-100.0 * t);
}

// Rodzaje metod numerycznych użytych do rozwiązania zadania
enum class Method
{
    ExplicitEuler,
    ImplicitEuler,
    Trapezoid
};

/**
 * Zwraca nazwę tekstową wybranej metody numerycznej
 *
 * @param m Wybrana metoda
 *
 * @return Napis opisujący metodę
 */
const char* method_name(const Method m)
{
    switch (m)
    {
        case Method::ExplicitEuler:
        {
            return "Euler jawny";
        }

        case Method::ImplicitEuler:
        {
            return "Euler pośredni";
        }

        case Method::Trapezoid:
        {
            return "Metoda trapezów";
        }
    }

    return "Nieznana metoda";
}

/**
 * Wykonuje jeden krok metody numerycznej dla zadanego punktu siatki
 *
 * @param m Wybrana metoda całkowania
 * @param t Aktualny czas
 * @param y Aktualna wartość przybliżenia rozwiązania
 * @param h Krok siatki czasowej
 *
 * @return Wartość rozwiązania po jednym kroku metody
 */
double step(const Method m, const double t, const double y, const double h)
{
    const double a_n   = a(t);
    const double t_np1 = t + h;
    const double a_np1 = a(t_np1);

    switch (m)
    {
        case Method::ExplicitEuler:
        {
            return y - h * a_n * ( y - 1.0 );
        }

        case Method::ImplicitEuler:
        {
            return ( y + h * a_np1 ) / ( 1.0 + h * a_np1 );
        }

        case Method::Trapezoid:
        {
            return ( y * ( 1.0 - 0.5 * h * a_n ) + 0.5 * h * ( a_n + a_np1 ) ) / ( 1.0 + 0.5 * h * a_np1 );
        }
    }

    return y;
}

/**
 * Oblicza maksymalny błąd bezwzględny na całym przedziale dla danej metody i kroku
 *
 * @param m Wybrana metoda numeryczna
 * @param h Krok siatki czasowej
 *
 * @return Maksymalna wartość bezwzględnego błędu na przedziale [0, T]
 */
double max_abs_error(const Method m, const double h)
{
    const std::size_t n = std::llround(T / h);

    double t       = 0.0;
    double y       = y_init;
    double max_err = 0.0;

    for (std::size_t i = 0; i <= n; ++i)
    {
        const double err = std::abs(y - exact(t));
        max_err          = std::max(max_err, err);

        if (i < n)
        {
            y = step(m, t, y, h);
            t += h;
        }
    }

    return max_err;
}

/**
 * Zapisuje do pliku CSV przebieg rozwiązania numerycznego oraz dokładnego
 *
 * @param file Ścieżka do pliku wyjściowego
 * @param m Wybrana metoda numeryczna
 * @param h Krok siatki czasowej
 */
void write_solution_csv(const std::string_view file, const Method m, double h)
{
    const std::size_t n = std::llround(T / h);

    std::ofstream out{ std::string(file) };
    out << std::setprecision(17) << std::fixed;
    out << "t,exact,numerical,error\n";

    double t       = 0.0;
    double y       = y_init;
    double max_err = 0.0;

    for (std::size_t i = 0; i <= n; ++i)
    {
        const double ye  = exact(t);
        const double err = std::abs(y - ye);
        max_err          = std::max(max_err, err);

        out << t << ',' << ye << ',' << y << ',' << err << '\n';

        if (i < n)
        {
            y = step(m, t, y, h);
            t += h;
        }
    }

    std::print("{}: h = {:.6g}, max|błąd| = {:.6e}\n", method_name(m), h, max_err);
}

/**
 * Zapisuje dane do wykresu zbieżności dla trzech metod numerycznych
 *
 * @param file Ścieżka do pliku wyjściowego
 */
void write_convergence_csv(const std::string_view file)
{
    std::ofstream out{ std::string(file) };
    out << std::setprecision(17) << std::scientific;
    out << "h,explicit_euler,implicit_euler,trapezoid\n";

    // Generuje kolejne kroki jako potęgi dwójki od 2^-6 do 2^-20
    for (int k = 6; k <= 20; ++k)
    {
        const double h  = std::ldexp(1.0, -k);
        const double e1 = max_abs_error(Method::ExplicitEuler, h);
        const double e2 = max_abs_error(Method::ImplicitEuler, h);
        const double e3 = max_abs_error(Method::Trapezoid, h);

        out << h << ',' << e1 << ',' << e2 << ',' << e3 << '\n';
    }
}

[[nodiscard]] bool run_python_script(const std::string_view script_path)
{
    const std::string command = "PYTHONWARNINGS=ignore ../.venv/bin/python " + std::string(script_path);
    return std::system(command.c_str()) == 0;
}

int main()
{
    std::print("Rozwiązuję równanie na przedziale [0, 1].\n");

    write_solution_csv("../output/csv/explicit_stable.csv", Method::ExplicitEuler, 0.01);
    write_solution_csv("../output/csv/explicit_border.csv", Method::ExplicitEuler, 0.02);
    write_solution_csv("../output/csv/explicit_unstable.csv", Method::ExplicitEuler, 0.03);

    write_solution_csv("../output/csv/implicit_euler.csv", Method::ImplicitEuler, 0.02);
    write_solution_csv("../output/csv/trapezoid.csv", Method::Trapezoid, 0.02);

    write_convergence_csv("../output/csv/convergence.csv");

    std::print("Gotowe pliki CSV zapisane w katalogu bieżącym.\n");

    if (!run_python_script("../create_graphs.py"))
    {
        std::print("Failed to run create_graphs.py!\n");
    }
}
