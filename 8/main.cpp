#include <cmath>
#include <hprint>
#include <vector>
#include <fstream>

/**
 * Funkcja testowa f(x) = cos(x), dla której liczona jest pochodna numerycznie
 *
 * @tparam T Typ zmiennych rzeczywistych, np. double albo long double
 * @param x Punkt, w którym obliczana jest wartość funkcji
 *
 * @return Wartość cos(x) w typie T
 */
template <class T>
[[nodiscard]] constexpr T f(T x)
{
    return std::cos(x);
}

/**
 * Pochodna funkcji f(x) = cos(x)
 *
 * @tparam T Typ zmiennych rzeczywistych
 * @param x Punkt, w którym obliczana jest dokładna pochodna
 *
 * @return Dokładna wartość f'(x) = -sin(x)
 */
template <class T>
[[nodiscard]] constexpr T exact_d1(T x)
{
    return -std::sin(x);
}

/**
 * Dwupunktowe przybliżenie jednostronne do przodu
 * Korzysta z wartości funkcji w punktach x oraz x + h
 * Wzór: f'(x) ≈ (f(x+h) - f(x)) / h
 *
 * Jest to przybliżenie rzędu O(h), więc po zmniejszeniu kroku h
 * błąd powinien maleć w przybliżeniu liniowo
 *
 * @tparam T Typ zmiennych rzeczywistych
 * @param x Punkt, w którym liczymy pochodną
 * @param h Krok siatki
 *
 * @return Przybliżenie pochodnej metodą dwupunktową jednostronną do przodu
 */
template <class T>
[[nodiscard]] constexpr T d1_forward_2(T x, T h)
{
    return ( f(x + h) - f(x) ) / h;
}


/**
 * Trzypunktowe przybliżenie jednostronne do przodu.
 * Korzysta z wartości funkcji w punktach x, x+h oraz x+2h.
 * Wzór: f'(x) ≈ (-3f(x) + 4f(x + h) - f(x + 2h)) / (2h)
 *
 * Jest to przybliżenie rzędu O(h^2), więc dokładniejsze od wzoru dwupunktowego.
 *
 * @tparam T Typ zmiennych rzeczywistych
 * @param x Punkt, w którym liczymy pochodną
 * @param h Krok siatki
 *
 * @return Przybliżenie pochodnej metodą trzypunktową jednostronną do przodu
 */
template <class T>
[[nodiscard]] constexpr T d1_forward_3(T x, T h)
{
    return ( -T(3) * f(x) + T(4) * f(x + h) - f(x + T(2) * h) ) / ( T(2) * h );
}

/**
 * Dwupunktowe przybliżenie jednostronne do tyłu.
 * Korzysta z wartości funkcji w punktach x oraz x-h.
 * Wzór: f'(x) ≈ (f(x) - f(x - h)) / h
 *
 * Jest to przybliżenie rzędu O(h), analogiczne do wzoru do przodu,
 * ale wykorzystujące punkty po lewej stronie x
 *
 * @tparam T Typ zmiennych rzeczywistych
 * @param x Punkt, w którym liczymy pochodną
 * @param h Krok siatki
 *
 * @return Przybliżenie pochodnej metodą dwupunktową jednostronną do tyłu
 */
template <class T>
[[nodiscard]] constexpr T d1_backward_2(T x, T h)
{
    return ( f(x) - f(x - h) ) / h;
}

/**
 * Trzypunktowe przybliżenie jednostronne do tyłu
 * Korzysta z wartości funkcji w punktach x, x-h oraz x-2h
 * Wzór: f'(x) ≈ (3f(x) - 4f(x - h) + f(x - 2h)) / (2h)
 *
 * Jest to przybliżenie rzędu O(h^2), więc dokładniejsze od wersji dwupunktowej
 *
 * @tparam T Typ zmiennych rzeczywistych
 * @param x Punkt, w którym liczymy pochodną
 * @param h Krok siatki
 *
 * @return Przybliżenie pochodnej metodą trzypunktową jednostronną do tyłu
 */
template <class T>
[[nodiscard]] constexpr T d1_backward_3(T x, T h)
{
    return ( T(3) * f(x) - T(4) * f(x - h) + f(x - T(2) * h) ) / ( T(2) * h );
}

/**
 * Dwupunktowe przybliżenie centralne.
 * Korzysta z wartości funkcji po obu stronach punktu x: x-h i x+h.
 * Wzór: f'(x) ≈ (f(x+h) - f(x-h)) / (2h)
 *
 * Jest to przybliżenie rzędu O(h^2), a więc tak samo dokładne asymptotycznie
 * jak wzory trzypunktowe, ale z mniejszym błędem stałej w wielu sytuacjach
 *
 * @tparam T Typ zmiennych rzeczywistych
 * @param x Punkt, w którym liczymy pochodną
 * @param h Krok siatki
 *
 * @return Przybliżenie pochodnej metodą centralną
 */
template <class T>
[[nodiscard]] constexpr T d1_central_2(T x, T h)
{
    return ( f(x + h) - f(x - h) ) / ( T(2) * h );
}

/**
 * Opis jednej metody przybliżania pochodnej
 * Przechowuje nazwę oraz wskaźnik do funkcji realizującej obliczenia
 *
 * @tparam T Typ zmiennych rzeczywistych
 */
template <class T>
struct Approximation
{
    std::string_view name;        ///< nazwa metody do wypisywania w tabeli i wykresach
    T (*             eval)(T, T); ///< wskaźnik do funkcji liczącej aproksymację pochodnej
};

/**
 * Jeden wiersz wyników dla konkretnego kroku h
 * Zawiera sam krok oraz błędy bezwzględne dla wszystkich metod
 *
 * @tparam T Typ zmiennych rzeczywistych
 */
template <class T>
struct Row
{
    T                h{};      ///< krok siatki
    std::array<T, 5> errors{}; ///< błędy bezwzględne dla 5 metod
};

/**
 * Pomocnicze rzutowanie do double
 * Dla logarytmów i dobrego wypisania
 *
 * @tparam T Typ wejściowy
 * @param x Wartość do zamiany
 *
 * @return X rzutowane na double
 */
template <class T>
[[nodiscard]] double to_double(T x)
{
    return static_cast<double>(x);
}

/**
 * Dla jednego punktu x i kolejnych kroków h buduje tabelę wyników
 * Dla każdego h:
 *  - liczy dokładną pochodną
 *  - liczy pięć przybliżeń różnicowych
 *  - zapisuje błędy bezwzględne |approx - exact|
 *
 * Dzięki temu można później:
 *  - wypisać tabelę na konsolę
 *  - dopasować nachylenie na wykresie log-log
 *  - wyeksportować dane do CSV
 *
 * @tparam T Typ zmiennych rzeczywistych
 * @param x Punkt, w którym badamy pochodną
 * @param n_steps Liczba kolejnych kroków siatki
 *
 * @return Wektor wierszy z wynikami dla kolejnych kroków h
 */
template <class T>
[[nodiscard]] std::vector<Row<T>> compute_rows(T x, std::size_t n_steps)
{
    const std::array<Approximation<T>, 5> methods{ { { "forward 2-point", &d1_forward_2<T> }, { "forward 3-point", &d1_forward_3<T> }, { "backward 2-point", &d1_backward_2<T> }, { "backward 3-point", &d1_backward_3<T> }, { "central 2-point", &d1_central_2<T> } } };

    std::vector<Row<T>> rows;
    rows.reserve(n_steps);

    // Startowy krok wybrany tak, żeby w punktach 0, pi/4 i pi/2
    // wszystkie używane wzory dało się bezpiecznie policzyć
    T h = std::numbers::pi_v<T> / T(8);

    for (std::size_t i = 0; i < n_steps; ++i)
    {
        Row<T> row{};
        row.h = h;

        // Dokładna wartość pochodnej w punkcie x
        const T exact = exact_d1<T>(x);

        // Błędy dla każdej metody
        for (std::size_t m = 0; m < methods.size(); ++m)
        {
            const T approx = methods[m].eval(x, h);
            row.errors[m]  = std::abs(approx - exact);
        }

        rows.push_back(row);

        // Krok zmniejszany o połowę, żeby obserwować zbieżność i wpływ błędów maszynowych
        h /= T(2);
    }

    return rows;
}

/**
 * Szacuje doświadczalny rząd dokładności na podstawie wykresu log10(|błąd|) względem log10(h)
 * W praktyce wykonujemy regresję liniową na punktach, dla których błąd jeszcze maleje jak w reżimie błędu obcięcia
 * Gdy pojawia się dominacja błędów maszynowych dalsze punkty nie są dobre do wyznaczania rzędu
 *
 * @tparam T Typ zmiennych rzeczywistych
 * @param rows Tabela wyników dla kolejnych kroków
 * @param method_index Indeks metody, dla której liczymy nachylenie
 *
 * @return Nachylenie prostej w skali log-log, czyli doświadczalny rząd dokładności
 */
template <class T>
[[nodiscard]] double estimate_slope_loglog(const std::vector<Row<T>>& rows, std::size_t method_index)
{
    std::vector<std::pair<double, double>> pts;
    pts.reserve(rows.size());

    // Wybieramy tylko punkty sensowne do logarytmowania
    // co znaczy: h > 0, błąd > 0 oraz oba są skończone
    for (const auto& r : rows)
    {
        const double h = to_double(r.h);
        if (const double e = to_double(r.errors[method_index]); h > 0.0 && e > 0.0 && std::isfinite(h) && std::isfinite(e))
        {
            pts.emplace_back(std::log10(h), std::log10(e));
        }
    }

    if (pts.size() < 2)
    {
        return std::numeric_limits<double>::quiet_NaN();
    }

    // Brany pod uwagę jest tylko ten fragment, dla którego błąd jeszcze maleje
    // Gdy zaczyna rosnąć, najczęściej wchodzi już wpływ błędów maszynowych -> std::size_t end = pts.size();
    std::size_t end = pts.size();

    for (std::size_t i = 1; i < rows.size(); ++i)
    {
        if (to_double(rows[i].errors[method_index]) > to_double(rows[i - 1].errors[method_index]))
        {
            end = i;
            break;
        }
    }

    // Jeśli punkty są bardzo krótkie, brane jest przynajmniej kilka pierwszych
    if (end < 3)
    {
        end = std::min<std::size_t>(pts.size(), 6);
    }

    // Regresja liniowa y = a x + b, gdzie x = log10(h), y = log10(|błędu|)
    double      sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0;
    std::size_t n  = 0;

    for (std::size_t i = 0; i < end; ++i)
    {
        const auto [x, y] = pts[i];
        sx                += x;
        sy                += y;
        sxx               += x * x;
        sxy               += x * y;
        ++n;
    }

    const auto   dn          = static_cast<double>(n);
    const double denominator = dn * sxx - sx * sx;
    if (std::abs(denominator) < 1e-30)
    {
        return std::numeric_limits<double>::quiet_NaN();
    }

    return ( dn * sxy - sx * sy ) / denominator;
}

/**
 * Zapisuje dane do pliku CSV
 * Każdy wiersz zawiera: krok h, log10(h) i błędy dla wszystkich metod
 *
 * @tparam T Typ zmiennych rzeczywistych
 * @param filename Nazwa pliku wyjściowego bez katalogu
 * @param rows Wszystkie wiersze wyników do zapisania
 */
template <class T>
void write_csv(const std::string_view filename, const std::vector<Row<T>>& rows)
{
    std::ofstream out{ "../data/" + std::string(filename) };
    out << "h,log10_h,forward_2,forward_3,backward_2,backward_3,central_2\n";

    for (const auto& r : rows)
    {
        out << to_double(r.h) << ',' << std::log10(to_double(r.h)) << ',' << to_double(r.errors[0]) << ',' << to_double(r.errors[1]) << ',' << to_double(r.errors[2]) << ',' << to_double(r.errors[3]) << ',' << to_double(r.errors[4]) << '\n';
    }
}

/**
 * Wykonuje pełne badanie dla jednego punktu x i jednego typu liczbowego T
 * - liczy tabelę błędów
 * - wypisuje wyniki na konsolę
 * - wyznacza doświadczalne rzędy dokładności
 * - zapisuje CSV do folderu data/
 *
 * @tparam T Typ zmiennych rzeczywistych
 * @param type_name Nazwa typu, używana w nazwach plików i nagłówkach
 * @param point_name Etykieta punktu, np. "0", "pi_over_4", "pi_over_2"
 * @param x Punkt, w którym badamy pochodną
 * @param n_steps Liczba kroków h
 */
template <class T>
void run_study_for_point(std::string_view type_name, std::string_view point_name, T x, const std::size_t n_steps)
{
    const auto rows = compute_rows<T>(x, n_steps);

    hp::printf("{} | x = {}", type_name, point_name);
    std::print("exact derivative = {:.20e}\n", to_double(exact_d1<T>(x)));
    std::print("{:<22} {:<22} {:<22} {:<22} {:<22} {:<22}\n", "h", "f2 err", "f3 err", "b2 err", "b3 err", "c2 err");

    for (const auto& r : rows)
    {
        std::print("{:<14.16e} {:>20.16e} {:>20.16e} {:>20.16e} {:>20.16e} {:>20.16e}\n", to_double(r.h), to_double(r.errors[0]), to_double(r.errors[1]), to_double(r.errors[2]), to_double(r.errors[3]), to_double(r.errors[4]));
    }

    constexpr std::array<std::string_view, 5> names{ "forward 2-point", "forward 3-point", "backward 2-point", "backward 3-point", "central 2-point" };

    std::print("\nEstimated orders from log10|error| vs log10 h:\n");
    for (std::size_t m = 0; m < names.size(); ++m)
    {
        const double slope = estimate_slope_loglog(rows, m);
        std::print("{:<18} p ≈ {:.16f}\n", names[m], slope);
    }

    const std::string file = std::string(type_name) + "_" + std::string(point_name) + ".csv";
    write_csv(file, rows);
    std::print("\nCSV written: {}\n\n", file);
}

/**
 * Uruchamia badanie dla wszystkich trzech punktów testowych
 * dla jednego typu zmiennych rzeczywistych
 *
 * @tparam T Typ zmiennych rzeczywistych
 * @param type_name Nazwa typu, np. "double" lub "long_double"
 * @param n_steps Liczba kolejnych kroków h
 */
template <class T>
void run_all_for_type(const std::string_view type_name, const std::size_t n_steps)
{
    const T pi = std::numbers::pi_v<T>;

    // Punkt lewy -> używa wzorów jednostronnych do przodu
    run_study_for_point<T>(type_name, "0", T(0), n_steps);

    // Punkt środkowy -> wzór centralny
    run_study_for_point<T>(type_name, "pi_over_4", pi / T(4), n_steps);

    // Punkt prawy -> używa wzorów jednostronnych do tyłu
    run_study_for_point<T>(type_name, "pi_over_2", pi / T(2), n_steps);
}

[[nodiscard]] bool run_python_script(const std::string_view script_path)
{
    const std::string command = "PYTHONWARNINGS=ignore ../.venv/bin/python " + std::string(script_path);
    return std::system(command.c_str()) == 0;
}

int main()
{
    constexpr std::size_t n_steps = 30;

    std::print("Numerical differentiation of f(x)=cos(x)\n");
    std::print("Approximations: 2-point and 3-point one-sided, plus central 2-point.\n");
    std::print("Target interval: [0, pi/2]\n");

    run_all_for_type<double>("double", n_steps);
    run_all_for_type<long double>("long_double", n_steps);

    if (!run_python_script("../data/scripts/create_sheet.py"))
    {
        std::print("Failed to run create_sheet.py!\n");
    }

    if (!run_python_script("../data/scripts/create_graphs.py"))
    {
        std::print("Failed to run create_graphs.py!\n");
    }

    /**
     * Błąd totalny: O(h^p) + O(ϵ / h)
     * Skoro dzielone jest przez h, to jeśli jest za małe to jest utrata precyzji
     * I wtedy błąd strasznie rośnie
     */
}
