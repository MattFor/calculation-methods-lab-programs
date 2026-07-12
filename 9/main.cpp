#include <print>
#include <cmath>
#include <vector>
#include <fstream>

// Przechowuje punkty siatki oraz odpowiadające im wartości rozwiązania
struct RozwiazanieSiatki
{
    std::vector<double> x;
    std::vector<double> u;
};

// Przechowuje cztery wektory potrzebne do metody Thomasa
// dolna przekątna, główna przekątna, górna przekątna oraz prawa strona
struct MacierzTrojdiagonalna
{
    std::vector<double> dolna;
    std::vector<double> glowna;
    std::vector<double> gorna;
    std::vector<double> prawa;
};

// Rozwiązuje układ trójdiagonalny algorytmem Thomasa
std::vector<double> rozwiaz_thomasa(const MacierzTrojdiagonalna& macierz)
{
    // Liczba równań w układzie
    const std::size_t   n = macierz.glowna.size();

    // Wektory pomocnicze używane w eliminacji i podstawianiu wstecznym
    std::vector gorna_po_mod(n, 0.0);
    std::vector prawa_po_mod(n, 0.0);
    std::vector wynik(n, 0.0);

    // Inicjalizacja pierwszego kroku eliminacji
    gorna_po_mod[0] = macierz.gorna[0] / macierz.glowna[0];
    prawa_po_mod[0] = macierz.prawa[0] / macierz.glowna[0];

    // Eliminacja w przód dla kolejnych wierszy macierzy
    for (std::size_t i = 1; i < n; ++i)
    {
        // Oblicza nowy mianownik po eliminacji
        const double mianownik = macierz.glowna[i] - macierz.dolna[i] * gorna_po_mod[i - 1];

        // Wyznacza zmodyfikowaną górną przekątną jeśli nie jesteśmy na końcu układu
        if (i + 1 < n)
        {
            gorna_po_mod[i] = macierz.gorna[i] / mianownik;
        }

        // Wyznacza zmodyfikowaną prawą stronę
        prawa_po_mod[i] = ( macierz.prawa[i] - macierz.dolna[i] * prawa_po_mod[i - 1] ) / mianownik;
    }

    // Podstawianie wsteczne od ostatniego równania
    wynik[n - 1] = prawa_po_mod[n - 1];
    for (std::size_t i = n - 1; i-- > 0;)
    {
        wynik[i] = prawa_po_mod[i] - gorna_po_mod[i] * wynik[i + 1];
    }

    return wynik;
}

// Zwraca wartość wyrazu wymuszającego w równaniu różniczkowym
double wymuszenie(const double x)
{
    return 0.5 * x * x * x;
}

// Zwraca rozwiązanie dokładne używane do porównania z wynikami numerycznymi
double rozwiazanie_dokladne(const double x)
{
    const double pierwiastek_5 = std::sqrt(5.0);
    const double r1            = -1.0 + pierwiastek_5;
    const double r2            = -1.0 - pierwiastek_5;

    // Część szczególna rozwiązania analitycznego
    const double up = ( ( ( ( x / 8.0 ) + 3.0 / 16.0 ) * x ) + 3.0 / 8.0 ) * x + 9.0 / 32.0;

    // Stałe wyznaczone z warunków brzegowych
    constexpr double A  = 55.0 / 32.0;
    constexpr double B  = -95.0 / 32.0;
    const double     e1 = std::exp(r1);
    const double     e2 = std::exp(r2);
    const double     c1 = ( B - A * e2 ) / ( e1 - e2 );
    const double     c2 = A - c1;

    // Składa rozwiązanie ogólne z części jednorodnej i szczególnej
    return c1 * std::exp(r1 * x) + c2 * std::exp(r2 * x) + up;
}

// Rozwiązuje zagadnienie brzegowe metodą różnic skończonych
RozwiazanieSiatki roznice_skonczone(const std::size_t liczba_odcinkow)
{
    // Krok siatki oraz liczba punktów wewnętrznych
    const double      krok                = 1.0 / static_cast<double>(liczba_odcinkow);
    const std::size_t liczba_wewnetrznych = liczba_odcinkow - 1;

    // Przygotowanie układu trójdiagonalnego dla punktów wewnętrznych
    MacierzTrojdiagonalna macierz;
    macierz.dolna.assign(liczba_wewnetrznych, 0.0);
    macierz.glowna.assign(liczba_wewnetrznych, 0.0);
    macierz.gorna.assign(liczba_wewnetrznych, 0.0);
    macierz.prawa.assign(liczba_wewnetrznych, 0.0);

    // Współczynniki dyskretyzacji pochodnych na siatce jednorodnej
    const double odwrotnosc_kroku  = 1.0 / krok;
    const double odwrotnosc_kroku2 = odwrotnosc_kroku * odwrotnosc_kroku;
    const double a                 = odwrotnosc_kroku2 - odwrotnosc_kroku;
    const double b                 = -2.0 * odwrotnosc_kroku2 - 4.0;
    const double c                 = odwrotnosc_kroku2 + odwrotnosc_kroku;

    // Wypełnia układ dla wszystkich punktów wewnętrznych
    for (std::size_t i = 0; i < liczba_wewnetrznych; ++i)
    {
        const double x    = static_cast<double>(i + 1) * krok;
        macierz.dolna[i]  = a;
        macierz.glowna[i] = b;
        macierz.gorna[i]  = c;
        macierz.prawa[i]  = -wymuszenie(x);
    }

    // Uwzględnia warunki brzegowe w pierwszym i ostatnim równaniu
    macierz.prawa[0]                       -= a * 2.0;
    macierz.prawa[liczba_wewnetrznych - 1] -= c * -2.0;

    // Rozwiązuje powstały układ liniowy algorytmem Thomasa
    const std::vector<double> rozwiazanie_wewnetrzne = rozwiaz_thomasa(macierz);

    // Tworzy pełne rozwiązanie na całej siatce
    RozwiazanieSiatki rozwiazanie;
    rozwiazanie.x.resize(liczba_odcinkow + 1);
    rozwiazanie.u.resize(liczba_odcinkow + 1);

    // Wpisuje punkty skrajne oraz wartości wewnętrzne
    for (std::size_t i = 0; i <= liczba_odcinkow; ++i)
    {
        rozwiazanie.x[i] = static_cast<double>(i) * krok;

        if (i == 0)
        {
            rozwiazanie.u[i] = 2.0;
        }
        else if (i == liczba_odcinkow)
        {
            rozwiazanie.u[i] = -2.0;
        }
        else
        {
            rozwiazanie.u[i] = rozwiazanie_wewnetrzne[i - 1];
        }
    }

    return rozwiazanie;
}

// Przechowuje stan układu równań pierwszego rzędu dla metody strzałów
struct Stan
{
    double u;
    double v;
};

// Zwraca prawą stronę układu dla zadania niejednorodnego
Stan prawa_strona_niejednorodna(const double x, const Stan& stan)
{
    return { stan.v, -2.0 * stan.v + 4.0 * stan.u - wymuszenie(x) };
}

// Zwraca prawą stronę układu jednorodnego używanego w metodzie strzałów
Stan prawa_strona_jednorodna(double /*x*/, const Stan& stan)
{
    return { stan.v, -2.0 * stan.v + 4.0 * stan.u };
}

// Wykonuje jeden krok klasycznego schematu Rungego-Kutty czwartego rzędu
template <class PrawaStrona>
Stan krok_rk4(double x, const double krok, const Stan& stan, PrawaStrona prawa_strona)
{
    // Oblicza kolejne przybliżenia nachylenia
    const Stan k1 = prawa_strona(x, stan);
    const Stan s2{ stan.u + 0.5 * krok * k1.u, stan.v + 0.5 * krok * k1.v };
    const Stan k2 = prawa_strona(x + 0.5 * krok, s2);
    const Stan s3{ stan.u + 0.5 * krok * k2.u, stan.v + 0.5 * krok * k2.v };
    const Stan k3 = prawa_strona(x + 0.5 * krok, s3);
    const Stan s4{ stan.u + krok * k3.u, stan.v + krok * k3.v };
    const Stan k4 = prawa_strona(x + krok, s4);

    // Zwraca nowy stan po jednym kroku RK4
    return { stan.u + ( krok / 6.0 ) * ( k1.u + 2.0 * k2.u + 2.0 * k3.u + k4.u ), stan.v + ( krok / 6.0 ) * ( k1.v + 2.0 * k2.v + 2.0 * k3.v + k4.v ) };
}

// Rozwiązuje zagadnienie brzegowe metodą strzałów
RozwiazanieSiatki metoda_strzalow(const std::size_t liczba_odcinkow)
{
    // Krok siatki
    const double krok = 1.0 / static_cast<double>(liczba_odcinkow);

    // Rozwiązanie podstawowe dla układu niejednorodnego
    std::vector<Stan> rozwiazanie_podstawowe(liczba_odcinkow + 1);

    // Rozwiązanie jednorodne potrzebne do dopasowania warunku brzegowego
    std::vector<Stan> rozwiazanie_jednorodne(liczba_odcinkow + 1);

    // Warunki początkowe dla obu przebiegów
    rozwiazanie_podstawowe[0] = { 2.0, 0.0 };
    rozwiazanie_jednorodne[0] = { 0.0, 1.0 };

    // Całkuje oba układy krok po kroku
    for (std::size_t i = 0; i < liczba_odcinkow; ++i)
    {
        const double x                = static_cast<double>(i) * krok;
        rozwiazanie_podstawowe[i + 1] = krok_rk4(x, krok, rozwiazanie_podstawowe[i], prawa_strona_niejednorodna);
        rozwiazanie_jednorodne[i + 1] = krok_rk4(x, krok, rozwiazanie_jednorodne[i], prawa_strona_jednorodna);
    }

    // Wyznacza współczynnik sklejający oba rozwiązania tak, aby spełnić warunek brzegowy w x = 1
    const double wspolczynnik = ( -2.0 - rozwiazanie_podstawowe.back().u ) / rozwiazanie_jednorodne.back().u;

    // Składa końcowe rozwiązanie na całej siatce
    RozwiazanieSiatki rozwiazanie;
    rozwiazanie.x.resize(liczba_odcinkow + 1);
    rozwiazanie.u.resize(liczba_odcinkow + 1);

    // Łączy rozwiązanie podstawowe z jednorodnym
    for (std::size_t i = 0; i <= liczba_odcinkow; ++i)
    {
        rozwiazanie.x[i] = static_cast<double>(i) * krok;
        rozwiazanie.u[i] = rozwiazanie_podstawowe[i].u + wspolczynnik * rozwiazanie_jednorodne[i].u;
    }

    return rozwiazanie;
}

// Oblicza maksymalny błąd bezwzględny względem rozwiązania dokładnego
double maksymalny_blad_bezwzgledny(const RozwiazanieSiatki& rozwiazanie)
{
    double maksimum = 0.0;
    for (std::size_t i = 0; i < rozwiazanie.x.size(); ++i)
    {
        maksimum = std::max(maksimum, std::abs(rozwiazanie.u[i] - rozwiazanie_dokladne(rozwiazanie.x[i])));
    }

    return maksimum;
}

// Zapisuje dane do pliku CSV w formacie używanym przez skrypt Pythona
void zapisz_csv(const std::string& nazwa_pliku, const RozwiazanieSiatki& fd, const RozwiazanieSiatki& strz, const std::vector<std::pair<double, std::pair<double, double>>>& zbieznosc)
{
    std::ofstream plik(nazwa_pliku);
    plik << "typ,n,h,x,fd,shoot,exact,err_fd,err_sh\n";

    // Liczba odcinków i krok siatki dla wykresu porównawczego
    const std::size_t n = fd.x.size() - 1;
    const double      h = 1.0 / static_cast<double>(n);

    // Zapisuje rozwiązanie dla jednego wybranego kroku siatki
    for (std::size_t i = 0; i < fd.x.size(); ++i)
    {
        plik << "solution," << n << ',' << h << ',' << fd.x[i] << ',' << fd.u[i] << ',' << strz.u[i] << ',' << rozwiazanie_dokladne(fd.x[i]) << ',' << std::numeric_limits<double>::quiet_NaN() << ',' << std::numeric_limits<double>::quiet_NaN() << '\n';
    }

    // Zapisuje błędy dla kolejnych siatek do wykresu zbieżności
    for (const auto& [krok_siatki, bledy] : zbieznosc)
    {
        plik << "convergence," << static_cast<std::size_t>(1.0 / krok_siatki) << ',' << krok_siatki << ',' << std::numeric_limits<double>::quiet_NaN() << ',' << std::numeric_limits<double>::quiet_NaN() << ',' << std::numeric_limits<double>::quiet_NaN() << ',' << std::numeric_limits<double>::quiet_NaN() << ',' << bledy.first << ',' << bledy.second << '\n';
    }
}

// Uruchamia skrypt Pythona odpowiedzialny za tworzenie wykresów
[[nodiscard]] bool run_python_script(const std::string_view sciezka_skryptu)
{
    const std::string polecenie = "PYTHONWARNINGS=ignore ../.venv/bin/python " + std::string(sciezka_skryptu);
    return std::system(polecenie.c_str()) == 0;
}

// Główna funkcja programu
int main()
{
    std::print("Równanie: u'' + 2u' - 4u + x^3/2 = 0,  u(0)=2, u(1)=-2\n\n");

    // Nazwa pliku CSV, do którego trafiają wszystkie dane
    const std::string nazwa_pliku = "../data.csv";

    // Liczba odcinków użyta do wykresu porównawczego
    constexpr std::size_t liczba_odcinkow_dla_wykresu = 64;
    const auto            rozwiazanie_fd              = roznice_skonczone(liczba_odcinkow_dla_wykresu);
    const auto            rozwiazanie_strz            = metoda_strzalow(liczba_odcinkow_dla_wykresu);

    // Wektor przechowujący wyniki zbieżności dla kolejnych siatek
    std::vector<std::pair<double, std::pair<double, double>>> zbieznosc;

    // Wypisuje błąd dla wybranej siatki testowej
    std::print("Porównanie na siatce n={} (h=1/{}):\n", liczba_odcinkow_dla_wykresu, liczba_odcinkow_dla_wykresu);
    std::print("  max |błąd| FD   = {:.16e}\n", maksymalny_blad_bezwzgledny(rozwiazanie_fd));
    std::print("  max |błąd| strz = {:.16e}\n\n", maksymalny_blad_bezwzgledny(rozwiazanie_strz));

    // Wypisuje nagłówek tabeli zbieżności
    std::print("Tabela zbieżności (wartości do wykresu log-log):\n");
    std::print("{:>8} {:>14} {:>18} {:>18}\n", "n", "h", "err_FD", "err_STRZ");

    // Liczy błędy dla coraz drobniejszych siatek
    for (std::size_t n = 8; n <= 65536; n *= 2)
    {
        const double h         = 1.0 / static_cast<double>(n);
        const auto   fd        = roznice_skonczone(n);
        const auto   strz      = metoda_strzalow(n);
        const double blad_fd   = maksymalny_blad_bezwzgledny(fd);
        const double blad_strz = maksymalny_blad_bezwzgledny(strz);

        zbieznosc.emplace_back(h, std::pair{ blad_fd, blad_strz });
        std::print("{:8} {:14.6e} {:18.10e} {:18.10e}\n", n, h, blad_fd, blad_strz);
    }

    // Zapisuje wszystkie dane do jednego pliku CSV
    zapisz_csv(nazwa_pliku, rozwiazanie_fd, rozwiazanie_strz, zbieznosc);

    std::print("\nPlik CSV do Pythona: {}\n", nazwa_pliku);

    // Uruchamia skrypt do wygenerowania wykresów
    if (!run_python_script("../create_graphs.py"))
    {
        std::print("Nie udało się uruchomić create_graphs.py\n");
    }
}