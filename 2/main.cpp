#include <print>
#include <cmath>
#include <fstream>
#include <iomanip>

using ld = long double;

ld f_direct(const ld x)
{
    // Bezpośrednia implementacja wzoru z zadania
    // f(x) = x^3 / ( 6 (sinh(x) - x) )
    // sinh(x) pochodzi z biblioteki standardowej <cmath>
    // ale tutaj jest wersja ld -> sinhl()

    return x * x * x / ( 6.0L * ( sinhl(x) - x ) );
}

// Pochodne:
// sinh'(x) = cosh(x)
// sinh''(x) = sinh(x)
// sinh'''(x) = cosh(x)
// dla zera sinh(0) = 0, cosh(0) = 1
// Więc zostają tylko nieparzyste potęgi
// sinh(x) = x + x^3/3! + x^5/5! + ...
// sinh(x) - x = x^3/3! + x^5/5! + ...
// f(x) = x^3/6(sinh(x) - x)
// f(x) = x^3/6(x^3/3! + x^5/5! + ...)
// sinh(x) - x = x^3(1/6 + x^2/120 + ...) -> x^3 przed nawias
// f(x) = x^3/ (6 * x^3 * (1/6 + x^2/120 + ...)) -> podstawienie sinh(x) - x
// f(x) = 1/6 * (1/6 + x^2/120 + ...)
// f(x) = 1/(1 + x^2/20 + ...)   [6 * 1/120 = 1/20]
// Ostateczny szereg to
// f(x) = 1/(1 + x^2/20 + x^4/840 + ...)
ld f_taylor(const ld x)
{
    // Dla x = 0 wynik to dokładnie 1 (lim x->0 f(x) = 1)
    if (x == 0.0L)
    {
        return 1.0L;
    }

    // Obliczam x^2 raz i zapisuje, żeby nie liczyć tego wielokrotnie w pętli
    const ld x2 = x * x;

    // Aktualny wyraz szeregu -> pierwszy wyraz szeregu = 1
    ld wyraz = 1.0L;

    // Suma wszystkich wyrazów szeregu
    ld suma = 1.0L;

    // Dodawanie kolejnych wyrazów rozwiniecia aż do tego aż zbiegnie
    for (int k = 1; k < 10000; ++k)
    {
        // Rekurencyjne wyznaczanie kolejnego wyrazu szeregu
        wyraz *= x2 / ( ( 2 * k + 2 ) * ( 2 * k + 3 ) );

        // Nawet jeśli wyraz nie ma znaczenia do sumy pędzimy dalej bóg dał computing power a my go używamy

        // Dodanie nowego wyrazu do sumy szeregu
        suma += wyraz;
    }

    // Funkcja f(x) jest odwrotnością tej sumy
    return 1.0L / suma;
}

// Automatycznie wybiera najlepszą metodę obliczeń
ld f(const ld x)
{
    // Dla bardzo małych x pojawia się problem numeryczny sinh(x) ≈ x
    // więc w wyrażeniu "sinh(x) - x" odejmujemy dwie prawie równe liczby co powoduje utrate cyfr znaczących

    if (fabsl(x) < 1.0L) //1e-5L)
    {
        // Dla małych x używamy rozwinięcia Taylora
        return f_taylor(x);
    }

    // Dla większych x można bezpiecznie użyć wzoru bezpośredniego
    return f_direct(x);
}

int main()
{
    // log10(x) x dokładna_wartość
    std::ifstream in("../data.txt");
    std::ofstream out("../output.txt");

    // Do pliku zapisywane w formacie scientific (z e) z precyzją maksymalną danego typu
    out << std::scientific << std::setprecision(std::numeric_limits<ld>::max_digits10);

    std::string line;
    while (std::getline(in, line))
    {
        std::stringstream ss(line);

        ld x;
        ld logx;  // log10(x)
        ld exact; // Dokładna wartość funkcji

        // Jeśli nie w formacie który chce to ignoruje
        if (!( ss >> logx >> x >> exact ))
        {
            continue;
        }

        // Bezpośrednie podstawienie do wzoru f(x) = x^3 / (6 (sinh(x) - x))
        ld fd = f_direct(x);

        // Rozwinięcie funkcji w szereg Taylora stosowane głównie dla małych wartości x unika odejmowania liczb prawie równych
        ld ft = f_taylor(x);

        // Automatycznie wybiera najlepszy sposób obliczeń
        // Taylor dla małych x, metoda bezpośrednia dla większych x
        ld fs = f(x);

        // Błędy względne |wyliczone-dokładne| / dokładne
        ld err_d = fabsl(( fd - exact ) / exact);
        ld err_t = fabsl(( ft - exact ) / exact);
        ld err_s = fabsl(( fs - exact ) / exact);
        // Usuwanie -INF gdzie błąd jest taki niski, że nie ma znaczenia, żeby graf wyglądał fajnie
        // if (ld err_s = fabsl(( fs - exact ) / exact); err_d > 0 && err_t > 0 && err_s > 0)
        // {
        out << logx << " " << log10l(err_d) << " " << log10l(err_t) << " " << log10l(err_s) << "\n";
        // }
    }

    system("gnuplot -p -e \"\
        set xlabel 'log10(x)';\
        set ylabel 'log10(|blad wzgledny|)';\
        set grid;\
        plot '../output.txt' using 1:2 with lines lc rgb '#ff00ff' lw 2 title 'Bezposrednia',\
        '../output.txt' using 1:3 with lines lc rgb '#0066ff' lw 2 title 'Taylor',\
        '../output.txt' using 1:4 with lines lc rgb '#00aa00' lw 3 title 'Auto zmiana przy 1e-0L'\"");
}